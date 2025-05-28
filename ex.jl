# SciML Tools
using OrdinaryDiffEq, ModelingToolkit, DataDrivenDiffEq, SciMLSensitivity, DataDrivenSparse, SymbolicRegression
using Optimization, OptimizationOptimisers, OptimizationOptimJL, LineSearches, Bumper, MLJ, LinearAlgebra, Statistics
using ComponentArrays, Lux, Zygote, Plots, StableRNGs, SymbolicUtils, Latexify

gr()

# Set a random seed for reproducible behaviour
rng = StableRNG(1111)

function enzymes!(du, u, p, t)
    a, b, c, d0, d1, d2, s1, s2, k = p
    du[1] = a - b * u[1] #mRNA
    du[2] = c * u[1] - d0 * u[2] - s1 * u[2] + k * u[3]  # non-phosphoryleated protein
    du[3] = s1 * u[2] - k * u[3] - d1 * u[3] - s2 * u[3] + k * u[4]   # single phosphorylated protein
    du[4] = s2 * u[3] - k * u[4] - d2 * u[4]   # double phosphorylated protein
end

# Define the experimental parameter
tspan = (0.0, 7.0)

#u0 = 5.0f0 * rand(rng, 2)
u0 = [0.0, 0.0, 0.0, 0.0]
#p_ = [1, 1, 1, 1, 1, 1, 1, 1, 1]
p_ = [1.2, 1, 0.6, 1, 1.3, 1, 0.7, 1.1, 1]

prob = ODEProblem(enzymes!, u0, tspan, p_)
solution = solve(prob, Vern7(), abstol = 1e-12, reltol = 1e-12, saveat = 0.25)

# Add noise in terms of the mean
X = Array(solution)
t = solution.t

x̄ = mean(X, dims = 4)

# Impact on symbolic regression, 5e-3 standard
noise_magnitude = 5e-3

# Only allow positive noise+values
Xₙ = X .+ (noise_magnitude * x̄) .* randn(rng, eltype(X), size(X))

pl_solution = plot(solution, alpha = 0.75, color = :black, label = ["True Data" nothing])
scatter!(t, transpose(Xₙ), color = :red, label = ["Noisy Data" nothing])

# Multilayer FeedForward DNN
#const U = Lux.Chain(Lux.Dense(2, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 2))
rbf(x) = exp.(-(x .^ 2))

# Using one hidden layer per equation +1. Using number of nodes in hidden layers scale with complexity of missing physics. Lookup Jans code

const U = Lux.Chain(
    Lux.Dense(4, 8, rbf),
    Lux.Dense(8, 8, rbf),
    Lux.Dense(8, 8, rbf),
    Lux.Dense(8, 8, rbf),
    Lux.Dense(8, 8, rbf),
    Lux.Dense(8, 8, rbf),
    Lux.Dense(8, 4))

# Get the initial parameters and state variables of the model
p, st = Lux.setup(rng, U)
const _st = st

# Hybrid model
function ude_dynamics!(du, u, p, t, p_true)
    û = U(u, p, _st)[1] # Network prediction
    du[1] =  û[1] - p_true[2]*u[1]
    du[2] =  û[2] + p_true[3]*u[1] - p_true[4]*u[2]
    du[3] =  û[3] - p_true[5]*u[3]
    du[4] =  û[4] - p_true[6]*u[4]
end

# Closure with the known parameter
nn_dynamics!(du, u, p, t) = ude_dynamics!(du, u, p, t, p_)
# Define the problem
prob_nn = ODEProblem(nn_dynamics!, Xₙ[:, 1], tspan, p)

function nn_predict(θ, X = Xₙ[:, 1], T = t)
    _prob = remake(prob_nn, u0 = X, tspan = (T[1], T[end]), p = θ)
    Array(solve(_prob, Vern7(), saveat = T,
        abstol = 1e-6, reltol = 1e-6,
        sensealg = QuadratureAdjoint(autojacvec = ReverseDiffVJP(true))))
end

function loss(θ)
    X̂ = nn_predict(θ)
    mean(abs2, Xₙ .- X̂)
end

losses = Float64[]

callback = function (state, l)
    push!(losses, l)
    if length(losses) % 100 == 0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentVector{Float64}(p))

res1 = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(), callback = callback, maxiters = 6000)
println("Training loss (ADAM) after $(length(losses)) iterations: $(losses[end])")

optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(
    optprob2, LBFGS(linesearch = BackTracking()), callback = callback, maxiters = 1000)
println("Final training loss (LBFGS) after $(length(losses)) iterations: $(losses[end])")

# Rename the best candidate
p_trained = res2.u

# Plot the losses
pl_losses = plot(1:6000, losses[1:6000], yaxis = :log10, xaxis = :log10,
    xlabel = "Iterations", ylabel = "Loss", label = "ADAM", color = :blue)
plot!(6001:length(losses), losses[6001:end], yaxis = :log10, xaxis = :log10,
    xlabel = "Iterations", ylabel = "Loss", label = "LBFGS", color = :red)

## Analysis of the trained network
# Plot the data and the approximation
ts = first(solution.t):(mean(diff(solution.t)) / 2):last(solution.t)
X̂ = nn_predict(p_trained, Xₙ[:, 1], ts)
# Trained on noisy data vs real solution
pl_trajectory = plot(ts, transpose(X̂), xlabel = "t", ylabel = "x(t), y(t)", color = :red,
    label = ["UDE Approximation" nothing])
scatter!(solution.t, transpose(Xₙ), color = :black, label = ["Measurements" nothing])

# Ideal unknown interactions of the predictor

#VectP1 = [p_[1] for i in 1:length(X̂[2, :])]

Ȳ = [(p_[1] )';(- p_[7] * X̂[2, :] .+ p_[9] * X̂[3, :])';(p_[7] * X̂[2, :] .- p_[9] * X̂[3, :] - p_[8] * X̂[3, :] .+ p_[9] * X̂[4, :])';(p_[8] * X̂[3, :] .- p_[9] * X̂[4, :])']


# Neural network guess
Ŷ = U(X̂, p_trained, st)[1]

pl_reconstruction = plot(ts, transpose(Ŷ), xlabel = "t", ylabel = "U(x,y)", color = :red,
    label = ["UDE Approximation" nothing])
plot!(ts, transpose(Ȳ), color = :black, legend =:right, label = ["True Interaction" nothing])

# Plot the error
pl_reconstruction_error = plot(ts, norm.(eachcol(Ȳ - Ŷ)), yaxis = :log, xlabel = "t",
    ylabel = "L2-Error", label = nothing, color = :red)
pl_missing = plot(pl_reconstruction, pl_reconstruction_error, layout = (2, 1))

pl_overall = plot(pl_trajectory, pl_missing)

X1x = (u1=X̂[1,:], u2=X̂[2,:], u3=X̂[3,:], u4=X̂[4,:])
Y1x = (y1 = Ŷ[1,:], y2 = Ŷ[2,:], y3 = Ŷ[3,:], y4 = Ŷ[4,:])

model = MultitargetSRRegressor(
    maxsize=10,
    #niterations=100,
    #ncycles_per_iteration=5000,
    batching=true, # noisy data
    #elementwise_loss=L1DistLoss(),
    constraints=[(*)=>(1, 1)],
    #populations=50,
    #complexity_of_constants=3,
    #model_selection = "accuracy",
    turbo = true,
    should_simplify = true,
    print_precision = 2,
    bumper = true,
    save_to_file = false,
    #precision=64,
    #parsimony=0.01,
    #cluster_manager="slurm",
    #procs=100,
    #niterations=10000000,
    #ncycles_per_iteration=10000,
    binary_operators=[+, -, *]
    #adaptive_parsimony_scaling=1000
    #complexity_of_operators = {"*": 2}
    #unary_operators=[exp],
)

mach = machine(model, X1x, Y1x)
fit!(mach)
rep = report(mach)

# Generic plotting function for SR discovered equations
plot_titles = ["True Eq = 1.2", "True Eq = -0.7u2 + u3", "True Eq = 0.7u2 -2.1 u3 + u4 ", "True Eq = 1.1u3 - u4"]

pl_SR = function (report, machine, X_data, Y_data, titles)
    # Loop over the discovered equations
    for i in range(1, length(report.equations))
        index_array = [1,1,1,1]
        scplot = scatter(xlabel = "Truth", ylabel = "SR Eq Prediction", title = titles[i], legendfont = (10), legend = :best, size = (800,800))
        for j in range(1, length(report.equations[i]))
            index_array[i] = j
            eqn = string(node_to_symbolic(report.equations[i][j]))
            formatted_eqn = replace(eqn, "x"=>"u")
            formatted_dat = "C="*string(report.complexities[i][j])*", S="*string(round(report.scores[i][j], digits=2))*", L="*string(report.losses[i][j])
            formatted_res = latexify(formatted_eqn, fmt = x -> round(x, sigdigits=2))*"\n"*latexify(formatted_dat,fmt = x -> round(x, sigdigits=2))
            #; report.losses[i][j])#, complexity=report.complexities[i][j], score=report.scores[i][j])
            #println(extra)
            #println(i,"\t",j,"\t",round(report.scores[i][j]/report.losses[i][j]/report.complexities[i][j],digits=3),"\t\t\t\t",eqn)
            #if report.losses[i][j] > 1e-3    
            SR_predictions = MLJ.predict(machine, (data=X_data, idx=index_array))
            scatter!(Y_data[i], SR_predictions[i], label=formatted_res)
            #end
        end
        display(scplot)
    end
    return
end

pl_SR(rep, mach, X1x, Y1x, plot_titles)