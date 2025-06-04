#We need data for x(t) (prey) and y(t) (predator). 
#We can generate this with DifferentialEquations.jl.
#We need x˙(t) and y˙(t).
#Step 1: Generate/Load Data
using DifferentialEquations 
using DataDrivenDiffEq 
using Plots
using ModelingToolkit 
using SymbolicRegression
using DataDrivenSparse
using StableRNGs
using LinearAlgebra
using OrdinaryDiffEq


# For @variables
# using LinearAlgebra 
# For sparse regression, if done manually

# Generate Lotka-Volterra data (from Day 3)
function lotka_volterra_sindy!(du, u, p, t) # Renamed
    a, β, δ, g = p
    x_s, y_s = u # Renamed state vars
    du[1] = a*x_s - β*x_s*y_s
    du[2] = δ*x_s*y_s - g*y_s
end
u0_lv_sindy = [2.0, 1.0]
p_lv_true_sindy = [1.5, 1.0, 0.8, 1.0] # α, β, δ, γ
tspan_data_sindy = (0.0, 20.0)
dt_save_sindy = 0.1
prob_data_sindy = ODEProblem(lotka_volterra_sindy!, u0_lv_sindy, tspan_data_sindy, p_lv_true_sindy)
data_sol_sindy = solve(prob_data_sindy, Tsit5(), saveat=dt_save_sindy)

X_sindy = Array(data_sol_sindy)
t_sindy = data_sol_sindy.t
println("SInDy Data generated. Size of X_sindy: ", size(X_sindy))

#Step 2: Estimate Derivatives (DataDrivenDiffEq can help)
#DataDrivenDiffEq.jl can often handle this internally or provide tools.
#For this example, we can also extract true derivatives if the generating model is known (for testing).
X_derivs_sindy = similar(X_sindy)
for i in 1:length(t_sindy)
    # Call the true ODE function to get derivatives at each data point
    # Note: lotka_volterra_sindy! modifies its first argument in-place
    # We need a temporary array for du for each call.
    du_temp = similar(u0_lv_sindy)
    lotka_volterra_sindy!(du_temp, X_sindy[:, i], p_lv_true_sindy, t_sindy[i])
    X_derivs_sindy[:, i] = du_temp
end
println("Derivatives for SInDy estimated/retrieved.")

# Step 3, 4 & 5: Using DataDrivenDiffEq.jl for SInDy
# Create a DataDrivenProblem
# We provide X (states), t (time), and DX (derivatives)

prob_dd_sindy = ContinuousDataDrivenProblem(X_sindy, t_sindy, DX = X_derivs_sindy)

# Define a basis of candidate functions using ModelingToolkit variables
@variables (u_sindy(t_sindy[1]))[1:size(X_sindy, 1)] # Define symbolic variables, u_sindy[1] is x, u_sindy[2] is y

# Create a basis with polynomial terms up to degree 2.
# This will include u_sindy[1], u_sindy[2], u_sindy[1]^2, #u_sindy[2]^2, u_sindy[1]*u_sindy[2]
# and a constant term (1) by default.
basis_sindy = Basis(polynomial_basis(u_sindy, 2), u_sindy)
println("SInDy Basis created with ", length(basis_sindy), " terms.")

# Choose a sparse optimizer (STLSQ - Sequentially Thresholded Least Squares)
# The threshold is a hyperparameter and may need tuning.
# It can be a single value or a vector of values to try.
# kommentar

optimizer_sindy = STLSQ(exp10.(-10:0.5:0)) # Example: thresholds from 1e-5 to 1

println("Attempting SInDy with DataDrivenDiffEq.jl...")
# Perform the sparse regression to find the model.
# `maxiter` might be needed for convergence.
# `denoise` and `normalize` can be helpful for real data.
#sindy_options = DataDrivenCommonOptions(maxiter = 10000, denoise = false, normalize = true, verbose=false)
#sindy_options = DataDrivenCommonOptions(data_processing = sampler, digits = 1)
sindy_result = DataDrivenDiffEq.solve(prob_dd_sindy, basis_sindy, optimizer_sindy)

#discovered_system = result(sindy_result) 
# Get the identified system structure

system = get_basis(sindy_result)
params = get_parameter_map(system)

println("\nDiscovered System Equations (SInDy):")
# Print the equations from the `discovered_system`
# This depends on how `DataDrivenDiffEq` formats its output.
# print_equations(discovered_system, show_parameter = true)
println(system) # A basic print

# We can also get the coefficient matrix Xi
Xi_discovered = params #coefficients(system)
println("\nDiscovered Coefficient Matrix (Xi):")
display(Xi_discovered) 
# Rows = basis functions, Columns = state variables

