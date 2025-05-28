using OrdinaryDiffEq, Plots
# Define the Logistic model
default_K = 10.0
function logistic!(du, u, p, t)
    a, K = p
    du[1] = a * u[1] * (1 - u[1]/K)
end

# Initial conditions and parameters
u0 = [1.0]
p = [1.5, default_K]  # a, K
tspan = (0.0, 10.0)

# Create an ODEProblem
prob = ODEProblem(logistic!, u0, tspan, p)

# Solve the problem
sol = solve(prob)

# Plot the results
plot(sol, label=["Population"], xlabel="Time", ylabel="Population", title="Logistic Growth Model")
display(plot(sol, xlabel="Time", ylabel="Population", title="Logistic Growth Model"))
savefig("logistic_growth.png")