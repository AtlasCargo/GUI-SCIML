using OrdinaryDiffEq, Plots
# Define the Lotka-Volterra model
function lotkavolterra!(du, u, p, t)
    x, y = u
    a, b, c, d = p 
    du[1] = a*x-b*x*y #Prey
    du[2] = c*x*y- d*y # Predator
end

# Initial conditions and parameters
u0 = [1.0 , 3.0]
p = [1.5, 1.0, 0.5, 1.0]
tspan = (0.0, 10.0)

# Create an ODEProblem
prob = ODEProblem(lotkavolterra!, u0, tspan, p)

# Solver the problem
sol = solve(prob)

# Plot the results
plt1 = plot(sol, label=["Prey" "Predator"], xlabel="Time", ylabel="Population", title="Lotka-Volterra Model")
plt2 = plot(sol, vars=(1,2), xlabel="Prey", ylabel="Predator", title="Phase Plane")
# Combine both plots side by side
combined = plot(plt1, plt2, layout=(1,2), size=(900,400))
display(combined)
savefig(combined, "lotkavolterra_combined.png")