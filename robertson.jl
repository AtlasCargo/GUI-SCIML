# Robertson's problem (stiff ODE system)
function robertson!(du, u, p, t)
    # Hardcoded parameter values
    k1 = 0.04
    k2 = 1e4
    k3 = 3e7

    y1, y2, y3 = u

    du[1] = -k1*y1 + k2*y2*y3
    du[2] =  k1*y1 - k2*y2*y3 - k3*y2^2
    du[3] =  k3*y2^2
end

# ICs and integration time span
u0_rob = [1.0, 0.0, 0.0]
tspan_rob = (0.0, 1e5)

# Define ODE problem (parameters hardcoded, so `p` not needed)
prob_rob_stiff = ODEProblem(robertson!, u0_rob, tspan_rob)

# Attempt solving with a non-stiff solver (Tsit5) — usually slow or fails
# Uncomment the following lines if you want to try the non-stiff solver:
# println("Solving Robertson with Tsit5 (non-stiff solver)...")
# @time sol_rob_nonstiff = solve(prob_rob_stiff, Tsit5())

# Solve using a stiff solver
println("Solving Robertson with Rodas5 (stiff solver)...")
@time sol_rob_stiff = solve(prob_rob_stiff, Rodas5())
println("Solution with Rodas5 successful.")

# Plot a shorter timespan to show early dynamics clearly
plot(sol_rob_stiff, tspan=(0.0, 1.0),
     xlabel="Time", ylabel="Concentration",
     title="Robertson's Kinetics",
     label=["y₁" "y₂" "y₃"],
     legend=:right)
