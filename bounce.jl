using DifferentialEquations, Plots

function free_fall!(du, u, p_ff, t) # Renamed p
    g = p_ff[1] # Gravity
    du[1] = u[2]       # dy/dt = v
    du[2] = -g         # dv/dt = -g
end

y0_ball = 10.0  # Initial height, renamed
v0_ball = 0.0   # Initial velocity, renamed
u0_ball_sys = [y0_ball, v0_ball] # Renamed
tspan_ball = (0.0, 10.0)
p_ball_params = [9.81] # Gravity, renamed

# Condition: ball hits the ground (height u[1] is zero)
function condition_ground(u, t, integrator)
    return u[1] # Trigger when height is zero
end

# Affect: reverse velocity (bounce)
function affect_bounce!(integrator)
    e = 0.8 # Coefficient of restitution
    integrator.u[2] = -e * integrator.u[2] # Reverse and dampen velocity
end

# Create the callback
bounce_cb = ContinuousCallback(condition_ground, affect_bounce!)

prob_ball = ODEProblem(free_fall!, u0_ball_sys, tspan_ball, p_ball_params)
sol_ball = solve(prob_ball, Tsit5(), callback=bounce_cb, dtmax=0.01) # dtmax for smoother plot

display(plot(sol_ball, vars=(0,1), xlabel="Time (s)", ylabel="Height (m)", label="Height",title="Bouncing Ball with Callback", lw=2))
plot!(sol_ball, vars=(0,2), label="Velocity", ls=:dash)
ylims!(-15, 15) # Adjust y-limits to see velocity
# savefig("bouncing_ball_plot.png")
println("Bouncing ball plot would be generated here.")
