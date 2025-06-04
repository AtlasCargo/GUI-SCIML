using Plots

# 1. Define the model function
# Takes parameters `p` (a vector, e.g., [theta_0, theta_1]) and independent variable `t`
function linear_model(p, t)
    theta_0, theta_1 = p # Unpack parameters
    return theta_0 + theta_1 * t
end

# 2. Some synthetic experimental data
# (In a real scenario, this would come from an experiment)
t_data = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
y_data = [0.9, 3.1, 5.2, 6.8, 9.1, 11.3] # y = 2*t + 1 with some noise

# 3. Define the Sum of Squared Errors (SSE) objective function
# Takes current parameters `p_current` and the data (t_data, y_data)
function sse_objective(p_current, t_obs, y_obs)
    model_predictions = linear_model.(Ref(p_current), t_obs) # Broadcasting model over t_obs
                                                            # Ref(p_current) keeps p_current as a single entity for broadcasting
    residuals = y_obs .- model_predictions
    sse = sum(residuals.^2)
    return sse
end

# Test with some initial guess for parameters p = [theta_0, theta_1]
initial_guess_p = [0.5, 1.5] # An initial guess for [theta_0, theta_1]
current_sse = sse_objective(initial_guess_p, t_data, y_data)
println("SSE with initial guess ", initial_guess_p, ": ", current_sse)

# True parameters (for comparison, we wouldn't know this in reality)
true_p = [1.0, 2.0]
sse_true_params = sse_objective(true_p, t_data, y_data)
println("SSE with 'true' parameters ", true_p, ": ", sse_true_params) # Should be lower

# Assume theta_0 is fixed at 1.0, we only estimate theta_1
function linear_model_theta1_only(theta_1_val, t)
    theta_0_fixed = 1.0
    return theta_0_fixed + theta_1_val * t
end

function sse_objective_theta1_only(theta_1_val, t_obs, y_obs)
    model_predictions = linear_model_theta1_only.(theta_1_val, t_obs)
    residuals = y_obs .- model_predictions
    sse = sum(residuals.^2)
    return sse
end

# Range of theta_1 values to explore
theta_1_range = 0.5:0.05:3.5 # Explore values for theta_1
sse_values = [] # To store SSE for each theta_1

for th1_val in theta_1_range
    current_sse_th1 = sse_objective_theta1_only(th1_val, t_data, y_data)
    push!(sse_values, current_sse_th1)
end

# Plot the objective function landscape
plll=plot(theta_1_range, sse_values,
     xlabel="Parameter theta_1",
     ylabel="SSE (Objective Function Value)",
     title="Objective Function Landscape for theta_1",
     label="SSE vs. theta_1",
     lw=2, legend=:topright)
vline!([2.0], label="Approx. True theta_1 (for our synthetic data)", ls=:dash, color=:red)
display(plll)
# In a real case, we wouldn't know the "true" value to plot with vline!
# savefig("sse_landscape_1D.png")
println("1D SSE landscape plot generated (conceptually).")