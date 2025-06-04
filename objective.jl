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