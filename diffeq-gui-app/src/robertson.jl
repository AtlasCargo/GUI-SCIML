module RobertsonODE

using OrdinaryDiffEq

export robertson!, get_problem, default_u0, default_params, default_tspan, param_labels, param_explanations

function robertson!(du, u, p, t)
    k1 = 0.04
    k2 = 1e4
    k3 = 3e7
    # Defensive: always resize du to length 3
    if length(du) != 3
        resize!(du, 3)
    end
    y1 = length(u) >= 1 ? u[1] : 0.0
    y2 = length(u) >= 2 ? u[2] : 0.0
    y3 = length(u) >= 3 ? u[3] : 0.0
    du[1] = -k1*y1 + k2*y2*y3
    du[2] =  k1*y1 - k2*y2*y3 - k3*y2^2
    du[3] =  k3*y2^2
end

default_u0 = [1.0, 0.0, 0.0]
default_params = [0.0]  # dummy, not used
default_tspan = (0.0, 1e5)
param_labels = ["(dummy)"]
param_explanations = ["Parameters are hardcoded in the ODE (k1, k2, k3)"]

function get_problem(u0::Vector{<:Real}=default_u0, p::Vector{<:Real}=default_params, tspan::Tuple=default_tspan)
    # Defensive: always pass a copy of u0 with length 3
    u0_fixed = copy(u0)
    if length(u0_fixed) < 3
        append!(u0_fixed, zeros(3 - length(u0_fixed)))
    elseif length(u0_fixed) > 3
        u0_fixed = u0_fixed[1:3]
    end
    ODEProblem(robertson!, u0_fixed, tspan, p)
end

end # module
