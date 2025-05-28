module LogisticODE

using OrdinaryDiffEq

export logistic!, get_problem, default_params, default_u0, default_tspan, param_labels, param_explanations

default_u0 = [1.0]
default_params = [1.5, 10.0]  # a, K
default_tspan = (0.0, 10.0)
param_labels = ["a", "K"]
param_explanations = [
    "a: growth rate",
    "K: carrying capacity"
]

function logistic!(du, u, p, t)
    a, K = p
    du[1] = a * u[1] * (1 - u[1]/K)
end

function get_problem(u0::Vector{<:Real}=default_u0, p::Vector{<:Real}=default_params, tspan::Tuple=default_tspan)
    ODEProblem(logistic!, u0, tspan, p)
end

end # module
