module LotkaVolterraODE

using OrdinaryDiffEq

export lotkavolterra!, get_problem, default_params, default_u0, default_tspan, param_labels, param_explanations

default_u0 = [1.0, 3.0]
default_params = [1.5, 1.0, 0.5, 1.0]  # a, b, c, d
default_tspan = (0.0, 10.0)
param_labels = ["a", "b", "c", "d"]
param_explanations = [
    "a: prey growth rate",
    "b: predation rate",
    "c: predator reproduction rate",
    "d: predator death rate"
]

function lotkavolterra!(du, u, p, t)
    x, y = u
    a, b, c, d = p
    du[1] = a*x - b*x*y
    du[2] = c*x*y - d*y
end

function get_problem(u0::Vector{<:Real}=default_u0, p::Vector{<:Real}=default_params, tspan::Tuple=default_tspan)
    ODEProblem(lotkavolterra!, u0, tspan, p)
end

end # module
