# filepath: /diffeq-gui-app/diffeq-gui-app/src/enzymes.jl
module Enzymes

using OrdinaryDiffEq, StableRNGs

export enzymes!, run_simulation, default_params, default_u0, default_tspan

# ODE system
default_u0 = [0.0, 0.0, 0.0, 0.0]
default_params = [1.2, 1, 0.6, 1, 1.3, 1, 0.7, 1.1, 1]
default_tspan = (0.0, 7.0)

function enzymes!(du, u, p, t)
    a, b, c, d0, d1, d2, s1, s2, k = p
    du[1] = a - b * u[1]
    du[2] = c * u[1] - d0 * u[2] - s1 * u[2] + k * u[3]
    du[3] = s1 * u[2] - k * u[3] - d1 * u[3] - s2 * u[3] + k * u[4]
    du[4] = s2 * u[3] - k * u[4] - d2 * u[4]
end

function run_simulation(u0::Vector{<:Real}=default_u0, p::Vector{<:Real}=default_params, tspan::Tuple=default_tspan)
    prob = ODEProblem(enzymes!, u0, tspan, p)
    solve(prob, Vern7(), abstol=1e-12, reltol=1e-12, saveat=0.25)
end

end # module