using GLMakie
include("enzymes.jl")
using .Enzymes

# Helper to create sliders for each parameter
default_labels = ["a", "b", "c", "d0", "d1", "d2", "s1", "s2", "k"]

function gui_app()
    fig = Figure(size = (1000, 700))

    nparams = length(Enzymes.default_params)
    # --- LAYOUT FIX: Place sliders in a column to the left, plot in the center, explanations to the right ---
    # Create a new grid layout for better control
    gl = fig[1, 1] = GridLayout()
    # Place the axis (plot) in the center column, spanning all rows
    ax = GLMakie.Axis(gl[1:nparams, 2]; xlabel = "t", ylabel = "u")

    # Sliders for each parameter and their labels
    sliders = Vector{Slider}(undef, nparams)
    for i in 1:nparams
        sliders[i] = Slider(gl[i, 1], range = 0.0:0.05:2.0, startvalue = Enzymes.default_params[i])
        Label(gl[i, 1], default_labels[i], tellwidth=false)
    end

    # Button to update
    btn = Button(gl[nparams+1, 1], label = "Update Simulation")

    # Observable for parameters
    params = Observable(copy(Enzymes.default_params))

    # Update params when sliders move
    for (i, s) in enumerate(sliders)
        on(s.value) do v
            params[][i] = v
            notify(params)
        end
    end

    # Initial plot
    sol = Enzymes.run_simulation(Enzymes.default_u0, params[], Enzymes.default_tspan)
    for i in 1:length(sol.u[1])
        lines!(ax, sol.t, [u[i] for u in sol.u], label = "u$i")
    end

    # Update plot on button press
    on(btn.clicks) do _
        sol = Enzymes.run_simulation(Enzymes.default_u0, params[], Enzymes.default_tspan)
        ax.plots[] = () # clear
        for i in 1:length(sol.u[1])
            lines!(ax, sol.t, [u[i] for u in sol.u], label = "u$i")
        end
        autolimits!(ax)
    end

    # Add parameter explanations
    param_explanations = [
        "a: mRNA production rate",
        "b: mRNA degradation rate",
        "c: protein translation rate",
        "d0: non-phosphorylated protein degradation rate",
        "d1: single phosphorylated protein degradation rate",
        "d2: double phosphorylated protein degradation rate",
        "s1: phosphorylation rate 1",
        "s2: phosphorylation rate 2",
        "k: dephosphorylation rate"
    ]
    # Place explanations in the third column
    for (i, txt) in enumerate(param_explanations)
        Label(gl[i, 3], txt)
    end

    # Set column widths for better proportions
    colsize!(gl, 1, Relative(0.18))
    colsize!(gl, 2, Relative(0.62))
    colsize!(gl, 3, Relative(0.2))

    # Make the plot box square regardless of window size
    ax.aspect = DataAspect()
    ax.width = 600
    ax.height = 600
    # Pad y-limits to match x-limits for a square data box
    GLMakie.xlims!(ax, Enzymes.default_tspan)
    ypad = Enzymes.default_tspan[2] - Enzymes.default_tspan[1]
    ycenter = 0.5 * (maximum(vcat(sol.u...)) + minimum(vcat(sol.u...)))
    GLMakie.ylims!(ax, ycenter - ypad/2, ycenter + ypad/2)

    fig
end

gui_app()