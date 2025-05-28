using GLMakie
using GLMakie: Menu, Button, Slider, Label, DataAspect, colsize!, GridLayout, xlims!, ylims!
include("enzymes.jl")
include("logistic.jl")
include("lotkavolterra.jl")
using .Enzymes
using .LogisticODE
using .LotkaVolterraODE

dict_odes = Dict(
    "Enzymes" => (
        get_problem = (u0, p, tspan) -> Enzymes.get_problem(u0, p, tspan),
        default_u0 = Enzymes.default_u0,
        default_params = Enzymes.default_params,
        default_tspan = Enzymes.default_tspan,
        param_labels = ["a", "b", "c", "d0", "d1", "d2", "s1", "s2", "k"],
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
    ),
    "Logistic" => (
        get_problem = (u0, p, tspan) -> LogisticODE.get_problem(u0, p, tspan),
        default_u0 = LogisticODE.default_u0,
        default_params = LogisticODE.default_params,
        default_tspan = LogisticODE.default_tspan,
        param_labels = LogisticODE.param_labels,
        param_explanations = LogisticODE.param_explanations
    ),
    "Lotka-Volterra" => (
        get_problem = (u0, p, tspan) -> LotkaVolterraODE.get_problem(u0, p, tspan),
        default_u0 = LotkaVolterraODE.default_u0,
        default_params = LotkaVolterraODE.default_params,
        default_tspan = LotkaVolterraODE.default_tspan,
        param_labels = LotkaVolterraODE.param_labels,
        param_explanations = LotkaVolterraODE.param_explanations
    )
)

# Helper to create sliders for each parameter
function gui_app()
    fig = Figure(size = (1000, 700))
    main_grid = fig[1, 1] = GridLayout()
    # Left column: controls (dropdown, sliders, button, explanations)
    controls_grid = main_grid[1, 1] = GridLayout()
    # Center: plot
    plot_grid = main_grid[1, 2] = GridLayout()
    # Right: explanations
    explanation_grid = main_grid[1, 3] = GridLayout()

    ode_names = collect(keys(dict_odes))
    ode_select = Menu(controls_grid[1, 1], options = ode_names, width = 200)

    # Find the maximum number of parameters among all ODEs
    maxparams = maximum(length(ode.default_params) for ode in values(dict_odes))
    nparams = length(dict_odes[ode_names[1]].default_params)
    ax = GLMakie.Axis(plot_grid[1, 1]; xlabel = "t", ylabel = "u")

    # Sliders and labels (allocate for maxparams)
    sliders = Vector{Slider}(undef, maxparams)
    slider_labels = Vector{Label}(undef, maxparams)
    for i in 1:maxparams
        # Use default value and label from the first ODE, or fallback if not enough
        val = i <= nparams ? dict_odes[ode_names[1]].default_params[i] : 1.0
        lab = i <= nparams ? dict_odes[ode_names[1]].param_labels[i] : "param$i"
        sliders[i] = Slider(controls_grid[i+1, 1], range = 0.0:0.05:2.0, startvalue = val)
        slider_labels[i] = Label(controls_grid[i+1, 2], lab, tellwidth=false)
        # Hide unused by default
        if i > nparams
            try
                sliders[i].parent.visible[] = false
            catch
            end
            try
                slider_labels[i].parent.visible[] = false
            catch
            end
        end
    end
    btn = Button(controls_grid[maxparams+2, 1], label = "Update Simulation")
    params = Observable(copy(dict_odes[ode_names[1]].default_params))

    # Explanations (allocate for maxparams)
    param_explanations = dict_odes[ode_names[1]].param_explanations
    explanation_labels = Vector{Label}(undef, maxparams)
    for i in 1:maxparams
        txt = i <= nparams ? param_explanations[i] : ""
        explanation_labels[i] = Label(explanation_grid[i, 1], txt)
        if i > nparams
            try
                explanation_labels[i].parent.visible[] = false
            catch
            end
        end
    end

    # Set column widths
    colsize!(main_grid, 1, Relative(0.18))
    colsize!(main_grid, 2, Relative(0.62))
    colsize!(main_grid, 3, Relative(0.2))

    # Make the plot box square regardless of window size
    ax.aspect = DataAspect()
    ax.width = 600
    ax.height = 600

    # Initial plot
    sol = dict_odes[ode_names[1]].get_problem(params[], params[], dict_odes[ode_names[1]].default_tspan)
    sol = solve(sol)
    # Clear previous lines (Makie 0.22: use empty!(ax) to clear axis)
    empty!(ax)
    if isa(sol.u[1], AbstractArray)
        for i in 1:length(sol.u[1])
            lines!(ax, sol.t, [u[i] for u in sol.u], label = "u$i")
        end
    else
        lines!(ax, sol.t, sol.u, label = "u")
    end
    GLMakie.xlims!(ax, dict_odes[ode_names[1]].default_tspan)
    ypad = dict_odes[ode_names[1]].default_tspan[2] - dict_odes[ode_names[1]].default_tspan[1]
    ycenter = 0.5 * (maximum(vcat(sol.u...)) + minimum(vcat(sol.u...)))
    GLMakie.ylims!(ax, ycenter - ypad/2, ycenter + ypad/2)

    # Update plot on button press
    on(btn.clicks) do _
        sol = dict_odes[ode_select.selection[]].get_problem(params[], params[], dict_odes[ode_select.selection[]].default_tspan)
        sol = solve(sol)
        empty!(ax)
        if isa(sol.u[1], AbstractArray)
            for i in 1:length(sol.u[1])
                lines!(ax, sol.t, [u[i] for u in sol.u], label = "u$i")
            end
        else
            lines!(ax, sol.t, sol.u, label = "u")
        end
        GLMakie.xlims!(ax, dict_odes[ode_select.selection[]].default_tspan)
        ypad = dict_odes[ode_select.selection[]].default_tspan[2] - dict_odes[ode_select.selection[]].default_tspan[1]
        ycenter = 0.5 * (maximum(vcat(sol.u...)) + minimum(vcat(sol.u...)))
        GLMakie.ylims!(ax, ycenter - ypad/2, ycenter + ypad/2)
    end

    # Update everything when ODE system is changed
    on(ode_select.selection) do name
        ode = dict_odes[name]
        # Update sliders and labels
        for i in 1:length(ode.default_params)
            sliders[i].value[] = ode.default_params[i]
            slider_labels[i].text[] = ode.param_labels[i]
            explanation_labels[i].text[] = ode.param_explanations[i]
        end
        # Hide unused sliders/labels if system has fewer params
        for i in (length(ode.default_params)+1):maxparams
            if i <= length(sliders) && i <= length(slider_labels) && i <= length(explanation_labels)
                try
                    if !isnothing(sliders[i]) && hasproperty(sliders[i], :parent) && !isnothing(sliders[i].parent)
                        sliders[i].parent.visible[] = false
                    end
                catch; end
                try
                    if !isnothing(slider_labels[i]) && hasproperty(slider_labels[i], :parent) && !isnothing(slider_labels[i].parent)
                        slider_labels[i].parent.visible[] = false
                    end
                catch; end
                try
                    if !isnothing(explanation_labels[i]) && hasproperty(explanation_labels[i], :parent) && !isnothing(explanation_labels[i].parent)
                        explanation_labels[i].parent.visible[] = false
                    end
                catch; end
            end
        end
        # Show used sliders/labels
        for i in 1:length(ode.default_params)
            if i <= length(sliders) && i <= length(slider_labels) && i <= length(explanation_labels)
                try
                    if !isnothing(sliders[i]) && hasproperty(sliders[i], :parent) && !isnothing(sliders[i].parent)
                        sliders[i].parent.visible[] = true
                    end
                catch; end
                try
                    if !isnothing(slider_labels[i]) && hasproperty(slider_labels[i], :parent) && !isnothing(slider_labels[i].parent)
                        slider_labels[i].parent.visible[] = true
                    end
                catch; end
                try
                    if !isnothing(explanation_labels[i]) && hasproperty(explanation_labels[i], :parent) && !isnothing(explanation_labels[i].parent)
                        explanation_labels[i].parent.visible[] = true
                    end
                catch; end
            end
        end
        params[] = copy(ode.default_params)
        notify(params)
        # Update plot
        sol = ode.get_problem(params[], params[], ode.default_tspan)
        sol = solve(sol)
        empty!(ax)
        # Defensive: handle empty solution and avoid passing ODESolution to lines!
        if length(sol.u) > 0
            if isa(sol.u[1], AbstractArray)
                try
                    arr = reduce(hcat, sol.u)
                    for i in 1:size(arr, 1)
                        yvals = arr[i, :]
                        if length(sol.t) == length(yvals) && eltype(yvals) <: Number
                            lines!(ax, sol.t, collect(yvals), label = "u$i")
                        end
                    end
                catch err
                    @warn "Error plotting ODE solution (matrix/array case)" exception=err
                end
            elseif eltype(sol.u) <: Number
                try
                    if length(sol.t) == length(sol.u)
                        lines!(ax, sol.t, sol.u, label = "u")
                    end
                catch err
                    @warn "Error plotting ODE solution (vector case)" exception=err
                end
            else
                @warn "Unrecognized ODE solution type for plotting" typeof_sol_u=typeof(sol.u)
            end
        end
        GLMakie.xlims!(ax, ode.default_tspan)
        ypad = ode.default_tspan[2] - ode.default_tspan[1]
        ycenter = 0.5 * (maximum(vcat(sol.u...)) + minimum(vcat(sol.u...)))
        GLMakie.ylims!(ax, ycenter - ypad/2, ycenter + ypad/2)
    end

    fig
end

gui_app()