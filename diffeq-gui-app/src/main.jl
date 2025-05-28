# filepath: /diffeq-gui-app/diffeq-gui-app/src/main.jl

using Gtk
using OrdinaryDiffEq
include("gui.jl")
include("enzymes.jl")

function main()
    # Initialize the GUI
    win = create_main_window()
    
    # Set up the differential equation parameters
    p_ = [1.2, 1, 0.6, 1, 1.3, 1, 0.7, 1.1, 1]
    u0 = [0.0, 0.0, 0.0, 0.0]
    tspan = (0.0, 7.0)

    # Create the ODE problem
    prob = ODEProblem(enzymes!, u0, tspan, p_)

    # Connect GUI elements to update the parameters and solve the ODE
    connect_gui_signals(win, prob)

    # Show the window
    Gtk.showall(win)

    # Start the main event loop
    Gtk.main()
end

main()