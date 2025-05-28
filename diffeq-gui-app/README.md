# Differential Equation GUI Application

This project provides a graphical user interface (GUI) for interacting with a differential equation system. Users can adjust parameters in real-time and visualize the effects on the system's behavior.

## Project Structure

- `src/main.jl`: Entry point for the application. Initializes the GUI and sets up the main event loop.
- `src/gui.jl`: Contains definitions for GUI components, including sliders and buttons for parameter adjustments.
- `src/enzymes.jl`: Defines the differential equation system and includes functions to update the system based on user input.
- `src/utils.jl`: Utility functions for data processing, plotting, and other helper methods.
- `Project.toml`: Configuration file specifying project dependencies and their versions.
- `Manifest.toml`: Complete list of all dependencies and their specific versions for reproducibility.

## Setup Instructions

1. Ensure you have Julia installed on your machine.
2. Clone the repository:
   ```
   git clone <repository-url>
   ```
3. Navigate to the project directory:
   ```
   cd diffeq-gui-app
   ```
4. Activate the project environment:
   ```
   julia --project=.
   ```
5. Install the required dependencies:
   ```
   using Pkg
   Pkg.instantiate()
   ```

## Usage Guidelines

To run the application, execute the following command in the Julia REPL:

```julia
include("src/main.jl")
```

This will launch the GUI, allowing you to interact with the differential equation system. Use the provided sliders and buttons to adjust parameters and observe the changes in real-time.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any suggestions or improvements.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.