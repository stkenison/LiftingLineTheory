# LiftingLineTheory

## Introduction
This is a Python script that utilizes lifting line theory to calculate the aerodynamic characteristics of a wing. It was developed as part of the coursework for MAE 5500 at Utah State University. The program reads wing configuration parameters from an input JSON file, performs numerical calculations based on lifting line theory, and outputs aerodynamic coefficients such as lift coefficient, induced drag coefficient, rolling moment coefficient, and yawing moment coefficient.

## Features
- Computes aerodynamic properties of a wing using lifting line theory.
- Reads input parameters from a structured JSON file.
- Supports various wing configurations including aspect ratio, taper ratio, washout distribution, and aileron effects.
- Outputs key aerodynamic coefficients and writes results to an output text file.
- Includes visualization capabilities for lift distribution and wing planform.

## Dependencies
This program requires Python 3.11 (or a similar version) and the following Python libraries:

- numpy
- matplotlib
- json

These dependencies can be installed using pip:
```sh
pip install numpy matplotlib
```

## Installation
Clone the repository from GitHub:
```sh
git clone https://github.com/yourusername/LiftingLineTheory.git
cd LiftingLineTheory
```

## Running the Program
1. Ensure that Python and the required dependencies are installed.
2. Prepare the `input.json` file with the required wing configuration parameters.
3. Run the `main.py` file using the command:
   ```sh
   python main.py
   ```
4. The program will generate an `output.txt` file containing the calculated aerodynamic coefficients.
5. If enabled in the input file, the program will also generate graphical plots of the wing planform and lift distribution.

## Input File Format
The program reads an input JSON file (`input.json`) containing the wing configuration parameters. Below is an example structure:
```json
{
    "wing": {
        "planform": {
            "aspect_ratio": 6.0,
            "taper_ratio": 0.5
        },
        "nodes_per_semispan": 20,
        "airfoil_lift_slope": 5.7,
        "washout": {
            "distribution": "linear",
            "amount[deg]": 2.0
        },
        "aileron": {
            "begin[z/b]": 0.3,
            "end[z/b]": 0.7,
            "begin[cf/c]": 0.2,
            "end[cf/c]": 0.3,
            "hinge_efficiency": 0.8
        }
    },
    "condition": {
        "alpha_root[deg]": 5.0,
        "aileron_deflection[deg]": 10.0
    },
    "view": {
        "planform": true
    }
}
```

## Output
After execution, the program produces:
- `output.txt` containing the calculated aerodynamic coefficients and intermediate matrix values.
- Printed results of aerodynamic coefficients in the terminal.
- Optional plots of the wing planform and lift distribution if enabled.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.

## Contributing
If you would like to contribute to this project, feel free to fork the repository and submit a pull request with your improvements.

## Author
Spencer Kenison

For any questions or issues, feel free to open an issue on GitHub.

