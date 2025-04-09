
<h1 align="center">
  <br>
  <img src="https://github.com/cld-rostock/Christensen_FailureIndex/raw/main/img/plot.png" alt="Plot" width="240">
  <br>
  <br>
  <b>ChriFT</b>
  <br>
  <sub><sup>Christensen Failure Theory</sup></sub>
  <br>
</h1>

<p align="center">
  Implementation of a linear failure index and corresponding ductility number in accordance with the Christensen criterion for isotropic materials.
</p>


<p align="center">
  <a href="https://github.com/CLD-Rostock/ChriFT/issues">Report a bug</a> · 
  <a href="https://github.com/CLD-rostock/ChriFt/">Cite the software</a> · 
  <a href="https://www.cld.uni-rostock.de">Read the paper</a>
  <br>
  <br>
  <br>
</p>


<!-- TABLE OF CONTENTS -->
## Contents
1. [About the project](#about-the-project)
2. [Choose the right implementation](#choose-the-right-implementatioin)
3. [Python package](#python-package)
4. [Python GUI Plugin](#python-gui-plugin)
5. [Fortran Subroutine](#fortran-subroutine)
6. [Contact](#contact)


<!-- ABOUT THE PROJECT -->
## About the project


💡 
This project implements a linear failure index for the Christensen criterion — a stress-based failure model for isotropic materials ranging from fully brittle to fully ductile, including the classical limit cases of maximum normal stress and von Mises.

📐 The index is based on projecting the stress state onto the Christensen failure surface at constant angles in spherical coordinates (principal stress space).

🚀 Made available in this repository in 3 implementations:

1. 🧩 ABAQUS CAE postprocessing plugin
2. 🧾 ABAQUS Fortran UMAT subroutine
3. 🐍 Python module for stress analysis

📄 The method is described in the paper (under review). If you use this code, please cite:
```
Hach, M., Radtke, A., Weissgraeber, P. (2025). A linear failure index for the Christensen criterion
```
when using this code.

<img src="https://github.com/CLD-Rostock/Christensen_FailureIndex/raw/main/img/visualization-algorithm.png" alt="algorithm-visualization" width="200"/>


For more background info on the Christensen failure criterion please refer to the following literature:

- [The Theory of Materials Failure](https://www.failurecriteria.com/newbook.html)
- [Failure Mechanics—Part I: The Coordination Between Elasticity Theory and Failure Theory for all Isotropic Materials](https://doi.org/10.1115/1.4027753)
- [Failure Mechanics—Part II: The Central and Decisive Role of Graphene in Defining the Elastic and Failure Properties for all Isotropic Materials ](https://doi.org/10.1115/1.4028407)


<!-- CHOOSE THE RIGHT IMPLEMENTATION -->
## Choose the right implementation

Three different implementations are available


  1. Using the Python GUI Plugin for Simulia Abaqus (tested for V2021)
  - Use for linear elements
  - Use for non linearelastic material behavior
  - Use for beginner users
  2. Using the Fortran Subroutine for Simulia Abaqus (tested for V2021)
  - Use for nonlinear elements
  - A Fortran Compiler linked to Abaqus required
  3. Using the Python package for any list of stress states
  - Use for any other FE Codes or analytical problems

<!-- PYTHON GUI PLUGIN -->
## Python GUI Plugin

### Preparation

1. Download and store the folder **Plugin** in the **abaqus_plugins** folder of your local Abaqus installation
2. Restart Abaqus and the Plugin should be visible in the Plug-ins drop-down menu


### Usage
1. Run the simulation
2. Open the Plugin
3. Choose the .odb file
4. Choose the Material, for which you want to calculate the failure index


<!-- FORTRAN SUBROUTINE -->
## Fortran Subroutine

### Preparation

1. Check, that your abaqus installation is connected to a fortran compiler by executing Abaqus Verification
2. Download and store the [subroutine](https://github.com/CLD-Rostock/Christensen_FailureIndex/blob/main/Subroutine/christensen_subroutine.for)
3. During simulation setup, define the material as

| #  | Property               |
|----|------------------------|
| 1  | Young's Modulus        |
| 2  | Poisson's Ratio        |
| 3  | Tensile Strength       |
| 4  | Compressive Strength   |

  4. Activate SDV's in a Field Output Request
  5. Select the subroutine in the job definition
  6. Run the job
  7. Once the job is completed, the failure index can be viewed as SDV1 and the Failure Number as SDV2


<!-- PYTHON PACKAGE -->
## Python package

This package implements the Christensen failure criterion, providing methods to compute the failure index and failure number based on a given stress state and material properties. It supports 3D, plane strain, and plane stress conditions.

### Features

- Compute failure index based on Christensen’s theory
- Determine failure number to quantify material damage
- Principal stress calculation & polar coordinate transformation
- Handles both tensile (T) and compressive (C) strength asymmetry

### Usage
```python
from Christensen_FailureIndex import Christensen_class

# Define stress tensor components & material properties
stress_tensor = [100, 50, 30, 20, 10, 5]  # 3D stress state
christensen_inst = Christensen_class(stress_tensor, T=100, C=300)

# Compute failure index and failure number
result = christensen_inst.calc_main()

# Access data option a
print(f"Failure Index: {christensen_inst.failure_index}")
print(f"Failure Number: {christensen_inst.Fn}")

# Access data option b
print(f"Failure Index: {result[0]}")
print(f"Failure Number: {result[1]}")
```
### Requirements

- Python 3.x
- NumPy

## Contact
E-mail: mathis.hach@uni-rostock.de - Web: https://www.cld.uni-rostock.de
