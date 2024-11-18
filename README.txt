
# Vortex Panel Project

## Introduction

"VortexPanelProject.py" is a python script written by Spencer Kenison
to estimate the coefficients of lift, leading edge moment, and quarter-chord 
moment for an airfoil using the vortex panel method. This code was written 
for MAE 5500 at Utah State University.

## Dependencies

The script was written using Python 3.11 interpreter but any similar 
version should function normally.

The script requires the following libraries to be installed:
- numpy
- matplotlib.pyplot
using PIP or other similar package-management system.

## Running program

Run the "VortexPanelProject.py" file using a Python interpreter. Can be run
in terminal using "python CFD_Project3.py" command or from IDE such
as VSCode. The python program will estimate the coefficients of lift, 
leading edge moment, and quarter-chord moment for an airfoil using the 
vortex-panel method. 

The code references the input.json file to know the desired input
airfoil geometry file, freestream velocity, and angle(s) of attack for analysis.
The program can take an array of angle(s) of attack, but currently 
only supports one freestream velocity. 
