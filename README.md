# Two-phase flow using Multi-Point Flux Approximation/Finite Volume Method (2ph-MPFA)

## Description
This repository is a collection of implementations using cell centered finite volume methods, 
specifically Multi-Point Flux Approximation, in order to discretize the incompressible two-phase 
flow equation (two-dimensional, only).

To approximate the solutions at each time step, we can use:
- a classical Newton iterative scheme (a sequential implicity);
- IMplicity Pressure Explicit Saturation scheme (IMPES);
- Streamline approach.
- modified-IMPES. 

Examples of usage can be found in:
- bl2-D.m  
- fivespot.m
- simpleFracture.m

For more information on the MPFA discretization, we refer to:
   I. Aavatsmark: An introduction to multipoint flux approximations for quadrilateral grids, Comput. Geosci. 6(3-4) 405-432.

## Requirements
The modules are compatible with the MATLAB, and it is assumed that all folders are in the MATLAB path. The code has been 
tested with Matlab R2018a and R2019a.

## Cite
If you use 2ph-MPFA, please cite: