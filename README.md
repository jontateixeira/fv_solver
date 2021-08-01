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

For more information on the MPFA discretization, we refer to:
 * I. Aavatsmark: An introduction to multipoint flux approximations for quadrilateral grids, Comput. Geosci. 6(3-4) 405-432.
 * Gao, Z. and Wu, J. (2011), A linearity-preserving cell-centered scheme for the heterogeneous and anisotropic diffusion equations on general meshes. Int. J. Numer. Meth. Fluids, 67: 2157-2183. https://doi.org/10.1002/fld.2496
 * F.R.L. Contreras, P.R.M. Lyra, M.R.A. Souza, D.K.E. Carvalho (2016), A cell-centered multipoint flux approximation method with a diamond stencil coupled with a higher order finite volume method for the simulation of oil–water displacements in heterogeneous and anisotropic petroleum reservoirs. Computers & Fluids, Volume 127, Pages 1-16.
 * J. C. Teixeira, L. J. N. Guimarães, D. K. E. Carvalho (2021), Streamline-based simulation in highly heterogeneous and anisotropic petroleum reservoirs using a non-orthodox MPFA method and an adaptive timestep strategy with unstructured meshes, Jour. Petr. Sci. and Eng., Volume 201.

## Requirements
The modules are compatible with the MATLAB, and it is assumed that all folders are in the MATLAB path. The code has been 
tested with Matlab R2018a and R2019a.

## Cite
If you use 2ph-MPFA, please cite: