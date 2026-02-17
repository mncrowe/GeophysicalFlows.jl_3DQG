# 3D Quasi-Geostrophy Solver

This repository contains scripts which implement the 3D quasi-geostrophic equations in the FourierFlows.jl framework. This code may eventually be converted to a module and merged with GeophysicalFlows.jl. The method uses spectral collocation with Chebyshev points in the vertical direction in order to allow top and bottom boundaries.

This code is currently set up to evolve an SQG modon on a beta plane using [`QGDipoles.jl`](https://github.com/mncrowe/QGDipoles.jl).

## Files:

- `3D_QG.jl`: implements the 3D QG model
- `SQG_modon.jl`: defines functions to run a simulation of an SQG modon on a beta plane in 3D QG
- `run.jl`: runs SQG modon simulations by calling the functions from `SQG_modon.jl`
