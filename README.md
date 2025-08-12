# EuropeanPG-V2



## Model Description

This repository contains code for performing non-normality analysis of the European power grid, along with the associated dataset. For a deeper explanation of the methods and theoretical background, see Ref. [1]. The workflow is built on open-source MATLAB tools: MATPOWER [2] for steady-state power flow calculations. We recommend MATLAB R2024a to ensure smooth execution. 


## Usage

- `Cal_nonnormal.m` : This script performs the non-normality analysis of the European power grid. It first calls `build_model.m` to construct the linearized system representation, then calculates the non-normality index by comparing the spectral properties of the system matrix and its symmetric part. The script proceeds to extract the first k dominant modes, compute and sum the modal reactivity for each bus, and finally generates a plot showing the cumulative reactivity distribution in ascending bus order.
-
-
- fault initialization, parameter setup, and dynamic simulation procedures.
- `build_model.m` :
- `EUR_V2.mat` : Power grid data, provided in a `mpc` format compatible with MATPOWER [4], with power flows pre-adjusted to match the generation profiles shown in Figures 2d–f.


Matpower 6.0 (https://matpower.org/download/) is required for the power flow calculation.

## Reference

 1. Y. Wang, A. N. Montanari, and A. E. Motter, “Nonnormal Frequency Dynamics under High-Renewable Penetration: A Case Study of the Iberian Blackout,” *submitted for publication*.
