# EuropeanPG-V2



## Model Description

The European power-grid model incorporates **7,343** transmission lines, **3,809** buses, and **1,089** generators, including **470** renewable units representing available wind and solar installations. The grid topology is obtained from the published dataset in Ref. [2], while information on wind and solar installations is curated from Ref. [3]. The renewable generation data for Portugal, Spain, and France (collectively referred to as PSF region) reflect installations and capacities as of June 2025.

For power system analysis, we utilize widely adopted open-source tools, including the power flow solver MATPOWER [4] and the time-domain dynamic simulation framework described in Ref. [2]. Both tools are implemented in MATLAB, and we recommend using version R2024a for compatibility. See the reference below for more details.


## Usage

- `PSF_renewable.mat` : Power grid data, provided in a `mpc` format compatible with MATPOWER [4], with power flows pre-adjusted to match the generation profiles shown in Figures 2d–f.
- `Dynamic_analysis.m` : Time-domain simulation code. This script includes fault initialization, parameter setup, and dynamic simulation procedures.

Matpower 6.0 (https://matpower.org/download/) is required for the power flow calculation.

## Reference

 1. Y. Wang, A. N. Montanari, and A. E. Motter, “Nonnormal Frequency Dynamics under High-Renewable Penetration: A Case Study of the Iberian Blackout,” *submitted for publication*.
