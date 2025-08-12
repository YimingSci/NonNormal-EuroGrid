# NonNormal EuroGrid



## Repository Description

This repository contains code for performing non-normality analysis of the European power grid, along with the associated dataset. For a deeper explanation of the methods and theoretical background, see Ref. [1]. The workflow is built on open-source MATLAB tools: MATPOWER [2] for steady-state power flow calculations. We recommend MATLAB R2024a to ensure smooth execution. 


## Energy Mix
Using data from the respective transmission system operators (TSOs), Fig. S1a–c illustrates the electricity generation patterns for Portugal [3], Spain [3], and France [4] on April 28, 2025. A blackout occurred at approximately 12:32 p.m. (UTC+01:00). To replicate these conditions in our simulations, the file EUR_2025.mat was prepared with generator outputs and power flows adjusted so that the modeled generation mix in the three countries matches the actual system state on that day.

<img width="2685" height="2090" alt="Energy_mix" src="https://github.com/user-attachments/assets/8b0e53cf-fe2f-4099-b6e8-032ab95a96f3" />

**Fig. S2:**
*(a,b,c) Generation mix of Portugal, Spain, and France on April 28, 2025 as reported by the TSO. Portugal and Spain exhibit similar structures, with solar and wind power dominating their generation portfolios. In contrast, nuclear power constitutes the dominant source in France.
(d,e,f) Adjusted generation mix used in the grid model for the corresponding countries.*



## Usage

- `Cal_nonnormal.m` : **This script performs the non-normality analysis of the European power grid.** It first calls `Build_model.m` to construct the linearized system representation, then calculates the non-normality index by comparing the spectral properties of the system matrix and its symmetric part. The script proceeds to extract the first k dominant modes, compute and sum the modal reactivity for each bus, and finally generates a plot showing the cumulative reactivity distribution in ascending bus order.
- `Build_model.m` : **Constructs the linearized state-space representation of the European power grid from the provided dataset.** It identifies generator and load buses, computes the network Laplacian from branch parameters, and rearranges it to separate generator and load dynamics. Using machine inertias, primary control gains, and load frequency coefficients, the function assembles the system matrices and returns the extended system matrix A_ext along with the total number of buses N_bus. This model serves as the foundation for subsequent non-normality and modal analysis.
- `EUR_2025.mat` : Power grid data, provided in a `mpc` format compatible with MATPOWER [2].


Matpower 6.0 (https://matpower.org/download/) is required for the power flow calculation.


## License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".


## Reference

 1. Y. Wang, A. N. Montanari, and A. E. Motter, “Nonnormal Frequency Dynamics under High-Renewable Penetration: A Case Study of the Iberian Blackout,” *submitted for publication*.
 2. R. D. Zimmerman, C. E. Murillo-Sánchez, and R. J. Thomas, “MATPOWER: Steady-State Operations, Planning, and Analysis Tools for Power Systems Research and Education,” *IEEE Trans. Power Syst.*, vol. 26, no. 1, pp. 12–19, 2011.   [https://doi.org/10.1109/TPWRS.2010.2051168](https://doi.org/10.1109/TPWRS.2010.2051168)
 3. ENTSO-E, "Iberian blackout on 28 April 2025", June 2025.   [https://www.entsoe.eu/publications/blackout/28-april-2025-iberian-blackout/](https://www.entsoe.eu/publications/blackout/28-april-2025-iberian-blackout/)
 4. RTE France, "eco2mix – Power generation by energy source".   [https://www.rte-france.com/en/eco2mix/power-generation-energy-source](https://www.rte-france.com/en/eco2mix/power-generation-energy-source) (accessed Jul. 3, 2025)
