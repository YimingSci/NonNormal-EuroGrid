# Non-normality analysis of the European power grid



This repository contains code for performing non-normality analysis of the European power grid, along with the associated dataset. The data file `EUR_2025.mat` provides a  2025 model of the European power grid, including data on the grid transmission network, generator and load geographical locations (with a high share of renewable sources), and generation mix at the time of the blackout. For details on the methods and theoretical background, see Ref. [1]. The workflow is built on the open-source MATLAB toolbox MATPOWER [2] for steady-state power flow calculations. All codes were tested on MATLAB R2024a. 


## Iberian Blackout on April 28, 2025
Using data from the respective transmission system operators (TSOs), Fig. S1(a–c) illustrates the electricity generation patterns for Portugal [3], Spain [3], and France [4] on April 28, 2025. A blackout occurred at approximately 12:32 p.m. (UTC+01:00). To replicate these conditions in our simulations, the file `EUR_2025.mat` was prepared with generator outputs and power flows adjusted so that the modeled generation mix in the three countries matches the actual system state on that day.

Fig. S1(d) summarizes the initial disturbance setup, where more than ten individual power-loss events are grouped into six representative buses, each named after a nearby major city. The Bus ID corresponds to the simulated disturbance node and can be matched to the bus numbering in `EUR_2025.mat`.

<img width="1000" alt="图片2" src="https://github.com/user-attachments/assets/150a4c68-7639-4d7e-b3e4-4f67d656fe49" />


**Fig. S1:**
(**a–c**) TSO-reported generation mixes for Portugal, Spain, and France on April 28, 2025. Portugal and Spain share a similar profile, with solar and wind power forming the bulk of their generation. In contrast, France relies primarily on nuclear power. (**d**) Locations and associated power losses of the initial disturbances.



## Usage

- `EUR_2025.mat` : **European power grid data.** The data file is provided in `mpc` format, which is compatible with MATPOWER [2].
- `Cal_nonnormal.m` : **Non-normality analysis of the European power grid.** This code first calls `Build_model.m` to construct the linearized system representation, then calculates the non-normality index by comparing the spectral properties of the system matrix and its symmetric part. The script proceeds to extract the first K dominant modes, compute the modal non-normality for each bus, and finally generates a plot showing the cumulative non-normality distribution.
- `Build_model.m` : **Linearized state-space representation of the European power grid.** This code identifies generator and load buses, computes the network Laplacian from branch parameters, and rearranges it into separate generator and load dynamics. Using the generator inertias, primary control gains, and load frequency coefficients, the function assembles the system matrices and returns the extended system matrix `A_ext` along with the total number of buses `N_bus`. This model serves as the foundation for the non-normality analysis.
- `Simulation.m ` : **Simulates the dynamic frequency response of the European power grid under multiple generator power disturbances.** The script loads the network model, identifies generator and load buses, and applies specified active power deviations to selected generators to emulate fault events. Using generator inertias, damping coefficients, and load frequency sensitivities, it integrates the system’s swing equations over time (via the radau5 solver) to track generator frequencies and rotor angles. This enables analysis of system stability under transient conditions.


Matpower 6.0 (https://matpower.org/download/) is required for the power flow calculation.

## Dependency

The following codes were used for simulations presented in Ref. [5]. If you use them
in any future work, please provide proper credit by citing Ref. [5].

- `NRsolver.m` : Solves power flow equations using the Newton-Raphson method.
- `Radau5.m`   : Integrates stiff differential equations using the Radau IIA method.

These scripts are included here with attribution for reproducibility and completeness.

## License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".


## Reference

 1. Y. Wang, A. N. Montanari, and A. E. Motter, “Nonnormal Frequency Dynamics under High-Renewable Penetration: A Case Study of the Iberian Blackout,” *submitted for publication*.
 2. R. D. Zimmerman, C. E. Murillo-Sánchez, and R. J. Thomas, “MATPOWER: Steady-State Operations, Planning, and Analysis Tools for Power Systems Research and Education,” *IEEE Trans. Power Syst.*, vol. 26, no. 1, pp. 12–19, 2011.   [https://doi.org/10.1109/TPWRS.2010.2051168](https://doi.org/10.1109/TPWRS.2010.2051168)
 3. ENTSO-E, "Iberian blackout on 28 April 2025", June 2025.   [https://www.entsoe.eu/publications/blackout/28-april-2025-iberian-blackout/](https://www.entsoe.eu/publications/blackout/28-april-2025-iberian-blackout/)
 4. RTE France, "eco2mix – Power generation by energy source".   [https://www.rte-france.com/en/eco2mix/power-generation-energy-source](https://www.rte-france.com/en/eco2mix/power-generation-energy-source) (accessed Jul. 3, 2025)
 5. M. Tyloo, L. Pagnier, and P. Jacquod, “The key player problem in complex oscillator networks and electric power grids: Resistance centralities identify local vulnerabilities,” *Science Advances*, vol. 5, no. 11, eaaw8359, 2019. [https://doi.org/10.1126/sciadv.aaw8359](https://doi.org/10.1126/sciadv.aaw8359)
