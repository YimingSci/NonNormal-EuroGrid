%% European Power Grid Dynamic Fault Response Simulation
%
% This script simulates the transient frequency and rotor angle dynamics of
% the European transmission grid in response to multiple generator power
% disturbances. The network model is loaded from `EUR_2025.mat` and contains
% bus, branch, and generator parameters for the high-voltage system.
%
% The workflow includes:
%   - Loading the European grid model and identifying generator/load buses
%   - Performing a power flow calculation to determine initial operating points
%   - Applying specified active power disturbances to selected generators
%   - Integrating the swing equations using a Radau5 implicit solver
%   - Storing generator frequency and rotor angle trajectories for analysis
%
% This framework supports:
%   - Studying system stability under multiple simultaneous faults
%   - Assessing damping and inertia effects on transient response
%   - Guiding control design for frequency stability enhancement
%
% For detailed application, methodology, and case studies, please refer to:
%   * "Nonnormal Frequency Dynamics under High-Renewable Penetration:
%      A Case Study of the Iberian Blackout"
%     Y. Wang, A. N. Montanari, and A. E. Motter
%     (submitted for publication)
%
% Please cite the above article if you use this code in your research.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025  Y. Wang, A. N. Montanari & A. E. Motter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
%
% Last modified by Y. Wang on 2025-08-12



clear;
load('EUR_2025.mat');

%% Basic settings
Sb = 100;                 % System base power (MVA)
N_bus  = length(pant.bus);
N_line = length(pant.branch);

dt  = 1E-3;               % Time step
Ndt = 5000;               % Number of simulation steps

%% Fault event settings (generator IDs and power changes in MW)
id_fault1 = 3335; dP_real1 = 355;
id_fault2 = 3159; dP_real2 = 917;
id_fault3 = 3324; dP_real3 = 550;
id_fault4 = 3036; dP_real4 =  23;
id_fault5 = 3704; dP_real5 =  34;
id_fault6 = 2537; dP_real6 = 37.5;

%% Load and generation setup
L = pant.bus(:,3) / Sb;          % Per-unit loads
G = zeros(N_bus,1);              % Per-unit generation

is_producing = pant.gen(:,2) > 0;
id_gen  = pant.gen(is_producing,1);
id_load = setdiff(1:length(L), id_gen)';

N_gen  = length(id_gen);
N_load = length(id_load);

G(id_gen) = pant.gen(is_producing,2) / Sb;
P = -L + G;                      % Net injections
P = P - mean(P);                 % Zero-mean adjustment

%% Power flow calculation (using Newton-Raphson)
A = sparse([pant.branch(:,1); pant.branch(:,2)], ...
           [1:N_line, 1:N_line], ...
           [ones(N_line,1); -ones(N_line,1)]);
b = -1i ./ pant.branch(:,4);     % -j/X
Ybus = conj(A * sparse(1:N_line,1:N_line,b) * A');
Q = zeros(N_bus,1);
V = ones(N_bus,1);
theta = zeros(N_bus,1);
[~, theta, ~, ~] = NRsolver(Ybus, V, theta, -P, Q, [], 9, 1E-11, 1000);

P_gen  = P(id_gen);
P_load = P(id_load);

%% Build edge matrix for network topology
edges = zeros(N_line,2);
line_start = pant.branch(:,1);
line_end   = pant.branch(:,2);

for i = 1:N_gen
    edges(line_start == id_gen(i),1) = i;
    edges(line_end   == id_gen(i),2) = i;
end
for i = 1:N_load
    edges(line_start == id_load(i),1) = i + N_gen;
    edges(line_end   == id_load(i),2) = i + N_gen;
end

line_susceptance = 1 ./ pant.branch(:,4);

%% Apply fault power disturbances (in per-unit)
faults = [id_fault1, dP_real1;
          id_fault2, dP_real2;
          id_fault3, dP_real3;
          id_fault4, dP_real4;
          id_fault5, dP_real5;
          id_fault6, dP_real6];

for f = 1:size(faults,1)
    dP = -faults(f,2) / Sb;
    row_idx = find(id_gen == faults(f,1));
    P_gen(row_idx) = P_gen(row_idx) + dP;
end

%% Original dynamic parameters
M_gen_orig = pant.gen_inertia(is_producing);
D_gen_orig = pant.gen_prim_ctrl(is_producing) + pant.load_freq_coef(id_gen);
D_load     = pant.load_freq_coef(id_load);

%% Initial states
omega_gen_init = zeros(N_gen,1);
theta_gen_init = theta(id_gen);
theta_load_init = theta(id_load);

%% Incidence matrix and graph structure
incidence_mat = sparse([edges(:,1); edges(:,2)], ...
                       [1:N_line, 1:N_line], ...
                       [ones(N_line,1); -ones(N_line,1)]);

G_graph = graph(pant.branch(:,1), pant.branch(:,2));
A_adj = adjacency(G_graph);
deg = sum(A_adj, 2);
D_inv = spdiags(1./deg, 0, N_bus, N_bus);
W = D_inv * A_adj;

%% Node sizes for visualization
sizes = ones(N_bus,1) * 15;
sizes(id_gen) = 30;

%% Adjusted inertia and damping
M_gen = M_gen_orig;
D_gen = D_gen_orig;

%% Simulation initialization
omega_gen = omega_gen_init;
theta_gen = theta_gen_init;
theta_load = theta_load_init;

m = 10;  % Downsampling factor for storing results
omega_t = zeros(N_gen, floor(Ndt/m));
delta_t = zeros(N_gen, floor(Ndt/m));
k = 1;

%% Time integration using Radau5 solver
for i = 1:Ndt
    y = radau5(omega_gen, theta_gen, theta_load, ...
               M_gen, D_gen, D_load, P_gen, P_load, ...
               incidence_mat, line_susceptance, dt, 14, 1E-6);

    if mod(i, 200) == 0
        fprintf('Completed %d steps...\n', i);
    end

    if mod(i, m) == 0
        omega_t(:,k) = y(1:N_gen);
        delta_t(:,k) = y(N_gen+1:2*N_gen);
        k = k + 1;
    end

    omega_gen  = y(1:N_gen);
    theta_gen  = y(N_gen+1:2*N_gen);
    theta_load = y(2*N_gen+1:end);
end
