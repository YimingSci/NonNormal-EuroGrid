clear; clc;
%% European Power Grid Non-Normality Analysis
%
% This code implements an analysis of the non-normality of the European
% transmission grid, based on a linearized model derived from the dataset
% `EUR_2025.mat`. The study focuses on high-renewable penetration scenarios,
% with emphasis on Iberian system dynamics.
%
% The workflow includes:
%   - Building the extended system matrix from the aggregated "pantagruel" model
%   - Constructing the symmetric (Hermitian) component of the system matrix
%   - Quantifying non-normality via the difference between Hermitian and
%     non-Hermitian dominant eigenvalues
%   - Extracting the first k dominant modes and computing their cumulative reactivity
%   - Sorting and visualizing reactivity across buses to identify the most
%     sensitive nodes under varying inertia conditions
%
% This framework supports:
%   - Stability assessment in high-renewable scenarios
%   - Evaluation of inertia-related vulnerability in interconnected grids
%   - Informed design of damping and control strategies
%
% For detailed application, analysis, and conclusions, please refer to:
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

load('EUR_2025.mat');

%% Build linearized model
pantagruel = pant;
[A_ext, N_bus] = Build_model(pantagruel);
H = (A_ext + A_ext') / 2;

%% Non-normality difference
lambda_max_H = eigs(H, 1, 'la');
lambda_max_A = eigs(A_ext, 1, 'largestreal');
lambda_delta = lambda_max_H - lambda_max_A;

%% First k modes and cumulative reactivity
kk_list = 1:20;
[V_all, D_all] = eigs(H, max(kk_list), 'la');
lambda_all = diag(D_all);
SumReactivity = zeros(N_bus,1);
for kk = kk_list
    v_tmp = V_all(:, kk) / norm(V_all(:, kk));
    v_re_total = abs(lambda_all(kk)) * lambda_delta * v_tmp;
    SumReactivity = SumReactivity + log10((v_re_total(1:N_bus)).^2);
end

%% Sort and plot
[vals_sorted, ~] = sort(SumReactivity, 'ascend');
figure(10); clf; hold on;
plot(1:N_bus, vals_sorted, '-', 'LineWidth', 2);
xlabel('Sorted Bus Index');
ylabel('Cumulative Reactivity');
title('Total Modal Reactivity under Varying Iberian Inertia');
grid on; hold off;
