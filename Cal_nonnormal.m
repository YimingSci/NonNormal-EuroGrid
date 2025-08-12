clear; clc;
load('EUR_V2.mat');

%% Build linearized model
pantagruel = pant;
[A_ext, N_bus] = build_model(pantagruel);
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
