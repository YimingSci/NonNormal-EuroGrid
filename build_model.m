function [A_ext, N_bus] = Build_model(pantagruel)
    % System base power
    Sb = pantagruel.baseMVA;
    N_bus  = length(pantagruel.bus);
    N_line = length(pantagruel.branch);

    %% === Load & Generator identification ===
    L = pantagruel.bus(:,3) / Sb;      
    G = zeros(N_bus,1);               

    is_producing = pantagruel.gen(:,2) > 0;
    id_gen  = pantagruel.gen(is_producing,1);
    id_load = setdiff(1:N_bus, id_gen)';

    G(id_gen) = pantagruel.gen(is_producing,2) / Sb;
    P = -L + G;
    P = P - mean(P);

    %% === Power flow solution ===
    result = pantagruel;
    A = sparse( ...
        [pantagruel.branch(:,1); pantagruel.branch(:,2)], ...
        [1:N_line, 1:N_line], ...
        [ones(N_line,1); -ones(N_line,1)] ...
    );

    V = ones(N_bus,1);
    B = 1 ./ result.branch(:,4);
    branch_weight = B .* (V(result.branch(:,1)) .* V(result.branch(:,2)));

    L_full = A * spdiags(branch_weight, 0, N_line, N_line) * A';
    L_full = full(L_full);

    %% Rearrange Laplacian matrix
    gen_idx  = id_gen(:,1);
    load_idx = id_load(:,1);
    n_gen  = length(gen_idx);
    n_load = length(load_idx);

    new_order = [gen_idx; load_idx];
    L_new = L_full(new_order, new_order);

    % Machine and load parameters
    M  = result.gen_inertia(is_producing);
    Dg = result.gen_prim_ctrl(is_producing) + result.load_freq_coef(id_gen);
    Dl = result.load_freq_coef(id_load);

    rr = sum(L_new, 2) - diag(L_new);

    %% Build system matrices
    n_dimension = n_gen*2 + n_load;

    % --- Matrix B1 ---
    B1 = zeros(n_dimension, n_dimension);
    B1(1:n_gen,1:n_gen) = diag(Dg);

    L13 = L_new(1:n_gen, 1:n_gen);
    B1(1:n_gen, n_gen+n_load+1:end) = L13 - diag(diag(L13)) + (-diag(rr(1:n_gen)));

    L21 = L_new(n_gen+1:n_gen+n_load, 1:n_gen);
    B1(n_gen+1:n_gen+n_load, 1:n_gen) = L21;

    L22 = L_new(n_gen+1:n_gen+n_load, n_gen+1:n_gen+n_load);
    B1(n_gen+1:n_gen+n_load, n_gen+1:n_gen+n_load) = ...
        L22 - diag(diag(L22)) + (-diag(rr(n_gen+1:end)));

    L31 = eye(n_gen);
    B1(n_gen+n_load+1:end, 1:n_gen) = L31;

    % --- Matrix B3 ---
    B31 = zeros(n_dimension, n_load);
    B31(1:n_gen,1:n_load) = L21';

    R   = -(L22 - diag(diag(L22)) + (-diag(rr(n_gen+1:end))));
    B32 = inv(R);

    B33 = zeros(n_load, n_dimension);
    B33(:, n_gen+1:n_gen+n_load) = -diag(Dl);
    B33(:, n_gen+n_load+1:end)  = L21;

    B3 = B31 * B32 * B33;

    %% Build B0
    B0 = zeros(n_dimension, n_dimension);
    B0(1:n_gen, 1:n_gen) = diag(M);
    B0(n_gen+1:n_gen+n_load, n_gen+1:n_gen+n_load) = diag(Dl);
    B0(n_gen+n_load+1:end, n_gen+n_load+1:end) = eye(n_gen);

    %% Final system matrix
    A1 = inv(B0) * B1;
    A2 = inv(B0) * B3;

    A_ext = A1 + A2;
end

