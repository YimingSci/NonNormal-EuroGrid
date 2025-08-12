function y=Radau5(omega_gen, theta_gen, theta_load, M_gen, D_gen, D_load, P_gen, P_load, incidence_mat, line_susceptance, dt, maxiter, tol)
%=========================================================
%   For information on Radau methods, see, for instance,
%   E. Hairer and G. Wanner, Stiff differential equations solved by Radau
%   methods, J. Comput. Appl. Math. 11(1-2): 93-111 (1999)
%=========================================================

    % butcher tableau for radau 5 
    a = [(88-7*sqrt(6))/360, (296-169*sqrt(6))/1800, (-2+3*sqrt(6))/225
        (296+169*sqrt(6))/1800, (88+7*sqrt(6))/360, (-5-3*sqrt(6))/225
       (16-sqrt(6))/36, (16+sqrt(6))/36, 1/9];
    Ns = 3;
    
    % diagonalisation of A
    [Tr, lambda, Tl] = eig(inv(a));
    Tl = Tl';
    lambda = diag(lambda);
    d2 = diag(Tl*Tr);
    Tr(:,1) = Tr(:,1)/d2(1);
    Tr(:,2) = Tr(:,2)/d2(2);
    Tr(:,3) = Tr(:,3)/d2(3);
    
    y0 = [omega_gen; theta_gen; theta_load];
    
    N_line = length(line_susceptance);
    Ap = abs(incidence_mat);
    B = sparse(1:N_line, 1:N_line, line_susceptance);
    N_gen = length(omega_gen);
    N_bus = N_gen + length(theta_load);
    N_var = N_gen + N_bus;
    
    M1 = kron(Tr,sparse(1:N_var, 1:N_var,ones(N_var,1)));
    M2 = kron(diag(lambda)*Tl,sparse(1:N_var, 1:N_var,ones(N_var,1)));
    
    mat1 = incidence_mat(1:N_gen,:)*B;
    mat2 = incidence_mat(N_gen+1:end,:)*B;
    dtheta = incidence_mat'*y0(N_gen+1:end);
    
    % Jacobian matrix
    M = sparse(1:N_var, 1:N_var, [M_gen; ones(N_gen,1); D_load]);
    J0 = Ap * sparse(1:N_line, 1:N_line, line_susceptance .* cos(dtheta)) * Ap';
    J0 = J0 - sparse(1:N_bus, 1:N_bus,2*diag(J0));
    
    J = [sparse(1:N_gen, 1:N_gen, -D_gen) J0(1:N_gen,:);
        sparse(1:N_gen,1:N_gen,ones(N_gen,1)) sparse(N_gen,N_bus);
        sparse(N_bus-N_gen,N_gen) J0(N_gen+1:end,:)];

    P = [P_gen; zeros(N_gen,1); P_load];
   
    Y = repmat(y0,Ns,1);
    F = zeros(length(y0)*Ns, 1);
    W = zeros(length(y0)*Ns, 1);
    dW = zeros(length(y0)*Ns, 1);
    iter = 0;
    error = 2*tol;
    % Newton method
	while(error > tol && iter < maxiter)
        for s = 1:Ns
            F(N_var*(s-1)+1:N_var*s) = M*(Y(N_var*(s-1)+1:N_var*s)-y0);
            for o = 1:Ns
                dtheta = incidence_mat'*Y(N_var*(o-1)+N_gen+1:N_var*o);
                F(N_var*(s-1)+1:N_var*s) = F(N_var*(s-1)+1:N_var*s) - dt*a(s,o)*(P + [-D_gen.*Y(N_var*(o-1)+1:N_var*(o-1)+N_gen) - mat1*sin(dtheta); Y(N_var*(o-1)+1:N_var*(o-1)+N_gen); -mat2*sin(dtheta)]);
            end
        end
        F2 = -M2*F;
        for s = 1:Ns
            dW(N_var*(s-1)+1:N_var*s) = (lambda(s)*M - dt*J) \ F2(N_var*(s-1)+1:N_var*s);
        end
        W = W + dW;
        error = max(abs(M1*dW));
        Y = M1*W + repmat(y0,Ns,1);
        iter = iter + 1;
	end
    if(iter == maxiter)
        disp(['Max iteration reached, error ' num2str(error)])
    end
    y = Y((Ns-1)*N_var+1:end); 
end
