function [V,theta,iter,J] = NRsolver(Ybus, V, theta, P0, Q0, id_PQ, id_slack, epsilon, maxiter)
%=========================================================
%   For information on solving the power flow equations with 
%	Newton-Raphson, see, for instance,
%	V. Vittal and A. Bergen, Power systems analysis,
%   Prentice Hall, 1999.
%=========================================================
    if(nargin < 9)
        maxiter = 50;
    end
    n = length(Ybus);
    error  = 2*epsilon;
    iter = 0;  
    id = [1:id_slack-1 id_slack+1:n];
    while(error > epsilon && iter < maxiter)
        Vc = V.*exp(1i*theta);
        S = Vc.*(conj(Ybus)*conj(Vc));
        dPQ = [real(S(id))-P0(id); imag(S(id_PQ))-Q0(id_PQ)];
        temp1 = -1i*sparse(1:n,1:n,Vc)*conj(Ybus) * sparse(1:n,1:n,conj(Vc)) + 1i*sparse(1:n,1:n,Vc.*conj(Ybus*Vc));
        temp2 = sparse(1:n,1:n,Vc)*conj(Ybus) * sparse(1:n,1:n,exp(-1i*theta)) + sparse(1:n,1:n,exp(1i*theta).*conj(Ybus*Vc));
        J = [real(temp1(id,id)) real(temp2(id,id_PQ)); imag(temp1(id_PQ,id)) imag(temp2(id_PQ,id_PQ))];
        x = J\dPQ;
        theta(id) = theta(id) - x(1:n-1);
        if(~isempty(id_PQ))
            V(id_PQ) = V(id_PQ) - x(n:end);
        end
        error = max(abs(dPQ));
        iter = iter + 1;
    end
	if(iter == maxiter)
        disp(['Max iteration reached, error ' num2str(error)])
    end
end
