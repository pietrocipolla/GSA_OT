function [X,etas]=transportLDN(C,w,v,lam)
% TRANSPORTLDN approximately solve a transport problem
%   using the gradient descend algorithm for a linear diagonal network
% 
% Haoyue Wang, Promit Ghosal, Rahul Mazumder:
% Linear programming using diagonal linear networks, arXiv 2023

    scale = 2;
    max_iter=8000; 
    bnd = 1e-2; % 16*eps;
    % reparametrized GD
    L = 2;
    U = exp(-C/(2*lam)); X = U.^2;
    %etas=nan(max_iter,1);
    for t=1:max_iter
        g = sum(X, 2) - w;
        h = sum(X, 1) - v';
		err = max( max(abs(g)), max(abs(h)) );
		%if(t<=32 || bitand(t,3)==0)
         nrm = max(X,[],'all'); % max(X(:))
         if err<bnd*nrm || ~isfinite(nrm), break; end
	     eta = scale*min(0.25/err, 0.2/(L*nrm));
        % etas(t)=eta;
        %end
        %U = U - eta * (g.*U + U.*h) + eta^2 * (g*h) .* U;
        U = (1 - eta*((g+h))).*U; 
		%U = (1 - eta*((g+h) - eta*(g*h))).* U;
        X = U.^2;
    end
% complete last step
    eta = scale*min(0.25/err, 0.2/(L*nrm));
 	U = (1 - eta*((g+h))).*U;
    X = U.^2;	
end