function B=mysqrtm(A)
% MYSQRTM robust symmetric matrix square root
 [U,S,V]=svd(A,'econ');
 B=U*diag(sqrt(diag(S)))*V';
end