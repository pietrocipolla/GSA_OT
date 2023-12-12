%% 
% Implementation of our Greenkhorn algorithm for matrix scaling.
% Specifically, each iteration, choose the row or column that gives
% maximum "gain" (w.r.t. our gain function), and then normalize P to have
% total mass 1. See paper for further details.
%
% Input parameters:
%  -- A:  positive square matrix to project onto U_{r,c} transport polytope (dims: nxn)
%  -- r:  desired row sums (marginals)         (dims: nx1)
%  -- c:  desired column sums (marginals)      (dims: 1xn)
%  -- updates: number of row/col updates
%  -- compute_otvals: flag whether to compute otvals (slow but used in some plots)
%  -- C:  cost matrix for OT
%  -- ot_downsampling: how often to compute current iterate's OT value (this step is slow)
% 
% Output:
%  -- P:   final scaled matrix
%  -- err: sum of row and column violations at each iteration
%  -- ot:  values of optimal transport of matrix iterates

function [P, err, otvals] = greenkhorn(A,r,c,updates,compute_otvals,C,ot_downsampling)
%GREENKHORN  Matrix scaling via greedy Sinkhorn algorithm

% addpath(genpath('algorithms/'));
P = A;

% compute full row and column marginals once
r_P = sum(P,2);
c_P = sum(P,1);

% compute gains for each row and column
r_gain = r_P - r + r.*log(r./r_P);
c_gain = c_P - c + c.*log(c./c_P);

err = zeros(updates+1,1);
err(1) = norm(r_P-r,1)+norm(c_P-c,1);

if compute_otvals
    % initialize OT
    otvals = zeros(updates+1,1);
    otvals(1) = frobinnerproduct(round_transpoly(P,r,c),C);
end

for t=1:updates
    % find row or column with maximum gain
    [r_gain_max, i] = max(r_gain);
    [c_gain_max, j] = max(c_gain);
    
    if r_gain_max > c_gain_max        
        % update row i
        scaling = r(i)/r_P(i);
        old_row = P(i,:);
        new_row = old_row*scaling;
        P(i,:)  = new_row;
        
        % renormalize (can also be done implicitly if one wants to optimize)
        P = P/sum(sum(P));
        
        % compute full row and column marginals
        r_P = sum(P,2);
        c_P = sum(P,1);
        
        % compute gains for each row and column
        r_gain = r_P - r + r.*log(r./r_P);
        c_gain = c_P - c + c.*log(c./c_P);

        % % tricks to speed up computation if we are not renormalizing
        % % matrix each time
        %         % update row and column sums in O(n) time
        %         r_P(i)  = r(i);
        %         c_P     = c_P - old_row + new_row;
        %
        %         % update row and column gains in O(n) time
        %         r_gain(i) = 0;
        %         c_gain    = c_P - c + c.*log(c./c_P);
        
        err(t+1) = norm(r_P-r,1)+norm(c_P-c,1);        
    else
        % update column j
        scaling = c(j)/c_P(j);
        old_col = P(:,j);
        new_col = old_col*scaling;
        P(:,j)  = new_col;
        
        % renormalize (can also be done implicitly if one wants to optimize)
        P = P/sum(sum(P));
        
        % compute full row and column marginals
        r_P = sum(P,2);
        c_P = sum(P,1);
        
        % compute gains for each row and column
        r_gain = r_P - r + r.*log(r./r_P);
        c_gain = c_P - c + c.*log(c./c_P);
        
        % % tricks to speed up computation if we are not renormalizing
        % % matrix each time
        %       % update row and column sums in O(n) time
        %       c_P(j)  = c(j);
        %       r_P     = r_P - old_col + new_col;
        %
        %       % update row and column gains in O(n) time
        %       c_gain(j) = 0;
        %       r_gain = r_P - r + r.*log(r./r_P);
        
        err(t+1) = norm(r_P-r,1)+norm(c_P-c,1);
    end
    
    if compute_otvals && (mod(t,ot_downsampling) == 0)
        otvals(t+1) = frobinnerproduct(round_transpoly(P,r,c),C);
    end
end
P=round_transpoly(P,r,c); % project
end
% make file self-contained
function d=frobinnerproduct(A,B)
% FROBINNERPRODUCT Frobenius product of matrices
 d=A(:)'*B(:); % originally sum(sum(A.*B))
end

function A = round_transpoly( X,r,c )
% Implementation of our algorithm for rounding a matrix onto the U_{r,c}
% transport polytope. See Section 2 of the paper for details.
A=X;
n=size(A,1);
r_A = sum(A,2);
for i=1:n
    scaling = min(1,r(i)/r_A(i));
    A(i,:)=scaling*A(i,:);
end

c_A = sum(A,1);
for j=1:n
    scaling = min(1,c(j)/c_A(j));
    A(:,j)=scaling*A(:,j);
end

r_A = sum(A,2);
c_A = sum(A,1);
err_r = r_A - r;
err_c = c_A - c;

A = A + err_r*err_c/sum(abs(err_r));

end