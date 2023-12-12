function [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options)

% sinkhorn_log - stabilized sinkhorn over log domain with acceleration
%
%   [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);
%
%   mu and nu are marginals.
%   c is cost
%   epsilon is regularization
%   coupling is 
%       gamma = exp( (-c+u*ones(1,N(2))+ones(N(1),1)*v')/epsilon );
%
%   options.niter is the number of iterations.
%   options.tau is an avering step. 
%       - tau=0 is usual sinkhorn
%       - tau<0 produces extrapolation and can usually accelerate.
%
%   options.rho controls the amount of mass variation. Large value of rho
%   impose strong constraint on mass conservation. rho=Inf (default)
%   corresponds to the usual OT (balanced setting). 
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;
niter = getoptions(options, 'niter', 1000);
tau  = getoptions(options, 'tau', -.5);
verb = getoptions(options, 'verb', 1);
rho = getoptions(options, 'rho', Inf);
nskip = getoptions(options, 'nskip', 1);

lambda = rho/(rho+epsilon);
if rho==Inf
    lambda=1;
end

N = [size(mu,1) size(nu,1)];
H1 = ones(N(1),1);
H2 = ones(N(2),1);

ave = @(tau, u,u1)tau*u+(1-tau)*u1;


%lse = @(A)log(sum(exp(A),2));
lse = @(A)logsumexp(A,2); % stable version (EP)
%M = @(u,v)(-c+u*H2'+H1*v')/epsilon;
M = @(u,v)(-c+u+v')/epsilon; % implicit binary singleton expansion
% kullback divergence
H = @(p)-sum( p(:).*(log(p(:)+1e-20)-1) );
KL  = @(h,p)sum( h(:).*log( h(:)./p(:) ) - h(:)+p(:) );
KLd = @(u,p)sum( p(:).*(exp(-u(:))-1) );
dotp = @(x,y)sum(x(:).*y(:));

err = [];
u = zeros(N(1),1); 
v = zeros(N(2),1);

Muv = M(u,v); 
for i=1:niter
    if verb==1
        progressbar(i,niter);
    end
    u1 = u;
    u = ave(tau, u, ...
        lambda*epsilon*log(mu) - lambda*epsilon*lse( Muv ) + lambda*u );
    v = ave(tau, v, ...
        lambda*epsilon*log(nu) - lambda*epsilon*lse( M(u,v)' ) + lambda*v );
    
    Muv = M(u,v);
    if i==niter ||  mod(i-1,nskip)==0 
        % coupling
    gamma = exp(Muv);
    if rho==inf % marginal violation
        Wprimal(i) = dotp(c,gamma) - epsilon*H(gamma);
        Wdual(i) = dotp(u,mu) + dotp(v,nu) ...
            - epsilon*sum( gamma(:) );
        err(i,1) = norm( sum(gamma,2)-mu );
    else % difference with previous iterate
        Wprimal(i) = dotp(c,gamma) - epsilon*H(gamma) ...
            + rho*KL(sum(gamma,2),mu) ...
            + rho*KL(sum(gamma,1),nu);
        Wdual(i) = -rho*KLd(u/rho,mu) - rho*KLd(v/rho,nu) ...
            - epsilon*sum( gamma(:) );
        err(i,1) = norm(u(:)-u1(:), 1);
    end
    end
end

end

%%
function v = getoptions(options, name, v, mendatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mendatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mendatory
    error(['You have to provide options.' name '.']);
end 
end

%%
function progressbar(n,N,w)

% progressbar - display a progress bar
%
%    progressbar(n,N,w);
%
% displays the progress of n out of N.
% n should start at 1.
% w is the width of the bar (default w=20).
%
%   Copyright (c) Gabriel PeyrŽ 2006

if nargin<3
    w = 20;
end

% progress char
cprog = '.';
cprog1 = '*';
% begining char
cbeg = '[';
% ending char
cend = ']';

p = min( floor(n/N*(w+1)), w);

global pprev;
if isempty(pprev)
    pprev = -1;
end

if not(p==pprev)
    ps = repmat(cprog, [1 w]);
    ps(1:p) = cprog1;
    ps = [cbeg ps cend];
    if n>1
        % clear previous string
        fprintf( repmat('\b', [1 length(ps)]) );
    end
    fprintf(ps);
end
pprev = p;
if n==N
    fprintf('\n');
end
end