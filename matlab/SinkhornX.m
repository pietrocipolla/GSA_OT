function [T,Tl,Ts]=SinkhornX(lambda,x,y,p,q)
% SINKHORNX Sinkhorn-Wasserstein with Richardson Extrapolation.
% [TX,T1,T2]=SINKHORNX(L,X,Y,P,Q)
% computes the squared Euclidean Wasserstein distance estimate
% given discrete probability masses P at X and Q at Y
% using Entropy parameter L (and sqrt(2)L).
% TX is the extrapolated Wsserstein distance from T1 / Sinkhorn with L
% and T2 / Sinkhorn with sqrt(2)L
% 
% from Chizat et al. 2020

 verbose = true;
 [m,k]=size(y);
 [n,l]=size(x);
 
 if(nargin<4)
  p=ones(n,1)/n; % uniform masses
  q=ones(m,1)/m;
 end
 
 c=-x*y'+(sum(x.^2,2)+sum(y.^2,2)');

 sigma=sqrt(2)*lambda;
 vl=zeros(n,1);
 vs=zeros(n,1);

 maxiter=200;
 converged=false;
 iter=0;
 while ~converged
  ul=-lambda*logsumexp( (vl'-c)/lambda,2,q); 
  vl=-lambda*logsumexp( (ul-c)/lambda,1,p)';
  us=-sigma*logsumexp( (vs'-c)/sigma,2,q);
  vs=-sigma*logsumexp( (us-c)/sigma,1,p)';
  iter=iter+1;
  
  if(verbose)), disp(sprintf('inter %d',iter)); end
  converged= ~all(isfinite(ul)) || iter==maxiter;
 end
% disp(['Error bounds ' num2str(max (abs(c(:))).^2./[lambda sigma]/iter)])
 Tl=(ul'*p+vl'*q);
 Ts=(us'*p+vs'*q);
 T=2*Tl-Ts;
end