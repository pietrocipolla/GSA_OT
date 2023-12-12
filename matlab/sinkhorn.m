function dM=sinkhorn(M,lambda,r,c)
% Algorithm 1 Computation of d_M(r; c) using 
% Sinkhorn-Knopp's fixed point iteration
%
% Inputs: M       dxd cost matrix
%         LAMBDA  Lagrange multiplyer (>>0)
%         R       d vector containing histogram
%         C       dxm matrix containing multiple histograms
	
% Cuturi 2013

if(isempty(lambda)), lambda=1/median(M(:)); end
I=(r>0); r=r(I); M=M(I,:); K=exp(-lambda*M);
x=ones(size(r,1),size(c,2))/length(r);
converged=false;
iter=0;
while ~converged 
 xold=x;
 x=diag(1./r)*K*(c.*(1./(K'*(1./x))));
 iter=iter+1;
 converged=all(abs(x-xold)<max(max(x))*50*eps) || iter==500;
end
u=1./x; v=c.*(1./(K'*u));
dM=sum(u.*((K.*M)*v));
end

function testsink
%%
n=3000;
a=randn(n,1);
b=randn(n,1)*2+1;
W22=mean((sort(a)-sort(b)).^2)

M=(a-b').^2;
W22s=sinkhorn(M,[],ones(n,1)/n,ones(n,1)/n)

%%
end

function testsink2
%%
ishigami
n=1024;
x=trafo(sobolpoints(n,k));
y=model(x);
D=(y-y').^2;
p=ones(n,1)/n;
xt=trafo(zeros(1,k));
M=12;
q=[];
for m=1:M
  m
  xu=trafo(ones(1,k)*m/M);
  % select between xt and xu
  ii=x<=xu & x>xt;
  q=[q,ii./sum(ii)];
  xt=xu;
end
%%
tic
d=1;
fail=false
W0=zeros(1,k*M); % nans are replaced by 0s, might not be wanted
while ~fail
 W1=sinkhorn(D,d,p,q);
 ii=isnan(W1); 
 disp(sprintf('Used %g as multiplyer, failure %d out of %d',d,sum(ii),k*M));
 W0(~ii)=W1(~ii);
d=d*sqrt(2);
fail=all(ii);
end
toc
W=reshape(W0,k,M);  
plot(1:M,W)
mean(sqrt(W),2)
%%
end