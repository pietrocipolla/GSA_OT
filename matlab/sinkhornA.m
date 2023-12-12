function [OT,retargs]=sinkhornA(C,lambda,x,y,alfa,beta)
% SINKHORNA Optimal Transport via Altschuler Sinkhorn.
%
% [W,R]=SINKHORNA(C,L,X,Y,P,Q) returns the Euclidean Wasserstein 
%   distance (squared) in W and further results in R.
%   Use fieldnames(R) for a list of available data
%
%   Inputs: C   cost function c(x,y) (default squared Euclidean), N x M matrix
%	          L   regularization
%           X,Y     sample locations (sizes N and M)
%           P,Q     weight vectors (sizes N and M) (uniform defaults)

% Near-linear time approximation algorithms for optimal transport
% via Sinkhorn iteration
% Jason Altschuler, Jonathan Weed, and Philippe Rigollety
 maxiter=10;
 miniter=0;
 maxerr=1e-5;
 verbose=true; %false; % 
 
% assume: 
% sum(p)==sum(q)==1
% size(c)==[n,n], all(all(c>=0)), size(p)==[n,1],size(q)==[n,1]
 n=size(x,1);
 m=size(y,1);
 if(nargin<5)
  alfa=ones(n,1)/n; % uniform masses
  beta=ones(m,1)/m;
 end
 % Euclidean norm squared
 if isempty(C)
   C=@(x,y)sum(x.^2,2)+sum(y.^2,2)'-2*x*y';
 end
 cxy=C(x,y);nc=max(abs(cxy(:)));
 
 Delta=max(sqrt(cxy(:))); % diameter
 
 if(isempty(lambda))
  eta   = 2*log(n*m)/maxerr;lambda=1/eta;
  disp(['SINKHORN: lambda set to ' num2str(lambda) ]);
 else
  eta=1/lambda;
 end

 errrel = maxerr/(8*nc);
 
 isodd=@(k)bitand(k,1);
 rowsum=@(a)sum(a,2);
 colsum=@(a)sum(a,1)';
 distA=@(a,r,c)sqrt(mean((rowsum(a)-r).^2))+sqrt(mean((colsum(a)-c).^2));
 
 converged=false;
 
 A=exp(-eta*cxy);K=A; % save for return 
 A=A/sum(abs(A(:)));
 
 x=zeros(n,1);
 y=zeros(m,1);
 iter=0;
 while ~converged
   iter = iter+1;
   if(isodd(iter))
    x=x+log(alfa./rowsum(A));
    x(isnan(x))=1;
  else
    y=y+log(beta./colsum(A));
    y(isnan(y))=1;
  end
  A=exp(y)'.*A.*exp(x);
  err=distA(A,alfa,beta);
  
  if verbose   
  disp(sprintf('iter %d, error %g',iter,err)); 
  end
  converged=(iter>miniter && err<errrel) || iter==maxiter;
 end
 A=A.*min(1,alfa./rowsum(A));
 A=min(1,beta./colsum(A))'.*A;
 erow=alfa-sum(A,2);ecol=beta-sum(A,1)';
 P=A-erow*ecol'/sum(abs(erow));
 
 OT=P(:)'*cxy(:)
 T2=eta*(x'*alfa+y'*beta)
 
 if(nargout>1)
  Q=beta'.*P.*alfa;
  W2prime=(Q(:)'*cxy(:));
  SinkhornPrimal=W2prime+lambda*alfa'*(P.*(log(P+1e-20)-1))*beta; 
  SinkhornDual=T2-lambda*alfa'*P*beta; % sum(P(:)); 
   retargs = struct('SinkhornDivergence',OT,'W2Dual',T2,...
      'W2Primal',W2prime,'Kernel',K,'Lambda',lambda,...
      'NumIterations',iter','LastError',err,...
      'Cost',cxy,'Coupling',Q,...
      'SinkhornPrimal',SinkhornPrimal,'SinkhornDual',SinkhornDual);
 end
end
function y=expz(x)
% EXPZ exp with expz(-infty)=0
% B. Schmitzer
 ii=isinf(x) & x<0;
 y=exp(x);
 y(ii)=0;
end


function testsink
%%
n=100;
a=randn(n,1);
b=randn(n,1)*2+1;
W22=mean((sort(a)-sort(b)).^2)

W22s=sinkhorn3([],[],a,b,ones(n,1)/n,ones(n,1)/n)

%%
end

function testsink2
%%
ishigami
n=1024;
x=trafo(sobolpoints(n,k));
y=model(x);
D=(y-y').^2/2;
p=ones(n,1)/n;
M=10;
xt=trafo(zeros(1,k));
for m=1:M
  m
  xu=trafo(ones(1,k)*m/M);
  % select between xt and xu
  ii=x<=xu & x>xt;
  q=ii./sum(ii);
  for i=1:k
    W(i,m)=sinkhorn2(D,.1,p,q(:,i)); %no light speed here
  end
  %
 % tic;
 % W(:,m)=sinkhorn(D*2,4,p,q);toc
  xt=xu;
end
plot(1:M,W)
mean(sqrt(W),2)
%%
end