function [OT,retargs]=sinkhorn5(C,lambda,x,y,alfa,beta)
% Compute Optimal Transport via symmetric Sinkhorn algorithm
%
% Inputs: C       cost function c(x,y) (default squared Euclidean)
%	        lambda  regularization
%         x,y     n,m sample location
%         p,q     n,m weight vectors (uniform defaults)	
% partly taken from Jean Feydy 2020

 maxiter=1000;
 miniter=5;
 maxerr=1e-5;
 verbose=false; % true;
 
% assume: 
% sum(p)==sum(q)==1
% size(c)==[n,n], all(all(c>=0)), size(p)==[n,1],size(q)==[n,1]
 n=size(x,1);
 m=size(y,1);
 if(nargin<5)
  alfa=ones(n,1)/n; % uniform masses
  beta=ones(m,1)/m;
 end
 % recomputing may be faster than storing
 if isempty(C)
   C=@(x,y)(sum(x.^2,2)+sum(y.^2,2)'-2*x*y')/2;
 end
 cxy=C(x,y);%nc=norm(cxy)
 cxx=C(x,x);
 cyy=C(y,y);
 
 Delta=max(sqrt(cxy(:))) % diameter
 
 if(isempty(lambda))
  alfai=min([alfa;beta]); 
  if(alfai>0)
   lambda=-1/10/log(alfai);
  else
   lambda=max(cxy(:));
  end
  disp(['SINKHORN: lambda set to ' num2str(lambda) ]);
 end

 f = cxy*beta;
 g = (alfa'*cxy)';
 ff = cxx*alfa;
 gg = cyy*beta;
 
% v0=zeros(n,1);
% u0=zeros(n,1);

 converged=false;
 Tlold=0;
 Tl=inf();
 iter=0;
 while ~converged
    f_ =.5*(f-lambda*log( exp((g'-cxy)/lambda)*beta));
    g  =.5*(g-lambda*log( alfa'*exp((f-cxy)/lambda))');
    f  = f_; % fake parallel execution
    ff =.5*(ff-lambda*log( exp((ff'-cxx)/lambda)*alfa));
    gg =.5*(gg-lambda*log( beta'*exp((gg-cyy)/lambda))');
  iter=iter+1;
  Tlold=Tl;
  Tl=2*((f-ff)'*alfa+(g-gg)'*beta);
  err=abs(Tl-Tlold);
  if verbose,   
    disp(sprintf('error %g, iter %d',err, iter)); 
  end
  % relative error dependent on cost range
  converged=(iter>miniter && err<maxerr*Delta) || isnan(Tl) || iter==maxiter;
 end
 if(isnan(Tl))
     OT=Tlold;
     T2=0;
 else
     OT=2*((f-ff)'*alfa+(g-gg)'*beta);
     T2=2*(f'*alfa+g'*beta); % without divergence correction
 end
 if(nargout>1)
 %% factor 2 missing everywhere
  K=exp(-cxy/lambda);
  P=beta'.*exp((f+g'-cxy)/lambda).*alfa;
  W2prime=P(:)'*cxy(:);
  SinkhornPrimal=W2prime+lambda*sum(P(:).*(log(P(:)+1e-20)-1)); 
  SinkhornDual=T2-lambda*sum(P(:)); 
   retargs = struct('SinkhornDivergence',OT,'W2Dual',T2,
      'W2Primal',W2prime,'Kernel',K,'Lambda',lambda,...
      'PotentialF',f,'PotentialG',g,'PotentialFF',ff,'PotentialGG',gg,...
      'NumIterations',iter','LastError',err,...
      'Cost',cxy,'Coupling',P,...
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