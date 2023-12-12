function [W22,retargs]=sinkfast2(l,x,y)
% SINKFAST symmetric Sinkhorn Knight iteration on EXP(-C/L) with C=||X-Y||^2.
% [W,R]=SINKFAST(L,X,Y) returns the Euclidean Wasserstein 
%   distance (squared) in W and further results in R.
%   Use fieldnames(R) for a list of available data.

% written by elmar.plischke@tu-clausthal.de
verbose=true; % false; % 
meanshift=false;
n=size(x,1); % assume uniform weights
m=size(y,1);
maxerr=1e-3;maxiter=1000;
if(meanshift)
mx=mean(x);my=mean(y); % center 
x=x-mx;y=y-my;
else
mx=0;my=0;
end
Cneg=2*x*y'-sum(x.^2,2)-sum(y.^2,2)';

% entries in K will be wiped out by exp(- * /l)
lminr=-min(max(Cneg,[],2))/log(realmax); %row wipeout
lminc=-min(max(Cneg,[],1))/log(realmax); %column wipeout
rowcritical=lminr>lminc;
if(verbose)
 disp(['Critical lambdas ' num2str(lminr) ' (row) ' num2str(lminc) ' (col)']);
end
if isempty(l), l=max(lminr,lminc)*1.11; end

K=exp(Cneg/l);%clear('Cneg');
uv=ones(n+m,1);
err=inf;iter=0;

% old for comparison
%u=ones(n,1); 
%Ktu=K'*u;
 Ktu=K'*uv(1:n);
while(iter<2 || isfinite(err) && err>maxerr && iter <maxiter)
 vu=uv;
 uv=1./[n*K*uv(n+1:end); m*Ktu]; % m*K'*uv(1:n)];
% uv(isinf(uv))=1e10;  % avoid 0*inf = nan
  Ktu=K'*uv(1:n);
 if mod(iter,2)==0
   % old for comparison
%   v=(1/m)./(K'*u); % == second component in even iters
 
   err=sum(abs(vu(n+1:end).*Ktu-1/m));
% else
%   u=(1/n)./(K*v); % == first component in odd iters
    
 end

 iter=iter+1;
if(verbose)
 disp(sprintf('iteration %d, error %f stddev ratio %g\n',iter,err, std(uv./vu)));
end
end
%uv=(uv+1./[n*K*uv(n+1:end); m*K'*uv(1:n)])/2;
uv=(uv+vu)/2; % or a weighted mean?
% sqrt <P,C>
% l*(mean(log(u))+mean(log(v)))
% scalings
u=uv(1:n);v=uv(n+1:end);
% potentials
f=l*log(u);g=l*log(v);
meanSqDiff=sum((mx-my).^2);
% cost (dual)
W22=meanSqDiff+(mean(f)+mean(g));

P=v'.*K.*u; % = Gamma exp((f+g'-C)/l)
% but sum(P,1), sum(P,2) != 1 
SinkhornDual=W22-l*mean(P(:)); % Peyre
% cost (primal) <P,C>
W2prime=meanSqDiff-(P(:)'*Cneg(:))/(m*n);
SinkhornPrimal=W2prime+l*P(:)'*(log(P(:)+1e-20)-1)/(m*n); % Peyre

% put all intermediates and results in a struct
 retargs = struct('W2Dual',W22,'W2Primal',W2prime,'Kernel',K,...
           'Cost',-Cneg,'Coupling',P,'ScalingU',u,'ScalingV',v,...
           'PotentialF',f,'PotentialG',g,'MeanSqShift',meanSqDiff,...
           'NumIterations',iter','LastError',err,...
           'Lambda',l,'LambdasCritical',[lminr,lminc],...
           'SinkhornPrimal',SinkhornPrimal,'SinkhornDual',SinkhornDual);
end

function testsink
%%
n=100;
a=randn(n,1);
b=randn(n,1)*2+1;
W22a=(mean(a)-mean(b)).^2+(std(a)-std(b)).^2
W22=mean((sort(a)-sort(b)).^2)

[W22s,R]=sinkfast2([],a,b);W22s

%%
end