function [W2,p,x,A]=wasseropt2(mu,nu)
% wasserstein via optimal transport / linear programming

% costs = distance between mu and nu
C=(sum(mu.^2,2)+sum(nu.^2,2)'-2*mu*nu')';

m=size(mu,1);n=size(nu,1)
% Jarre Stoer: Constraints for double stochastic matrix
A=[kron(eye(m),ones(1,n));kron(ones(1,m),eye(n))];
try
% octave: lower bound is zeros  
[x,fopt]=glpk(C(:),A,ones(n+m,1),[],[],repmat('S',1,n+m),...
        repmat('C',1,n*m));
catch
%matlab: equality constraints
%lb=zeros(1,length(C));
%linprog(f,A,b,Aeq,beq,lb,ub,options) 
%,'Algorithm','dual-simplex'
[x,fopt]=linprog(C(:),[],[],A,ones(n+m,1),zeros(n*m,1),[]);
end
W2=fopt/n;
p=(reshape(x,m,n))'*(1:m)'; % permutation
end
function testwasseropt
%%
clc
close all
clear all
N=5;
%mu=randn(N,1);
%nu=randn(N,1)*2+1;
mu=xlsread('DataforExperiment.xlsx','A1:A5')
nu=xlsread('DataforExperiment.xlsx','B1:B5')

W2ana=(1-0)^2+(1-2)^2;
W2=mean((sort(mu)-sort(nu)).^2);
tic
[W21p,p,x,A]=wasseropt2(mu,nu);
W21p
A
CC=xlsread('DataforExperiment.xlsx','C7:G12')

toc
%%
end

function wasser_on_datasaurus
  %%
  nogfx=1;
  datasaurus
  d=size(xx,2);
  W=zeros(d,d);
  for k1=1:d-1
    for k2=k1+1:d
     [W2,p]=wasseropt2([xx(:,k1),yy(:,k1)],[xx(:,k2),yy(:,k2)]);
     plot([xx(:,k1),xx(p,k2)]',[yy(:,k1),yy(p,k2)]','*-k');
     axis([0 100 0 100]);
     W(k1,k2)=W2;
     pause(.5);
    end
  end
  %%
end
