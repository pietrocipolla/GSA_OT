function [W2,p]=wasseropt(mu,nu)
% wasserstein via optimal transport / linear programming

% costs = distance between mu and nu
C=(sum(mu.^2,2)+sum(nu.^2,2)'-2*mu*nu')';

m=size(mu,1);n=size(nu,1);
% Jarre Stoer: Constraints for double stochastic matrix
A=[kron(eye(m),ones(1,n));kron(ones(1,m),eye(n))];
try
% octave: lower bound is zeros  
[x,fopt]=glpk(C(:),A,ones(n+m,1),[],[],repmat('S',1,n+m),...
        repmat('C',1,n*m));
catch
%matlab: equality constraints
[x,fopt]=linprog(C(:),[],[],A,ones(n+m,1),zeros(n*m,1),'Algorithm','dual-simplex');
end
W2=fopt/n;
p=(reshape(x,m,n))'*(1:m)'; % permutation
end
function testwasseropt
%%
mu=randn(50,1);
nu=randn(50,1)*2+1;
W2ana=(1-0)^2+(1-2)^2;
W2=mean((sort(mu)-sort(nu)).^2);
W2lp=wasseropt(mu,nu)
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
     [W2,p]=wasseropt([xx(:,k1),yy(:,k1)],[xx(:,k2),yy(:,k2)]);
     plot([xx(:,k1),xx(p,k2)]',[yy(:,k1),yy(p,k2)]','*-k');
     axis([0 100 0 100]);
     W(k1,k2)=W2;
     pause(.5);
    end
  end
  %%
end
