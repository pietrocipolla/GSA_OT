function [R,S,Es]=sinkhornranksort(x,a,m)
% SINKHORNRANKSORT Continuous ranking and sorting
% Cuturi et al. 2019

epsl= 5e-3;
eta=  1e-3;
h=@abs;      % for W1
% h=@(x)x.^2 % for W2
y=linspace(1/(2*m),1-1/(2*m),m)';
b=ones(m,1)/m;

% rescale and squash x
x=x-mean(x,1);x=x./std(x,1);
%x=1./(1+exp(-x)); % (1+tanh(x))/2; % arctan or logistic map
x=.5+x/2.96;
[n,s]=size(x);
%% assert [n,s]==size(a), all(sum(a))==1
if any(sum(a,1)~=1), error('a must be probability mass vector.'); end
if any([n s]~=size(a)), error ('x and a: dimension mismatch.'); end

softMin=@(M,s)-s*log(sum(exp(-M/s),2));

elb=epsl*log(b);

%if(s>1), error('Sorry, unimplemented'); end
R=zeros(n,s);
S=zeros(m,s);
% parfor

for t=1:s
 C=h(x(:,t)-y');
 ela=epsl*log(a(:,t));

 alfa=zeros(n,1);
 beta=zeros(m,1);
%
converged=false;
iter=0;
D=C-alfa;

while ~converged
 beta = elb+softMin(D'-beta,epsl)+beta;
 alfa = ela+softMin(D-beta',epsl)+alfa;
 D = C-alfa; 
 
 E = exp((beta-D')/epsl);
 err=mean(abs(sum(E,2)-b));
% 
 iter=iter+1;
% if(~mod(iter,10))
%  [iter, err, ...
%  ... once divided by epsilon, their exponentials sum to one...
%  sum(sum(exp(-(D'-beta)/epsl))) ]
% end
 converged = err<eta;
end
iter
 R(:,t)=(E'*cumsum(b))./a(:,t);
 S(:,t)=(E*x(:,t))./b;
 Es{t}=E;
end 
end

%%
function titanic
% test sinkhornranking for sensitivity
n=8192*2; %1024;
ishigami
x=trafo(sobolpoints(n,k));
%x=trafo(rand(n,k));
y=model(x);
%%
plt=0;
for M=[8,16,32,64,128,256]
Y=zeros(M,k);
    [R,S,E]=sinkhornranksort(x,ones(n,k)/n,M);
for i=1:k
    %subplot(2,k,i);  plot(R,y,'.');
    %subplot(2,k,k+i);plot(1:M,M*E*y,'.');
    Y(:,i)=M*E{i}*y;
end
plt=plt+1;subplot(2,3,plt);plot(1:M,Y); title(['M=' num2str(M)]);
a=axis;a(1)=1;a(2)=M;axis(a);drawnow
mean((Y-mean(y)).^2)/var(y)
sum(S.*(Y-mean(Y)).^2)./(var(y).*sum(S))
end
 
%%
% same idea with sinkfast? expand to higher dimensions?
%M=128;
%z=trafo(mhalton(M,k)); Mingled Halton
%z=trafo(sobolpoints(M,k));

M=257; % uniform design
V=mod((1:M)'*[1,26,99,113],M);z=trafo((V+.5)/M);
Y=zeros(M,2^k-1);
plt=0;
%s=sqrt(2^k);
 str='1234'
for i=1:2^k-1 % 1:2^k  ... 
  alfa=logical(bitget(i,1:k))
[W22,u,v,K]=sinkfast([],x(:,alfa),z(:,alfa));
P=u.*K.*v'; % transposed
%R=(P*linspace(1/M,1,M)')*n;
%subplot(2^k,ceil(s),round(s));
Y(:,i)=M*P'*y;
if(sum(alfa)==1)
 plt=plt+1;
 subplot(k/2,k+1,plt);
 plot(z(:,alfa),Y(:,i),'.'); title(str(alfa));
elseif(sum(alfa)==2)
 plt=plt+1;
 subplot(k/2,k+1,plt);
 zz=z(:,alfa);plot3(zz(:,1),zz(:,2),Y(:,i),'.');title(str(alfa));
end
end
tau=mean((Y-mean(Y)).^2)/var(y)
%% Fast Yates Trafo
fyt([0,tau]',[1,0;-1,1])'
%%
end
    
