function [W,D,E]=wassgrad(x,C,M,l)
% WASSGRAD Entropic Wasserstein Sensitivity with a gradient method.
% [W,D,E]=WASSGRAD(X,C,M,L)
% INPUTS
% X input matrix (N x K), 
% C cost matrix (N x N) or output vector (N x I, I~=N)
% M partition size
% L scaling factor for entropy, negative for relative to max(cost)
% OUTPUTS
% W primal separation matrix (M x K)
% D dual separation matrix ( M x K)
% E convergence info (1,2: converged, -1: diverged)
 [n,k]=size(x);
 if(size(C,2)~=n)
 % assume vector/multivariate output and no cost matrix
  C=-2*(C*C')+(sum(C.^2,2)+sum(C.^2,2)'); % Euclidean squared
 end
 if nargin<3 || isempty(M)
  M=ceil(n^(1/3));
 end
 if(nargin<4 || isempty(l))
     l=max(C(:))*.002; %0.05 % *.0015 exploits floating point range.
 end
 if(l<0), l=max(C(:))*(-l); end % interpret negative as relative 
 ms=round(linspace(0,n,M+1));
 [~,ix]=sort(x);

 K=exp(-C/l); % 
 W=zeros(M,k);
 D=zeros(M,k); % for duals
 E=zeros(M,k); % for convergence / error
 for m=1:M
     for i=1:k
      ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
      % Cohen Madry Tsipras Vladu 
      % MM=@(x,y)K(:,ii).*exp(x-y');
      % g=@(x,y)sum(MM(x,y),'all')-(mean(x)-mean(y));
      % rM=@(M)sum(M,2);cM=@(M)sum(M,1)';
      % gradM=@(M)[rM(M)-1/n; -(cM(M)-1/nn)];
      % dg=@(x,y)gradM(MM(x,y));
      % hessM=@(M)[diag(rM(M)),-M;-M',diag(cM(M))];
      % d2g=@(x,y)hessM(M(x,y));
       
      u=zeros(n,1);v=zeros(nn,1);
      P=K(:,ii).*exp(u-v');
      muv=(mean(u)-mean(v));
      z=sum(P,'all')-muv;
      while(true)
       Ju=sum(P,2)-1/n;
       Jv=-(sum(P,1)'-1/nn);
       h= 1/sqrt(Ju'*Ju+Jv'*Jv); % 1
       % avoid step size getting too large (all J entries <.001)
       if(h>1e6/(n+nn)),E(m,i)=1;break; end
       z1=sum(P.*exp(h*(Jv'-Ju)),'all')-muv-h*(mean(Jv)-mean(Ju));
       z2=z;
       while(z1>=z)
        h=h*.95;  % h too large
        z1=sum(P.*exp(h*(Jv'-Ju)),'all')-muv-h*(mean(Jv)-mean(Ju));
       end
       while(z1<z && z1<z2)
        z2=z1;h2=h;
        h=h*1.5;
        z1=sum(P.*exp(h*(Jv'-Ju)),'all')-muv-h*(mean(Jv)-mean(Ju));
       end
       % quadratic line search
       h0=0;c0=z;h1=h;
       c1=(z1-c0)/(h1-h0);
       c2=((z2-c0)/(h2-h0)-c1)/(h2-h1);
       h0=.5*((h0+h1)-c1/c2);
       % 
       u=u-h0*Ju;v=v-h0*Jv;
       P=K(:,ii).*exp(u-v');
       muv=(mean(u)-mean(v));
       if(~isfinite(muv)), E(m,i)=-1; break; end
       z=sum(P,'all')-muv;
       if(abs(z-c0)<=1e-4), E(m,i)=2; break; end % no change
      end
      % entropic OT, maybe corrected for Wasserstein?
      W(m,i)=sum(P.*C(:,ii),'all');
      D(m,i)=muv;
     end
 end
end
%%
function testishi
%%
 ishigami
 n=1024;
 x=trafo(sobolpoints(n,k));
 y=model(x);
 l=bigwassersteintest(x,y,16,2^6+2);
 [m,d]=wassgrad(x,y,16);
%%
end
%%
function oldexperiments
%%
 b=[1,5,2,8,2]';
 a=[3,8,1,6]';
 C=[3 4 6 8 9; 2 2 4 5 5; 2 2 2 3 2; 3 3 2 4 2];
 
for l=2:10 %2:5:50
 K=exp(-C*l); % lambda=1/l
 
%Sinkhorn 
 u=ones(size(K,1),1);
for i=1:(50+l)
 v=b./(K'*u);
 u=a./(K*v);
end
P=diag(u)*K*diag(v);
% optimal 61

ff=1/l*log(u./a);
gg=1/l*log(v./b);
[P(:)'*C(:), a'*ff+b'*gg]
end

%%
% Allen-Zhu, Li, Oliveira, Wigderson
f=@(x)a'*log(sum(K.*exp(x'),2))-b'*x;
% derivative of logsumexp is softmax
df=@(x)((a'*K).*exp(x'))'./(sum(K,1).*exp(x'))'-b;

% Cohen Madry Tsipras Vladu
% sum(...,'all')  
g=@(x,y)sum(sum(K.*exp(x-y')))-(a'*x-b'*y);
 
M=@(x,y)K.*exp(x-y');
rM=@(M)sum(M,2);cM=@(M)sum(M,1)';

gradM=@(M)[rM(M)-a; -(cM(M)-b)];
dg=@(x,y)gradM(M(x,y));
hessM=@(M)[diag(rM(M)),-M;-M',diag(cM(M))];
d2g=@(x,y)hessM(M(x,y));

%optimal values from Sinkhorn
x=log(u);y=-log(v);
g(x,y), dg(x,y)
%%
x=zeros(size(a));y=zeros(size(b));
%for l=1:3
h=.5;
z=g(x,y);
n=size(a,1);
m=size(b,1);
for i=1:600
 J=dg(x,y);
 p=x-h*J(1:n);
 q=y-h*J(n+1:end);
 r=g(p,q);
 if abs(r-z)<10*eps, break; end
 while(r>z)% if while
  h=h*.9; % h=h*.1 -> full brake
  p=x-h*J(1:n);
  q=y-h*J(n+1:end);
  r=g(p,q);
  %h=h*.95;
  %p=x-h*J(1:n);
  %q=y-h*J(n+1:end);
  %r=g(p,q);
 %else
 end
  h=h*1.2; % 1.5
 %end
 z=r; 
 if(i>400), [i z,h], end
 mp=max(p);
 x=p-mp; % fix first coordinate
 y=q-mp;
 
end
%end
%% second ordner method with correction for non-strict convexity
y=-log(min(C,[],1))';x=log(min(C,[],2));
xy=[x;y];
for i=1:10
xy=xy-(1/(n+m)+d2g(xy(1:n),xy(n+1:end)))\dg(xy(1:n),xy(n+1:end));
g(xy(1:n),xy(n+1:end))
end

%% line search?
g(x,y)
J=dg(x,y)
H=d2g(x,y);
h=-(J'*J)/(J'*H*J);
x=x+h*J(1:n);y=y+h*J(n+1:end);
g(x,y)
%%
end