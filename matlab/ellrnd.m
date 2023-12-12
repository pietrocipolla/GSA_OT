function z=ellrnd(n,mu,Sigma,g)
% Random numbers from elliptic distribution
if nargin==0, testell(); end

d=size(mu,2);
if nargin<3 || isempty(Sigma), Sigma=eye(d); end
if nargin<4 || isempty(g), g=@(u)u; end
 
%u=rand(n,d+1);
u=sobolpoints(n+1,d+1);u(1,:)=[]; % [0.5 0.5] -> [0,0] -> norm 0 -> inv: NaN 
% without stats toolbox
norminv_=@(z)-sqrt(2)*erfinv(1-2*z);
x=norminv_(u(:,1:d)); x=x./sqrt(sum(x.^2,2)); % points on a sphere
r=g(u(:,end)); % radius

z=mu+(r.*x)*sqrtm(Sigma); % chol(Sigma); %

end

function testell
%% 
 try
 pkg load signal
 pkg load statistics
 end

% x=ellrnd(200,[-1,1],[3,-1;-1,1]);
 x=ellrnd(160,[-1,1],[3,-1;-1,1],@(u)betainv(u,4,1));
 y=ellrnd(240,[1,-1],6*[1,-1;-1,3],@(u)betainv(u,4,1));
 
 plot(x(:,1),x(:,2),'*',y(:,1),y(:,2),'+');
 
 [W,R]=sinkfast(.008,x,y);W 
 Cx=cov(x);
 Cy=cov(y);Ry=sqrtm(Cy);
 W2=sum((mean(x)-mean(y)).^2)+trace(Cx+Cy-2*sqrtm(Ry*Cx*Ry)) 
 
 n=length(x);m=length(y);mn=gcd(n,m);
 % 
 C=-2*x*y'+sum(x.^2,2)+sum(y.^2,2)'; % distance squared
 tic
 P=transsimp(ones(n,1)*m/mn,ones(m,1)*n/mn,C); % enforce integer solution
 toc
 W3=(P(:)'*C(:))/lcm(n,m)
 
 hold on
 for i=1:n
   p=find(P(i,:));
   if(~isempty(p))
    plot([x(i,1)*ones(length(p),1),y(p,1)]',...
    [x(i,2)*ones(length(p),1),y(p,2)]','g-');
  end
 end
 hold off
%% 
end