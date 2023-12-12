% Ishigami 2D with modified parameters
trafo = @(u)(u-.5)*2*pi;
model = @(x)[sin(x(:,1)).*(1.0+0.1*x(:,3).^4)+7.0*(sin(x(:,2))).^2, ...
             sin(x(:,1)).*(2.0-0.3*x(:,3).^4)+5.0*(sin(x(:,2))).^2];
k=3;
n=4000;
x=trafo(rand(n,k));
y=model(x);

% some methods take their time ...
%l1=bigwassersteintest(x,y,32)

% just one
l2=bigwassersteintest(x,y,32,{'bures'})
drawnow
pause(2)
% the sinkhorns
l3=bigwassersteintest(x,y,32,{'sinkhorn','gradient','sinkmem'})

mean(l3.Wsinkhorn)/2
% in contrast, kernel based stuff
m1=mmd2si(x,y,32,'dot') % sum of variance, should be smaller than OT/2
% same as
s1=sum(cosi(x,y,-8))

q1=mmd2si(x,y,-32,'dot') % quick mode? is broken

m2=mmd2si(x,y,32,'ed') % Gini / Energy Distance
m3=mmd2si(x,y,32,'exp')