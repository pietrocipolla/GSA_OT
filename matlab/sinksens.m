function [c,results]=sinksens(x,y,M,entropy)
%% SINKSENS Optimal Transport for Sensitivity via Sinkhorn Algorithm
% SINKSENS(X,Y,M,P) returns the Wasserstein22 sensitivity and separations and
%  given the NxD matrix X as input and the NxL matrix Y as output.
% M denotes the number of partitions in X, P the entropic penalty

% written by elmar.plischke@tu-clausthal.de

chatty=true; % console output 
gfx=true; % graphics output
results=struct();

if(nargin<3 || isempty(M))
    M=32; % default parttition size (maybe modify for small sample size?)
end
if(nargin<4 || isempty(entropy))
    entropy=0; % default entropic parameter
end
 
[n,k]=size(x);

 if(chatty)
  disp('cost matrix construction');
  tic
 end
 C=-2*(y*y')+(sum(y.^2,2)+sum(y.^2,2)'); % Euclidean squared
 if(chatty), toc, end

ms=linspace(0,1,M+1); %linspace(0,n,M+1)
%[xs,ix]=sort(x);
xr=empcdf(x); % allow for ties

%disp('Sinkhorn');
W=zeros(M,k);
if(entropy==0)
    entropy=-0.001;
	if(chatty)
     fprintf('Entropic penalty %g\n', entropy);
	end
end
tic
for m=1:M
    if(chatty), m, end
    ii=xr>ms(m) & xr<=ms(m+1);
    N(m,:)=sum(ii); % realizations per partition
    for i=1:k
        [S,R]=sinkfastCost(entropy,C(:,ii(:,i))); 
        W(m,i)=S;
     end
end
t=toc;
if(chatty), t,end
c=sum(W.*N)/n;
results.separation=W;
results.binsize=N;
results.time=t;

if(gfx)
plot(1:M,W);title('Sinkhorn');drawnow
end
end