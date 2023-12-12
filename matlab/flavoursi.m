function CMi=flavoursi(x,y,MP,gfx)
% the flavour method for CvM
% inspired by Gamboa, Klein, Lagnoux
[n,k]=size(x);
ys=sort(y);
%%
if((nargin<=2)||isempty(MP))
P=40;
M=10;
else
P=MP(2);
M=MP(1);
end
yy=[];
ps=linspace(1/2/P,1-1/2/P,P);
for p=ps % quantiles
    yy=[yy,y<ys(ceil(n*p))];
%    yy=[yy,y.*(y<ys(ceil(n*p/P)))];
end
try
 [ViT,V]=cosi(x,yy,M,'Unscaled',true); % non-scaled, multi-output
catch
 % someone stole signal toolbox licenses ...
 [ViT,V]=easi(x,1.0*yy,M,'Unscaled',true); 
end
CMi=6*mean(ViT);       % normalization
%%
if(nargin>3)
clf; plot(ps,6*ViT);
xlabel('Cutoff Output Quantile');
ylabel('Reliability Variance');
title([ gfx ': Flavour Method for Gini/CvM']);
legend(strcat('x_{',num2str((1:k)'),'}'))
end
%%
end
