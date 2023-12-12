%%
% create as in MultivariateGaussianNormal
%%
M=16;
[n,k]=size(X);
[nn,l]=size(Y);
ms=round(linspace(0,n,M+1));
[~,ix]=sort(X);
R=30;
for r=1:R
[~,ir]=sort(rand(size(X))); % random permutations
W=zeros(M,k);
for m=1:M
    m
     mc=ms(m+1)-ms(m);
     
    for i=1:k
    % conditional output
     yc=Y(ix(ms(m)+1:ms(m+1),i),:);
    % random pick
     yr=Y(ir(ms(m)+1:ms(m+1),i),:);    
     
     W(m,i)=sinkfast(1,yc,yr);
    end
end
plot(1:M,W);hold on;
set(gca,'ColorOrderIndex',1);
end
hold off