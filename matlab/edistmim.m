function [e,E]=edistmim(x,y,M)
% Moment indep. measure using Cramer/van Mises energy distance
[n,k]=size(x);
if(nargin<3), M=32; end
%[yr,ind]=sort(y);yr(ind)=1:n; % output ranks
[~,ind]=sort(x);

Ms=round(linspace(0,n,M+1));
E=zeros(k,M);
%F=zeros(k,M);
vec=@(M)M(:);
%B=abs(bsxfun(@minus,y,y'));
% Multivariate output: Kernel is Euclidean 
y2=sum(y.^2,2);
B=sqrt(max( (y2+y2')-2*(y*y'),0));
EB=mean(vec(B));
%EBR=(n^2-1)/(3*n^2);
for m=1:M
    I=ind((Ms(m)+1):Ms(m+1),:);
    for i=1:k
        nc=sum(I(:,i));
        EA=mean(vec(B(:,I(:,i))));
        EC=mean(vec(B(I(:,i),I(:,i))));
%        E(i,m)=2*EA-EB-EC;
         E(i,m)=2*EA-EB*n/(n-1)-EC*nc/(nc-1); % debiasing
%        JR=yr(I(:,i));
%        EAR=.5+1/(2*n^2*length(JR))*sum(2*JR.*(JR-1)-(2*JR-1)*n);
%        %EAR=1/(2*n^2*length(JR))*sum((JR-1).*JR+(n-JR).*(n-JR+1));
%        CR=abs(bsxfun(@minus,JR,JR'))/n;ECR=mean(CR(:));
%        F(i,m)=2*EAR-EBR-ECR;
    end
end
%
%plot((Ms(1:end-1)+diff(Ms)/2)/n,E/2);
%
e=sum(bsxfun(@times,E,diff(Ms)/n),2)'/2;
%f=sum(bsxfun(@times,F,diff(Ms)/n),2)'/2;
end
