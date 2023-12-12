function [ xc ] = empcdf( x )
%EMPCDF Computes the empirical cumulative distribution function
%   C=EMPCDF(X) assigns ranks (including ties)

[n,k]=size(x);
xc=zeros(n,k);
for j=1:k
 [xs,indj]=sort(x(:,j));
 xr=(1:n)';						 % ranks	
 tie_loc=[find(diff(xs)==0);n+2]; % stopper
 tie_next=diff(tie_loc);
 maxt=numel(tie_loc);
 i=1;while(i < maxt)
	run=tie_loc(i);len=1;
	while(tie_next(i)==1), i=i+1;len=len+1; end
	xr( run : run + len) = run+len/2;
	i=i+1;
 end	
 xc(indj,j)=(xr)/n;
end
end

