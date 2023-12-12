function [X,Y,tau,sigma]=wasserstein_quick(X,Y)
% Compute the Wasserstein distance
% Turn Puccetti 2017 into quicksort

n=size(X,1);
sigma=1:n;
tau=1:n;
X=X';
Y=Y';
converged=false;
% 
stack=zeros(ceil(n/2),2);
top=0;
% push initial values
top=top+1;stack(top,:)=[1 n];
while(top>0)
 % pop
 p=stack(top,1);r=stack(top,2);top=top-1;
 % start partition
 %randomization also should swap x?
 %s=unidrand(r-p+1)-1+p;
 %ss=sigma(s);sigma(s)=sigma(r);sigma(r)=ss;
 %ys=Ys(:,s); Ys(:,s)=Ys(:,r);Ys(:,r)=ys;
 % pivot element (last entry)  
  xr=X(:,r); yr=Y(:,r);
	
  i=p-1;
  for j=p:r-1
   if( xr-X(:,j))'*(yr-Y(:,j))<0 % Puccetti comparison
    i=i+1;
%    if(i==r), disp(sprintf('i==r in (i,j,p,r)=%d %d %d %d\n',i,j,p,r)), end
    if(i~=j), 
% swap i and j
    si=sigma(i);sigma(i)=sigma(j);sigma(j)=si;
    yi=Y(:,i); Y(:,i)=Y(:,j);Y(:,j)=yi;
% transitivity step: also exchange x
    ti=tau(i);tau(i)=tau(j);tau(j)=ti;
    xi=X(:,i); X(:,i)=X(:,j);X(:,j)=xi;
   end 
  end
 end
 q=i+1;
 if(q~=r), 
% swap q and r (antisymmetry step: just y)
  sq=sigma(q);sigma(q)=sigma(r);sigma(r)=sq;
  yq=Y(:,q); Y(:,q)=Y(:,r);Y(:,r)=yq;
 end
 % end of partition
 if(q-1>p), top=top+1;stack(top,:)=[p q-1]; end
 if(q+1<r), top=top+1;stack(top,:)=[q+1 r]; end
end
X=X';
Y=Y';
end


% Cormen et al.: Algorithms
function quicksort(A,p,r)
 if p<r
  q=partition(A,p,r)
  quicksort(A,p,q-1)
  quicksort(A,q+1,r)
 end
end

function j=partition(A,p,r)
 x=A(r);
 i=p-1;
 for j=p:r-1
  if(A(j)<=x)
   i=i+1;
   A([i j])=A([j i]);
  end
 end
 j=i+1;
 A([j r])=A([r j]);
end
 
%% iterative version using stack
function A=quicksort_iter(A,lo,hi)

stack=zeros(ceil((hi-lo)/2),2); % ceil((hi-lo)/2)
top=0;
% push initial values
top=top+1;stack(top,:)=[lo hi];
while(top>0)
% pop
p=stack(top,1);r=stack(top,2);top=top-1;
% partition
 x=A(r);
 i=p-1;
 for j=p:r-1
  if(A(j)<=x)
   i=i+1;
   A([i j])=A([j i]);
  end
 end
 q=i+1;
 A([q r])=A([r q]);
 % end of partition
 if(q-1>p), top=top+1;stack(top,:)=[p q-1]; end
 if(q+1<r), top=top+1;stack(top,:)=[q+1 r]; end
end
end
 

