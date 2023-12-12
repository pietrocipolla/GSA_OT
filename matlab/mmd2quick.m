% investigate broken quick mode
M=16;

 ranks=zeros(n,k);
 [~,xi]=sort(x);
 for i=1:k
     ranks(xi(:,i),i)=(1:n)'; 
 end
 ri=ceil(ranks*M/n); % binned ranks
 means=zeros(k,M);
 for j=1:k
  means(j,:)=accumarray(ri(:,j),1)./n; % bin count
 end

 maxRank=max(ri,[],'all');

seps=zeros(k,maxRank);
Ki=zeros(1,k);

kern=@(y,z)sum((y-z).^2,2);

   ramp=(1:n)';
   for j=1:k
    % an ugly kludge to process matrices
    mm = accumarray(ri(:,j),ramp,[], @(z)Ekernxx(y,y(z,:),kern));
    seps(j,:)=mm;
    Ki(j)=means(j,:)*mm;
   end

function m=Ekernquick(x,y,kern)
 % Gretton et al. quick (linear) version
 % if x and y are indep. and of same size
 n=size(x,1);
 l=size(y,1);
 % if sizes don't match, use random subset selection
 if(n<l)
     y=y(randperm(l,n),:); % select n out of l
 elseif(n>l)
     x=x(randperm(n,l),:);
     n=l;
 end
 m=0;
 for i=1:n-1
     m=m+kern(x(i,:),x(i+1))+kern(y(i,:),y(i+1)) ...
        -kern(x(i,:),y(i+1))-kern(y(i,:),x(i+1));
 end
 % close the loop
 m=m+kern(x(n,:),x(1))+kern(y(n,:),y(1)) ...
        -kern(x(n,:),y(1))-kern(y(n,:),x(1));
 m=m/n; 
 % in contrast to the even/odd version in Gretton et.al. this 
 % moving average version is 1-dependent which might pose problems in
 % bootstrapping
end
 
function m=Ekernxx(x,y,kern)
 n=size(x,1);
 l=size(y,1);
 % if sizes don't match, use random subset selection
 if(n<l)
     y=y(randperm(l,n),:); % select n out of l
 elseif(n>l)
     x=x(randperm(n,l),:);
     n=l;
 end
 m=0;
 for i=1:n
     for j=1:n 
     m=m+kern(x(i,:),x(j))+kern(y(i,:),y(j)) ...
        -kern(x(i,:),y(j))-kern(y(i,:),x(j));
     end
 end
 m=m/n^2; 
 % in contrast to the even/odd version in Gretton et.al. this 
 % moving average version is 1-dependent which might pose problems in
 % bootstrapping
end