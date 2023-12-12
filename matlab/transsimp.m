function [Popt,fopt]=transsimp(a,b,C)
% Transportation simplex as of Luenberger Yu, Jarre Stoer

[n,m]=size(C);
verbose= false; % true; % 
initfeas='Vogel'; % 'NW2'; %  'NW'; % 'minRow'; % 'minCol'; % 'Nozicka'; % 
assert(size(a,1)==n && size(b,1)==m,'mass and const matrix dimension mismatch')
assert(sum(a)==sum(b), 'No feasible solution for the transport problem.');

% preallocate working arrays
iB=zeros(n+m-1,2); % indices of basis
sP=zeros(n+m-1,1);   % sparse solution
sC=zeros(n+m-1,1);   % sparse costs (unused)
u=nan(n,1);v=nan(m,1); % simplex multipliers
known=false(m+n-1,2); % for backsubstitution
iCC=zeros(length(iB),1); % indices for cycle of change / stepping stones

%% Step 1. Compute an initial basic feasible solution using the Northwest Corner
% Rule or some other method.
switch initfeas
case 'NW'
%
i=1;j=1;
supply=a(i);demand=b(j);

for(cnt=1:(n+m-1))
 iB(cnt,:)=[i,j];
  if(supply>demand) 
   sP(cnt)=demand; %P(i,j)=demand; 
   supply=supply-demand;
   if(j<m), j=j+1; end, demand=b(j); 
  else 
   sP(cnt)=supply; %P(i,j)=supply; 
   demand=demand-supply;
   if(i<n), i=i+1; end, supply=a(i); 
  end
end
case 'NW2'
% Northwestern Rule, alternative
i=1;j=1;
supply=b(j);demand=a(i);

for(cnt=1:(n+m-1))
 iB(cnt,:)=[i,j];
  if(supply>demand) 
   sP(cnt)=demand; 
   supply=supply-demand;
   if(i<n), i=i+1; end, demand=a(i); 
  else 
   sP(cnt)=supply; 
   demand=demand-supply;
   if(j<m), j=j+1; end, supply=b(j); 
  end
end
case 'minRow'
% Dempe/Schreier: Wähle die ungestrichene Zeile mit dem kleinsten Index.
%   Wähle in dieser Zeile ein ungestrichenes Feld mit dem kleinsten
%   Kostenkoeffizienten.
%
as=a;bs=b;
n0=1:n;
m0=1:m;

for(cnt=1:(n+m-1))
 i=n0(1); % first unprocessed index
 [~,k]=min(C(i,m0));
 j=m0(k);
 aa=as(i);bb=bs(j);
 if(aa<bb)
  iB(cnt,:)=[i,j];
  sP(cnt)=aa;
  bs(j)=bb-aa;
  as(i)=0;
  n0(n0==i)=[];
 else
  iB(cnt,:)=[i,j];
  sP(cnt)=bb;
  bs(j)=0;
  as(i)=aa-bb;
  m0(m0==j)=[];
 end
end

case 'minCol'
%
as=a;bs=b;
n0=1:n;
m0=1:m;

for(cnt=1:(n+m-1))
 j=m0(1); % first unprocessed index
 [~,k]=min(C(n0,j));
 i=n0(k);
 aa=as(i);bb=bs(j);
 if(aa<=bb)
  iB(cnt,:)=[i,j];
  sP(cnt)=aa;
  bs(j)=bb-aa;
  as(i)=0;
  n0(n0==i)=[];
 else
  iB(cnt,:)=[i,j];
  sP(cnt)=bb;
  bs(j)=0;
  as(i)=aa-bb;
  m0(m0==j)=[];
 end
end
%{
case 'Nozicka'
% Nozicka et al. 1972 
[~,iI]=sort(C(:));
as=a;bs=b;
bcnt=0; % base entry
dcnt=0; % degenerated base entry
top=n+m;
l=1; 
% basic feasible solution has size m+n-1
while(any(as)) % || any(bs) )
  assert(l<=m*n,'undistributed mass at end of loop')
  assert(bcnt<top,'base vector overflow')
  assert(sum(as)==sum(bs),'loss of mass')
  i=1+mod(iI(l)-1,n); % vector to matrix coordinates
  j=1+(iI(l)-i)/n;
  l=l+1;
  
  aa=as(i);
  bb=bs(j);
%  if(aa==0 || bb==0), continue; end % nothing to redistribute
  d=aa-bb;
%  sP(l)=min(aa,bb); %abs(d);
%  iB(l,:)=[i,j];
  if(d>0)
   as(i)=d;
   bs(j)=0;
   if(bb==0) % degenerate point: conditionally add to top
    dcnt=dcnt+1;
    if(top>bcnt+dcnt) iB(top-dcnt,:)=[i,j]; end
   else %% aa>bb>0   
    bcnt=bcnt+1;
    sP(bcnt)=bb; % add to base
    iB(bcnt,:)=[i,j];    
   end
  elseif(d<0)
   as(i)=0;
   bs(j)=-d;
   if(aa==0)
    dcnt=dcnt+1;
    if(top>bcnt+dcnt) iB(top-dcnt,:)=[i,j]; end
   else %% aa>bb>0   
    bcnt=bcnt+1;
    sP(bcnt)=aa; %abs(d);
    iB(bcnt,:)=[i,j];    
   end  
  else
   as(i)=0;
   bs(j)=0;
   if(aa==0)
    dcnt=dcnt+1;
    if(top>bcnt+dcnt) iB(top-dcnt,:)=[i,j]; end
   else %% aa>bb>0   
    bcnt=bcnt+1;
    sP(bcnt)=aa; 
    iB(bcnt,:)=[i,j];    
   end  
  end
end
%}
%{
case 'Vogel'
%% Vogel Approximation Method
 as=a;bs=b;c=C;
 bcnt=0; % # base entries
 rowsUsed=0;
 colsUsed=0;
 l=1; 
 while rowsUsed+colsUsed<n+m-1 % (any(as))
  [rmax,j]=max(diff(mink(c,2,1))); % minimal and minimal-but-one
  [cmax,i]=max(diff(mink(c,2,2)'));
  if(rmax>=cmax)
   [~,i]=min(c(:,j));
  else
   [~,j]=min(c(i,:));
  end

  aa=as(i);bb=bs(j);
  d=aa-bb;
  if d>0
   as(i) = d;
   bs(j) = 0;
 
    bcnt=bcnt+1;
    sP(bcnt)=bb; 
    iB(bcnt,:)=[i,j];    

    c(:,j)= nan; % uncomparable
    colsUsed=colsUsed+1;
%   visitedj(j)=-inf;
  else %if d<0
    as(i) = 0;
    bs(j) =-d;

    bcnt=bcnt+1;
    sP(bcnt)=aa;
    iB(bcnt,:)=[i,j];    

    c(i,:)= nan;
    rowsUsed=rowsUsed+1;
%   visitedi(i)=-inf;
%  else
%   as(i) = 0;
%   bs(j) = 0;
% 
%    bcnt=bcnt+1;
%    sP(bcnt)=bb; 
%    iB(bcnt,:)=[i,j];    
%
%   c(i,:)= inf;
%   c(:,j)= inf;
  end
 end
% Dempe/Schreier: Existiert nur noch eine ungestrichene Zeile oder eine ungestrichene
% Spalte, dann gehe zur Nordwesteckenregel über.
 disp('debug me');
%} 
case 'Vogel'
%% Vogel Approximation Method, second try
 as=a;bs=b;%c=C;
 bcnt=0; % # base entries
 rowsUsed=0;
 colsUsed=0;
 n0=(1:n)';
 m0=(1:m)';

 while rowsUsed<n-1 && colsUsed<m-1 % (any(as))
  c=C(n0,m0);
  [rmax,l]=max(diff(mink(c,2,1)));   % minimal and minimal-but-one
  [cmax,k]=max(diff(mink(c,2,2)')); 
  if(rmax>=cmax)
   [~,k]=min(c(:,l));
  else
   [~,l]=min(c(k,:));
  end
  i=n0(k);j=m0(l);
  
  aa=as(i);bb=bs(j);
  d=aa-bb;
  if d>0
   as(i) = d;
   bs(j) = 0;
 
    bcnt=bcnt+1;
    sP(bcnt)=bb; 
    iB(bcnt,:)=[i,j];    

    m0(m0==j)= [];
    colsUsed=colsUsed+1;
  else %if d<0
    as(i) = 0;
    bs(j) =-d;

    bcnt=bcnt+1;
    sP(bcnt)=aa;
    iB(bcnt,:)=[i,j];    

    n0(n0==i)= [];
    rowsUsed=rowsUsed+1;
  end
 end
% Dempe/Schreier: Existiert nur noch eine ungestrichene Zeile oder eine ungestrichene
% Spalte, dann gehe zur Nordwesteckenregel über.
 % either n0 or m0 will be a scalar
 iB(bcnt+1:end,:)=[n0*ones(size(m0)),m0*ones(size(n0))];
end % select method
if(verbose)
disp('Initial feasible solution found.')
[iB,sP]'
%spy(sparse(iB(:,1),iB(:,2),sP))
end
%% Step 2. Compute the simplex multipliers and the relative cost coefficients. 
  u(:)=nan;v(:)=nan;
while (true)
  % Simplex multipliers
  % ********************************************* 
  % Step 1. Assign an arbitrary value to any one of the multipliers.
  u(:)=nan;v(:)=nan;
%{  
  u(1)=0; % v(m)=0;
  
  while(any(u~=u) || any(v~=v)) % implicit isnan()
  % Step 2. Scan the rows and columns of the array until a circled element cij is found
  % such that either ui or vj (but not both) has already been determined.
   for ij=iB' % flipud(iB)'
     c=C(ij(1),ij(2));ui=u(ij(1));vj=v(ij(2));
     % Step 3. Compute the undetermined ui or vj from the equation ci j = ui + vj. 
     if(isnan(ui) && ~isnan(vj)) 
      u(ij(1))=c-vj;
     elseif(~isnan(ui) && isnan(vj))
      v(ij(2))=c-ui;
     end
   end
  end
%}
% Jarre / Stoer 
% Wir können ohne Einschränkung u1 = 0 wählen.
% Für die Nachbarn Dl von S1, diemit S1 durch eine Kante in G(J ) verbunden sind, folgen
% aus u1 = 0 und u1 +v ell = c1ell die Werte vell . Für deren Nachbarn Sk (in G(J )) sind dann
% wiederum die Werte uk durch uk +vell = ckelll eindeutig gegeben. So lassen sich sukzessive
% alle ui und v j bestimmen.
  known(:)=false;
  k=1;
  u(iB(k,1))=0;

  k=iB(:,1)==iB(k,1); 
  known(k,1)=true; 
  for l=iB(k,:)'; v(l(2))=C(l(1),l(2)); end  
  
  while(~all(any(known,2)))
   k=any(iB(:,2)==iB(k,2)',2)  & ~known(:,1);
   known(k,2)=true;
   for l=iB(k,:)'; u(l(1))=C(l(1),l(2))-v(l(2)); end  
   
   k=any(iB(:,1)==iB(k,1)',2)  & ~known(:,2);
   known(k,1)=true;
   for l=iB(k,:)'; v(l(2))=C(l(1),l(2))-u(l(1)); end   
   
  end
  if(verbose)
  disp('Simplex multipliers');
  u',v'
  end
  % **********************************************
  R=C-u-v'; % relative costs
  % enforce exact zeros for all indices in iB
  for i=1:m+n-1; R(iB(i,1),iB(i,2))=0; end  
% If all relative cost coefficients are nonnegative, stop; the solution is optimal. 
  if all(R>-10*eps), break; end
%% Step 3. Select a nonbasic variable corresponding to a negative cost coefficient to
% enter the basis (usually the one corresponding to the most negative cost coefficient).
[cr,ij]=min(R(:));i=1+mod(ij-1,n);j=1+(ij-i)/n;
 if(cr>=0), error('oops, negative costs expected'); end
 if(verbose)
  disp('New coordinate'); [i,j]
 end
 % Compute the cycle of change

%{
%% last working copy ...
B=iB;
ii=1:(n+m-1); % true(n+m-1,1); %
n0=1:n;n0(i)=[];
m0=1:m;m0(j)=[];
while(true)
 hh=B(:,1)==n0;
 vv=B(:,2)==m0;
% remove all dead-end elements
 hleaf=sum(hh)==1;
 vleaf=sum(vv)==1; 
 dd=any([hh(:,hleaf),vv(:,vleaf)],2); 
 if(dd==0), break; end % any([])==0
 B(dd,:)=[];
 ii(dd)=[];
 n0(hleaf)=[];
 m0(vleaf)=[];
end
%}

%% improvements: 
% work with blacklist (0) instead of deletion ([]) -> No (61s vs 54s)
% only test those indices which pairmates have been blacklisted (48 with [], 49 with 0)
iC=[iB;i,j]; % copy once, close the loop with candidate index
ii=1:(n+m-1); % true(n+m-1,1); %
n0=1:n;n0(i)=[];
m0=1:m;m0(j)=[];
% remove all dead-end elements
while(true)
 if(~isempty(n0))
  hh=iC(:,1)==n0;
  hleaf=sum(hh)==1;
  hind=any(hh(:,hleaf),2);
  m1=iC(hind,2)'; 
 else
  hind=[];
  m1=[];
 end
 if(~isempty(m0))
  vv=iC(:,2)==m0;
  vleaf=sum(vv)==1;
  vind=any(vv(:,vleaf),2);
  n1=iC(vind,1)'; 
 else
  vind=[];
  n1=[];
 end
 
 dd = any([hind, vind],2);
 if(isempty(dd) || all(dd==0)), break; end % any([],2)==0, but any( [ [],[] ],2)==[] ?
 iC(dd,:)=0;
 ii(dd)=0;
 n0=n1;
 m0=m1;
end
ii=ii(ii>0);
B=iC(ii,:); 

if(verbose)
disp('Cycle index');
[i,j;B]' % the cycle index
end
ci=i;cj=j;

iCC(:)=0;
first=true;
while first || cj~=j
% column scanning
 ff=B(:,2)==cj;  if(~first), ff(jj)=false; end
 jj=find(ff);
% just one entry, use it
 iCC(ii(jj))=-1;
  ci=B(jj,1);
 % 
 first=false;
  %row scanning
 ff=B(:,1)==ci;  ff(jj)=false;
 jj=find(ff);
 if(isempty(jj)), break; end
 iCC(ii(jj))=1;
 cj=B(jj,2);
end

%% Set theta equal to the smallest basic variable
% with a minus assigned to it. 
kk=find(iCC<0);
thta=inf;
for l=1:length(kk)
 %q=P(iB(kk(l),1),iB(kk(l),2));
 q=sP(kk(l));
 if(thta>q)
  thta=q; k=kk(l);
 end
end
if verbose
thta
end
%Update the solution.
if(thta~=0.0), sP=sP+thta*iCC; end
%for l=ii
% P(iB(l,1),iB(l,2))=P(iB(l,1),iB(l,2))+thta*iCC(l);
%end
%P(i,j)=thta;
%P(iB(k,1),iB(k,2))  % should be zero
%P
sP(k)=thta;
iB(k,:)=[i,j];
 % debugging
% Popt=zeros(n,m);
% for i=1:(n+m-1)
%  Popt(iB(i,1),iB(i,2))=sP(i);
% end
% assert(all(sum(Popt)==b'),'corrupt row sum'); 
% assert(all(sum(Popt,2)==a),'corrupt column sum');
end % while
% construct full matrix
Popt=zeros(n,m);
for i=1:(n+m-1)
Popt(iB(i,1),iB(i,2))=sP(i);
end
fopt=C(:)'*Popt(:);
end % function

function xxtest
%%
 b=[1,5,2,8,2]';
 a=[3,8,1,6]';
 C=[3 4 6 8 9; 2 2 4 5 5; 2 2 2 3 2; 3 3 2 4 2];
 transsimp(a,b,C)
%%
end