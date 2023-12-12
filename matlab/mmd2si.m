function [Ki,Si,seps,means]=mmd2si(x,y,M,kern)
% MMD2SI Maximum Mean Descrepancy Sensitivity using Kernels
%  KI=MMD2SI(X,Y,M,KERN) computes the sensitivity indices using M partitions of
%  input X.
%  M<0 enables quick mode (Gretton et al. MMD2_linear estimator).
%  M<1 enables normalization.
%  KERN defaults to Energy Distance. See code for the implemented kernels.
%  For standardizing use Y./std(Y)
%
% written by elmar.plischke@tu-clausthal.de

% QMC seems to be 1/sqrt(n) off for Si -- 'dot'/var(y)
% QMC and quick mode fails

% mmd2si(x,y/std(y),1./M,'dot') == easi/cosi/...
% mmd2si(*,'ed') == ginimim() first row
% mmd2si(x,uniformer(y),M,'ed')*6 == ginimim 2nd row / flavoursi
 [nn,d]=size(y);
 if(any(isnan(y),'all')), error('Refusing to work with NaNs.'); end
 if nargin<4 || isempty(kern)
     if d>1
      kern=@(z,y)-.5*sqrt(sum((z-y).^2,2));
     else
      kern=@(z,y)-.5*abs(z-y);
     end
 end
 if ischar(kern), kern=fetchKern(kern); end
 [n,k]=size(x);
 assert(n==nn,'In- and Output sizes mismatch.');
 usequick=false; 
 if(M<0)
     M=-M; usequick=true; 
 end
 % and more parameter misuse
 usescale=false;
 if(M<1)
     M=round(1/M); usescale=true; 
 end

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
 if(nargout>1 && d==1)
 % Correlation Ratio estimator
  Ey=mean(y);
  Si=zeros(1,k);
  for j=1:k
   Si(j)=means(j,:)*(accumarray(ri(:,j),y,[],@mean)-Ey).^2;
  end
  Si=Si/var(y)*(M/(M-1)); % biased estimator
 else
  Si=nan(k,1); % hi-dim Si ?
 end

 maxRank=max(ri,[],'all');
 if(~usequick)
 %{
 my2=EkernGramian(y,kern);

 if(d==1)
  for i=1:k
   mycy = accumarray(ri(:,i),y,[], @(z)EkernMix(y,z,kern));
   mcy2 = accumarray(ri(:,i),y,[], @(z)EkernGramian(z,kern));
   Ki(i)= means(i,:)*(my2+mcy2-2*mycy);
   if(nargout>2)
    seps(i,:)=my2+mcy2-2*mycy;
   end
  end
 else
  ramp=(1:n)';
  for i=1:k
   mycy = accumarray(ri(:,i),ramp,[], @(z)EkernMix(y,y(z,:),kern));
   mcy2 = accumarray(ri(:,i),ramp,[], @(z)EkernGramian(y(z,:),kern));
   Ki(i)= means(i,:)*(my2+mcy2-2*mycy);
   if(nargout>2)
    seps(i,:)=my2+mcy2-2*mycy;
   end
  end
 end
 %}
 %% insert the Eq (47) estimator here
  E=zeros(n,k+2); 
  dbase=cell(maxRank,k); % for saving classes/bins
  if(n>9999 || d>20)
      stat=101;%0
      h=[];%h=waitbar(stat,'Chewing on kernels'); 
  else 
      h=[]; 
      stat=101; % status 101% keeps waitbar from opening
  end

  for i=1:n
     stat1=round(i/n*100);
     if(stat1>stat), stat=stat1; waitbar(stat/100,h); end

     Ky=kern(y(i,:),y);
     E(i,k+2)=Ky(i); % save diagonal entry
     Ky(i)=0; % Ky(i) might be inf or nan, depending on kernel
     E(i,k+1)=sum(Ky)/(n-1); % unconditional, E[K(Y*,Yi)] without diagonal entry
	 
     for j=1:k
         m=ri(i,j);
         ii=dbase{m,j}; % fetch partition class
         if isempty(ii) % if not available
          ii=find(ri(:,j)==m); % compute
          dbase{m,j}=ii;
         end
		  % conditional, E[K(Yj,Yi) | Xi,Xj from the same class,i~=j ]
         % E(i,j)=(sum(Ky(ii))-Ky(i))/(length(ii)-1);
         E(i,j)=(sum(Ky(ii)))/(length(ii)-1);
     end
  end
  if(~isempty(h)), delete(h); end % waitbar
  meanE=mean(E);
  if(~usescale)
  Ki=meanE(1:k)-meanE(end-1); % paper has wrong sign?
  else
  Ki=( meanE(1:k)-meanE(end-1) )./( meanE(end)-meanE(end-1) ) ; % in [0,1]?
 % Kelley debiasing -- not working
 % Ki=(n*Ki-M)/(n-M);
  end
  seps=E;
 else
  % use quick version with subset selection
  seps=zeros(k,maxRank);
  Ki=zeros(1,k);
  if(d==1)
   for j=1:k
      mm = accumarray(ri(:,j),y,[],@(z)Ekernquick(y,z,kern));
      seps(j,:)=mm;
      Ki(j)=means(j,:)*mm;
   end
  else
   ramp=(1:n)';
   for j=1:k
    % an ugly kludge to process matrices
    mm = accumarray(ri(:,j),ramp,[], @(z)Ekernquick(y,y(z,:),kern));
    seps(j,:)=mm;
    Ki(j)=means(j,:)*mm;
   end
  end
 end
% if(nargout>=5)
%     % test debiased version
%     Ki2=Eq47estimatorlight(y,ri,kern);
% end
end 
%{
function my=EkernGramian(y,kern)
% for symmetric kernel
% 
 m=size(y,1);
 my=0;
 for i=2:m
  my=my+sum(kern(y(i,:),y(1:i-1,:)));
 end
 my=2*my/(m*(m-1));
%
%{
 m=size(y,1);
 my=0;
 for i=1:m
  my=my+sum(kern(y(i,:),y(1:i,:)));
 end
 my=2*my/(m*(m+1));
%}
%{
 m=size(y,1);
 my=0;
 % diagonal
 md=kern(y(1,:),y(1,:));
 for i=2:m
  % diagonal
  md=md+kern(y(i,:),y(i,:));
  %off-diagonal
  my=my+sum(kern(y(i,:),y(1:(i-1),:)));
 end
 my=(2*my+md)/m^2;
%}
end

function mxy=EkernMix(x,y,kern)
% for independent data
n=size(x,1);
m=size(y,1);
 mxy=0;
 if n<m
  for i=1:n
   mxy=mxy+sum(kern(x(i,:),y));
  end
 else
  for i=1:m
   mxy=mxy+sum(kern(y(i,:),x));
  end
 end
 mxy=mxy/(n*m);
end
%}
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

function k=fetchKern(str)
  switch(lower(str))
  case { 'ed','energy','l2'}
   k= @(y,z)-.5*sqrt(sum((z-y).^2,2)); % energy distance
  case {'product','variance','inner' ,'linear','dot'}
   k=@(y,z)sum(y.*z,2); % variance based
  case {'exp','hsic','gauss','rbd'}
   % k=@(y,z)exp(-.5*sqrt(sum((z-y).^2,2))); % exponential / HSIC
   k=@(y,z)exp(-.5*sum((z-y).^2,2)); 
  case 'cauchy'
   k=@(y,z)prod(1./(1+(z-y).^2),2); % Cauchy
  case 'laplace'
   k=@(y,z)exp(-sum(abs(z-y),2)); % Laplace
  case {'linf','max'}
   k=@(y,z)-max(abs(y-z),[],2); % L infinity
  case {'l1','taxicab','manhattan'}
   k=@(y,z)-sum(abs(y-z),2); % L infinity
  case 'matern32' % Rasmussen Williams: characteristic kernels
   matern32=@(r)(1+sqrt(3)*r).*exp(-sqrt(3)*r);
   k=@(y,z)matern32( sqrt(sum((z-y).^2,2)) );
  case 'matern52'
   matern52=@(r)(1+sqrt(5)*r+5/3*r.^2).*exp(-sqrt(5)*r);
   k=@(y,z)matern52( sqrt(sum((z-y).^2,2)) );
      otherwise
   error('unknown kernel');
  end
end

function res=testme
%%
 try
     if isempty(gcp('nocreate'))
         parpool('threads'); % maybe not for timing
     end
 catch
 end
%%
 n=1024; % 2^14; % 8192;
 k=4;
 R=5;
 res=[];
 for r=1:R
  fprintf('Run %d of %d\n',r,R);
  u=rand(n,k);
  % u=sobolpoints(n,k);
  x=[1,1,1,0]+norminv(u)*chol([1 .5 .5 0; .5 1 .5 0; .5 .5 1 0; 0 0 0 1]);
  y=x*[4 -2 1 0;  2 5 -1 0]'; % [4 -2 1 0]'; %
%

  Ms=2.^(2:8);
  Vy=sum(var(y),2);
  for M=Ms
    M
    [Ei,Si,~,~,Ei2]=mmd2si(x,y,M);
    [W2,W22sep]=bwsi(x,y,M);
    res(end+1,:)=[r,M,Ei2,Ei,W2,mean(W22sep)/(2*Vy)]; % Si ->Ei2
  end
 end
 %%
 clf
 for i=1:R
  subplot(2,2,1);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+k+(1:k)),'o-')
  if(i==1)
   grid on; set(gca,'XTick',[10,10^2])
   title('Energy distance');xlabel('Partition size');
   try
    legend(strcat('x_',num2str((1:k)')),'AutoUpdate','off','NumColumns',2)
   catch
    legend(strcat('x_',num2str((1:k)')));
   end
   hold on
  end
  subplot(2,2,2);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+(1:k)),'o-')
  if(i==1)
   grid on; set(gca,'XTick',[10,10^2])
   title('Sobol sensitivity');xlabel('Partition size');
   hold on
  end
  subplot(2,2,3);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+2*k+(1:k)),'o-')
  if(i==1)
   grid on; set(gca,'XTick',[10,10^2])
   title('Wasserstein2');xlabel('Partition size');
   hold on
  end
  subplot(2,2,4);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+3*k+(1:k)),'o-')
  if(i==1)
   grid on; set(gca,'XTick',[10,10^2])
   title('Wasserstein2 squared normalized');xlabel('Partition size');
   hold on
  end
 end
 %%
end
%%
%{
function L=Eq47estimatorlight(Y,partition,kern)
% 4/n(n-1) sum i\neq j K(Yi,Yj)-4/n sum h=1 n 1/d_i sum_j in Xh(i) K(Xi,Xj)
 n=size(Y,1); % sample size
 d=size(partition,2); % imput dimension
 E=zeros(n,d+1);
 dbase=cell(max(partition,[],'all'),d); % for saving classes/bins
 for i=1:n
     Ky=kern(Y(i,:),Y);
     E(i,d+1)=(sum(Ky)-Ky(i))/(n-1);
     for j=1:d
         m=partition(i,j);
         ii=dbase{m,j};
         if isempty(ii)
          ii=find(partition(:,j)==m);
          dbase{m,j}=ii;
         end
         E(i,j)=(sum(Ky(ii))-Ky(i))/(length(ii)-1);
         %E(i,j)=sum(Ky(ii))/length(ii);
     end
 end
 L=mean(E(:,1:d))-mean(E(:,end)); % wrong sign?
end
%}
function testforbias
%%
ishigami
for R=1:1
n=2000;
x=trafo(rand(n,k));
y=model(x);
ms=2:2:400;
K=zeros(length(ms),k);
L=zeros(length(ms),k);
% hack M<1 for  normalization
for i=1:length(ms), [K(i,:)]=mmd2si(x,y,1/ms(i),'ed'); 
% L(i,:)=susi(x,y,ms(i));
end
subplot(1,2,1)
set(gca,'ColorOrderIndex',1);
plot(ms,sqrt(max(K,0)),'-'); hold on
%plot(ms,L,'--');
%set(gca,'ColorOrderIndex',1);
%plot(ms,(n*L-ms')./(n-ms'),':')
% it's variance-based, so we know the analytical values
if(R==1)
    legend('x_1','x_2','x_3','x_4','AutoUpdate','off','NumColumns',2)
title(['MMD^2 sensitivity, sample size ' num2str(n)]);
ylabel('Fraction of variance');
xlabel('Partition size')
    a=axis;
 plot(a(1:2),[.3139,.3139;.4428,.4428; 0 0]','-.k')
end
%hold off

subplot(1,2,2)
set(gca,'ColorOrderIndex',1);
plot(ms,L,'-'); hold on
if(R==1)
 legend('x_1','x_2','x_3','x_4','AutoUpdate','off','NumColumns',2)
title(['Pearson CR, sample size ' num2str(n)]);
ylabel('Fraction of variance');
xlabel('Partition size')
 a=axis;
 plot(a(1:2),[.3139,.3139;.4428,.4428; 0 0]','-.k')
end

end
hold off

%%
end

%%
function chhortest
%%
% Julien Chhor on fitting with singular kernels + MMD2

ishigami
n=2000;
x=trafo(rand(n,k));y=model(x);

kernel=@(r) (r<1).*(r.^(-1/3)); % power between (0 and dimension/2)
%Scale=[100,10,1,.8,.5,.2,.1,.08,.05,.02,.01,.008,.005,.002,.001];
Scale=logspace(-1.5,0,45);
for i=1:length(Scale)
    Ki(i,:)=mmd2si(x,y*Scale(i),12,@(y,z)kernel(abs(y-z')));
end
semilogx(Scale,Ki);
title('output scale factor')
legend('x_1','x_2','x_3','dummy');
ylabel('MMD2 sensitivity')
title('Ishigami with singular kernel')
xlabel('output scale factor')
%%
end