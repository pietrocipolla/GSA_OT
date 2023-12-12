% Oh, Tannenbaum et al. on Kernel Wasserstein Bures
n=2000;
y=-.5+4*rand(n,2); % not elliptically contoured, but nevertheless ...
z=randn(n,2)*[1,-.6;0,.8];

plot(y(:,1),y(:,2),'o',z(:,1),z(:,2),'+');

advective0=mean((mean(y)-mean(z)).^2)
Cy=cov(y);Ry=sqrtm(Cy);
Cz=cov(z);Rz=sqrtm(Cz);
diffusivea0=trace(Cy+Cz)
diffusiveb0=2*trace(Ry*Rz) % shortcut, trace(RyCzRy)^1/2 = trace(RyRz)  

Wasserstein0=advective0+diffusivea0-diffusiveb0

Kyy=y*y';
Kzz=z*z';
Kyz=y*z';

advective1=.5*(mean(Kyy,'all')+mean(Kzz,'all')-2*mean(Kyz,'all'))

J=eye(n,n)-ones(n,n)/n; % J^2 =J J'= J !

diffusivea1=trace(Kyy*J)/(n-1)+trace(Kzz*J)/(n-1)
diffusivea1_alt=trace(Kyy)/n-mean(Kyy,'all')+trace(Kzz)/n-mean(Kzz,'all')
%diffusive1a_alt=sum(Kyy.*J,'all')/n+sum(Kzz.*J,'all')/n
tic
diffusiveb1=2*real(trace(sqrtm(Kyz'*J*Kyz)))/n % need pos def. for the trick above
toc,tic
diffusiveb1_alt=2*mean(svd(J*Kyz)) % 2*mean(svd(Kyz-mean(Kyz)))
toc
Wasserstein1=advective1+diffusivea1-diffusiveb1

%%
dCov2=trace(Kzz*J*Kyy)/n
dCov2_alt=sum(Kzz.*(J*Kyy),'all')/n
dCov= sum(Kyy.*Kzz,'all')/n-mean(Kyy,1)*mean(Kzz,2)
%%

Wasserstein2 = mean(diag(Kyy))+mean(diag(Kzz)) ...
  -1/2*(mean(Kyy,'all')+mean(Kzz,'all'))-mean(Kyz,'all')...
  -2*mean(svd(Kyz-mean(Kyz)))
