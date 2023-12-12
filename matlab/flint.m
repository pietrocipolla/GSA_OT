function [Wi,Si,Sy]=flint(x,y,M)
% FLINT Approximate Wasserstein sensitivity using DCT.
 gfx='bla';
[n,k]=size(x);
[~,ix]=sort(x);
Ey=mean(y);
Sy=std(y);
%ys=y(ix)-Ey;
yy=(y-Ey)/Sy;  % standardized
ys=yy(ix);

  xx=linspace(0,1,n); % gfx

%  ym = diff(ys).^2; % Chatterjee estimator
%  yc = diff(ys.^2).^2;
%  Wi=mean( sqrt(ym+(sqrt(yc+1))))*Sy;
%  Si(1,:)= 1-mean(ym)/2; 
%% Now, orthonormal regression/projection
 if(nargin<3) || isempty(M)
  M=8;
 end

  allcoef=dct([ys,ys.^2]); 
  selcoef=zeros(n,2*k);selcoef(1:(M+1),:)=allcoef(1:(M+1),:);
  ym=idct(selcoef);

  if(~isempty(gfx))
    subplot(1,3,1);
   % plot(1:n,Ey+ys(:,i),'.',1:n,Ey+ym(:,i),'LineWidth',3);
   plot(xx,Ey+Sy*ym(:,1:k),[0,1],Ey*[1,1],'k--','LineWidth',3);
   title('Deviation from mean')
    subplot(1,3,2);
    plot(xx,Sy*sqrt(ym(:,k+(1:k))),[0,1],Sy*[1,1],'k--','LineWidth',3);
   title('Deviation from stddev')
   subplot(1,3,3);
    plot(Ey+Sy*ym(:,1:k),Sy*sqrt(max(ym(:,k+(1:k)),0)),'-',Ey,Sy,'k*','LineWidth',2,'MarkerSize',5);
    hold on;
    phi=linspace(-pi,pi)';circ=[cos(phi),sin(phi)];
    axis('equal');a=axis;
    for i=1:8
        r2=Sy*(i/8)^2;
        plot(Ey+r2*circ(:,1),Sy+r2*circ(:,2),'k:');
        xlabel('mean');ylabel('stddev');
    end
    axis(a);
    hold off
    
  end
  Si=mean(ym(:,1:k).^2);
  Wi=Sy*sqrt(mean(ym(:,1:k).^2+(1-sqrt(max(ym(:,k+(1:k)),0))).^2));
 end 

function testflint
%% Feuerstein ist kein wasserstein
  
%%  study a Gaussian linear model
  model=@(x)x*[4;-2;1];
  trafo=@(u)norminv(u);
%%
  R=20;
  opts1=struct('GaussQuadrature',true,'RandomSource',@sobolpoints);
  opts2=struct('GaussQuadrature',false,'RandomSource',@rand,'ShowEval',false');
  L=7;
  for i=1:L
      % Gauss quadrature / QMC
   l=mydoubleloop(3,[4*2^i,16],model,trafo,[],opts1);
   Ws1(i,:)=l(7,:);
   for r=1:R
   m=mydoubleloop(3,[4*2^i,16],model,trafo,[],opts2);
   Ws2(r,i,:)=m(7,:);
   end
  end
  subplot(4,1,1)
  boxplot([Ws2(:,:,1),Ws2(:,:,2),Ws2(:,:,3)]); hold on
  plot(1:3*L,Ws1(:),'o'); 
  a=axis;plot([7.5,7.5],a(3:4),':k',[14.5,14.5],a(3:4),':k');hold off

  title('Double-loop design. Rings: Gaussian quadrature + QMC, boxes: Riemann integration + MC, 20 replicates');
  set(gca,'XTickLabel',{'392','784','1568','3136','6272','12544','25088',...
      '392','784','1568','3136','6272','12544','25088',...
      '392','784','1568','3136','6272','12544','25088'},'XTickLabelRotation',-45)
  %%
  k=3;
  for i=1:L
      n=20*2^i;
      for r=1:R
      x=trafo(rand(n,k));
     % x=trafo(sobolpoints(n,k)); % need scramble
      y=model(x);
      l=wassersi(x,y);
         Ws3(r,i,:)=l.W2;
         m=wassermim(x,y);
         Ws4(r,i,:)=m(1,:);
         Ws5(r,i,:)=flint(x,y);
   end
  end
  subplot(4,1,2);
  boxplot([Ws3(:,:,1),Ws3(:,:,2),Ws3(:,:,3)]);
  hold on;a=axis;plot([7.5,7.5],a(3:4),':k',[14.5,14.5],a(3:4),':k');hold off
  title('Given data binning, resampled conditional. MC, 20 replicates');
    set(gca,'XTickLabel',{'40','80','160','320','640','1280','2560',...
      '40','80','160','320','640','1280','2560',...
      '40','80','160','320','640','1280','2560'},'XTickLabelRotation',-45)

    subplot(4,1,3);
  boxplot([Ws4(:,:,1),Ws4(:,:,2),Ws4(:,:,3)]);
  hold on;a=axis;plot([7.5,7.5],a(3:4),':k',[14.5,14.5],a(3:4),':k');hold off
  title('Given data binning, inverse quantile. MC, 20 replicates');
    set(gca,'XTickLabel',{'40','80','160','320','640','1280','2560',...
      '40','80','160','320','640','1280','2560',...
      '40','80','160','320','640','1280','2560'},'XTickLabelRotation',-45)
    subplot(4,1,4);
  boxplot([Ws5(:,:,1),Ws5(:,:,2),Ws5(:,:,3)]);
  hold on;a=axis;plot([7.5,7.5],a(3:4),':k',[14.5,14.5],a(3:4),':k');hold off
    title('Given data, lower bound from regression curves. MC, 20 replicates');
    set(gca,'XTickLabel',{'40','80','160','320','640','1280','2560',...
      '40','80','160','320','640','1280','2560',...
      '40','80','160','320','640','1280','2560'},'XTickLabelRotation',-45)
%%   
  ishigami
  n=8192;
  x=trafo(sobolpoints(n,k));
  y=model(x);
  
  [Wi,Si]=flint(x,y);
  W=wassersi(x,y,16)
%%
end