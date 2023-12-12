function [ S,Wi ] = wassermim( x,y,M,trafo,gfx)
%WASSERMIM Moment independent measures utilizing Wasserstein metrics.
%   S=WASSERMIM(X,Y) returns a matrix of Wasserstein W1, W2, and Winf
%   sensitivity measures using linear interpolation.
%   S=WASSERMIM(X,Y,M) specifies the partition size M(1) and the 
%   number of boundary points to drop M(2).
%   S=WASSERMIM(X,Y,[],TRAFO,GFX) specifies an input transformation and
%   a graphical output.

% written by elmar.plischke@tu-clausthal.de
[n,k]=size(x);
[ys,yi]=sort(y);
showSeparation=1;
usePP=false;

if(nargin<3)||isempty(M), M=[20,3]; end
if(nargin<5), gfx=''; end
if(nargin<4||isempty(trafo))
    [~,ix]=sort(x); 
    for j=1:k
        x(ix(:,j),j)=(1:n)/n;
    end
    trafo=@(u)u;
end
Wi=zeros(M(1),k); % Wasserstein W1
W2=zeros(M(1),k); % Wasserstein W2^2
Winf=zeros(M(1),k); % Wasserstein W_oo

W=zeros(M(1),k); % weights per segment

if(~isempty(gfx))
 cols=jet(M(1));
 L=sqrt(k);
 if(ceil(L)*floor(L)>=k), myround=@floor; else myround=@ceil; end
end

for m=1:M(1)
    % binary singleton extension
 I=bsxfun(@gt,x(yi,:),trafo((m-1)/M(1)*ones(1,k)))...
  & bsxfun(@le,x(yi,:),trafo(m/M(1)*ones(1,k)));
 w=sum(I,1);
 W(m,:)=w/n;
 C=bsxfun(@rdivide,cumsum(I,1),w); % conditional cdfs
 iC=interp1(linspace(0,1,n),ys,C,'spline'); % Inverse cdfs
 Diffs=bsxfun(@minus,iC((1+M(end)):(end-M(end)),:),...
     ys((1+M(end)):(end-M(end))));
 Wi(m,:)= mytrapz(C((1+M(2)):end-M(2),:),abs(Diffs)); %mean(abs(Diffs));
 W2(m,:)= sqrt(mytrapz(C((1+M(2)):end-M(2),:),Diffs.^2)); %sqrt(mean(Diffs.^2));
 Winf(m,:)= max(abs(Diffs));
%Wi(m,:)= mean(abs(iC(:,(1+M(end)):(end-M(end)))-ys(:,(1+M(end)):(end-M(end)))))
if~isempty(gfx)
for j=1:k
     subplot(myround(L)+showSeparation,ceil(L),j);
	 if(usePP)
     % PP view
      plot(linspace(0,1,n),C(:,j),'Color',cols(m,:));
	 else
     % QQ view
      plot(iC(:,j),ys,'Color',cols(m,:));
	 end
     if(m==1)
        xlabel('y unconditional');
        ylabel('y conditional');
        title(gfx)
        hold on;
     end
     
end
end
end
if~isempty(gfx)
for j=1:k,subplot(myround(L)+showSeparation,ceil(L),j);
	if(usePP)
    % PP view
     plot([0,1],[0,1],'k--');hold off;
    else
	% QQ view
     plot([ys(1),ys(end)],[ys(1),ys(end)],'k--');hold off;
	end
end
end
S(1,:)=mean(W2); %sum(W2.*W);
S(2,:)=mean(Wi); % sum(Wi.*W);
S(3,:)=mean(Winf); % sum(Winf.*W);

if~isempty(gfx)&&(showSeparation)
    LL=myround(L)+showSeparation;LR=3*LL;
    t=linspace(0,1,M(1));
    subplot(LL,3,LR-1);plot(t,Wi);title('W_1 separation');
    subplot(LL,3,LR-2);plot(t,W2);title('W_2 separation');
    subplot(LL,3,LR);plot(t,Winf);title('W_\infty separation');
end
end
%%
function I=mytrapz(x,y)
%MYTRAPZ trapezoidal rule with matrix input
% octave does not choke on using trapz for this.
I=.5*sum(diff(x).*(y(1:end-1,:)+y(2:end,:)));
end
%%
function testwasser
%%
 ishigami
 n=2^16; %8192;
 x=trafo(sobolpoints(n,k));
 y=model(x);
 %%
 Ws=wassermim(x,y,[64 1]);
 %%
 [Di,seps]=mydoubleloop(k,[1024 32],model,trafo,[],struct('GaussQuadrature',false));
 
 t=linspace(0,1,32);
  %%
 subplot(1,3,1)
 hold on;set(gca,'ColorOrderIndex',1);
 plot(t,[seps{1}(6,:);seps{2}(6,:);seps{3}(6,:);seps{4}(6,:)],'--');
 subplot(1,3,2)
 hold on;set(gca,'ColorOrderIndex',1);
 plot(t,[seps{1}(7,:);seps{2}(7,:);seps{3}(7,:);seps{4}(7,:)],'--');
  subplot(1,3,3)
 hold on;set(gca,'ColorOrderIndex',1);
 plot(t,[seps{1}(8,:);seps{2}(8,:);seps{3}(8,:);seps{4}(8,:)],'--');
%% 
end

