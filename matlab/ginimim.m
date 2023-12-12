function [ g,s,t,p,q,Gstar,S,K,P,Q,W ] = ginimim( x,y,M,gfx,trafo)
%GINIMIM Moment independent measures utilizing the CDF.
%   G=GINIMIM(X,Y) returns the Gini/Energy Distance Measure
%   [G,S,K,P,Q]=GINIMIM(...) S Smirnov, K Kuiper, P CR ECV, Q CR CDF

% written by elmar.plischke@tu-clausthal.de
[n,k]=size(x);
[ys,yi]=sort(y);
if(nargin<3), M=20; end
if(nargin<4), gfx=''; end
if(nargin<5)
    [~,ix]=sort(x); 
    for j=1:k
        x(ix(:,j),j)=(1:n)/n;
    end
    trafo=@(u)u;
end
G=zeros(M,k); % Gini
Gstar=zeros(M,k); % Gini
K=zeros(M,k); % Kuiper
S=zeros(M,k); % Smirnov
P=zeros(M,k); % Pearson Correlation Ratio: Expectation Cond. Var.
Q=zeros(M,k); % Pearson from cdf
W=zeros(M,k); % weights per segment
if(~isempty(gfx))
    cols=hsv(M);
end   
for m=1:M
    % binary singleton extension
 I=bsxfun(@gt,x(yi,:),trafo((m-1)/M*ones(1,k)))...
  & bsxfun(@le,x(yi,:),trafo(m/M*ones(1,k)));
 w=sum(I,1);
 W(m,:)=w/n;
 C=bsxfun(@rdivide,cumsum(I,1),w); % conditional cdfs
 D=bsxfun(@minus,C,linspace(1/n,1,n)');     % cdf differences
if(~isempty(gfx))
    for j=1:k
        subplot(1,k,j)
        plot(ys,-D(:,j),'.','Color',cols(m,:),'MarkerSize',5);hold on;
    end
end
 G(m,:)=D(1:end-1,:)'.^2*diff(ys);        % Gini measure
 Gstar(m,:)=mean(D(1:end-1,:)'.^2,2);       % transformation invariant Gini
 S(m,:)=max(abs(D));
 K(m,:)=max(D)-min(D);
 for j=1:k
  P(m,j)=var(ys(I(:,j)));
 end
 %Q(m,:)=(diff([ zeros(1,k); D])'*ys).^2;
 Q(m,:)=(D(1:end-1,:)'*diff(ys)).^2;
end
g=[sum(G.*W);6*sum(Gstar.*W)]; % normalization by Gamboa Klein Lagnoux
t=sum(K.*W);
s=sum(S.*W);
Vy=var(y);
q=sum(Q.*W)/Vy;
p=1-sum(P.*W)./Vy;
if(~isempty(gfx))
    for j=1:k
        subplot(1,k,j)
    title([gfx ', given x_{' num2str(j) '}']);
    xlabel('Output');
    ylabel('\Delta CDF');
    hold off
end
end
end
%%
function testgini
%%
 product21
 n=2^18;
 
 x=trafo(sobolpoints(n,k));
 y=model(x);
 ns=round(logspace(log10(512),log10(n),37));
 gs=[];bs=[];es=[];fs=[];ts=[];
 for N=ns
    N
 [gamma,beta,kappa,eta,iota]=ginimim(x(1:N,:),y(1:N),ceil(N^(1/3)),'',trafo);
 gs=[gs;gamma];
 bs=[bs;beta];ts=[ts;kappa];es=[es;eta];fs=[fs;iota];
 end
 
%%
subplot(1,5,1);plot(ns,gs);   
subplot(1,5,2);plot(ns,bs);  
subplot(1,5,3);plot(ns,ts); 
subplot(1,5,4);plot(ns,es);
subplot(1,5,5);plot(ns,fs);
%% 
end

