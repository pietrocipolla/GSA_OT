function [Si,V,yhats,coeffs]=cosi(x,y,M,varargin)
%% COSI Calculation of sensitivity indices from given data.
%     SI = COSI(X,Y) returns the sensitivity indices for input arguments X 
%     and output arguments Y (per line) based upon a discrete cosine 
%     transformation.
%     SI = COSI(X,Y,M) specifies the max. harmonic cutoff.
%     SI = COSI(X,Y,M,'Gfx Title') provides a figure.
%     [Vi,V]=COSI(X,Y,M) with M<0 returns absolute contributions

%      Use -coefs(2,:)/std(y)/sqrt(n) for direction of change
%
%     References: 
%      E. Plischke, "An Effective Algorithm for Computing 
%       Global Sensitivity Indices (EASI)",
%       Reliability Engineering & Systems Safety, 95(4), 354-360, 2010
%      E. Plischke, "How to compute variance-based sensitivity 
%       indicators with your spreadsheet software", 
%       Environmental Modelling & Software, In Press, 2012

%%
%     Written by Elmar Plischke, elmar.plischke@tu-clausthal.de

 [n,k]=size(x);
 [nn,kk]=size(y);
 if nn~=n, error('Input/output sizes mismatch!'), end
% default [graphics] options
opts=struct('GfxTitle','',...
            'GfxCols',min(k,2),...
            'Labels','',...
            'OutputLabel','y',...
            'ShowScatter',true,...
			'Unscaled',false,...
			'PlotOpts',{ {'-','LineWidth',2}},...
            'ScatterData',[]); 
%
if(nargin>=4) && ~isempty(varargin)
	if isstruct(varargin{1})
        opts_in=varargin{1};
    elseif length(varargin)==1
        opts_in=struct('GfxTitle',varargin{1});
      else
        opts_in=struct(varargin{:});
    end
    members=fieldnames(opts);
    for i=1:length(members)
            o=members{i};
            if isfield(opts_in,o)
                opts.(o)=opts_in.(o);
            end
    end
end
%%

 GfxCols=opts.GfxCols;
 
 [xr,index]=sort(x);

 if kk==1
% sort output
    yr=y(index);
 else
    yr=zeros(n,k*kk);
    for i=1:kk
        z=y(:,i);
        yr(:,(i-1)*k+(1:k))=z(index);
    end
 end
 
 %% frequency selection
if (nargin==2) || (isempty(M))
  M=max(ceil(sqrt(n)),3);
 fprintf('COSI: Using %d coefficients.\n',M);
end 

if(M<0), M=-M; opts.Unscaled=true; end
% consider M terms
d=zeros(1,n);
d(1+(1:M))=1;

%% Compute transformation
allcoeff=dct(yr);

% transformation is orthogonal, so by Parseval's Theorem
V = sum(allcoeff(2:end,:).^2);
if(M==0)
% estimate approximation error
% or mean
Si=1-median(n*cumsum(allcoeff(end:-1:2,:).^2)./ ((1:(n-1))'*V));
else
Vi= sum(allcoeff(1+(1:M),:).^2);
if(~opts.Unscaled)
Si= Vi./V;
else
Si=Vi/(n-1);
V=V/(n-1);
end
end
%% Return prediction: Could be merged with Gfx (assumes kk=1)
if(nargout>=3)
 yhats=zeros(n,k);
 coeffs=allcoeff(1:(M+1),:);
 for i=1:k
  yhat=zeros(n,1);
  yhat(1:(M+1))=allcoeff(1:(M+1),i);
  yhats(:,i)=idct(yhat);
 end    
end
%% Gfx (assumes kk=1)
if ~isempty(opts.GfxTitle)
 for i=1:k
  if(k>1), subplot(ceil(k/GfxCols),GfxCols,i); end
  yhat=zeros(n,1);
  yhat(1:(M+1))=allcoeff(1:(M+1),i);
  if(opts.ShowScatter)
    plot(x(index(:,i),i),y(index(:,i)),'.',x(index(:,i),i),idct(yhat),opts.PlotOpts{:});
  else
    plot(x(index(:,i),i),idct(yhat),opts.PlotOpts{:});
  end
  if(~isempty(opts.Labels))
    xlabel(opts.Labels{i});
  else
    xlabel(['x_{',num2str(i),'}']);
  end
  ylabel(opts.OutputLabel)
  title(opts.GfxTitle);
 end
end
%%%
 if kk>1, Si=reshape(Si',k,kk)'; V=reshape(V',k,kk);end
 return
end