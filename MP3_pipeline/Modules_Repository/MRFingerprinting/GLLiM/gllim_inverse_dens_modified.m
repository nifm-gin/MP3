function [x_dens,psi,theta] = gllim_inverse_dens_modified(y,theta,chi,x_samples,verb)

%%%%%%%%%% Inverse Conditional Density from Gllim Parameters %%%%%%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Calculate Gaussian mixture parameters of the inverse
% conditional density p(x|y;theta) in space Y using gllim parameters
% theta, where y may be a set of T observations with missing components
% specified by chi. Evaluate the density at points x_samples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - y (DxT)               % Input observations to map
% - theta  (struct)       % Gllim model parameters
%   - theta.c (LxK)       % Gaussian means of X's prior
%   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
%   - theta.pi (1xK)      % Gaussian weights of X's prior
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
% - chi (DxT)             % Indicate non-missing observations (def ones(y))
% - x_samples (LxN)       % Points where to evaluate p(x|y;theta) (def [])
% - verb {0,1,2}          % Verbosity (def 1)
%%%% Output %%%%
% - x_dens (1xN)          % Density p(x|y;theta) evaluated at x_samples
% - psi (struct)          % Gaussian mixture parameters of p(x|y;theta)
%   - psi.S (LxLxK)       % Gaussian covariance matrices
%   - psi.mu (LxK)        % Gaussian means
%   - psi.alpha (1xK)     % Gaussian weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,T]=size(y);
[L,K]=size(theta.c);
% =========================Default input values============================
if(~exist('chi','var')||isempty(chi));chi=ones(size(y));end;
if(~exist('x_samples','var'));x_samples=[];end;
if(~exist('verb','var')||isempty(verb));verb=1;end;

% ======================Inverse density parameters=========================
if(verb>=1);fprintf(1,'Compute inverse conditional density parameters\n');end;
% pre-precomputations:
log2piL=L*log(2*pi);
sqrtChibar=sqrt(sum(chi,2)); % Dx1
% Parameters to estimate:
mu=NaN(L,K);         % Conditional means
S=NaN(L,L,K);        % Conditional covariance matrices
logalpha=zeros(1,K); % Conditional log-weights log(p(Z=k|y;theta))
pXyz=zeros(size(x_samples,2),K); % Probability p(x|y,Z=k;theta)
% Add by Me
for k=1:K
    if(verb>=2), fprintf(1,'  k=%d ',k); end
    if(verb>=2), fprintf(1,'AbcG '); end
    Ak=reshape(theta.A(:,:,k),D,L); % DxL
    bk=reshape(theta.b(:,k),D,1); % Dx1
    Sigmak=reshape(theta.Sigma(:,:,k),D,D); %DxD
    ck=reshape(theta.c(:,k),L,1); % Lx1
    Gammak=reshape(theta.Gamma(:,:,k),L,L); % LxL 
    weighted_Ak = bsxfun(@times,Ak,sqrtChibar); % DxL
    
    if(verb>=2), fprintf(1,'invS '); end  
    invSk=weighted_Ak'/Sigmak*weighted_Ak+inv(Gammak); % LxL
    S(:,:,k)=inv(invSk);
    
    if(verb>=2), fprintf(1,'mu '); end 
    diff=chi.*bsxfun(@minus,y,bk); % DxT
    mu(:,k)=invSk\(Ak'/Sigmak*sum(diff,2)+Gammak\ck);    
    muk=reshape(mu(:,k),L,1); % Lx1  
    
    if(verb>=2), fprintf(1,'logalpha '); end 
    logalpha(k)=log(det(invSk))+log(det(Gammak))+...
                ck'/Gammak*ck-muk'*invSk*muk;
    invSigmak=inv(Sigmak);
    for t=1:T
        if(verb>=2), fprintf(1,'t'); end 
        dat=(chi(:,t)==1); % Indexes of the D' non-missing data in y_t
        dchi=y(dat,t)-bk(dat); % D'x1
        iSchi=invSigmak(dat,dat); %D'xD'
        logalpha(k)=logalpha(k)+dchi'*iSchi*dchi;        
    end
    logalpha(k)=log(theta.pi(k))-logalpha(k)/2;
    
    if(~isempty(x_samples))
        if(verb>=2), fprintf(1,'pXyz '); end       
        diff=bsxfun(@minus,x_samples,muk);  %LxN
        pXyz(:,k)=(log(det(invSk))-log2piL-sum((diff'*invSk)'.*diff))/2;
    end

    if(verb>=2), fprintf(1,'\n'); end 
    
    invSigmaks2             = eye(L)+Gammak*Ak'/Sigmak*Ak;
    theta.A_s(:,:,k)        = invSigmaks2\Gammak*Ak'/Sigmak; 
    theta.b_s(:,:,k)        = invSigmaks2\(ck-Gammak*Ak'/Sigmak*bk);
    theta.c_s(:,:,k)        = Ak*ck+bk;
    theta.Sigma_s(:,:,k)    = invSigmaks2\Gammak;
    theta.Gamma_s(:,:,k)    = Sigmak+Ak*Gammak*Ak';
    
end
den=logsumexp(logalpha,2); % Tx1
logalpha=bsxfun(@minus,logalpha,den); % TxK Normalization
alpha=exp(logalpha); % TxK


if(~isempty(x_samples))
    x_dens=sum(bsxfun(@times,pXyz',alpha'));
else
    x_dens=[];
end

psi.S=S;
psi.mu=mu;
psi.alpha=alpha;
end
