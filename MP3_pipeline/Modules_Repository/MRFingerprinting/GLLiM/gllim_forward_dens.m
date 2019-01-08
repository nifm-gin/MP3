function [y_dens,psi] = gllim_forward_dens(x,theta,y_samples,verb)
%%%%%%%%%%% Forward Conditional Density from Gllim Parameters %%%%%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Calculate Gaussian mixture parameters of the forward
% conditional density p(y|x;theta) in space Y using a single observation x
% and gllim parameters theta. Evaluate density at points y_samples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - x (Lx1)               % Input observations to map
% - theta  (struct)       % Gllim model parameters
%   - theta.c (LxK)       % Gaussian means of X's prior
%   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
%   - theta.pi (1xK)      % Gaussian weights of X's prior
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
% - y_samples (DxN)       % Points where to evaluate p(y|x;theta) (def [])
% - verb {0,1,2}          % Verbosity (def 1)
%%%% Output %%%%
% - y_dens (1xN)          % Density p(y|x;theta) evaluated at y_samples
% - psi (struct)          % Gaussian mixture parameters of p(y|x;theta)
%   - psi.S (DxDxK)       % Gaussian covariance matrices
%   - psi.mu (DxK)        % Gaussian means
%   - psi.alpha (1xK)     % Gaussian weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=size(x,1);
[D,K]=size(theta.b);
% =========================Default input values============================
if(~exist(y_samples,'var'));y_samples=[];end;
if(~exist(verb,'var')||isempty(verb));verb=1;end;

% ======================Forward density parameters=========================
if(verb>=1);fprintf(1,'Compute foward conditional density parameters\n');end;
% Parameters to estimate:
mu=NaN(D,K);         % Conditional means
logalpha=zeros(1,K); % Conditional log-weights log(p(Z=k|x;theta))
pXyz=zeros(size(xgrid,1),K); % Probability p(y|x,Z=k;theta)
for k=1:K
    if(verb>=2), fprintf(1,'  k=%d ',k); end
    if(verb>=2), fprintf(1,'Ab '); end
    Ak=reshape(theta.A(:,:,k),D,L); % DxL
    bk=reshape(theta.b(:,k),D,1); % Dx1
    
    if(verb>=2), fprintf(1,'mu '); end 
    mu(:,k)=Ak*x+bk; % Dx1 
    
    if(verb>=2), fprintf(1,'logalpha '); end 
    logalpha(k)=log(theta.pi(k))+...
                loggausspdf(x,theta.c(:,k),theta.Gamma(:,:,k));

    if(~isempty(y_samples))
        if(verb>=2), fprintf(1,'pXyz '); end           
        pXyz(:,k)=logalpha(k)+...
                  loggausspdf(y_samples,mu(:,k),theta.Sigma(:,:,k)); %Nx1
    end

    if(verb>=2), fprintf(1,'\n'); end 
end
den=logsumexp(logalpha);
logalpha=logalpha-den; % Normalization
alpha=exp(logalpha);

if(~isempty(y_samples))
    y_dens=sum(bsxfun(@times,pXyz',alpha));
else
    y_dens=[];
end;

psi.S=theta.Sigma;
psi.mu=mu;
psi.alpha=alpha;

end
