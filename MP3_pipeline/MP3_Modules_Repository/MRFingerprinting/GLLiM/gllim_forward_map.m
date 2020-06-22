function [y_exp,alpha] = gllim_forward_map(x,theta,verb)
%%%%%%%%%%%%%%%%% Forward Mapping from Gllim Parameters %%%%%%%%%%%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Map N observations x using the forward conditional
% expectation E[y|x;theta] of the gllim model with parameters theta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - x (LxN)               % Input observations to map
% - theta  (struct)       % Gllim model parameters
%   - theta.c (LxK)       % Gaussian means of X's prior
%   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
%   - theta.pi (1xK)      % Gaussian weights of X's prior
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
% - verb {0,1,2}          % Verbosity (def 1)
%%%% Output %%%%
% - y_exp (DxN)           % Posterior mean estimates E[yn|xn;theta]
% - alpha (NxK)           % Soft assigments of points to transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L,N]=size(x);
[D,K]=size(theta.b);
% ======================Forward density parameters=========================
if(verb>=1);fprintf(1,'Compute K projections to Y space and weights\n');end;
% Parameters to estimate:
proj=NaN(D,N,K); % K projections to Y space
logalpha=zeros(N,K); % Conditional log-weights log(p(Z=k|x;theta))
for k=1:K
    if(verb>=2), fprintf(1,'  k=%d ',k); end
    if(verb>=2), fprintf(1,'Ab '); end
    Ak=reshape(theta.A(:,:,k),D,L); % DxL
    bk=reshape(theta.b(:,k),D,1); % Dx1
    
    if(verb>=2), fprintf(1,'projections '); end 
    proj(:,:,k)=bsxfun(@plus,Ak*x,bk); % DxN 
    
    if(verb>=2), fprintf(1,'logalpha '); end 
    logalpha(:,k)=log(theta.pi(k))+...
                    loggausspdf(x,theta.c(:,k),theta.Gamma(:,:,k)); % Nx1

    if(verb>=2), fprintf(1,'\n'); end 
end
den=logsumexp(logalpha,2); % 1x1xK
logalpha=bsxfun(@minus,logalpha,den); % NxK Normalization
alpha=exp(logalpha); % NxK

y_exp=reshape(sum(bsxfun(@times,reshape(alpha,[1,N,K]),proj),3),D,N); %DxN

end
