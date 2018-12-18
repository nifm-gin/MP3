function [x_exp,alpha] = gllim_inverse_map(y,theta,verb)
%%%%%%%%%%%%%%%%% Inverse Mapping from Gllim Parameters %%%%%%%%%%%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Map N observations y using the inverse conditional
% expectation E[x|y;theta] of the gllim model with parameters theta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - y (DxN)               % Input observations to map
% - theta  (struct)       % Gllim model parameters
%   - theta.c (LxK)       % Gaussian means of X's prior
%   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
%   - theta.pi (1xK)      % Gaussian weights of X's prior
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
% - verb {0,1,2}          % Verbosity (def 1)
%%%% Output %%%%
% - x_exp (LxN)           % Posterior mean estimates E[xn|yn;theta]
% - alpha (NxK)           % Weights of the posterior GMMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,N]=size(y);
[L,K]=size(theta.c);
% ======================Inverse density parameters=========================
if(verb>=1);fprintf(1,'Compute K projections to X space and weights\n');end;
% Parameters to estimate:
proj=NaN(L,N,K);     % K projections to X space
logalpha=zeros(N,K); % Conditional log-weights log(p(Z=k|y;theta))
for k=1:K
    if(verb>=2), fprintf(1,'  k=%d ',k); end
    if(verb>=2), fprintf(1,'AbcG '); end
    Ak=reshape(theta.A(:,:,k),D,L); % DxL
    bk=reshape(theta.b(:,k),D,1); % Dx1
    Sigmak=reshape(theta.Sigma(:,:,k),D,D); %DxD
    ck=reshape(theta.c(:,k),L,1); % Lx1
    Gammak=reshape(theta.Gamma(:,:,k),L,L); % LxL 
    
    if(verb>=2), fprintf(1,'cks '); end
    cks=Ak*ck+bk;
    
    if(verb>=2), fprintf(1,'Gks '); end
    Gammaks=Sigmak+Ak*Gammak*Ak';

    if(verb>=2), fprintf(1,'iSks '); end
    invSigmaks2=eye(L)+Gammak*Ak'/Sigmak*Ak; %Sigmaks=invSigmaks2\Gammak

    if(verb>=2), fprintf(1,'Aks '); end
    Aks=invSigmaks2\Gammak*Ak'/Sigmak;     

    if(verb>=2), fprintf(1,'bks '); end
    bks=invSigmaks2\(ck-Gammak*Ak'/Sigmak*bk);   
    
    if(verb>=2), fprintf(1,'projections '); end 
    proj(:,:,k)=bsxfun(@plus,Aks*y,bks); % LxN 
    
    if(verb>=2), fprintf(1,'logalpha '); end 
    logalpha(:,k)=log(theta.pi(k))+...
                  loggausspdf(y,cks,Gammaks); % Nx1

    if(verb>=2), fprintf(1,'\n'); end 
end
den=logsumexp(logalpha,2); % Nx1
logalpha=bsxfun(@minus,logalpha,den); % NxK Normalization
alpha=exp(logalpha); % NxK

x_exp=reshape(sum(bsxfun(@times,reshape(alpha,[1,N,K]),proj),3),L,N); %LxN

end
