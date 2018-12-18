function theta = gmmj2gllim(psi,L)
%%%%%%%%% Conversion of gllim parameters to joint gmm parameters %%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Convert joint GMM parameters psi to GLLiM parameters theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - psi (struct)            % Joint GMM parameters
%   - psi.mu (L+DxK)        % Gaussian means
%   - psi.Sigma (L+DxL+DxK) % Gaussian covariances
%   - psi.weight (1xK)      % Gaussian weights
% - L (int)                 % Low dimensional variable dimensionality (T)
%%%% Output %%%%
% - theta  (struct)       % GLLiM parameters
%   - theta.c (LxK)       % Gaussian means of T
%   - theta.Gamma (LxLxK) % Gaussian covariances of T
%   - theta.pi (1xK)      % Gaussian weights of T
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DL,K]=size(psi.mu);
D=DL-L;
theta.c=psi.mu(1:L,:);
theta.Gamma=psi.Sigma(1:L,1:L,:);
theta.pi=psi.weight; 
theta.A=zeros(D,L,K);    
if(L==0)
   theta.b=psi.mu;
   theta.Sigma=psi.Sigma;
else
   theta.b=zeros(D,K);    
   theta.Sigma=zeros(D,D,K);           
   for k=1:K
       muxk=reshape(psi.mu(1:L,k),L,1);
       muyk=reshape(psi.mu(L+1:L+D,k),D,1);       
       Sxxk=reshape(psi.Sigma(1:L,1:L,k),L,L);
       Sxyk=reshape(psi.Sigma(1:L,L+1:L+D,k),L,D);
       Syyk=reshape(psi.Sigma(L+1:L+D,L+1:L+D,k),D,D);
       Ak=Sxyk'/Sxxk;
       theta.A(:,:,k)=Ak;
       theta.b(:,k)=muyk-Ak*muxk;
       theta.Sigma(:,:,k)=Syyk-Ak*Sxyk;
   end
end
end