function psi = gllim2gmmj(theta)
%%%%%%%%% Conversion of joint gmm parameters to gllim parameters %%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Convert joint GMM parameters psi to GLLiM parameters theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - theta  (struct)       % GLLiM parameters
%   - theta.c (LxK)       % Gaussian means of T
%   - theta.Gamma (LxLxK) % Gaussian covariances of T
%   - theta.pi (1xK)      % Gaussian weights of T
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
%%%% Output %%%%
% - psi (struct)            % Joint GMM parameters
%   - psi.mu (L+DxK)        % Gaussian means
%   - psi.Sigma (L+DxL+DxK) % Gaussian covariances
%   - psi.weight (1xK)      % Gaussian weights
% - L (int)                 % Latent variable dimensionality (T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,K]=size(theta.b);
L=size(theta.c,1);
psi.weight=theta.pi;
psi.mu=zeros(L+D,K);
psi.Sigma=zeros(L+D,L+D,K);
for k=1:K
    Ak=reshape(theta.A(:,:,k),[D,L]);
    ck=theta.c(:,k);
    bk=theta.b(:,k);
    Gammak=reshape(theta.Gamma(:,:,k),[L,L]);
    Sigmak=reshape(theta.Sigma(:,:,k),[D,D]);
    psi.mu(:,k)=[ck;Ak*ck+bk];
    psi.Sigma(:,:,k)=[Gammak,Gammak*Ak';Ak*Gammak,Sigmak+Ak*Gammak*Ak'];
end
end