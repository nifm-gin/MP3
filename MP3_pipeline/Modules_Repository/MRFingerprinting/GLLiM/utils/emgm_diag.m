function R = emgm_diag(X, C, it)
% Perform ONE iteration of EM algorithm for fitting a Gaussian mixture
% model with diagonal covariance matrices, given initial centroids
%   X: D x N data matrix
%   C: center (D x K)
%   R: posterior probabilities p(X|Z=k) (N x K)
% Written by Antoine Deleforge and Diego Di Carlo (inria.fr)
%% initialization

% Compute binary assignments
K = size(C,2);
[D,N] = size(X);
[~,label] = max(bsxfun(@minus,C'*X,dot(C,C,1)'/2),[],1);
R = full(sparse(1:N,label,1,N,K,N));

% Diagonal Sigma_k estimation (M-step) and posterior update (E-step)
for k=1:K
    % M-step
    diff = bsxfun(@minus,X,C(:,K)).^2;
    Sigmak = max(sum(bsxfun(@times,diff,R(:,k)'),2)/sum(R(:,k)), 1e-08);
    % E-step
    logRho(:,k) = loggausspdf_diag(X,C(:,k),Sigmak);
end
logRho = bsxfun(@plus,logRho,log(1/K));
T = logsumexp(logRho,2);
logR = bsxfun(@minus,logRho,T);
R = exp(logR);