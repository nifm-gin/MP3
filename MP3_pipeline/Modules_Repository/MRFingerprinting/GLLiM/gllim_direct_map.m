function [y_exp,alpha] = gllim_direct_map(x,theta,verb)

[L, N]  = size(x);
%[L,K]   = size(theta.c);


% Parameters to estimate:
proj        = NaN(L, N, K);     % K projections to X space
logalpha    = zeros(N, K);      % Conditional log-weights log(p(Z=k|y;theta))

for k = 1:K
    if(verb>=2), fprintf(1,'  k=%d ', k); end
    if(verb>=2), fprintf(1,'AbcG '); end
    
    Ak      = reshape(theta.A(:,:,k),D,L);      % DxL
    bk      = reshape(theta.b(:,k),D,1);        % Dx1
    Sigmak  = reshape(theta.Sigma(:,:,k),D,D);  % DxD
    ck      = reshape(theta.c(:,k),L,1);        % Lx1
    Gammak  = reshape(theta.Gamma(:,:,k),L,L);  % LxL 
    
    if(verb>=2), fprintf(1,'projections '); end 
    proj(:,:,k)=bsxfun(@plus,Aks*y,bks);        % LxN 
    
end

end

