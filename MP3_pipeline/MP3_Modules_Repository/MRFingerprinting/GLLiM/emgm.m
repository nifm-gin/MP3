function [label, model, llh, R] = emgm(X, init, maxiter,verb)
% Perform EM algorithm for fitting the Gaussian mixture model.
%   X: d x n data matrix
%   init: k (1 x 1) or posteriors (n x k) or center (d x k)
% Written by Michael Chen (sth4nth@gmail.com)
%% initialization
if(verb>=1);fprintf('     EM for Gaussian mixture: running ... \n');end
R = initialization(X,init);
[~,label(1,:)] = max(R,[],2);
R = R(:,unique(label));

tol = 1e-14;
llh = -inf(1,maxiter);
converged = false;
t = 0;
while ~converged && t < maxiter
    t = t+1;
    if(verb>=1);fprintf('     Step %d\n',t);end;
    model = maximization(X,R);
    [R, llh(t+1)] = expectation(X,model);
   
    [~,label(:)] = max(R,[],2);
    u = unique(label);   % non-empty components
    if size(R,2) ~= size(u,2)
        R = R(:,u);   % remove empty components
    else
        converged = llh(t+1)-llh(t) < tol*abs(llh(t+1));
    end
end
if(verb>=1);
    if converged
        fprintf('     Converged in %d steps.\n',t);
    else
        fprintf('     Did not converge in %d steps.\n',maxiter);
    end
end

function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with a model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
    count=0;
    while (k ~= length(u) && count<20)
        count=count+1;
        k=length(u);
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
        [u,~,label] = unique(label);
    end
    k=length(u);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == n  % initialize with posteriors
    R = init;
elseif size(init,1) == d  % initialize with only centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('     ERROR: init is not valid.');
end

function [R, llh] = expectation(X, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(X,2);
k = size(mu,2);
logRho = zeros(n,k);

for i = 1:k
    logRho(:,i) = loggausspdf(X,mu(:,i),Sigma(:,:,i));
end
logRho = bsxfun(@plus,logRho,log(w));
T = logsumexp(logRho,2);
llh = sum(T)/n; % loglikelihood
logR = bsxfun(@minus,logRho,T);
R = exp(logR);


function model = maximization(X, R)
[d,n] = size(X);
k = size(R,2);
%detS=zeros(1,k);

nk = sum(R,1);
w = nk/n;
mu = bsxfun(@times, X*R, 1./nk);

Sigma = zeros(d,d,k);
sqrtR = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(:,i));
    Xo = bsxfun(@times,Xo,sqrtR(:,i)');
    Sigma(:,:,i) = Xo*Xo'/nk(i);
    % add a prior for numerical stability
    Sigma(:,:,i) = Sigma(:,:,i)+eye(d)*(1e-08);
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;