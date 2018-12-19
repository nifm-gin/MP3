function y = loggausspdf(X, mu, Sigma)
[d,n] = size(X);
X = bsxfun(@minus,X,mu); % dxn
[U,p]= chol(Sigma); %dxd
if p ~= 0
    fprintf(1,'SNPD! ');
    y=-Inf(1,n); % 1xn
    return;
end

Q = U'\X; % dxn
q = dot(Q,Q,1); % 1xn quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U))); % 1x1 normalization constant
y = -(c+q)/2; % 1xn

end