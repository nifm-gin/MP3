function Y = inv_digamma(X,niter)
  % INV_DIGAMMA    Inverse of the digamma function.
  %
  % inv_digamma(y) returns x such that digamma(x) = y.

  % never need more than 5 iterations
  if nargin == 1
    niter = 3;
  end

  M = double(X >= -2.22);
  Y = M .*(exp(X) + 0.5) + (1-M) .* -1./(X-psi(1));

  % Newton iteration to solve digamma(x)-y = 0
  for iter=1:niter)
    Y = Y - (psi(Y)-X)./psi(1,Y);
  end

end
