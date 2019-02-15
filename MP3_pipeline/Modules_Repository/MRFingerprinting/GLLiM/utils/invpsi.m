function Y=invpsi(X)
% Y = INVPSI(X)
%
% Inverse digamma (psi) function.  The digamma function is the
% derivative of the log gamma function.  This calculates the value
% Y > 0 for a value X such that digamma(Y) = X.
%
% This algorithm is from Paul Fackler:
% http://www4.ncsu.edu/~pfackler/
%
  L = 1;
  Y = exp(X);
  while L > 10e-8
    Y = Y + L*sign(X-psi(Y));
    L = L / 2;
  end
