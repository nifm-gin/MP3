function [y, c, q] = loggausspdf_diag_lowk(X, mu, Adiag, C, B)
% LOG GAUS PDF when the Sigma matrix is in the following form
%	Sigma = A_diagonal + C_full*B_lowRank*C_full' with dimension DxD
%	the proposed algorithm uses the Woodbury matrix identity and
%	some trick-track with of the diagonality of A:
%
%	invSigma = inv(A + CBC') = invA - invA C inv(invB + C' invA C) C' invA
%            				 = invA - invA C invD C' invA
%
% X = DxN samples
% mu = DxN centroids/means
% 
% A = Dx1	diagonal of a diagonal matrix
% B = LwxLw low rank matric
% C = DxLw  full rank matrix

makesym =  @(A) triu(A) + (triu(A,1))';
B = makesym(B);

[d,n] = size(X);
Lw = size(B,1);

if (sum(Adiag == 0)>0) || all(Adiag < 10*eps)
    fprintf(1,'SNPD! ');
    y=-Inf(1,n); % 1xn
    return;
end

Xdiff = bsxfun(@minus,X,mu); % dxn difference

% WOODBURY MATRIX IDENTITY
% inversion of diagonal matrix
invA_diag = Adiag.^-1; % Dx1
D = inv(B) + bsxfun(@times,C,invA_diag)'*C; % LwxLw

% 1x1 NORMALIZATION CONSTANT
% logdetSigma = 2*sum(log(diag(chol(diag(Adiag) + C*B*C'))))
% matrix determinant lemma:
%	det(A + UWV') = det(invW + V'inAU)det(W)det(A)
%
logdetSigma = sum(log(Adiag)) + log(det(B)*det(D));

c = d*log(2*pi)+logdetSigma;

bar = zeros(1,n);
warning('off','all')

for i=1:n
    a = (Xdiff(:,i).*invA_diag)'*C;
    bar(i) = (a/D)*a';
end

b = sum(bsxfun(@times, invA_diag, Xdiff.^2),1); % 1xN
q = b - bar;

y = -(c+q)/2; % 1xn
end
