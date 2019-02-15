%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 		TEST SUITE
%

% run the following:
%	result = runtests('test_loggaus');

% data
N = 10000;
D = 2000;
Lw = 10;
alpha = 2;

X  = alpha*rand(D,N);			    % DxN
mu = alpha*rand(D,1);			    % Dx1
Sigma_diag = alpha*rand(D,1);     % Dx1
Sigma = Sigma_diag.*eye(D); % DxD

% exprexted values
y = loggausspdf(X, mu, Sigma);
% actual value
y_dia = loggausspdf_diag(X, mu, Sigma_diag);

assertMatrixWithAbsTol(y,     y, 'The results are not the same')
assertMatrixWithAbsTol(y_dia, y, 'The results are not the same')


makesym =  @(A) triu(A) + (triu(A,1))';
A_diagonal = alpha*rand(D,1);	% Dx1
C_full = alpha*rand(D,Lw);   	% DxLw
B_lowRank = makesym(cov(alpha*rand(D,Lw))); % must be simmetrical
CBCt = makesym(C_full*B_lowRank*C_full');   % DxD
Sigma = makesym(diag(A_diagonal) + CBCt); % DxD

assertIsCovariance(B_lowRank, 'Non positive semi-definite matrix.')
%assertIsCovariance(CBCt, 	  'Non positive semi-definite matrix.')
assertIsCovariance((Sigma),   'Non positive semi-definite matrix.')
assertSameDimension(Sigma, eye(D), 'Dimension are not the same')

[y, c, q] = loggausspdf(X, mu, Sigma);
[y_low, c_low, q_low] = loggausspdf_diag_lowk(X, mu, A_diagonal, C_full, B_lowRank);

assertMatrixWithAbsTol(c_low, c, 'The results are not the same')
assertMatrixWithAbsTol(q_low, q, 'The results are not the same')
assertMatrixWithAbsTol(y_low, y, 'The results are not the same')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%  	S O m E  H e L P E R S
%

function assertMatrixWithAbsTol(actVal,expVal,varargin)
% Helper function to assert equality within an absolute tolerance.
% Takes two values and an optional message and compares
% them within an absolute tolerance of 1e-6.
	tol = 1e-6;
	tf = any(abs(actVal-expVal) <= tol);
	assert(tf, varargin{:});
end

function assertSameDimension(actVal,expVal,varargin)
% Takes two values and an optional message and compares
	tf = any(abs(size(actVal)-size(expVal)) == 0);
	assert(tf, varargin{:});
end

function assertIsCovariance(actVal,varargin)
	symmetry = issymmetric(actVal);
	[~,D]=eig(actVal);
	eigenvalues = diag(D);
	eigenvalues(abs(eigenvalues)<eps) = 0;
	tf = all(eigenvalues>0) & symmetry;
	assert(tf, varargin{:});
end