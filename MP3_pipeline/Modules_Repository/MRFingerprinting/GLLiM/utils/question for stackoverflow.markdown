Matlab broadcast matrix moltipliation, woodbury and mahalanobis distance

Dear all,  
I need to compute the Mahalanobis distance distances given the point `Xdiff` and the covariance matrix `Sigma` in order to compute a multivariate Gaussian distribution. That is:

    (Xdiff)' * inv(Sigma) * (Xdiff)

The problem is that I need to do it for N sample, that is `Xdiff` is a `DxN` matrix and `Sigma` is a `DxD` matrix.
In my setup `N > 10e3` and `D>1e3`.


__baseline:__ I got the right results with this very naive code:
    
    N = 10000;
    D = 1000;
    
    Xdiff = rand(D,N);
    makesym =  @(A) triu(A) + (triu(A,1))'; % make it symmetrical
    Sigma = makesym(rand(D,D));

    res = zeros(N,1)
    invS = inv(Sigma)
    for i=1:n
        res(i) = (Xdiff(:,i)'/Sigma)*Xdiff(:,i);
    end

__searching for improvement:__ In my application, the covariance matrix `Sigma` is such that `Sigma = A + B*C*B'`. Where A is a diagonal matrix of dimension `DxD`, B is a full matrix of dimension `DxL` and C is a low-rank 'covariance' matrix of dimension `LxL`, where `L << D`.  
The [Woodbury matrix identity](https://en.wikipedia.org/wiki/Woodbury_matrix_identity) provide a way to compute inverse of `Sigma`, that is:

    invSigma = inv(A + B*C*B) 
             = inv(A) - inv(A)*B*inv(inv(C) + B'*inv(A)*C)*B'*inv(A)

Thanks to it, I would like to perform the above computation without creating/storing and compute an inverse on a `DxD` matrix.  
So I came up with the following code:

    N = 10000;
    D = 1000;
    L = 3

    Xdiff = rand(D,N);
    Adiag = rand(D,1);
    Bfull = rand(D,Lw);
    Clrnk = makesym(rand(D,L)); % must be simmetrical

    invAdiag = Adiag.^1;
    E = Bfull/invAdiag*Cfull';
    F = invAdiag * invAdiag' .* E;
    G = diag(invA_diag) - F;
    q1 = zeros(1,n);
    for i=1:n
       q1(i) = Xdiff(:,i)'*G*Xdiff(:,i);
    end