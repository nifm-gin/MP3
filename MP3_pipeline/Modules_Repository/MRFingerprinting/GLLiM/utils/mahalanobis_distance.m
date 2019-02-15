function y = mahalanobis_distance(X,mu,A,type,C,B)

switch type
    case 'full' % Full Matrix
        [~,n] = size(X);

        [U,p]= chol(A); %dxd
        if p ~= 0
            fprintf(1,'SNPD! ');
            y=-Inf(1,n); % 1xn
            return;
        end

        X = bsxfun(@minus,X,mu); % dxn

        Q = U'\X; % dxn
        y = dot(Q,Q,1); % 1xn quadratic term (M distance)

    case 'diag' % Diagonal Matrix
        [~,n] = size(X);

        if (sum(A == 0)>0)
            fprintf(1,'SNPD! ');
            y=-Inf(1,n); % 1xn
            return;
        end

        Xdiff = abs(bsxfun(@minus,X,mu)).^2; % dxn Squared difference
        Xdiff = bsxfun(@rdivide,Xdiff,A); %dxn

%         y = d*log(2*pi)+sum(log(Sigma_diag)); % 1x1 normalization constant

        y = sum(Xdiff,1); % 1*n
        
    case 'diag_lowr' % Diagonal Matrix + Low Rank Matrix
        if nargin ~= 6
            error('not enought inputs') 
        end
        
        [d,n] = size(X);
        Lw = size(B,1);

        if (sum(A == 0)>0)
            fprintf(1,'SNPD! ');
            y=-Inf(1,n); % 1xn
            return;
        end

        Xdiff = bsxfun(@minus,X,mu); % dxn difference

        % WOODBURY MATRIX IDENTITY
        % inversion of diagonal matrix
        invA_diag = 1./A; % Dx1
        D = inv(B) + bsxfun(@times,C,invA_diag)'*C; % LwxLw

        % 1x1 NORMALIZATION CONSTANT
        % logdetSigma = sum(log(A)) + log(det(eye(Lw) + C'*diag(invA_diag)*(C*B)));
%         logdetSigma = sum(log(A)) + log( ...
%                                         det(eye(Lw) ...
%                                         + bsxfun(@times, C, invA_diag)'*(C*B)));
%         c = d*log(2*pi)+logdetSigma;

        bar = zeros(1,n);
        for i=1:n
            a = (Xdiff(:,i).*invA_diag)'*C;
            bar(i) = (a/D)*a';
        end

        b = sum(bsxfun(@times, invA_diag, Xdiff.^2),1); % 1xN
        y = b - bar;
    otherwise
     error('not enough arguments')
end
end
