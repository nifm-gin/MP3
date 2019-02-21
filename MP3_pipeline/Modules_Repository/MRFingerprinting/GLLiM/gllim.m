function [theta,r,LLf] = gllim(t,y,in_K,varargin)
%%%%%%%% General EM Algorithm for Gaussian Locally Linear Mapping %%%%%%%%%
%%% Author: Antoine Deleforge (April 2013) - antoine.deleforge@inria.fr %%%
% Description: Compute maximum likelihood parameters theta and posterior
% probabilities r=p(z_n=k|x_n,y_n;theta) of a gllim model with constraints
% cstr using N associated observations t and y.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
%- t (LtxN)               % Training latent variables
%- y (DxN)                % Training observed variables
%- in_K (int)             % Initial number of components
% <Optional>
%- Lw (int)               % Dimensionality of hidden components (default 0)
%- maxiter (int)          % Maximum number of iterations (default 100)
%- in_theta (struct)      % Initial parameters (default [])
%                         | same structure as output theta
%- in_r (NxK)             % Initial assignments (default [])
%- cstr (struct)          % Constraints on parameters theta (default [],'')
%   - cstr.ct             % fixed value (LtxK) or ''=uncons.
%   - cstr.cw             % fixed value (LwxK) or ''=fixed to 0
%   - cstr.Gammat         % fixed value (LtxLtxK) or ''=uncons.
%                         | or {'','d','i'}{'','*','v'} (1)
%   - cstr.Gammaw         % fixed value (LwxLwxK) or ''=fixed to I
%   - cstr.pi             % fixed value (1xK) or ''=uncons. or '*'=equal
%   - cstr.A             % fixed value (DxL) or ''=uncons.
%   - cstr.b             % fixed value (DxK) or ''=uncons.
%   - cstr.Sigma         % fixed value (DxDxK) or ''=uncons.
%                         | or {'','d','i'}{'','*'} (1)
%- verb {0,1,2}           % Verbosity (default 1)
%%%% Output %%%%
%- theta  (struct)        % Estimated parameters (L=Lt+Lw)
%   - theta.c (LxK)       % Gaussian means of X
%   - theta.Gamma (LxLxK) % Gaussian covariances of X
%   - theta.pi (1xK)      % Gaussian weights of X
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
%- r (NxK)                % Posterior probabilities p(z_n=k|x_n,y_n;theta) 
%%% (1) 'd'=diag., 'i'=iso., '*'=equal for all k, 'v'=equal det. for all k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================Input Parameters Retrieval=========================
[Lw, maxiter, in_theta, in_r, cstr, verb] = ...
    process_options(varargin,'Lw',0,'maxiter',100,'in_theta',[],...
                             'in_r',[],'cstr',struct(),'verb',1);
Lt=size(t,1);
L=Lt+Lw;
[D,N]=size(y);                         
% ==========================Default Constraints============================
if(~isfield(cstr,'ct'));cstr.ct=[];end;
if(~isfield(cstr,'cw'));cstr.cw=[];end;
if(~isfield(cstr,'Gammat'));cstr.Gammat=[];end;
if(~isfield(cstr,'Gammaw'));cstr.Gammaw=[];end;
if(~isfield(cstr,'pi'));cstr.pi=[];end;
if(~isfield(cstr,'A'));cstr.A=[];end;
if(~isfield(cstr,'b'));cstr.b=[];end;
if(~isfield(cstr,'Sigma'));cstr.Sigma=[];end;

% ==========================EM initialization==============================
if(verb>=1), fprintf(1,'EM Initializations\n'); end
if(~isempty(in_theta))
    theta=in_theta;
    K=length(theta.pi);
    if(isempty(cstr.cw));
        cstr.cw=zeros(L,K); % Default value for cw
    end; 
    if(isempty(cstr.Gammaw));
        cstr.Gammaw=repmat(eye(Lw),[1,1,K]); % Default value for Gammaw
    end;
    [r,~,ec] = ExpectationZ(t,y,theta,verb); 
    [theta,cstr]=remove_empty_clusters(theta,cstr,ec); 
    [muw,Sw] = ExpectationW(t,y,theta,verb);    
    if(verb>=1);fprintf(1,'\n');end;    
else
    if(isempty(in_r))
        % Initialise posteriors with K-means + GMM on joint observed data
%         [~,C] = kmeans([t;y]',in_K);
%         [~, ~, ~, r] = emgm([t;y], C', 3, verb);
%         [~, ~, ~, r] = emgm([t;y], in_K, 3, verb); 
        [~, ~, ~, r] = emgm([t;y], in_K, 100, verb); 
%         [~,cluster_idx]=max(r,[],2);
%         fig=figure;clf(fig);    
%         scatter(t(1,:),t(2,:),200,cluster_idx','filled');        
        
%         weight=model.weight;
%         K=length(weight);
%         mutt=model.mu(1:Lt,:);
%         Sigmatt=model.Sigma(1:Lt,1:Lt,:);
%         normr=zeros(N,1);
%         r=zeros(N,K);
%         for k=1:K
%             r(:,k)=weight(k)*mvnpdf(t',mutt(:,k)',reshape(Sigmatt(:,:,k),Lt,Lt));
%             normr=normr+reshape(r(:,k),N,1);       
%         end
%         r=bsxfun(@rdivide,r,normr);        
%         fig=figure;clf(fig);
%         [~,classes]=max(r,[],2); % Nx1
%         scatter(y(1,:),y(2,:),200,classes','filled');
    end
    if(Lw==0)
        Sw=[];
        muw=[];
    else
        % Start by running an M-step without hidden variables (partial
        % theta), deduce Awk by local weighted PCA on residuals (complete
        % theta), deduce r, muw and Sw from E-steps on complete theta.
        theta = Maximization(t,y,r,[],[],cstr,verb);
        K=length(theta.pi);
        if(isempty(cstr.cw));
            cstr.cw=zeros(Lw,K); % Default value for cw
        end;
        theta.c=cat(1,theta.c,cstr.cw(:,1:K)); %LxK
        Gammaf=zeros(L,L,K);
        Gammaf(1:Lt,1:Lt,:)=theta.Gamma; %LtxLtxK
        if(isempty(cstr.Gammaw));
            cstr.Gammaw=repmat(eye(Lw),[1,1,K]); % Default value for Gammaw
        end;
        Gammaf(Lt+1:L,Lt+1:L,:)=cstr.Gammaw(:,:,1:K); %LwxLwxK    
        theta.Gamma=Gammaf;    
        % Initialize Awk with local weighted PCAs on residuals:
        Aw=zeros(D,Lw,K);
        for k=1:K
            rk_bar=sum(r(:,k));
            bk=reshape(theta.b(:,k),D,1);
            w=bsxfun(@minus,y,bk); %DxN
            if(Lt>0)
                Ak=reshape(theta.A(:,:,k),D,Lt);
                w=w-Ak*t;
            end
            w=bsxfun(@times,w,sqrt(r(:,k)'./rk_bar)); %DxN
            C=w*w'; % Residual weighted covariance matrix
            [U,Lambda]=eigs(C,Lw); % Weighted residual PCA U:DxLw
            % The residual variance is the discarded eigenvalues' mean           
            sigma2k=(trace(C)-trace(Lambda))./(D-Lw);
            theta.Sigma(:,:,k)=sigma2k*eye(D);
            Aw(:,:,k)=U*sqrt(Lambda-sigma2k*eye(Lw));
        end
        theta.A=cat(2,theta.A,Aw); %DxLxK
        [r,~,ec] = ExpectationZ(t,y,theta,verb);
        [theta,cstr]=remove_empty_clusters(theta,cstr,ec);        
        [muw,Sw] = ExpectationW(t,y,theta,verb);
        if(verb>=1);fprintf(1,'\n');end;        
    end
end

%===============================EM Iterations==============================
if(verb>=1), fprintf(1,'Running EM\n'); end
LL = -inf(maxiter,1);
iter = 0;
converged=false;
while (~converged && iter<maxiter)
    iter = iter + 1;
    
    if(verb>=1), fprintf(1,'      Iteration %d\n',iter); end;
%     fprintf(1,'mean r1=%g\n',mean(r(:,1)));      
    % =====================MAXIMIZATION STEP===========================
    theta = Maximization(t,y,r,muw,Sw,cstr,verb);
    
%     fprintf(1,'\nmean Sigma1=%g\n',mean(diag(theta.Sigma(:,:,1))));
%     fprintf(1,'mean Gamma1=%g\n',mean(diag(theta.Gamma(:,:,1))));
%     fprintf(1,'mean A1=%g\n',mean(mean(theta.A(:,:,1))));
%     fprintf(1,'mean c1=%g\n',mean(theta.c(:,1)));    
%     fprintf(1,'mean b1=%g\n',mean(theta.b(:,1)));        
%     fprintf(1,'pi1=%g\n',theta.pi(1));        
    
    % =====================EXPECTATION STEPS=========================== 
    [r,LL(iter),ec] = ExpectationZ(t,y,theta,verb);
    [theta,cstr]=remove_empty_clusters(theta,cstr,ec);
    
    [muw,Sw] = ExpectationW(t,y,theta,verb);    

    if(iter>=3)
            deltaLL_total=max(LL(1:iter))-min(LL(1:iter));
            deltaLL=LL(iter)-LL(iter-1);
            converged=(deltaLL < 0.001*deltaLL_total);
    end    
    
    if(verb>=1);fprintf(1,'\n');end;
end
%%% Final log-likelihood %%%%
LLf=LL(iter);


% =============================Final plots===============================
if(verb>=1);fprintf(1,'Converged in %d iterations\n',iter);end;
if(verb>=2)
    fig=figure;clf(fig);
    plot(LL);
    [~,cluster_idx]=max(r,[],2);
%     for d=1:D/2
%         fig=figure;clf(fig);    
%         scatter(y(2*d-1,:),y(2*d,:),200,cluster_idx');
%     end
%     fig=figure;clf(fig);    
%     scatter(t(1,:),t(2,:),200,cluster_idx','filled');
end
end

function  [r,LL,ec] = ExpectationZ(t,y,th,verb)
    if(verb>=1);fprintf(1,'  EZ'); end;
    if(verb>=3);fprintf(1,' k='); end;    
    [D,N]=size(y);
    K=length(th.pi);
    Lt=size(t,1);
    L=size(th.c,1);
    Lw=L-Lt;
    logr=NaN(N,K);
    for k=1:K
        if(verb>=3);fprintf(1,'%d,',k); end;          
        muyk=th.b(:,k); % Dx1
        covyk=reshape(th.Sigma(:,:,k),D,D); % DxD
        if(Lt>0);
            Atk=reshape(th.A(:,1:Lt,k),D,Lt); % DxLt     
            muyk=bsxfun(@plus,muyk,Atk*t); % DxN
        end
        if(Lw>0);
            Awk=reshape(th.A(:,Lt+1:L,k),D,Lw); % DxLw
            Gammawk=reshape(th.Gamma(Lt+1:L,Lt+1:L,k),Lw,Lw); % LwxLw
            cwk=th.c(Lt+1:L,k); % Lwx1
            covyk=covyk+Awk*Gammawk*Awk'; % DxD
            muyk=bsxfun(@plus,muyk,Awk*cwk); % DxN
        end
        logr(:,k) = log(th.pi(k))*ones(N,1);
        logr(:,k) = logr(:,k) + loggausspdf(y,muyk,covyk)';            
        if(Lt>0)
            logr(:,k) = logr(:,k)+...
                        loggausspdf(t,th.c(1:Lt,k),th.Gamma(1:Lt,1:Lt,k))';     
        end
    end
    lognormr=logsumexp(logr,2);
    LL=sum(lognormr);
    r=exp(bsxfun(@minus,logr,lognormr));
    
    % remove empty clusters
    ec=true(1,K); % false if component k is empty.
    for k=1:K
        if(sum(r(:,k))==0 || ~isfinite(sum(r(:,k))))
            ec(k)=false;
            if(verb>=1)
                fprintf(1,'     WARNING: CLASS %d HAS BEEN REMOVED\n',k);
            end
        end
    end
    if(sum(ec)==0)
        fprintf(1,'REINIT! ');
        [~, ~, ~, r] = emgm([t;y], K, 2, verb);
        ec=true(1,size(r,2));
    else
        r=r(:,ec);
    end
end

function [muw,Sw] = ExpectationW(t,y,th,verb)
    if(verb>=1);fprintf(1,'  EW'); end;
    if(verb>=3);fprintf(1,' k='); end;
    [D,N]=size(y);
    K=length(th.pi);
    Lt=size(t,1);
    L=size(th.c,1);
    Lw=L-Lt;
    if(Lw==0)
        muw=[];
        Sw=[];
        return;
    end
    Sw=zeros(Lw,Lw,K);
    muw=zeros(Lw,N,K);
    for k=1:K
        if(verb>=3);fprintf(1,'%d,',k); end;        
        Atk=reshape(th.A(:,1:Lt,k),D,Lt); %DxLt       
        Awk=reshape(th.A(:,Lt+1:L,k),D,Lw); %DxLw
        Sigmak=reshape(th.Sigma(:,:,k),D,D); %DxD
        Gammawk=reshape(th.Gamma(Lt+1:L,Lt+1:L,k),Lw,Lw);%LwxLw
        cwk=th.c(Lt+1:L,k); %Lwx1
        invSwk=eye(Lw)+Gammawk*Awk'/Sigmak*Awk; %LwxLw
        if(~isempty(t))
            Atkt=Atk*t;
        else
            Atkt=0;
        end
        muw(:,:,k)=invSwk\bsxfun(@plus,Gammawk*Awk'/Sigmak*...
                                bsxfun(@minus,y-Atkt,th.b(:,k)),cwk); %LwxN
        Sw(:,:,k)=invSwk\Gammawk;
    end
end

function  th = Maximization(t,y,r,muw,Sw,cstr,verb)
    if(verb>=1);fprintf(1,'  M'); end;
    if(verb>=3);fprintf(1,' k='); end;    
    K=size(r,2);
    [D,N]=size(y);
    Lt=size(t,1);
    Lw=size(muw,1);
    L=Lt+Lw;
    
    th.c=NaN(L,K);
    th.Gamma=zeros(L,L,K);
    if(Lw>0)
        th.c(Lt+1:L,:)=cstr.cw; % LwxK
        th.Gamma(Lt+1:L,Lt+1:L,:)=cstr.Gammaw; % LwxLwxK
    end
    th.pi=NaN(1,K);    
    th.A=NaN(D,L,K);
    th.b=NaN(D,K);
    th.Sigma=NaN(D,D,K);  
    
    rk_bar=zeros(1,K);
    for k=1:K
        if(verb>=3);fprintf(1,'%d:',k); end;    
        % Posteriors' sums
        rk=r(:,k)'; % 1xN         
        rk_bar(k)=sum(rk); % 1x1
        
        if(Lt>0)
            if(verb>=3);fprintf(1,'c'); end;
            % Compute optimal mean ctk  
            if(isempty(cstr.ct))
                th.c(1:Lt,k)=sum(bsxfun(@times,rk,t),2)./rk_bar(k); % Ltx1
            else
                th.c(1:Lt,k)=cstr.ct(:,k);
            end
            % Compute optimal covariance matrix Gammatk
            if(verb>=3);fprintf(1,'Gt'); end;
            diffGamma=bsxfun(@times,sqrt(rk),bsxfun(@minus,t,th.c(1:Lt,k))); % LtxN
            if(isempty(cstr.Gammat) || (length(cstr.Gammat)==1 && cstr.Gammat=='*'))
                %%%% Full Gammat
                th.Gamma(1:Lt,1:Lt,k)=diffGamma*diffGamma'./rk_bar(k); % DxD
                th.Gamma(1:Lt,1:Lt,k)=th.Gamma(1:Lt,1:Lt,k);                    
            elseif(~ischar(cstr.Gammat))
                %%%% Fixed Gammat
                th.Gamma(1:Lt,1:Lt,k)=cstr.Gammat(:,:,k);            
            elseif(cstr.Gammat(1)=='d' || cstr.Gammat(1)=='i')
                % Diagonal terms   
                gamma2=sum(diffGamma.^2,2)./rk_bar(k); %Ltx1
                if(cstr.Gammat(1)=='d')
                    %%% Diagonal Gamma
                    th.Gamma(1:Lt,1:Lt,k)=diag(gamma2); % LtxLt  
                else
                    %%% Isotropic Gamma
                    th.Gamma(1:Lt,1:Lt,k)=mean(gamma2)*eye(Lt); % LtxLt
                end
            elseif(cstr.Gammat(1)=='v')
                %%%% Full Gamma
                th.Gamma(1:Lt,1:Lt,k)=diffGamma*diffGamma'./rk_bar(k); % LtxLt
            else
                cstr.Gammat,
                error('  ERROR: invalid constraint on Gamma.');
            end
        end

        % Compute optimal weight pik
        th.pi(k)=rk_bar(k)/N; % 1x1

        if(Lw>0)
            x=[t;reshape(muw(:,:,k),Lw,N)]; % LxN
            Skx=[zeros(Lt),zeros(Lt,Lw);zeros(Lw,Lt),Sw(:,:,k)]; % LxL    
        else
            x=t; % LtxN
            Skx=zeros(Lt); %LtxLt
        end
        
        if(verb>=3);fprintf(1,'A'); end;
        if(isempty(cstr.b))
            % Compute weighted means of y and x
            yk_bar=sum(bsxfun(@times,y,rk),2)./rk_bar(k); % Dx1
            if(L>0);
                xk_bar=sum(bsxfun(@times,x,rk),2)./rk_bar(k); % Lx1
            else
                xk_bar=[];
            end;
        else
            yk_bar=cstr.b(:,k);
            xk_bar=zeros(L,1);
            th.b(:,k)=cstr.b(:,k);  
        end
        % Compute weighted, mean centered y and x
        weights=sqrt(rk); % 1xN        
        y_stark=bsxfun(@minus,y,yk_bar); % DxN
        y_stark=bsxfun(@times,weights,y_stark) ./ sqrt(rk_bar(k)); % DxN         
        if(L>0);
            x_stark=bsxfun(@minus,x,xk_bar); % LxN  
            x_stark=bsxfun(@times,weights,x_stark) ./ sqrt(rk_bar(k)); % LxN            
        else
            x_stark=[];
        end;
        
        % Robustly compute optimal transformation matrix Ak
        warning off MATLAB:nearlySingularMatrix;
        if(~all(all(Skx==0)))
            if(N>=L && det(Skx+x_stark*x_stark')>10^(-8))
                th.A(:,:,k)=y_stark*x_stark'/(Skx+x_stark*x_stark'); % DxL
            else
                th.A(:,:,k)=y_stark*x_stark'*pinv(Skx+x_stark*x_stark'); %DxL
            end
        elseif(~all(all(x_stark==0)))
            if(N>=L && det(x_stark*x_stark')>10^(-8))
               th.A(:,:,k)=y_stark*x_stark'/(x_stark*x_stark'); % DxL
            elseif(N<L && det(x_stark'*x_stark)>10^(-8))
               th.A(:,:,k)=y_stark/(x_stark'*x_stark)*x_stark'; % DxL
            else
                if(verb>=3);fprintf(1,'p'); end;
                th.A(:,:,k)=y_stark*pinv(x_stark);  % DxL
            end
        else
            % Correspond to null variance in cluster k or L=0:
            if(verb>=1 && L>0);fprintf(1,'null var\n');end;
            th.A(:,:,k)=0;  % DxL
        end

        if(verb>=3);fprintf(1,'b'); end;
        % Intermediate variable wk=y-Ak*x
        if(L>0)
            wk=y-reshape(th.A(:,:,k),D,L)*x; % DxN
        else
            wk=y;
        end

        % Compute optimal transformation vector bk
        if(isempty(cstr.b))
            th.b(:,k)=sum(bsxfun(@times,rk,wk),2)./rk_bar(k); % Dx1
        end

        if(verb>=3);fprintf(1,'S'); end;
        % Compute optimal covariance matrix Sigmak
        if(Lw>0)
            Awk=reshape(th.A(:,Lt+1:L,k),D,Lw);
            Swk=reshape(Sw(:,:,k),Lw,Lw);                
            ASAwk=Awk*Swk*Awk';
        else
            ASAwk=0;
        end
        diffSigma=bsxfun(@times,sqrt(rk),bsxfun(@minus,wk,th.b(:,k))); %DxN
        if(isempty(cstr.Sigma) || (length(cstr.Sigma)==1 && cstr.Sigma=='*'))
            %%%% Full Sigma  
            th.Sigma(:,:,k)=diffSigma*diffSigma'./rk_bar(k); % DxD
            th.Sigma(:,:,k)=th.Sigma(:,:,k)+ASAwk;                    
        elseif(~ischar(cstr.Sigma))
            %%%% Fixed Sigma
            th.Sigma=cstr.Sigma;
        elseif(cstr.Sigma(1)=='d' || cstr.Sigma(1)=='i')
            % Diagonal terms   
            sigma2=sum(diffSigma.^2,2)./rk_bar(k); %Dx1
            if(cstr.Sigma(1)=='d')
                %%% Diagonal Sigma
                th.Sigma(:,:,k)=diag(sigma2); % DxD
                th.Sigma(:,:,k)=th.Sigma(:,:,k)+diag(diag(ASAwk));                 
            else
                %%% Isotropic Sigma
                th.Sigma(:,:,k)=mean(sigma2)*eye(D); % DxD
                th.Sigma(:,:,k)=th.Sigma(:,:,k)+trace(ASAwk)/D*eye(D);
            end                             
        else
            cstr.Sigma,
            error('  ERROR: invalid constraint on Sigma.');
        end
        % Avoid numerical problems on covariances:
        if(verb>=3);fprintf(1,'n'); end;
        if(~isfinite(sum(sum(th.Gamma(1:Lt,1:Lt,k)))))
            th.Gamma(1:Lt,1:Lt,k)=0;
        end
        th.Gamma(1:Lt,1:Lt,k)=th.Gamma(1:Lt,1:Lt,k)+1e-8*eye(Lt);
        if(~isfinite(sum(sum(th.Sigma(:,:,k))))) 
            th.Sigma(:,:,k)=0;
        end
        th.Sigma(:,:,k)=th.Sigma(:,:,k)+1e-8*eye(D);
        if(verb>=3);fprintf(1,','); end;
    end
    if(verb>=3);fprintf(1,'end'); end;

    if(ischar(cstr.Sigma) && ~isempty(cstr.Sigma) && ...
                             cstr.Sigma(length(cstr.Sigma))=='*')
        %%% Equality constraint on Sigma
        th.Sigma=bsxfun(@times,th.Sigma,reshape(rk_bar,[1,1,K]));
        th.Sigma=bsxfun(@times,ones(size(th.Sigma)),sum(th.Sigma,3))./N;
    end
        
    if(cstr.Gammat=='v')
        %%% Equal volume constraint on Gamma
        detG=zeros(1,K);
        for k=1:K
            detG(k)=det(reshape(th.Gamma(1:Lt,1:Lt,k),Lt,Lt)); % 1x1
        end
        th.Gamma(1:Lt,1:Lt,:)=bsxfun(@rdivide,th.Gamma(1:Lt,1:Lt,:),...
                                              reshape(detG,1,1,K));
        th.Gamma(1:Lt,1:Lt,:)=sum(detG.^(1/Lt).*th.pi).*...
                              th.Gamma(1:Lt,1:Lt,:);
    end
    
    if(ischar(cstr.Gammat) && ~isempty(cstr.Gammat) && ...
                             cstr.Gammat(length(cstr.Gammat))=='*')
        %%% Equality constraint on Gammat
        th.Gamma(1:Lt,1:Lt,:)=bsxfun(@times,th.Gamma(1:Lt,1:Lt,:),...
                                           reshape(rk_bar,[1,1,K]));        
        th.Gamma(1:Lt,1:Lt,:)=bsxfun(@times,ones(Lt,Lt,K),...
                                          sum(th.Gamma(1:Lt,1:Lt,:),3))./N;
    end    
    
    if(~ischar(cstr.pi) || isempty(cstr.pi))
        if(~isempty(cstr.pi))
            th.pi=cstr.pi;
        end
    else
        if(cstr.pi(1)=='*')
            th.pi=1/K.*ones(1,K);
        else
            error('  ERROR: invalid constraint on pi.');            
        end
    end 
end

function [th,cstr]=remove_empty_clusters(th1,cstr1,ec)
    th=th1;
    cstr=cstr1;
    if(sum(ec)~=length(ec))
        if(~isempty(cstr.ct) && ~ischar(cstr.ct))
            cstr.ct=cstr.ct(:,ec);end;    
        if(~isempty(cstr.cw) && ~ischar(cstr.cw))
            cstr.cw=cstr.cw(:,ec);end;     
        if(~isempty(cstr.Gammat) && ~ischar(cstr.Gammat))
            cstr.Gammat=cstr.Gammat(:,:,ec);end;
        if(~isempty(cstr.Gammaw) && ~ischar(cstr.Gammaw))
            cstr.Gammaw=cstr.Gammaw(:,:,ec);end;    
        if(~isempty(cstr.pi) && ~ischar(cstr.pi))
            cstr.pi=cstr.pi(:,ec);end;       
        if(~isempty(cstr.A) && ~ischar(cstr.A))
            cstr.A=cstr.A(:,:,ec);end;  
        if(~isempty(cstr.b) && ~ischar(cstr.b))
            cstr.b=cstr.b(:,ec);end;     
        if(~isempty(cstr.Sigma) && ~ischar(cstr.Sigma))
            cstr.Sigma=cstr.Sigma(:,:,ec);end;        

        th.c=th.c(:,ec);
        th.Gamma=th.Gamma(:,:,ec);
        th.pi=th.pi(ec);
        th.A=th.A(:,:,ec);
        th.b=th.b(:,ec);
        th.Sigma=th.Sigma(:,:,ec); 
    end
end