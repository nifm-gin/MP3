function [theta, r, LLf, LL] = sllim(tapp, yapp, in_K, varargin)
  %%%%%%%% General EM Algorithm for Gaussian Locally Linear Mapping %%%%%%%%%
  %%% Author: Antoine Deleforge (April 2013) - antoine.deleforge@inria.fr %%%
  % Description: Compute maximum likelihood parameters theta and posterior
  % probabilities r=p(z_n=k|x_n,y_n;theta) of a sllim model with constraints
  % cstr using N associated observations t and y.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% Input %%%%
  %- t (LtxN) % Training latent variables
  %- y (DxN)  % Training observed variables
  %- in_K (int)  % Initial number of components
  % <Optional>
  %- Lw (int) % Dimensionality of hidden components (default 0)
  %- maxiter (int) % Maximum number of iterations (default 100)
  %- in_theta (struct)% Initial parameters (default [])
  %| same structure as output theta
  %- in_r (NxK)  % Initial assignments (default [])
  %- cstr (struct) % Constraints on parameters theta (default [],'')
  %- cstr.ct  % fixed value (LtxK) or ''=uncons.
  %- cstr.cw  % fixed value (LwxK) or ''=fixed to 0
  %- cstr.Gammat% fixed value (LtxLtxK) or ''=uncons.
  %| or {'','d','i'}{'','*','v'} (1)
  %- cstr.Gammaw% fixed value (LwxLwxK) or ''=fixed to I
  %- cstr.pi  % fixed value (1xK) or ''=uncons. or '*'=equal
  %- cstr.A  % fixed value (DxL) or ''=uncons.
  %- cstr.b  % fixed value (DxK) or ''=uncons.
  %- cstr.Sigma% fixed value (DxDxK) or ''=uncons.
  %| or {'','d','i'}{'','*'} (1)
  %- verb {0,1,2}% Verbosity (default 1)
  %%%% Output %%%%
  %- theta  (struct)  % Estimated parameters (L=Lt+Lw)
  %- theta.c (LxK) % Gaussian means of X
  %- theta.Gamma (LxLxK) % Gaussian covariances of X
  %- theta.A (DxLxK)  % Affine transformation matrices
  %- theta.b (DxK) % Affine transformation vectors
  %- theta.Sigma (DxDxK) % Error covariances
  %- theta.gamma (1xK)% Arellano-Valle and Bolfarine's Generalized t
  %-  phi  (struct)% Estimated parameters
  %- phi.pi (1xK)  % t weights of X
  %- phi.alpha (1xK)  % Arellano-Valle and Bolfarine's Generalized t
  %- r (NxK)  % Posterior probabilities p(z_n=k|x_n,y_n;theta)
  %- u (NxK)  % Posterior probabilities E[u_n|z_n=k,x_n,y_n;theta,phi)
  %%% (1) 'd'=diag., 'i'=iso., '*'=equal for all k, 'v'=equal det. for all k
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ======================Input Parameters Retrieval=========================
  [Lw, maxiter, in_theta, in_phi, in_r, cstr,verb] = ...
  process_options(varargin,'Lw',0,'maxiter',100,'in_theta',[],...
                  'in_r',[],'in_phi',[],'cstr',struct(),'verb',0);

  % ==========================Default Constraints============================
  if(~isfield(cstr,'ct'));cstr.ct=[];end
  if(~isfield(cstr,'cw'));cstr.cw=[];end
  if(~isfield(cstr,'Gammat'));cstr.Gammat=[];end
  if(~isfield(cstr,'Gammaw'));cstr.Gammaw=[];end
  if(~isfield(cstr,'pi'));cstr.pi=[];end
  if(~isfield(cstr,'A'));cstr.A=[];end
  if(~isfield(cstr,'b'));cstr.b=[];end
  if(~isfield(cstr,'Sigma'));cstr.Sigma=[];end

  % ==========================EM initialization==============================
  Lt = size(tapp,1);
  L  = Lt+Lw;
  [D,N] = size(yapp);

  if(verb>=1), fprintf(1,'EM Initializations\n'); end
  if(~isempty(in_theta) && ~isempty(in_phi))
    theta = in_theta;
    ph = in_phi;
    K = leng;
  end
  if(~isempty(in_theta) && ~isempty(in_phi))
      theta=in_theta;
      ph = in_phi;
      K=length(theta.pi);
      if(isempty(cstr.cw));
          cstr.cw=zeros(L,K); % Default value for cw
      end
      if(isempty(cstr.Gammaw));
          cstr.Gammaw=repmat(eye(Lw),[1,1,K]); % Default value for Gammaw
      end
      [r, ~, ec, ~, ~, ~] = expectation_ZU(tapp,yapp,theta,verb);
      [theta, ~, cstr] = remove_empty_clusters(theta,cstr,ec);
	  [muw, Sw] = ExpectationW(tapp,yapp,u,theta,verb);
      if(verb>=1);fprintf(1,'\n');end
  else
    if(isempty(in_r))
      % Initialise posteriors with K-means + GMM on joint observed data
      % Fast initialization for large datasets
      if verb>=1; fprintf('Running k-means... '); end
      if N>1000
          p = randperm(N);
          [~,C,~] = kmeans([tapp(:,p(1:1000));yapp(:,p(1:1000))],in_K);
      else
          [~,C,~] = kmeans([tapp;yapp],in_K);
      end
      fprintf('done.\n')
      [~, ~, ~, r] = emgm([tapp;yapp], C, 3, verb);
      if verb>=1
          fprintf('done.\n')
          fprintf('Running GMM... \n')
      end
    else
      r = in_r;
    end

    %Initializing parameters
    phi.pi = mean(r,2)'; % 1xK
    u = ones(N,in_K) ; % NxK
    K = size(r,2); % extracting K, the number of affine components
    phi.alpha = 20*ones(1,K); % 1xK
    mahalt=[];
    mahaly=[];

    if Lw == 0
        Sw=[];
        muw=[];
    else
        % Start by running an M-step without hidden variables (partial
        % theta), deduce Awk by local weighted PCA on residuals (complete
        % theta), deduce r, muw and Sw from E-steps on complete theta.
        [theta, phi] = Maximization_hybrid(tapp,yapp,r,u,phi,[],[],[],[],cstr,verb);

        K=length(phi.pi);
        if(isempty(cstr.cw))
            cstr.cw=zeros(Lw,K); % Default value for cw
        end
        theta.c=cat(1,theta.c,cstr.cw(:,1:K)); %LxK
        Gammaf=zeros(L,L,K);
        Gammaf(1:Lt,1:Lt,:)=theta.Gamma; %LtxLtxK
        if(isempty(cstr.Gammaw))
            cstr.Gammaw=repmat(eye(Lw),[1,1,K]); % Default value for Gammaw
        end
        Gammaf(Lt+1:L,Lt+1:L,:)=cstr.Gammaw(:,:,1:K); %LwxLwxK
        theta.Gamma=Gammaf;
        % Initialize Awk with local weighted PCAs on residuals:
        Aw=zeros(D,Lw,K);
        for k=1:K
            rk_bar=sum(r(:,k));
            bk=reshape(theta.b(:,k),D,1);
            w=bsxfun(@minus,yapp,bk); %DxN
            if(Lt>0)
                Ak=reshape(theta.A(:,:,k),D,Lt);
                w=w-Ak*tapp;
            end
            w=bsxfun(@times,w,sqrt(r(:,k)'./rk_bar)); %DxN
            C=w*w'; % Residual weighted covariance matrix
            [U,Lambda]=eigs(C,Lw); % Weighted residual PCA U:DxLw
            % The residual variance is the discarded eigenvalues' mean
            sigma2k=(trace(C)-trace(Lambda))./(D-Lw);
            theta.Sigma(:,k)=sigma2k*ones(D,1);
            Aw(:,:,k)=U*sqrt(Lambda-sigma2k*eye(Lw));
        end
        theta.A=cat(2,theta.A,Aw); %DxLxK
        [r, ~, ec, u, mahalt, mahaly] = ExpectationZU(tapp,yapp,theta,phi,verb);
        [theta, phi, cstr] = remove_empty_clusters(theta,phi,cstr,ec);
        [muw, Sw] = ExpectationW(tapp,yapp,theta,phi,verb);
        if(verb>=1) fprintf(1,'\n');  end
    end
  end

  %===============================EM Iterations==============================

  if(verb>=1), fprintf(1,'Running EM\n'); end
  LL = -inf(maxiter,1);
  iter = 0;
  converged=false;
  while (~converged && iter<maxiter)
      iter = iter + 1;
      if(verb>=1), fprintf(1,'      Iteration %d\n',iter); end

      % =====================MAXIMIZATION STEP===========================
      [theta, phi] = Maximization_hybrid(tapp,yapp,r,u,phi,muw,Sw,mahalt,mahaly,cstr,verb);

      % =====================EXPECTATION STEPS===========================
      [r, LL(iter), ec, u, mahalt, mahaly] = ExpectationZU(tapp,yapp,theta,phi,verb);

      [theta, phi, cstr] = remove_empty_clusters(theta,phi,cstr,ec);

      [muw, Sw] = ExpectationW(tapp,yapp,theta,phi,verb);

      if(iter>=3)
              deltaLL_total=max(LL(1:iter))-min(LL(1:iter));
              deltaLL=LL(iter)-LL(iter-1);
              converged=(deltaLL < 0.001*deltaLL_total);
      end

      if(verb>=1);fprintf(1,'\tLog-like: %f\n',LL(iter));end
  end

  %% Final log-likelihood %%%
  LLf=LL(iter);

  Sigma_full = zeros(D,D,K);
  for k=1:K
      Sigma_full(:,:,k) = diag(theta.Sigma(:,k));
  end
  theta.Sigma = Sigma_full;


  if isempty(cstr.Sigma)
      nbparSigma = D*(D+1)/2;
  elseif cstr.Sigma(1)=='d'
      nbparSigma = D;
  elseif cstr.Sigma(1)=='i'
      nbparSigma = 1;
  elseif cstr.Sigma(1) == '*'
      nbparSigma = D*(D+1)/(2*K);
  else
      error('  ERROR: invalid constraint on Sigma.');
  end

  if ~isempty(cstr.Gammat)
      if (cstr.Gammat(1) == 'i')
          nbparGamma = 1;
      elseif (cstr.Gammat(1) == 'd')
          nbparGamma = Lt;
      elseif (cstr.Gammat(1) == '')
          nbparGamma = Lt*(Lt+1)/2;
      else (cstr.Gammat(1) == '*')
          nbparGamma = Lt*(Lt+1)/(2*K);
      end
  else
    nbparGamma = Lt*(Lt+1)/2;
  end
  theta.phi = phi;
  theta.nbpar = (in_K-1) + in_K*(1 + D*L + D + Lt + nbparSigma + nbparGamma);
  % =============================Final plots===============================
  if(verb>=1);fprintf(1,'Converged in %d iterations\n',iter);end
end  % end sllim()

function [r, LL, ec, u, mahalt, mahaly] = ExpectationZU(tapp, yapp, th, ph, verb)
  if(verb>=1);  fprintf(1,'  EZ'); end
  if(verb>=3);  fprintf(1,'  k='); end

  [D,N] = size(yapp);
  K = length(ph.pi);
  Lt = size(tapp,1);
  L = size(th.c,1);
  Lw = L-Lt;
  logr = NaN(N,K);

  % Space allocation
  u = NaN(N,K);  % matrix of weights

  mahaly = zeros(N,K);
  mahalt = zeros(N,K);

  for k=1:K
    if(verb>=3);fprintf(1,'%d,',k); end
    alphak = ph.alpha(k);  % 1x1
    muyk=th.b(:,k); % Dx1
    covyk= th.Sigma(:,k); % Dx1

    if Lt>0
        Atk = th.A(:,1:Lt,k); %% DxLt
        muyk = Atk*tapp + muyk;% DxN
        Gammawk = [];
    end

    if Lw>0
      Awk=reshape(th.A(:,Lt+1:L,k),D,Lw); % DxLw
      Gammawk=reshape(th.Gamma(Lt+1:L,Lt+1:L,k),Lw,Lw); % LwxLw
      cwk=th.c(Lt+1:L,k); % Lwx1
%       AwkGammawkAwkt = Awk*Gammawk*Awk'; % DxD
      muyk= bsxfun(@plus,muyk,Awk*cwk); % DxN
    end

    if isempty(Gammawk)
        mahaly(:,k) = mahalanobis_distance(yapp,muyk,covyk,'diag');
    else
        mahaly(:,k) = mahalanobis_distance(yapp,muyk,covyk,'diag_lowr',Awk,Gammawk);
    end

    mahalt(:,k) = mahalanobis_distance(tapp,th.c(1:Lt,k,:),th.Gamma(1:Lt,1:Lt,k),'full');

    u(:,k) = (alphak + (D+Lt)/2)./(1 + 0.5*(mahalt(:,k) + mahaly(:,k)));

    %modif EP
    %[R,p] = chol(covyk);
    [R,p]= chol(diag(covyk)); %dxd
	if p ~= 0
	    fprintf(1,'SNPD! ');
	    sldR=-Inf; % 1xn
    else
      sldR=sum(log(diag(R)));
    end

    %%%% modif EP
    logr(:,k) = repmat(log(ph.pi(k)),N,1) ...
                + logtpdfL(sldR,mahalt(:,k),alphak, Lt) ...
                + logtpdfD(sum(log(diag(chol(th.Gamma(1:Lt,1:Lt,k))))),mahaly(:,k),alphak + Lt/2, (1+mahalt(:,k)/2), D);  %% log_rnk
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
     [~, ~, ~, r] = emgm([tapp;yapp], K, 2, verb);
     ec=true(1,size(r,2));
  else
     r=r(:,ec);
  end
end % end expectation_ZU()

function [muw,Sw] = ExpectationW(tapp,yapp,th,ph,verb)
  if(verb>=1);fprintf(1,'  EW'); end
  if(verb>=3);fprintf(1,' k='); end
  [D,N]=size(yapp);
  K=length(ph.pi);
  Lt=size(tapp,1);
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
    if(verb>=3);fprintf(1,'%d,',k); end
    Atk=reshape(th.A(:,1:Lt,k),D,Lt); %DxLt
    if (Lw>0)
      Awk=reshape(th.A(:,Lt+1:L,k),D,Lw); %DxLw
      Gammawk=reshape(th.Gamma(Lt+1:L,Lt+1:L,k),Lw,Lw);%LwxLw
      cwk=th.c(Lt+1:L,k); %Lwx1
      GammawkAwktinvSigmak = Gammawk*bsxfun(@rdivide, Awk, th.Sigma(:,k))'; % LwxD
      invSwk=eye(Lw)+GammawkAwktinvSigmak*Awk; %LwxLw
    end
    if(~isempty(tapp))
        Atkt=Atk*tapp;
    else
        Atkt=0;
    end
    if Lw==0
      muw=[];
      Sw=[];
    else
      muw(:,:,k)=invSwk\bsxfun(@plus, GammawkAwktinvSigmak *...
                              bsxfun(@minus,yapp-Atkt,th.b(:,k)),cwk); %LwxN
      Sw(:,:,k)=invSwk\Gammawk;
    end
  end
end

function [th, ph] = Maximization(tapp,yapp,r,u,phi,muw,Sw,mahalt,mahaly,cstr,verb)
    if(verb>=1);fprintf(1,'  M'); end
    if(verb>=3);fprintf(1,' k='); end
    K=size(r,2);
    [D,N]=size(y);
    Lt=size(t,1);
    Lw=size(muw,1);
    L=Lt+Lw;

    th.c = NaN(L,K);
    th.Gamma = zeros(L,L,K);
    if(Lw>0)
        th.c(Lt+1:L,:)=cstr.cw; % LwxK
        th.Gamma(Lt+1:L,Lt+1:L,:)=cstr.Gammaw; % LwxLwxK
    end
    th.A=NaN(D,L,K);
    th.b=NaN(D,K);
    th.Sigma=NaN(D,D,K); % DxDxK Sigma matrix

    ph.pi=NaN(1,K);

    rk_bar=zeros(1,K);
    for k=1:K
        if(verb>=3);fprintf(1,'%d:',k); end
        % Posteriors' sums
        rk=r(:,k)'; % 1xN
        rk_bar(k)=sum(rk); % 1x1

    		uk=u(:,k);  % 1xN
    		rk_tilde = rk.*uk; % 1xN
    		rk_bar_tilde = sum(rk_tilde); % 1xK

        if Lt>0
            if(verb>=3);fprintf(1,'c'); end
            % Compute optimal mean ctk
            if(isempty(cstr.ct))
                th.c(1:Lt,k)=sum(bsxfun(@times,rk_tilde,tapp),2)./rk_bar_tilde(k); % Ltx1
            else
                th.c(1:Lt,k)=cstr.ct(:,k);
            end

            % Compute optimal covariance matrix Gammatk
            if(verb>=3);fprintf(1,'Gt'); end
            diffGamma=bsxfun(@times,sqrt(rk),bsxfun(@minus,tapp,th.c(1:Lt,k))); % LtxN
            if( is.null(cstr.Gammat) || (length(cstr.Gammat)==1 && cstr.Gammat=='*')) % | ou ||?
               %%%% Full Gammat
               th.Gamma(1:Lt,1:Lt,k)=diffGamma*diffGamma'./rk_bar(k); % LtxLt
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
    		ph.pi(k) = rk_bar(k)/N; % 1x1

        if Lw > 0
          x=[tapp; reshape(muw(:,:,k),Lw,N)]; % LxN
          Skx=[zeros(Lt),zeros(Lt,Lw);zeros(Lw,Lt),Sw(:,:,k)]; % LxL
        else
          x=tapp; % LtxN
          Skx=zeros(Lt); %LtxLt
        end

        if(verb>=3);fprintf(1,'A'); end
        if(isempty(cstr.b))
            % Compute weighted means of y and x
            yk_bar=sum(bsxfun(@times,yapp,rk_tilde),2)./rk_bar_tilde(k); % Dx1
            if(L>0);
                xk_bar=sum(bsxfun(@times,x,rk_tilde),2)./rk_bar_tilde(k); % Lx1
            else
                xk_bar=[];
            end
        else
            yk_bar=cstr.b(:,k);
            xk_bar=zeros(L,1);
            th.b(:,k)=cstr.b(:,k);
        end

    		% Compute weighted, mean centered y and x
    		weights=sqrt(rk_tilde)/sqrt(rk_bar(k)); % 1xN
        y_stark=bsxfun(@minus,yapp,yk_bar); % DxN
        y_stark=bsxfun(@times,weights,y_stark); % DxN
        if(L>0)
          x_stark=bsxfun(@minus,x,xk_bar); % LxN
          x_stark=bsxfun(@times,weights,x_stark); % LxN
        else
          x_stark=[];
        end

        % Robustly compute optimal transformation matrix Ak
        if(~all(all(Skx==0)))
            if(N>=L && rcond(Skx+x_stark*x_stark')>1e-4)
                th.A(:,:,k)=y_stark*x_stark'/(Skx+x_stark*x_stark'); % DxL
            else
                th.A(:,:,k)=y_stark*x_stark'*pinv(Skx+x_stark*x_stark'); %DxL
            end
        elseif(~all(all(x_stark==0)))
            if(rcond(x_stark*x_stark')>1e-4)
               th.A(:,:,k)=y_stark*x_stark'/(x_stark*x_stark'); % DxL
            elseif(rcond(x_stark'*x_stark)>1e-4)
               th.A(:,:,k)=y_stark/(x_stark'*x_stark)*x_stark'; % DxL
            else
                if(verb>=3);fprintf(1,'p'); end
                th.A(:,:,k)=y_stark*pinv(x_stark);  % DxL
            end
        else
            % Correspond to null variance in cluster k or L=0:
            if(verb>=1 && L>0);fprintf(1,'null var\n');end
            th.A(:,:,k)=0;  % DxL
        end

        % Compute optimal transformation vector bk
        if(verb>=3);fprintf(1,'b'); end
        % Intermediate variable wk=y-Ak*x
        if(L>0)
            wk=yapp-reshape(th.A(:,:,k),D,L)*x; % DxN
        else
            wk=yapp;
        end
        if(isempty(cstr.b))
            th.b(:,k)=sum(bsxfun(@times,rk_tilde,wk),2)./rk_bar_tilde(k); % Dx1
        end

        % Compute optimal covariance matrix Sigmak
        if(verb>=3);fprintf(1,'S'); end
        if Lw>0
          Awk=reshape(th.A(:,Lt+1:L,k),D,Lw);
          Swk=reshape(Sw(:,:,k),Lw,Lw);
          AwkSwk = Awk*Swk; % D*Lw
          diag_ASAwk = dot(Awk,AwkSwk,2); % Dx1 diag(Awk*Swk*Awk')
        else
            ASAwk=0;
        end

        diffSigma=bsxfun(@times,sqrt(rk_tilde),bsxfun(@minus,wk,th.b(:,k))); %DxN

        if(isempty(cstr.Sigma) || (length(cstr.Sigma)==1 && cstr.Sigma=='*'))
            %%%% Full Sigma
            th.Sigma(:,:,k)=diffSigma*diffSigma'./rk_bar(k); % DxD
            th.Sigma(:,:,k)=th.Sigma(:,:,k)+ASAwk;
        elseif(~ischar(cstr.Sigma))
            %%%% Fixed Sigma
            th.Sigma=cstr.Sigma;
          elseif(cstr.Sigma(1)=='d' || cstr.Sigma(1)=='i')
              % Diagonal terms
              sigma2=sum(diffSigma.^2,2)./rk_bar(k); % Dx1
              %%% Diagonal Sigma
              th.Sigma(:,k)=sigma2 + diag_ASAwk; % Dx1
              %%% Isotropic Sigma
              if cstr.Sigma(1)=='i'
                  th.Sigma(:,k)=mean(th.Sigma(:,k))*ones(D,1);
              end
          else
              cstr.Sigma,
              error('  ERROR: invalid constraint on Sigma.');
          end

        % Avoid numerical problems on covariances:
        if(verb>=3);fprintf(1,'n'); end
        if(~isfinite(sum(sum(th.Gamma(1:Lt,1:Lt,k)))))
            th.Gamma(1:Lt,1:Lt,k)=0;
        end
        th.Gamma(1:Lt,1:Lt,k)=th.Gamma(1:Lt,1:Lt,k)+1e-8*eye(Lt);
        if(~isfinite(sum(sum(th.Sigma(:,k)))))
            th.Sigma(:,k)=0;
        end
        th.Sigma(:,k)=th.Sigma(:,k)+1e-8*ones(D,1);
        if(verb>=3);fprintf(1,','); end
    end
    if(verb>=3);fprintf(1,'end'); end

      if(ischar(cstr.Sigma) && ~isempty(cstr.Sigma) && ...
                               cstr.Sigma(length(cstr.Sigma))=='*')
          %%% Equality constraint on Sigma
          th.Sigma=bsxfun(@times,th.Sigma,reshape(rk_bar./N,[1,K]));
          th.Sigma=repmat(sum(th.Sigma,2),[1,K]);
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

      % Compute phi.alpha
      if(~isempty(mahalt) && ~isempty(mahaly))
        % potentielle erreur?passage difficile
        ph.alpha(k) ...
          = inv_psi(...
                (psi(...
                phi.alpha(k) + (D+Lt)/2) ...
                - (1/rk_bar(k)) ...
                .* sum( rk .* log(1 + (1/2) .* (mahaly(:,k) + mahalt(:,k)))) ...
             ));
        if verb>=3; fprintf(1, 'k = %d -> alpha = %f',k,ph.alpha(k)); end
      else
        ph.alpha = phi.alpha;
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

function [th, ph, cstr] = remove_empty_clusters(th1,phi,cstr1,ec)
    th=th1;
    ph=phi;
    cstr=cstr1;
    if(sum(ec) ~= length(ec))
      if(~isempty(cstr.ct) && ~ischar(cstr.ct))
          cstr.ct=cstr.ct(:,ec);end
      if(~isempty(cstr.cw) && ~ischar(cstr.cw))
          cstr.cw=cstr.cw(:,ec);end
      if(~isempty(cstr.Gammat) && ~ischar(cstr.Gammat))
          cstr.Gammat=cstr.Gammat(:,:,ec);end
      if(~isempty(cstr.Gammaw) && ~ischar(cstr.Gammaw))
          cstr.Gammaw=cstr.Gammaw(:,:,ec);end
      if(~isempty(cstr.pi) && ~ischar(cstr.pi))
          cstr.pi=cstr.pi(:,ec);end
      if(~isempty(cstr.A) && ~ischar(cstr.A))
          cstr.A=cstr.A(:,:,ec);end
      if(~isempty(cstr.b) && ~ischar(cstr.b))
          cstr.b=cstr.b(:,ec);end
      if(~isempty(cstr.Sigma) && ~ischar(cstr.Sigma))
          cstr.Sigma=cstr.Sigma(:,ec);end

      th.c=th.c(:,ec);
      th.Gamma=th.Gamma(:,:,ec);
      ph.pi=ph.pi(ec);
      th.A=th.A(:,:,ec);
      th.b=th.b(:,ec);
      th.Sigma=th.Sigma(:,ec);
    end
end
