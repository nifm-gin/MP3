function [th, ph] = Maximization_hybrid(tapp,yapp,r,u,phi,muw,Sw,mahalt,mahaly,cstr,verb)
    if(verb>=1);fprintf(1,'  M'); end
    if(verb>=3);fprintf(1,' k='); end
    K=size(r,2);
    [D,N]=size(yapp);
    Lt=size(tapp,1);
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
    th.Sigma=NaN(D,K); % DxK Sigma matrix

    ph.pi=NaN(1,K);
    ph.alpha=NaN(1,K);

    rk_bar=zeros(1,K);
    for k=1:K
        if(verb>=3);fprintf(1,'%d:',k); end
        % Posteriors' sums
        rk=r(:,k)'; % 1xN
        rk_bar(k)=sum(rk); % 1x1

        uk=u(:,k)';  % 1xN
        rk_tilde = rk.*uk; % 1xN
        rk_bar_tilde = sum(rk_tilde); % 1x1

        if Lt>0
            if(verb>=3);fprintf(1,'c'); end
            % Compute optimal mean ctk
            if(isempty(cstr.ct))
                th.c(1:Lt,k)=sum(bsxfun(@times,rk_tilde,tapp),2)./rk_bar_tilde; % Ltx1
            else
                th.c(1:Lt,k)=cstr.ct(:,k);
            end

            % Compute optimal covariance matrix Gammatk
            if(verb>=3);fprintf(1,'Gt'); end
            diffGamma=bsxfun(@times,sqrt(rk_tilde),bsxfun(@minus,tapp,th.c(1:Lt,k))); % LtxN
            if(isempty(cstr.Gammat) || (length(cstr.Gammat)==1 && cstr.Gammat=='*')) %# | ou ||?
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
            yk_bar=sum(bsxfun(@times,yapp,rk_tilde),2)./rk_bar_tilde; % Dx1
            if(L>0)
                xk_bar=sum(bsxfun(@times,x,rk_tilde),2)./rk_bar_tilde; % Lx1
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
            th.b(:,k)=sum(bsxfun(@times,rk_tilde,wk),2)./rk_bar_tilde; % Dx1
        end

        % Compute optimal covariance matrix Sigmak
        if(verb>=3);fprintf(1,'S'); end
        if Lw>0
          Awk=reshape(th.A(:,Lt+1:L,k),D,Lw);
          Swk=reshape(Sw(:,:,k),Lw,Lw);
          AwkSwk = Awk*Swk; % D*Lw
          diag_ASAwk = dot(Awk,AwkSwk,2); % Dx1 diag(Awk*Swk*Awk')
        else
            diag_ASAwk=0;
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

          % Compute phi.alpha
          if(~isempty(mahalt) && ~isempty(mahaly))
            % potentielle erreur?passage difficile
            ph.alpha(k) = ...
              inv_psi(...
                (psi(...
                  phi.alpha(k) + (D+Lt)/2) ...
                  - (1/rk_bar(k)) ...
                  .* sum( rk' .* log(1 + (1/2) .* (mahaly(:,k) + mahalt(:,k))))...
                )...
              )';
            if verb>=3; fprintf(1, 'k = %d -> alpha = %f',k,ph.alpha(k)); end
          else
            ph.alpha = phi.alpha;
          end
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
