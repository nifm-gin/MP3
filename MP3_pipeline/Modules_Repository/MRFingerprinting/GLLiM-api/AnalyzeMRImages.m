function [Estimation, Parameters] = AnalyzeMRImages(Sequences,Dico,Method,Parameters,References,Outliers,SNRMap)

if nargin < 3, error('Not enought input arguments'); end
if ~exist('Method','var'),      Method      = 'DBL'; end
if ~exist('Parameters','var'),  Parameters  = []; end
if ~exist('References','var'),  References  = []; end
if ~exist('Outliers','var'),    Outliers    = []; end
if ~exist('SNRMap','var'),      SNRMap      = []; end

% This parameters is only used to enable the parameter data normalisation
normalization = 1;

%This line is only required to stay compatible with old notation 
switch Method
    case 'ClassicMRF'
        Method = 'DBM';
    case 'RegressionMRF'
        Method = 'DBL';
end

if isempty(Sequences), Sequences = Dico{1}.MRSignals; end
% Can be a problem, remove following line if necessary - 20/02/2019
if isempty(Parameters),	Parameters = struct(); end 
if isstruct(Dico), Dico = {Dico}; end

switch length(size(Sequences))
    case 4
        [s1,s2,t,slices] = size(Sequences);
    case 3
        [s1,s2,t]   = size(Sequences);
        slices      = 1;
    case 2
        [s1,t]      = size(Sequences);
        s2          = 1;
        Sequences   = reshape(Sequences, s1,s2,t);
        slices      = 1;
        if ~isempty(References) && length(size(References))==2
            References = reshape(References, s1,s2,size(References,2));
        end
    otherwise
        error('Invalid Sequences argument size')
end


f = 1;
switch Method
    
    case 'DBM'
        
        tic
        for s = 1:slices
            %Estimation of parameters
            Estimation.GridSearch.Y(:,:,:,s) = ...
                reshape(EstimateParametersFromGrid(reshape(Sequences(:,:,:,s),s1*s2,t), Dico{f}.MRSignals, Dico{f}.Parameters.Par), s1,s2, []);
                        
            %Errors computation if a reference image is provided
            if ~isempty(References)
                [Estimation.GridSearch.Errors.Rmse(s,:), Estimation.GridSearch.Errors.Nrmse(s,:), Estimation.GridSearch.Errors.Mae(s,:), Estimation.GridSearch.Errors.Nmae(s,:)] = ...
                    EvaluateEstimation(reshape(References(:,:,:,s),s1*s2,size(References,3)), reshape(Estimation.GridSearch.Y(:,:,:,s),s1*s2,size(References,3)));
            end
        end
        Estimation.GridSearch.quantification_time = toc; 
        
        
    case 'DBL'
        
        tic
        if ~any(strcmp(fieldnames(Parameters),'theta'))
            
            %Normalize trainning data
            if normalization == 1
                Parameters.factors.Ymean	= nanmean(Dico{f}.Parameters.Par,1);
                Parameters.factors.Ystd     = nanstd(Dico{f}.Parameters.Par);
                Dico{f}.Parameters.Par      = (Dico{f}.Parameters.Par - Parameters.factors.Ymean) ./ Parameters.factors.Ystd;
                
                Parameters.factors.normalization = 1;
            else
                Parameters.factors.Ymean	= 0;
                Parameters.factors.Ystd     = 1;
                                
                Parameters.factors.normalization = 0;
            end
            
            Xtrain = Dico{f}.MRSignals;
            [~,Parameters] = EstimateParametersFromRegression(Xtrain, Xtrain, Dico{f}.Parameters.Par, Dico{f}.Parameters.Par, Parameters);
        
        else
            Dico{f}.MRSignals       = [];
            Dico{f}.Parameters.Par  = [];
        end
        Estimation.Regression.learning_time = toc; 
        
        tic
        for s = 1:slices
            
            %Estimation of parameters
            if ~isempty(SNRMap)
                warning('off')
                
                %Recompute Sigma_k of the GLLiM model (if SNRMap given)
                Yestim = nan(s1*s2,size(Parameters.theta.A,2));
                Cov = nan(size(Parameters.theta.A,2),size(Parameters.theta.A,2),s1*s2);
                Pik = nan(s1*s2,Parameters.K);
                
                Param_updated = Parameters;
                SNRMap(SNRMap < 2) = 2;
                SNRMap(SNRMap > 202) = 202;
                var_noise = (max(Sequences(:,:,:,s),[],3) ./ SNRMap(:,:,s)).^2;
                
                Nb_model = 40;
                var_bounds = sort(1 ./ [2:200/Nb_model:202 inf].^2);
                for m = 1:length(var_bounds)-1
                    if any(reshape( (var_noise >= var_bounds(m)) & (var_noise < var_bounds(m+1)) ,1,[]))
                        var_noise_ = var_bounds(m) + (var_bounds(m+1) - var_bounds(m))/2;
                        Param_updated.theta = updateSigma(Parameters.theta,var_noise_);

                        Sequences_in = Sequences(:,:,:,s) .* ( (var_noise >= var_bounds(m)) & (var_noise < var_bounds(m+1)) ) ;
                        Sequences_in(Sequences_in == 0) = nan;
                        [Yestim_in,~,Cov_in,~,Pik_in] = ...
                            EstimateParametersFromRegression(reshape(Sequences_in(:,:,:),s1*s2,t), Dico{f}.MRSignals, Dico{f}.Parameters.Par, [], Param_updated);

                        Yestim = nansum(cat(3,Yestim,Yestim_in),3);
                        Cov = nansum(cat(4,Cov,Cov_in),4);
                        %Pik = nansum(cat(3,Pik,Pik_in),3);
                        Pik = [];
                    end
                end
                warning('on')

%                 for i = 1:s1
%                     for j = 1:s2
%                         Param_updated.theta = updateSigma(Parameters.theta,var_noise(i,j));
% 
%                         [Yestim((i-1)*s2+j,:),~,Cov(:,:,(i-1)*s2+j),~,Pik((i-1)*s2+j,:)] = ...
%                             EstimateParametersFromRegression(reshape(Sequences(i,j,:,s),1,[]), Dico{f}.MRSignals, Dico{f}.Parameters.Par, [], Param_updated);
%                     end
%                 end
            else
                [Yestim,~,Cov,~,Pik] = ...
                    EstimateParametersFromRegression(reshape(Sequences(:,:,:,s),s1*s2,t), Dico{f}.MRSignals, Dico{f}.Parameters.Par, [], Parameters);
            end
            
            %Rescale
            if any(strcmp(fieldnames(Parameters),'factors')) && normalization == 1
                Yestim(:,1:end-Parameters.Lw) = (Yestim(:,1:end-Parameters.Lw) .* Parameters.factors.Ystd) + Parameters.factors.Ymean;
            end
            Estimation.Regression.Y(:,:,:,s)    = reshape(Yestim, s1,s2,[]);
            
            Cov  	= reshape(Cov,size(Cov,1),size(Cov,2),s1,s2);
            for ss = 1:s1
                for sss = 1:s2
                    if any(strcmp(fieldnames(Parameters),'factors')) && normalization == 1
                        Estimation.Regression.Cov(ss,sss,:,s)   = diag(Cov(:,:,ss,sss))';
                        Estimation.Regression.Cov(ss,sss,1:end-Parameters.Lw,s) = squeeze(Estimation.Regression.Cov(ss,sss,1:end-Parameters.Lw,s))' .* Parameters.factors.Ystd.^2;
                    else
                        Estimation.Regression.Cov(ss,sss,:,s)   = diag(Cov(:,:,ss,sss))';
                    end
                end
            end
            
            Estimation.Regression.Pik = Pik;
            
            %Remove outliers
            if ~isempty(Outliers)
                if size(Estimation.Regression.Y,3) == length(Outliers)
                    for o = 1:length(Outliers)
                        tmp = Estimation.Regression.Y(:,:,o,s);
                        tmp(tmp < Outliers{o}(1)) = nan;
                        tmp(tmp > Outliers{o}(2)) = nan;                        
                        Estimation.Regression.Y(:,:,o,s) = tmp;
                    end
                end
            end
            
            %Errors computation if a reference image is provided
            if ~isempty(References)
                [Estimation.Regression.Errors.Rmse(s,:), Estimation.Regression.Errors.Nrmse(s,:), Estimation.Regression.Errors.Mae(s,:), Estimation.Regression.Errors.Nmae(s,:)] = ...
                    EvaluateEstimation(reshape(References(:,:,:,s),s1*s2,size(References,3)), reshape(Estimation.Regression.Y(:,:,:,s),s1*s2,size(References,3)));
            end
        end
        
        Estimation.Regression.quantification_time = toc;
end

end %function
        