function [Estimation, Parameters] = AnalyzeMRImages(Sequences,Dico,Method,Parameters,References,Outliers,compute_ci_corr)

if nargin < 3, error('Not enought input arguments'); end
if ~exist('Method','var'),      Method = 'DBL'; end
if ~exist('Parameters','var'),  Parameters = []; end
if ~exist('References','var'),  References = []; end
if ~exist('Outliers','var'),    Outliers = []; end
if ~exist('compute_ci_corr','var'), compute_ci_corr = false; end

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
                reshape(EstimateParametersFromGrid(reshape(Sequences(:,:,:,s),s1*s2,t), abs(Dico{f}.MRSignals), Dico{f}.Parameters.Par), s1,s2, []);
                        
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
            
            Xtrain = abs(Dico{f}.MRSignals);
            [~,Parameters] = EstimateParametersFromRegression(Xtrain, Xtrain, Dico{f}.Parameters.Par, Dico{f}.Parameters.Par, Parameters);
        else
            Dico{f}.MRSignals       = [];
            Dico{f}.Parameters.Par  = [];
        end
        Estimation.Regression.learning_time = toc; 
        
        tic
        for s = 1:slices
            
            %Estimation of parameters
            [Yestim,~,Cov,~,Pik] = ...
                EstimateParametersFromRegression(reshape(Sequences(:,:,:,s),s1*s2,t), abs(Dico{f}.MRSignals), Dico{f}.Parameters.Par, [], Parameters);
            
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
            
            %CI coeff correction computation
            if compute_ci_corr == 1
                snr_test = logspace(1, 2.2, 20);
                
                idx = floor(size(Dico{1}.MRSignals,1) /11);
                idx = 1:idx:idx*11;
                
                for i = 1:length(idx)-1
                    for snr = 1:length(snr_test)
                        Xtest_noisy = AddNoise(Dico{1}.MRSignals(idx(i):idx(i+1),:), snr_test(snr));
                        Estim   = AnalyzeMRImages(Xtest_noisy, [], 'DBL', Parameters);

                        std(snr,:,i) = nanmean(squeeze(Estim.Regression.Cov),1).^.5;
                        err(snr,:,i) = EvaluateEstimation((Dico{1}.Parameters.Par(idx(i):idx(i+1),1:end-Parameters.Lw) .* Parameters.factors.Ystd) + Parameters.factors.Ymean, Estim.Regression.Y);
                    end
                end
                
                F = @(xdata,x)x(1)*exp(-x(2)./xdata) + x(3);
                ff = mean(err./std,3);
                for param = 1:size(Estim.Regression.Y,3)
                    a(param,:) = levenbergmarquardt(F, snr_test, ff(:,param,:)', [max(ff(param,:)) -1 min(ff(param,:))]);
                end
                
                Parameters.ci_correction.Func = F;
                Parameters.ci_correction.Var = a;
            end
        end
        Estimation.Regression.quantification_time = toc;
        
end
end
        