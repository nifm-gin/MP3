function [Estimation, Parameters] = AnalyzeMRImages_MarkovField(Sequences,Dico,Method,Parameters,References,Outliers)

if nargin < 3, error('Not enought input arguments'); end
if ~exist('Method','var'),      Method = 'RegressionMRF'; end
if ~exist('Parameters','var'),  Parameters = []; end
if ~exist('References','var'),  References = []; end
if ~exist('Outliers','var'),    Outliers = []; end

% This parameters is only used to enable the parameter data normalisation
normalization = 1;

% Can be a problem, remove this line if necessary - 20/02/2019
if isempty(Parameters), Parameters = struct(); end

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
    
    case 'ClassicMRF'
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
        
    case 'RegressionMRF'               
        
        if ~any(strcmp(fieldnames(Parameters),'theta'))
            
            %Normalize trainning data
            if normalization == 1
                Parameters.factors.Ymean	= nanmean(Dico{f}.Parameters.Par,1);
                Parameters.factors.Ystd     = nanstd(Dico{f}.Parameters.Par);
                Dico{f}.Parameters.Par      = (Dico{f}.Parameters.Par - Parameters.factors.Ymean) ./ Parameters.factors.Ystd;
                
%                 Parameters.factors.Xmean	= nanmean(Dico{f}.MRSignals,1);
%                 Parameters.factors.Xstd     = nanstd(Dico{f}.MRSignals);
%                 Dico{f}.MRSignals           = (Dico{f}.MRSignals - Parameters.factors.Xmean) ./ Parameters.factors.Xstd;
                
                Parameters.factors.normalization = 1;
            else
                Parameters.factors.Ymean	= 0;
                Parameters.factors.Ystd     = 1;
                
%                 Parameters.factors.Xmean	= 0;
%                 Parameters.factors.Xstd     = 1;
                
                Parameters.factors.normalization = 0;
            end
            
            Xtrain = abs(Dico{f}.MRSignals);
            [~,Parameters] = EstimateParametersFromRegression(Xtrain, Xtrain, Dico{f}.Parameters.Par, Dico{f}.Parameters.Par, Parameters);
        else
            Dico{f}.MRSignals       = [];
            Dico{f}.Parameters.Par  = [];
        end
        
        for s = 1:slices
            
            % patch   = [0 1 0;
            %            1 1 1;
            %            0 1 0];
            patch   = [1 zeros(1,s1-2) 1 1 1 zeros(1,s1-2) 1];

            MF      = eye(s1*s2);
            MF      = conv2(MF,patch,'same');
            
            %Normalize observations
%             if any(strcmp(fieldnames(Parameters),'factors')) && normalization == 1
%                 Sequences(:,:,:,s) = reshape((reshape(Sequences(:,:,:,s),s1*s2,t) - Parameters.factors.Xmean) ./ Parameters.factors.Xstd, s1,s2,[]);
%             end
            
            %Estimation of parameters
            [Yestim,~,Cov,~,Pik] = ...
                EstimateParametersFromRegression_MarkovField(reshape(Sequences(:,:,:,s),s1*s2,t), abs(Dico{f}.MRSignals), Dico{f}.Parameters.Par, [], Parameters, MF);
            
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
        
        
end
end
        