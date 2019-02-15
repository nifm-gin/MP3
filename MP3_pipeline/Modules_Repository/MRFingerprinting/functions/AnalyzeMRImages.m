function [Estimation, Parameters] = AnalyzeMRImages(Sequences,Dico,Method,Parameters,References,Outliers)

if nargin < 3, error('Not enought input arguments'); end
if ~exist('Method','var'),      Method = 'RegressionMRF'; end
if ~exist('References','var'),  References = []; end
if ~exist('Outliers','var'),    Outliers = []; end
if ~exist('Parameters','var'),  Parameters = struct();
elseif isempty(Parameters),     Parameters = struct(); end


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
            [~,~,tmp] = ...
                EstimateParametersFromGrid(reshape(Sequences(:,:,:,s),s1*s2,t), abs(Dico{f}.MRSignals), Dico{f}.Parameters.Par);
            Estimation.GridSearch.BestMatch(:,:,:,s) = reshape(tmp, s1,s2, []);
            %Errors computation if a reference image is provided
            if ~isempty(References)
                [~,Estimation.GridSearch.Errors(s,:)] = EvaluateEstimation(reshape(References(:,:,:,s),s1*s2,size(References,3)), reshape(Estimation.GridSearch.Y(:,:,:,s),s1*s2,size(References,3)));
            end
        end
        Parameters = [];
        
    case 'RegressionMRF'
                
        if ~any(strcmp(fieldnames(Parameters),'theta'))
            Xtrain = abs(Dico{f}.MRSignals);
            [~,Parameters] = EstimateParametersFromRegression(Xtrain, Xtrain, Dico{f}.Parameters.Par, Dico{f}.Parameters.Par, Parameters);
        else
            Dico{f}.MRSignals = [];
            Dico{f}.Parameters.Par = [];
        end
        for s = 1:slices
            %Estimation of parameters
            [Yestim,~,Cov] = ...
                EstimateParametersFromRegression(reshape(Sequences(:,:,:,s),s1*s2,t), abs(Dico{f}.MRSignals), Dico{f}.Parameters.Par, [], Parameters);
            Estimation.Regression.Y(:,:,:,s)    = reshape(Yestim, s1,s2,[]);
            Cov         = reshape(Cov,size(Cov,1),size(Cov,2),s1,s2);
            for ss = 1:s1
                for sss = 1:s2
                    Estimation.Regression.Cov(ss,sss,:,s) = diag(Cov(:,:,ss,sss))';
                end
            end
            
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
                [~,Estimation.Regression.Errors(s,:)] = EvaluateEstimation(reshape(References(:,:,:,s),s1*s2,size(References,3)), reshape(Estimation.Regression.Y(:,:,:,s),s1*s2,size(References,3)));
            end
        end
        
        
end
end
        