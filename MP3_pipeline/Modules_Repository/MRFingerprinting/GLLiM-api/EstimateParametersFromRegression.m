function [Yestim, Parameters, Cov, Kurt, Pik] = EstimateParametersFromRegression(Xtest, Xtrain, Ytrain, Ytest, Parameters)

if ~exist('Parameters','var')
    Parameters.maxiter     = 200;
    Parameters.cstr.Sigma  = 'd*';
    Parameters.cstr.Gammat = ''; 
    Parameters.cstr.Gammaw = '';
    
elseif isempty(Parameters)
    Parameters.maxiter     = 200;
    Parameters.cstr.Sigma  = 'd*';
    Parameters.cstr.Gammat = ''; 
    Parameters.cstr.Gammaw = '';
    
else
    if ~any(contains(fields(Parameters), 'maxiter')),   Parameters.maxiter = 200; end
    if ~any(contains(fields(Parameters), 'cstr')),      Parameters.cstr.Sigma  = 'd*';
                                                        Parameters.cstr.Gammat = '';
                                                        Parameters.cstr.Gammaw = ''; end
end
if ~exist('Ytest','var'), Ytest = []; end


% Find Lw using BIC
if ~any(contains(fields(Parameters), 'Lw'))
    Lw       = [0:10 12:2:16 20];
    if ~any(contains(fields(Parameters), 'K'))
        Parameters.Lw = FindOptimalLw(Xtrain, Ytrain, 30, Lw, Parameters.maxiter, Parameters.cstr, 0);
    else
        Parameters.Lw = FindOptimalK(Xtrain, Ytrain, Parameters.K, Lw, Parameters.maxiter, Parameters.cstr, 0);
    end
end
% Find K using BIC
if ~any(contains(fields(Parameters), 'K'))
    K       = 2:2:50;
    Parameters.K = FindOptimalK(Xtrain, Ytrain, K, Parameters.Lw, Parameters.maxiter, Parameters.cstr, 0);
end


if ~any(strcmp(fieldnames(Parameters),'theta'))
    % Estimations
    Parameters.theta = EstimateInverseFunction(Ytrain, Xtrain, Parameters.K, Parameters.Lw, Parameters.maxiter, Parameters.cstr, 0);
end

if isempty(Xtest)
    Yestim = [];
else
    [Yestim, Cov, Kurt, Pik] = EstimateParametersFromModel(Xtest, Parameters.theta, 0);
end

