function [f_estim, r, ll] = EstimateInverseFunction(Y, X, K, Lw, maxiter, cstr, verb)

narginchk(2, 7);
switch nargin
    case 2
        K = 5; Lw = 0; maxiter = 100; verb = 0; cstr.Sigma = 'i';
    case 3
        Lw = 0; maxiter = 100; verb = 0; cstr.Sigma = 'i';
    case 4 
        maxiter = 100; verb = 0; cstr.Sigma = 'i';
    case 5
        verb = 0; cstr.Sigma = 'i';
    case 6
        verb = 0;
end

addpath(genpath(fullfile(pwd, '..', 'the_GLLiM_toolbox')));

[f_estim, r, ll] = gllim(Y', X', K, 'Lw',Lw, 'maxiter',maxiter, 'verb',verb, 'cstr',cstr);