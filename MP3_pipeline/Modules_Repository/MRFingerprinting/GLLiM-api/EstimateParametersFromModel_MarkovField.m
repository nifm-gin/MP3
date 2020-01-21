function [Yestim, Cov, Kurt, Pik] = EstimateParametersFromModel_MarkovField(X, f_estim, MF, verb)

narginchk(2, 4);
if nargin == 2, verb = 0; end


m = ones(1,size(X,1));
beta = 10;
Yestim = gllim_inverse_map_mrf(X', f_estim, beta, MF, m, 200, verb);

Kurt = [];
Cov = nan(size(Yestim',2),size(Yestim',2),size(Yestim',1));
Pik = [];

Yestim = Yestim';