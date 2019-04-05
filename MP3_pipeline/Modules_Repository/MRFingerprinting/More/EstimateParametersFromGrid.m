function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb)

narginchk(3, 4)
if nargin == 3, verb = 0; end

% normalization
Xobs_normalized     = (1 ./ sum(Xobs.^2, 2) .^0.5) * ones(1, size(Xobs, 2))  .* Xobs;
Xgrid_normalized    = (1 ./ sum(Xgrid.^2, 2).^0.5) * ones(1, size(Xgrid, 2)) .* Xgrid;

% dot-product/scalar product comparison
score       = Xobs_normalized * Xgrid_normalized';
[~, idx]    = max(score, [], 2);

% R2 version
% SSoo = sum( (Xobs  - mean(Xobs,  2)).^2, 2);
% SSgg = sum( (Xgrid - mean(Xgrid, 2)).^2, 2);
% SSog = (Xobs - mean(Xobs, 2)) * (Xgrid - mean(Xgrid, 2))';
% 
% R2 = SSog.^2 ./ (SSoo * SSgg');
% [~, idx]    = max(R2, [], 2);

Yestim  = Ygrid(idx, :);
Xestim  = Xgrid(idx, :);