function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb, flagNorm)

narginchk(3, 5)
if nargin == 4, flagNorm = 1; end
if nargin == 3, verb = 0; flagNorm = 1; end

% normalization
switch flagNorm
    case 1
        Xobs_normalized     = (1 ./ sum(Xobs.^2, 2) .^0.5) * ones(1, size(Xobs, 2))  .* Xobs;
        Xgrid_normalized    = (1 ./ sum(Xgrid.^2, 2).^0.5) * ones(1, size(Xgrid, 2)) .* Xgrid;
        % dot-product/scalar product comparison
        score       = Xobs_normalized * Xgrid_normalized';
    case 0
        % dot-product/scalar product comparison
        score       = Xobs * Xgrid';
end

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