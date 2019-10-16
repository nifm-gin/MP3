function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb, flagNorm)

narginchk(3, 5)
if nargin == 4, flagNorm = 1; end
if nargin == 3, verb = 0; flagNorm = 1; end

if size(Xgrid,1) > 300000

    % normalization
    score = [];
    N = size(Xgrid,1);
    chunk = 1e5;
    I = floor(N/chunk);
    best = zeros(size(Xobs,1),1);
    idx = ones(size(Xobs,1),1);
    switch flagNorm
        case 1
            Xobs_normalized     = (1 ./ sum(Xobs.^2, 2) .^0.5) * ones(1, size(Xobs, 2))  .* Xobs;
            Xgrid_normalized    = (1 ./ sum(Xgrid.^2, 2).^0.5) * ones(1, size(Xgrid, 2)) .* Xgrid;
            % dot-product/scalar product comparison
            for i = 1:I
                Xgrid_chunk = Xgrid_normalized((i-1)*chunk +1 : i*chunk,:);
                score       = Xobs_normalized * Xgrid_chunk';
                %keep best result for this chunk
                [localBest, localIdx] = max(score, [], 2);
                idx(best<localBest) = localIdx(best<localBest)+(i-1)*chunk;
                best(best<localBest)=localBest(best<localBest);

            end
            if mod(N,chunk)~=0
                Xgrid_chunk =  Xgrid_normalized(i*chunk +1 : end,:);
                score = Xobs_normalized * Xgrid_chunk';
                %keep best result for this chunk
                [localBest, localIdx] = max(score, [], 2);
                idx(best<localBest) = localIdx(best<localBest)+(I-1)*chunk;
                best(best<localBest)=localBest(best<localBest);

            end

        case 0
            % dot-product/scalar product comparison
            for i = 1:I
                Xgrid_chunk = Xgrid((i-1)*chunk +1 : i*chunk,:);
                score       = Xobs * Xgrid_chunk';
                %keep best result for this chunk
                [localBest, localIdx] = max(score, [], 2);
                idx(best<localBest) = localIdx(best<localBest)+(i-1)*chunk;
                best(best<localBest)=localBest(best<localBest);

            end
            if mod(N,chunk)~=0
                Xgrid_chunk =  Xgrid(i*chunk +1 : end,:);
                score = Xobs * Xgrid_chunk';
                %keep best result for this chunk
                [localBest, localIdx] = max(score, [], 2);
                idx(best<localBest) = localIdx(best<localBest)+(I-1)*chunk;
                best(best<localBest)=localBest(best<localBest);

            end

    %         score       = Xobs * Xgrid';
    end
else
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
end



%[~, idx]    = max(score, [], 2);

% R2 version
% SSoo = sum( (Xobs  - mean(Xobs,  2)).^2, 2);
% SSgg = sum( (Xgrid - mean(Xgrid, 2)).^2, 2);
% SSog = (Xobs - mean(Xobs, 2)) * (Xgrid - mean(Xgrid, 2))';
% 
% R2 = SSog.^2 ./ (SSoo * SSgg');
% [~, idx]    = max(R2, [], 2);

Yestim  = Ygrid(idx, :);
Xestim  = Xgrid(idx, :);