function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb, flagNorm)

% MRF matching with norm options and chunking of the dictionary to
% avoid RAM overload


narginchk(3, 5)
if nargin == 4, flagNorm = 1; end
if nargin == 3, verb = 0; flagNorm = 1; end

if flagNorm
    Xobs = Xobs ./ vecnorm(Xobs,2,2);
    Xgrid = Xgrid ./ vecnorm(Xgrid,2,2);
else
    Xobs = Xobs;
    Xgrid = Xgrid;
end

if size(Xgrid,1) > 300000
    
    score = [];
    N = size(Xgrid,1);
    chunk = 5e4;
    I = floor(N/chunk);
    best = zeros(size(Xobs,1),1);
    idx = ones(size(Xobs,1),1);
    
    
    for i = 1:I
%         Xgrid_chunk = Xgrid((i-1)*chunk +1 : i*chunk,:);
%         scorre       = Xobs * Xgrid((i-1)*chunk +1 : i*chunk,:)';
        %keep best result for this chunk
%         [localBest, localIdx] = max(score, [], 2);
        [localBest, localIdx] = max(Xobs * Xgrid((i-1)*chunk +1 : i*chunk,:)', [], 2);
        idx(best<localBest) = localIdx(best<localBest)+(i-1)*chunk;
        best(best<localBest)=localBest(best<localBest);
        
    end
    if mod(N,chunk)~=0
%         Xgrid_chunk =  Xgrid(i*chunk +1 : end,:);
%         score = Xobs * Xgrid(i*chunk +1 : end,:)';
        %keep best result for this chunk
%         [localBest, localIdx] = max(score, [], 2);
        [localBest, localIdx] = max(Xobs*Xgrid(i*chunk +1 : end,:)', [], 2);
        idx(best<localBest) = localIdx(best<localBest)+(I)*chunk;
        best(best<localBest)=localBest(best<localBest);
        
    end
    
    
    %     switch flagNorm
    %         case 1
    % %             Xobs     = (1 ./ sum(Xobs.^2, 2) .^0.5) * ones(1, size(Xobs, 2))  .* Xobs;
    % %             Xgrid    = (1 ./ sum(Xgrid.^2, 2).^0.5) * ones(1, size(Xgrid, 2)) .* Xgrid;
    %             Xobs = Xobs ./ vecnorm(Xobs,2,2);
    %             Xgrid = Xgrid ./ vecnorm(Xgrid,2,2);
    %             % dot-product/scalar product comparison
    %             for i = 1:I
    %                 Xgrid_chunk = Xgrid((i-1)*chunk +1 : i*chunk,:);
    %                 score       = Xobs * Xgrid_chunk';
    %                 %keep best result for this chunk
    %                 [localBest, localIdx] = max(score, [], 2);
    %                 idx(best<localBest) = localIdx(best<localBest)+(i-1)*chunk;
    %                 best(best<localBest)=localBest(best<localBest);
    %
    %             end
    %             if mod(N,chunk)~=0
    %                 Xgrid_chunk =  Xgrid(i*chunk +1 : end,:);
    %                 score = Xobs * Xgrid_chunk';%     switch flagNorm
%         case 1
%             Xobs     = (1 ./ sum(Xobs.^2, 2) .^0.5) * ones(1, size(Xobs, 2))  .* Xobs;
%             Xgrid    = (1 ./ sum(Xgrid.^2, 2).^0.5) * ones(1, size(Xgrid, 2)) .* Xgrid;
%             % dot-product/scalar product comparison

    %                 %keep best result for this chunk
    %                 [localBest, localIdx] = max(score, [], 2);
    %                 idx(best<localBest) = localIdx(best<localBest)+(I-1)*chunk;
    %                 best(best<localBest)=localBest(best<localBest);
    %
    %             end
    %
    %         case 0
    %             % dot-product/scalar product comparison
    %             for i = 1:I
    %                 Xgrid_chunk = Xgrid((i-1)*chunk +1 : i*chunk,:);
    %                 score       = Xobs * Xgrid_chunk';
    %                 %keep best result for this chunk
    %                 [localBest, localIdx] = max(score, [], 2);
    %                 idx(best<localBest) = localIdx(best<localBest)+(i-1)*chunk;
    %                 best(best<localBest)=localBest(best<localBest);
    %
    %             end
    %             if mod(N,chunk)~=0
    %                 Xgrid_chunk =  Xgrid(i*chunk +1 : end,:);
    %                 score = Xobs * Xgrid_chunk';
    %                 %keep best result for this chunk
    %                 [localBest, localIdx] = max(score, [], 2);
    %                 idx(best<localBest) = localIdx(best<localBest)+(I-1)*chunk;
    %                 best(best<localBest)=localBest(best<localBest);
    %
    %             end
    %
    %     %         score       = Xobs * Xgrid';
    %     end
else
    % normalization
%     switch flagNorm
%         case 1
%             Xobs     = (1 ./ sum(Xobs.^2, 2) .^0.5) * ones(1, size(Xobs, 2))  .* Xobs;
%             Xgrid    = (1 ./ sum(Xgrid.^2, 2).^0.5) * ones(1, size(Xgrid, 2)) .* Xgrid;
%             % dot-product/scalar product comparison
%             score       = Xobs * Xgrid';
%         case 0
            % dot-product/scalar product comparison
%             score       = Xobs * Xgrid';
%     end
    
    [best, idx]    = max(Xobs * Xgrid', [], 2);
end



%[~, idx]    = max(score, [], 2);

% R2 version
% SSoo = sum( (Xobs  - mean(Xobs,  2)).^2, 2);
% SSgg = sum( (Xgrid - mean(Xgrid, 2)).^2, 2);
% SSog = (Xobs - mean(Xobs, 2)) * (Xgrid - mean(Xgrid, 2))';
%
% R2 = SSog.^2 ./ (SSoo * SSgg');
% [~, idx]    = max(R2, [], 2);

score = best;
Yestim  = Ygrid(idx, :);
Xestim  = Xgrid(idx, :);