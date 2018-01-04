function [aif,scores] = extraction_aif(SLICE,ROI, parameters)
% [aif,scores] = extraction_aif(SLICE,ROI)
%
% Select voxels for AIF
%
% INPUTS :
% SLICE : Slice for selecting voxels (3D : [Height,Width,Dynamics])
%
% OUTPUTS :
% aif : mean of the voxels selected
% score : score of each selected voxel and the number of warning and the
%         reasons
%
% 20/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.fr>)
% Last modified : 20/03/2013 (TP)
wmaxr = parameters.wmaxr; 
nb_vox_cand = parameters.nb_vox_cand;
nb_vox = parameters.nb_vox;


scores = cell(nb_vox,8);

[Hvox,Wvox,Dvox] = size(SLICE);
wmax = round(wmaxr * Dvox);

%%% Selection de la ROI
BINSEL = ROI;

%%% Enlever les voxels inf et Nan (pas forcement necessaire)%%%
BINSEL = BINSEL & ~any(isinf(SLICE),3) & ~any(isnan(SLICE),3);

%%% Enlever les voxel nuls %%%
nnul = all(SLICE,3);
BINSEL = BINSEL & nnul;

%%% Enlever les voxels 'bruits' %%%
MOY = mean(mean2(SLICE(repmat(BINSEL,[1 1 Dvox]))));
NOISY = mean(SLICE,3) < MOY;
BINSEL = BINSEL & ~NOISY;

%%% Calcul de la hauteur des pics %%%
BL = zeros(Hvox,Wvox);
MINVOX = zeros(Hvox,Wvox);
BL(BINSEL) = mean(reshape(SLICE(cat(3,false(Hvox,Wvox),repmat(BINSEL,[1 1 5]))),[],5),2);
MINVOX(BINSEL) = min(reshape(SLICE(repmat(BINSEL,[1 1 Dvox])),[],Dvox),[],2);
HP = BL - MINVOX;

%%% Calcul de la largeur des pics %%%
WP = Dvox*ones(Hvox,Wvox);
for i=findn(BINSEL).'
    if ~isempty(find(SLICE(i(1),i(2),:) <= (MINVOX(i(1),i(2))+HP(i(1),i(2))/2),1))
        WP(i(1),i(2)) = find(SLICE(i(1),i(2),:) <= (MINVOX(i(1),i(2))+HP(i(1),i(2))/2),1,'last') - find(SLICE(i(1),i(2),:) <= (MINVOX(i(1),i(2))+HP(i(1),i(2))/2),1);
    end
end
WP(WP == 0) = Dvox;

%%% Garder que les voxels non saturés %%%
SVOX = false(Hvox,Wvox);
% for i=findn(BINSEL).'
%     ECTVOX(i(1),i(2)) = std(SLICE(i(1),i(2),2:6),0,3);
%     SVOX(i(1),i(2)) = numel(find(SLICE(i(1),i(2),:) <= (MINVOX(i(1),i(2))+4*ECTVOX(i(1),i(2))))) > 2;
% end
for v=findn(BINSEL).'
    t1 = find(SLICE(v(1),v(2),:) <= min(SLICE(v(1),v(2),:),[],3)+4*std(SLICE(v(1),v(2),2:6),0,3),1,'first');
    t2 = find(SLICE(v(1),v(2),:) <= min(SLICE(v(1),v(2),:),[],3)+4*std(SLICE(v(1),v(2),2:6),0,3),1,'last');
    SVOX(v(1),v(2)) = any(diff(SLICE(v(1),v(2),t1:t2),2,3) < 0);
end
BINSEL = BINSEL & ~SVOX;

%%% Garder que les voxels dont la largeur est inférieur à une valeur %%%
BINSEL = BINSEL & WP < wmax & WP > 0;

TMPDATA = SLICE;
TMPDATA(~repmat(BINSEL,[1 1 Dvox])) = zeros;
TMPDATA = reshape(TMPDATA,[],Dvox);

%%% Trie des voxels, on recalcule la hauteur avant pour enlever les voxels
%%% dont la largeur est trop importante %%%
HP(~BINSEL) = 0;
TMPHP = HP(:);
[~,trie] = sort(TMPHP,'descend');

%%% Calcul du score %%%
[~,MININD] = min(SLICE,[],3);
score = zeros(nb_vox_cand,1);
for i=1:nb_vox_cand
    %%% Calcul du temps d'arrivee %%%
    t0 = BAT(TMPDATA(trie(i),:));
    initslop = HP(trie(i))/(MININD(trie(i))-t0);
    score(i) = (HP(trie(i)).*initslop) / (WP(trie(i)).*t0);
end

[~,trie_score] = sort(score,'descend');
aif = mean(TMPDATA(trie(trie_score(1:nb_vox)),:),1);

BINSEL = false(Hvox,Wvox);
for i=1:nb_vox
    BINSEL(trie(trie_score(i)))=true;
end

[I,J] = ind2sub([Hvox Wvox],trie(trie_score(1:nb_vox)));
for v=1:nb_vox
    warn = 0;
    scores(v,1:3) = {score(trie_score(v)) I(v) J(v)};
    
    %%% Test de la baseline pre-bolus
    ect_basepre = std(TMPDATA(trie(trie_score(v)),2:6),0,2);
    basepre = mean(TMPDATA(trie(trie_score(v)),2:6),2);
    if ect_basepre >= basepre/10
        warn = warn + 1;
        scores{v,4+warn} = 'La baseline pre-bolus du voxel est trop bruitee';
    end
    
    %%% Test de la baseline post-bolus
    ect_basepost = std(TMPDATA(trie(trie_score(v)),end-5:end),0,2);
    basepost = mean(TMPDATA(trie(trie_score(v)),end-5:end),2);
    if ect_basepost >= basepost/10
        warn = warn + 1;
        scores{v,4+warn} = 'La baseline post-bolus du voxel est trop bruitee';
    end
    
    %%% Test du point a t0
    t0 = BAT(TMPDATA(trie(trie_score(v)),:));
    if TMPDATA(trie(trie_score(v)),t0) >= 11*basepre/10
        warn = warn + 1;
        scores{v,4+warn} = 'La valeur a t0 du voxel est trop importante';
    end
    
    %%% Test de la longueur de la baseline pre-bolus
    if t0 < 8
        warn = warn + 1;
        scores{v,4+warn} = sprintf('La baseline du voxel est trop courte (T0 = %d dynamics)',t0);
    end
    scores{v,4} = warn;
end
end

function t0 = BAT(voxel)
% function t0 = BAT(voxel)
% Compute Bolus Arrival Time
%
% INPUT :
%
% OUTPUTS :
%
% 20/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.fr>)
% Last modified : 20/03/2013 (TP)

% Parameters of algorithm
window_size = 8;
th = 2.0;

D = numel(voxel);
moy = zeros(1,D-window_size);
ect = zeros(1,D-window_size);
for t = 1:D-window_size
    moy(t) = mean(voxel(t:t+window_size));
    ect(t) = std(voxel(t:t+window_size));
end
Tlog = voxel(window_size+1:D) < (moy - th.*ect);
[~,t0] = max(Tlog);
t0 = t0 + window_size - 1;
[~,ttp] = min(voxel);
if t0 >= ttp || all(~Tlog)
    t0 = 40;
end
end