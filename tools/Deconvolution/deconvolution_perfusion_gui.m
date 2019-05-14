function [CBV,CBF,MTT,TMAX,TTP,T0,CBVm,CBFm,MTTm,CBV_mask,CBF_mask,MTT_mask,volume_mask] = deconvolution_perfusion_gui(AIF,VOX,TR,TE)
% function [CBV,CBF,MTT,TMAX,TTP,T0,CBVm,CBFm,MTTm,CBV_mask,CBF_mask,MTT_mask] = deconvolution_perfusion_gui(AIF,VOX,TR,TE,handles)
% Compute hemodynamics parameters from persufion MRI by deconvolution
%
%
% INPUTS :
% AIF : AIF(t) arterial input function (raw signal)
% VOX : volume (raw signal) (4D : [Height,Width,Slices,Dynamics])
%
% TR : Repetition time
% TE : Echo time
%
% OUTPUTS :
% CBV : matrix of cerebral blood volume (in %) of the size of VOX volume
% CBF : matrix of cerebral blood flow (in %.min-1) of the size of VOX volume
% MTT : matrix of mean transit time (in s) of the size of VOX volume
% TMAX : matrix of delay of the residu of the size of VOX volume
% TTP : matrix of time to peak (in s) of the size of VOX volume
% T0 : matrix of bolus arrival time (s) of the size of VOX volume
% ###m : parameters with non computed voxels averaged with neighbours
% volume_mask : Mask of computed voxels
%
% 12/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.fr>)
% Last modified : 15/03/2013 (TP)

% TODO:
% ajouter le temps de calcul

%%% Global parameters
if ndims(VOX)==3
    VOX = permute(VOX,[1,2,4,3]);
end
[Hvox,Wvox,nb_sli,nb_dyn] = size(VOX);

% % remove neg voxels in case of the presence of negative values from the
% % realign procedure
VOX=abs(VOX);
VOX(isnan(VOX))=0;

%%% Deconvolution Parameters
zero_pad_fact = 2; % MUST be an integer !
OIth = 0.1;

%%% Compute mask for processing
% volume_mask = mask(VOX);
meansignal = mean(VOX,4);
volume_mask = meansignal>max(VOX(:))*0.01;

%%% Compute concentration and zero padding and T0 (on profite du calcul de
%%% deltaR2* qui a besoin de T0 pour le calculer maintenant)
[CVOX,T0] = DELTAR2star(VOX,TE,volume_mask);
AIF4d(1,1,1,:) = AIF;
Caif = squeeze(DELTAR2star(AIF4d,TE,true(1,1,1)));
Caif = Caif(:).'; % Force row structure
Caif = [Caif zeros(1,nb_dyn*(zero_pad_fact-1))];
CVOXpad = cat(4,CVOX,zeros(Hvox,Wvox,nb_sli,nb_dyn*(zero_pad_fact-1)));

%%% Compute SVD for Ca
Ca = toeplitz(Caif,[Caif(1) Caif(end:-1:2)]);
[U,S,V] = svd(Ca);

%%% Compute deconvolution of CVOXpad by Ca (see function below)
R = devonvolution_osvd(U,S,V,CVOXpad,nb_dyn,OIth,volume_mask);

%%% Compute hemodynamics parameters
T0 = TR.*T0;
CBV = 100.*trapz(CVOX,4)./trapz(Caif);
[RMAX,TMAX] = max(R,[],4);
TMAX(~volume_mask) = 0;
TMAX(volume_mask) = (TMAX(volume_mask)-(1/2)).*TR;
MTT = TR.*trapz(R,4)./RMAX;
MTT(isinf(MTT) | isnan(MTT)) = 0;
CBF = 60.*CBV./MTT;
CBF(isinf(CBF) | isnan(CBF)) = 0;
[~,TTP] = max(CVOX,[],4);
TTP(~volume_mask) = 0;
TTP = TTP.*TR;

%%% Remove inf voxels values
% The inf values can be due to the 'nice' mask which include noisy voxels
% which can contain zeros values on dynamics (and results in infinity
% values on concentration)
CBV(isinf(CBV) | isnan(CBV)) = 0;

%%% Compute masks of problematics voxels
CBV_mask = CBV <= 0;
CBV(CBV_mask) = 0;
MTT_mask = MTT < TR/2 | MTT > TR*nb_dyn;
MTT(MTT_mask) = 0;
CBF_mask = CBV_mask | MTT_mask;
CBF(CBF_mask) = 0;

%%% Compute mean of problematics voxels
CBVm = convn(CBV,ones(3,3,3)./26,'same');
CBVm(~volume_mask) = 0;
CBVm(~CBV_mask) = CBV(~CBV_mask);
MTTm = convn(MTT,ones(3,3,3)./26,'same');
MTTm(~volume_mask) = 0;
MTTm(~MTT_mask) = MTT(~MTT_mask);
CBFm = convn(CBF,ones(3,3,3)./26,'same');
CBFm(~volume_mask) = 0;
CBFm(~CBF_mask) = CBF(~CBF_mask);
end

function M = mask(VOX)
% function M = mask(VOX)
% Compute masks of brain volume
%
% INPUTS :
%         -
%
% OUTPUTS :
%
% 12/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.org>)
% Last modified : 12/03/2013 (TP)


[Hvox,Wvox,Svox,~] = size(VOX);
noise_lvl = SNR(VOX);
M = true(Hvox,Wvox,Svox);

%%% Voxels nuls
M(~all(VOX,4)) = false;

%%% Voxels dans le bruit
M(nanmean(VOX(:,:,:,2:6),4) < max(max(max(mean(VOX(:,:,:,2:6),4))/(noise_lvl*2)))) = false;
M(nanmean(VOX(:,end-10:end),2) < max(max(max(mean(VOX(:,end-10:end),2))/(noise_lvl*2)))) = false;

%%% Voxels qui 'bougent peu'
M(min(VOX,[],4) > (nanmean(VOX(:,:,:,2:6),4)-2*std(VOX(:,:,:,2:6),0,4))) = false;

%%% 'Smooth' to get a 'nice' mask
for sli = 1:Svox
    M(:,:,sli) = bwmorph(squeeze(M(:,:,sli)),'majority',20);
    M(:,:,sli) = imfill(squeeze(M(:,:,sli)),'holes');
end
end

function noise_level = SNR(VOX)
% function noise_level = SNR(VOX)
% Compute noise level (SNR) from MRI datas
Dvox = size(VOX,4);
M = all(VOX,4);
VOX = reshape(VOX(repmat(M,[1 1 1 Dvox])),[],Dvox);
ECT = nanstd(VOX(:,2:6),0,2);
VOX = reshape(VOX(repmat(ECT > 0,[1 Dvox])),[],Dvox);
noise_level = nanmedian((nanmean(VOX(:,2:6),2) - min(VOX,[],2))./std(VOX(:,2:6),0,2));

end

function [CONC,T0] = DELTAR2star(VOX,TE,mask_computation)
% function [CONC,T0] = deltaR2star(VOX,TE);
% Compute DELTAR2*
%
% INPUTS :
%
% OUTPUTS :
%
% 12/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.fr>)
% Last modified : 15/03/2013 (TP)

Dvox = size(VOX,4);
[T0,T0_mask] = BAT(VOX,mask_computation);
BL = sum(VOX .* T0_mask,4)./sum((VOX .* T0_mask) ~= 0,4);
BL(~mask_computation) = 0;
BL = repmat(BL,[1 1 1 Dvox]);
CONC = zeros(size(VOX));
CONC(repmat(mask_computation,[1 1 1 Dvox])) = -(1/TE).*(log(VOX(repmat(mask_computation,[1 1 1 Dvox]))./BL(repmat(mask_computation,[1 1 1 Dvox]))));
end

function R = devonvolution_osvd(U,S,V,CVOXpad,N,OIth,mask_computation)
% function R = devonvolution_osvd(U,S,V,CVOXpad,OIth,mask_computation)
% Devonvolution
%
% INPUTS :
% U,S,V : matrices from SVD of Ca matrix
% CVOXpad : zero-padded concentration voxels (total volume, 4D)
% N : Dynamics number (without zero-padding)
% OIth : oscillation threshold
% Mask_computation : mask of the volume (see mask function)
%
% OUTPUT :
% R : Residu function from deconvolution
%
% 13/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.fr>)
% Last modified : 15/03/2013 (TP)

[Hvox,Wvox,Svox,Dvox] = size(CVOXpad);
CVOXpad = CVOXpad.*repmat(mask_computation, [1 1 1 Dvox]);
data_in_vector = reshape(CVOXpad, [Hvox*Wvox*Svox Dvox]);
% R = zeros(Hvox,Wvox,Svox,N);
R = zeros(Hvox*Wvox*Svox,N);

Sv=diag(S);
% VT = numel(find(mask_computation));
% loopnb = 0;
tic
for vox = 1:size(data_in_vector,1) %  % findn(mask_computation).'
    if sum(data_in_vector(vox,:)) ~= 0 
        r = zeros(size(CVOXpad));
        
        %     if rem(loopnb,1000) == 0
        %         tstart = tic;
        %     end
        %     I = vox(1);
        %     J = vox(2);
        %     if Svox == 1
        %         sli = 1;
        %     else
        %         sli = vox(3);
        %     end
        Cvox = data_in_vector(vox,:)';
        %     Cvox = squeeze(CVOXpad(I,J,sli,:));
        
        th_prev = Dvox;
        th = 0;
        OI = Inf;
        first = true;
        while (OI > OIth || (abs(th_prev - th) > 1)) && th < Dvox
            
            %%% Calcul de l'inverse de S et filtrage %%%
            Sinv = diag([1./Sv(1:Dvox-th);zeros(th,1)]);
            
            %%% Inverse de Ca %%%
            Cainv = V*Sinv*(U');
            
            %%% Calcul du residu %%%
            r = Cainv*Cvox;
            
            %%% Calcul de OI %%%
            rOI = r(1:N);
            
            f = abs(rOI(3:N) - 2*rOI(2:N-1) + rOI(1:N-2));
            sum_oi = sum(f);
            
            if max(rOI) ~= 0
                OI = (1/N) * (1/max(rOI)) * sum_oi;
            else
                OI = Inf;
            end
            
            
            if ~first
                if OI > OIth
                    th_next = th + round(abs((th - th_prev))/2);
                else
                    th_next = th - round(abs((th - th_prev))/2);
                end
            else
                if OI > OIth
                    th_next = th + round(abs((th - th_prev))/2);
                else
                    th_next = th;
                end
            end
            
            th_prev = th;
            th = th_next;
            
            first = false;
        end
        %     R(I,J,sli,:) = r(1:N);
        R(vox,:) = r(1:N);
        
        %%% Estimating remaining computing time
        %     loopnb = loopnb + 1;
        %     if rem(loopnb,1000) == 0
        %         tend = toc(tstart);
        %         VR = VT - loopnb;
        %         ERCT = VR/1000*tend;
        %         ERCT_min = floor(ERCT/60);
        %         ERCT_sec = ERCT - 60*ERCT_min;
        %         tmpstr = sprintf('Computing slice number : %02d / %02d Estimated remaining time : %d min %02.0d sec',sli,Svox,ERCT_min,floor(ERCT_sec));
        %         set(handles.CPM_patient_selected, 'String', ['Current process:  ' tmpstr]);
        %         %         set(handles.infos,'String',tmpstr);
        %         drawnow expose
        %     end
    end
end
toc
R = reshape(R, [Hvox Wvox Svox N]);
end

function [T0,T0_mask] = BAT(VOX,mask_computation)
% function T0 = BAT(VOX,curve_type)
% Compute Bolus Arrival Time
%
% INPUT :
%
% OUTPUTS :
%
% 13/03/2013 (Thomas Perret : <thomas.perret@grenoble-inp.fr>)
% Last modified : 15/03/2013 (TP)

% Parameters of algorithm
window_size = 8;
th = 2.0;

[Hvox,Wvox,Svox,Dvox] = size(VOX);
T0_mask = false(Hvox,Wvox,Svox,Dvox);
moy = zeros(Hvox,Wvox,Svox,Dvox-window_size);
ect = zeros(Hvox,Wvox,Svox,Dvox-window_size);
for t = 1:Dvox-window_size
    moy(:,:,:,t) = mean(VOX(:,:,:,t:t+window_size),4);
    ect(:,:,:,t) = std(VOX(:,:,:,t:t+window_size),0,4);
end
Tlog = VOX(:,:,:,window_size+1:Dvox) < (moy - th.*ect);
[~,T0] = max(Tlog,[],4);
T0 = T0 + window_size - 1;
T0(~mask_computation | T0 == window_size) = 0;
[~,TTP] = min(VOX,[],4);

for vox = findn(mask_computation).'
    I = vox(1);
    J = vox(2);
    if numel(vox) == 3
        sli = vox(3);
    else
        sli = 1;
    end
    
    if T0(I,J,sli) < TTP(I,J,sli) && T0(I,J,sli) > window_size
        T0_mask(I,J,sli,2:T0(I,J,sli)-2) = true;
    else
        T0_mask(I,J,sli,2:window_size-2) = true;
    end
end
end
