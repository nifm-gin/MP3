% function [Maps, Err, T1t, Ctmap, Ctfit] = parametric_map_DCEpermeability(DCE, T10map, coord_AIF, dynsctime, TR)
function [vp_struc, kep_struc, Ktrans_struc, Ctfit_struc, Ctmap_struc, T1t_struc, Err_struc] = parametric_map_DCEpermeability(filename_DCE, filename_T1map,filename_DCEtardif, add_parameters)
% From Marine Beaumont's code
% Modified by Irene Tropres 2013-03-04
% And by B. Lemasson 2013-10-01

% DCE = [nx, ny, nslices, necho, ndyn]
% T10map = [nx, ny, nslices]
% coord_AIF = coordonnnes de l'AIF dns un tableau [x,y,no_cp]
% dynamic scan time (s)
% TR (s)
%
% renvoie 2 tableaux de donnees :
%   - les donnees de vp, kep, Ktrans dans Maps
%   - les donnees en erreur dans Err
%


%=========================================================================
tic

% load the DCE file
fid=fopen(filename_DCE ,'r');
if fid>0
    fclose(fid);
    DCEstruct = load(filename_DCE);
    DCE = single(squeeze(DCEstruct.uvascim.image.reco.data));
else
    warning_text = sprintf('##$ Can not calculate the DCE permeability maps because there is\n##$ Somthing wrong with the data\n##$DCE file=%s\n##$',...
        filename_DCE);
    msgbox(warning_text, 'DCE permeability maps warning') ;
    vp_struc = []; kep_struc = []; Ktrans_struc =[]; Ctfit_struc = []; Ctmap_struc  = []; T1t_struc  = [];Err_struc  = [];
    return
end

% load the T1map file
fid=fopen(filename_T1map ,'r');
if fid>0
    fclose(fid);
    T1mapstruct = load(filename_T1map);
    T1map = single(squeeze(T1mapstruct.uvascim.image.reco.data(:,:,1,:)));
    
else
    warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ Somthing wrong with the data\n##$T1map=%s\n##$',...
        filename_CBF);
    msgbox(warning_text, 'DCE permeability maps warning') ;
    vp_struc = []; kep_struc = []; Ktrans_struc =[]; Ctfit_struc = []; Ctmap_struc  = []; T1t_struc  = [];Err_struc  = [];
    return
end


% load the DCE tardif if selected
if ~strcmp(filename_DCEtardif, 'None')
    DCEtardif_option =1;
    fid=fopen(filename_DCE ,'r');
    if fid>0
        fclose(fid);
        DCErardifstruct = load(filename_DCEtardif);
        DCEtardif = single(squeeze(DCErardifstruct.uvascim.image.reco.data));
    else
        warning_text = sprintf('##$ Can not calculate the DCE permeability maps because there is\n##$ Somthing wrong with the data\n##$DCE file=%s\n##$',...
            filename_DCE);
        msgbox(warning_text, 'DCE permeability maps warning') ;
        vp_struc = []; kep_struc = []; Ktrans_struc =[]; Ctfit_struc = []; Ctmap_struc  = []; T1t_struc  = [];Err_struc  = [];
        return
    end
else
    DCEtardif_option =0;
end
% % remove neg voxels in case of the presence of negative values from the
% % realign procedure
DCE=abs(DCE);
DCE(isnan(DCE))=0;
DCEtardif=abs(DCEtardif);
DCE(isnan(DCEtardif))=0;

% if Realign Raw scan first? option is set to Yes :
% --> Convert data to nii
% --> Run realign and reslice pipeline in SPM_12
% --> Convert nii resliced to uvascim struct
if strcmp(add_parameters{:}(5), 'Yes')
    %    moved = load(Scan_to_realign_filename);
    [PATHSTR, NAME_scan_to_realign, ~] = fileparts(filename_DCE);
    NAME_scan_to_realign = [NAME_scan_to_realign '-for_spm'];
    
    
    NHdr=CreateNiiHdr_from_uvascim(filename_DCE, DCEstruct, NAME_scan_to_realign);
    data_4d = squeeze(DCEstruct.uvascim.image.reco.data(:,:,1,:,:));
    WriteNiiIm(fullfile(PATHSTR, [NAME_scan_to_realign '.nii']),NHdr, data_4d)
    
    for xx = 1:size(data_4d,4)
        if ~exist('Scan_to_realign_nii_filename', 'var')
            Scan_to_realign_nii_filename{1,:} =  fullfile(PATHSTR,[NAME_scan_to_realign '.nii,1']);
        else
            Scan_to_realign_nii_filename{size(Scan_to_realign_nii_filename,1)+1,:} =  fullfile(PATHSTR, [NAME_scan_to_realign '.nii,' num2str(xx)]); %#ok<AGROW>
        end
    end
     matlabbatch{1}.spm.spatial.realign.estwrite.data  = {Scan_to_realign_nii_filename};
    %%
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    jobs = repmat(matlabbatch, 1, 1);
    inputs = cell(0, 1);
    for crun = 1:1
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
    clear matlabbatch 
    
    date_str = date;
    if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
        movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
            [PATHSTR, filesep, NAME_scan_to_realign(1:end-8) '_SPM_realign.ps']);
    end
    
        % Load resliced data
        V =spm_vol(fullfile(PATHSTR, ['r', NAME_scan_to_realign, '.nii']));
        tmp = spm_read_vols(V);
        tmp = permute(tmp, [1 2 5 3 4]);
        tmp = permute(tmp, [2 1 3 4 5]);
        tmp = flip(tmp, 1);
        tmp = flip(tmp, 2);
        DCEstruct.uvascim.image.reco.data =tmp;
end


% figure; plot(squeeze(DCE(99,92,19,:)))
% Define constants
% if strcmp(DCE.texte(1:27),'# === DATA DESCRIPTION FILE')
%     dynsctime = strfind(DCE.texte,'Scan Duration [sec]');
% end
TR = DCEstruct.uvascim.image.acq.tr;
dynsctime = TR;


% initialisations
r1 = str2double(add_parameters{:}(1));
% 3.5 L/(mmol.s) at 3T in serum pour le Gd-Dota (ref : Sieber, Pharmaceutical
% and safety aspects of gadolinium-based contrast agents, EJHP
% Practice, vol 15, 2009/6
% idem Peter Reimer et al, Clinical MR Imaging, Third edition, p39
AIF_x = str2double(add_parameters{:}(2));
AIF_y = str2double(add_parameters{:}(3));
AIF_z = str2double(add_parameters{:}(4));

% hard coded for test
coord_AIF = [AIF_x AIF_y AIF_z];

if DCEtardif_option == 1 
  DCE = cat(4,DCE,DCEtardif);
end

nx = size(DCE,1);
ny = size(DCE,2);
nslices = size(DCE,3);
ndyn = size(DCE,4);
if DCEtardif_option == 1 
    tps = (0:dynsctime:(ndyn-1-size(DCEtardif,4))*dynsctime);
    tps = [tps tps(end)+str2double(add_parameters{:}(7)):dynsctime:tps(end)+str2double(add_parameters{:}(7))+size(DCEtardif,4)*dynsctime-1];
    tps = tps./ 60;
else
    
    tps = (0:dynsctime:(ndyn-1)*dynsctime) ./ 60; % en min
end
Maps = zeros(nx*ny,3,nslices,'single');
Err = zeros(nx*ny,3,nslices,'single');

% initialisation des parametres du fit
vp = 3;
Ktrans = 0.03;
kep = 0.1;


%% cartes T1(t) -> Ct(t)
%---------------------------------------------------------------------
% T1(t) = -TR / (ln(1-S(t)/S0))
% ou S0 = Sbaseline/(1-exp(-TR/T10))
% !!! faire une recherche sur le tps de l'arrivee du pdc ? !!!
Sbase    = nanmean(DCE(:,:,:,1:5),4); % jusqu'a quel point ? besoin de mettre abs(mean) ?
S0map    = Sbase ./ (1-exp(-TR./T1map)); % ok

S0map    = reshape(S0map,nx*ny*nslices,1);
T1map   = reshape(T1map,nx*ny*nslices,1);
DCE = reshape(DCE,nx*ny*nslices,ndyn,1);
T1t      = zeros(nx*ny*nslices,ndyn,'single');
Ctmap    = zeros(nx*ny*nslices,ndyn,'single');

for i=1:ndyn
    % T1(t) pour le tissu
    tmp = 1-(DCE(:,i)./S0map);
    %     tmp(tmp<0) = 0;
    %     tmp(~isfinite(tmp)) = 0;
    T1t(:,i) = -TR ./ log(tmp); % en s
    
    % Ct(t)
    Ctmap(:,i) = 1/r1 * (1./T1t(:,i) - 1./T1map); % en mmol/L
end
clear('tmp');

Ctmap = reshape(real(Ctmap),nx,ny,nslices,ndyn);
Cp    = squeeze(Ctmap(coord_AIF(1),coord_AIF(2),coord_AIF(3),:));

% figure; plot(Cp);
%% ajustement gamma des donnees -> vp, kep, Ktrans
% Ct(t) = Vp.Cp(t) + Ktrans.somme de 0 a t Cp(tau).exp(-(t-tau).kep) dtau
%---------------------------------------------------------------------
% calcul sur les pixels non-nuls de la carte T10 et non NaN de Ct
% ind = ~isfinite(Ct);
% Ct(ind) = 0;

% passage a 3D (versus 4D) pour ne pas trop charger la memoire
Ctmap  = reshape(Ctmap,nx*ny*nslices,ndyn);
Ctfit  = zeros(nx*ny*nslices,ndyn,'single');
Maps = zeros(nx*ny*nslices,3,'single');
Err = zeros(nx*ny*nslices,3,'single');
Err(:,:) = -1;
% initialize the structure of optimization parameters. We'll carry this
% around to exchange information across the recursive function calls
B_init = [vp*0.01 Ktrans kep]; % en % et min-1
b_num  = numel(B_init);
fit_fcn = @(xdat,BETA)parametric_map_gado_tissu(xdat,BETA,Cp);
pars       = struct(...
    'N',     ndyn,...                   % number of Y data points
    'b_num', b_num,...                  % number of parameters of the model
    'A',     zeros(b_num,'single'),...  % variance covariance matrix of partial derivatives
    'G',     zeros(b_num,1,'single'),...% gradient vector to update paramters
    'D',     zeros(b_num,1,'single'),...% vector of parameter std's
    'q',     [],...                     % ??? (some sort of normalized covariance matrix)
    'B',     B_init(:),...              % current parameter guess
    'INC',   1+eye(b_num)*0.001,...     % perturbation of B to estimate partial derivatives
    'L',     0.1*eye(b_num),...         % lambda
    'dL',    10,...                     % factorial increment or decrement of lambda
    'CONV',  0.001,...                   % convergence criterion
    'T1',    30,...                     % limit of number of iterations (and output status)
    'T2',    0,...                      % iteration counter
    'F',     fit_fcn,...            % model function
    'E',     [],...                     % error vector for current parameter estimate
    'SO',    [],...                     % sum of squares error of current parameter guess
    'verbose',0);                       % flag indicating whether progress messages should be printed

% h = waitbar(0,'1','Name','Computing Permeability maps...');

% npix = length(find(T10map>0));
% ppix = 0;

DCE = reshape(DCE,nx*ny*nslices,ndyn,1);

% ind = find(DCE(:,1)>max(DCE(:,1)) * 0.05 & T1map(:,1)>0);

threshold = max(DCE(:,1))* 0.05 ;
parfor j=1:size(Ctmap,1)
    if max(DCE(j,:))>threshold && T1map(j,1)>0
        tmp = sort(Ctmap(j,:));
        mintmp = mean(tmp(1:5));
        maxtmp = mean(tmp(end-5:end));
        if sum(isfinite(Ctmap(j,:))) == ndyn && (maxtmp-mintmp)/abs(mintmp)>3
            [BETA, D, ~] = levenbergmarquardt(pars,tps,reshape(Ctmap(j,:),ndyn,1));
            % passer par le levmar.c quand on aura une fct gado_tissu en .c (avec eq Ct(t) )
            %             BETA(1) = abs(BETA(1));
           
            % Ctfit is not use so far, so I comment it to save time BL
            %         Ctfit(ind(j),:) = parametric_map_gado_tissu(tps,BETA); % calcul des valeurs de la fonction fittee
            
            Maps(j,:) = BETA;
            Err(j,:)  = D;
            % D : error on the parameter estimates.
            % Numerical array of the same size as B_init
        end
    end
end

Maps = reshape(Maps,nx,ny,nslices,3);
Err = reshape(Err,nx,ny,nslices,3);
T1t = reshape(T1t,nx,ny,1,nslices,ndyn);
Ctmap = reshape(Ctmap,nx,ny,1,nslices,ndyn);
Ctfit = reshape(Ctfit,nx,ny,1,nslices,ndyn);


% save imformationsize(ttt
tempstruct_for_Maps_and_Err=T1mapstruct.uvascim.image; %3D image
tempstruct_for_T1t_Ctmap_and_Ctfit=DCEstruct.uvascim.image; % 4D image


% adjust stuctures
tempstruct_for_Maps_and_Err.reco.no_echoes = 1;
tempstruct_for_Maps_and_Err.texte = 'parametric_map_DCEpermeability.m function';

tempstruct_for_T1t_Ctmap_and_Ctfit.texte = 'parametric_map_DCEpermeability.m function';


ParamConfig=sprintf('##$QuantifMethod=Permeability maps from a DCE scan\n##$Relaxivity of the CA=%s\n##$DCE scan used=\n%s\n##$T1map scan used=\n%s\n##$AIF coordonates=\n%s %s %s\n',...
    add_parameters{:}{1},...
    filename_DCE,...
    filename_T1map,...
    add_parameters{:}{2},   add_parameters{:}{3},   add_parameters{:}{4});

tempstruct_for_Maps_and_Err.reco.paramQuantif = ParamConfig;
tempstruct_for_T1t_Ctmap_and_Ctfit.paramQuantif = ParamConfig;

%remove old data form the structure
tempstruct_for_Maps_and_Err.reco = rmfield(tempstruct_for_Maps_and_Err.reco, 'data');
tempstruct_for_T1t_Ctmap_and_Ctfit.reco= rmfield(tempstruct_for_T1t_Ctmap_and_Ctfit.reco, 'data');

%
vp_struc = tempstruct_for_Maps_and_Err;
kep_struc = tempstruct_for_Maps_and_Err;
Ktrans_struc = tempstruct_for_Maps_and_Err;
Err_struc = tempstruct_for_Maps_and_Err;

Ctfit_struc = tempstruct_for_T1t_Ctmap_and_Ctfit;
Ctmap_struc = tempstruct_for_T1t_Ctmap_and_Ctfit;
T1t_struc = tempstruct_for_T1t_Ctmap_and_Ctfit;


% vp
vp_struc.reco.data = Maps(:,:,:,1);
vp_struc.reco.data= permute(vp_struc.reco.data, [1 2 4 3]);
vp_struc.reco.unit = {'%'};
vp_struc.reco.echo_label = {'vp'};
% Ktrans_struc
Ktrans_struc.reco.data = Maps(:,:,:,2);
Ktrans_struc.reco.data= permute(Ktrans_struc.reco.data, [1 2 4 3]);
Ktrans_struc.reco.unit = {'min-1'};
Ktrans_struc.reco.echo_label = {'Ktrans'};
% kep
kep_struc.reco.data =  Maps(:,:,:,3);
kep_struc.reco.data= permute(kep_struc.reco.data, [1 2 4 3]);
kep_struc.reco.unit = {'min-1'};
kep_struc.reco.echo_label = {'Kep'};
% Ctfit_struc
Ctfit_struc.reco.data = Ctfit;
Ctfit_struc.reco.unit = {'mmol/l'};
Ctfit_struc.reco.echo_label = {'C(t) fittee'};
% Ctmap_struc
Ctmap_struc.reco.data = Ctmap;
Ctmap_struc.reco.unit = {'mmol/l'};
Ctmap_struc.reco.echo_label = {'Ct(t)'};
% T1t_struc
T1t_struc.reco.data = T1t;
T1t_struc.reco.unit = {'s'};
T1t_struc.reco.echo_label = {'T1(t)'};
% Err_struc
Err_struc.reco.data = Err;
Err_struc.reco.data= permute(Err_struc.reco.data, [1 2 4 3]);
Err_struc.reco.unit(1:1:3) = {''};
Err_struc.reco.echo_label = {'vp error' 'Ktrans error' 'kep error'};
toc
end

