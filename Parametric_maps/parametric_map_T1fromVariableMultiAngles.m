function T10map = parametric_map_T1fromVariableMultiAngles(FA_filenames, add_parameters)
% generate a T1 map a multi angle scans 
% this code come from the the T1_GE_MULTIPLE_ANGLES_UVASC_MODULE coded by
% I. Tropres and adapted by B. Lemasson
tic
% no additional parameter for now
seuil_du_fit=str2double(add_parameters{:}(1));

% load the FA files
NbAngles = size(FA_filenames,2);
for angles=1:NbAngles
    fid=fopen(FA_filenames{angles}{:} ,'r');
    if fid>0
        fclose(fid);
        FAuvasc = load(FA_filenames{angles}{:});
        FA(angles) = FAuvasc.uvascim.image;
    else
        warning_text = sprintf('##$ Can not calculate the T10 map because there is\n##$ something wrong with the data \n##$ Could not open the files');
        msgbox(warning_text, 'T10 map warning') ;
        T10map = [];
        return
    end
    temp_thickness(angles) = FA(angles).reco.thickness;
    temp_slice_nbr(angles) = FA(angles).reco.no_slices;
    temp_resolution(angles) = FA(angles).reco.no_samples;
end


% check data compatibility (slice thickness and slice number)
if  length(find(temp_thickness == FA(1).reco.thickness)) ~= numel(temp_thickness)
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ a slice thickness missmatch between the scans');
    msgbox(warning_text, 'T10map warning') ;
    T10map = [];
    return
end
if length(find(temp_resolution == FA(1).reco.no_samples)) ~= numel(temp_resolution)
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ a resolution mismatch between the scans');
    msgbox(warning_text, 'T10map warning') ;
    T10map = [];
    return
end

if length(find(temp_slice_nbr == FA(1).reco.no_slices)) ~= numel(temp_slice_nbr)
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ a slice number missmatch between the scans');
    msgbox(warning_text, 'T10map warning') ;
    T10map = [];
    return
end

%%
T10map = FA(1);
flip_angles = zeros(NbAngles,1,'single');
data_to_process = zeros(FA(1).reco.no_samples,FA(1).reco.no_views,FA(1).reco.no_slices,NbAngles);

%% Recuperation de la valeur des angles et mise en forme des donnees a fitter
%On poole les donnees, on recupere les angles, et on verifie les dimensions
for angles=1:NbAngles
    if strcmp(FA(angles).texte(1:27),'# === DATA DESCRIPTION FILE') % PHILIPS DATA
        liste = textscan(FA(angles).texte, '%s','delimiter', char(13));
        lgn = size(liste{1,1},1);
        parinfo = liste{1,1}{lgn-3,1}; % donne la derniere ligne des parametres
        champs = textscan(parinfo,'%s');
        flip_angles(angles) = str2num(champs{1,1}{36,1}); %#ok<*ST2NM>
    else % BRUKER DATA
        infosImpulsionExcitation = strsplit(scan_acqp('##$ExcPulse1=(',FA(angles).texte,0));
        flip_angles(angles) = str2double(infosImpulsionExcitation(3));
    end
    % Mise en forme donnees a fitter
    if numel(size((FA(angles).reco.data))) == 5
        data_to_process(:,:,:,angles) = squeeze(FA(angles).reco.data(:,:,1,:,1));
    else
        data_to_process(:,:,:,angles)= squeeze(FA(angles).reco.data(:,:,1,:));
    end
end


%% calculate T1map for all slices

%% if Realign Raw scan first? option is set to Yes :
% --> Convert data to nii to realign --- CHECK FOR MULTI ANGLE a faire !
% --> Run realign and reslice pipeline in SPM_12
% --> Convert nii resliced to uvascim struct
if strcmp(add_parameters{:}(4), 'Yes')
    Nifti_name = 'MutliAngles_RawData';
    tmp_structure.uvascim.image = FA(1);
    tmp_structure.uvascim.image.reco.data = data_to_process;
    [PATHSTR, ~, ~] = fileparts(FA_filenames{angles}{:});
    % create a Nifti file
    NHdr=CreateNiiHdr_from_uvascim(FA_filenames{1}{:}, tmp_structure, Nifti_name);
    WriteNii(fullfile(PATHSTR, [Nifti_name '.nii']), NHdr, data_to_process)
    clear tmp_structure
     
    matlabbatch{1}.spm.spatial.realign.estwrite.data = ...
            { fullfile(PATHSTR, [Nifti_name, '.nii,' num2str(1)])};
    for angles=2:NbAngles
        matlabbatch{1}.spm.spatial.realign.estwrite.data =  [matlabbatch{1}.spm.spatial.realign.estwrite.data' ...
            { fullfile(PATHSTR, [Nifti_name, '.nii,' num2str(angles)])}]';
    end
     matlabbatch{1}.spm.spatial.realign.estwrite.data = {matlabbatch{1}.spm.spatial.realign.estwrite.data};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    [SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
    jobs = repmat(matlabbatch, 1, 1);
    inputs = cell(0, 1);
    for crun = 1:1
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
    clear matlabbatch
   
    close(SPMinter)
    close(SPMgraph)
    
    % Load the realign data
    V =spm_vol(fullfile(PATHSTR, ['r', Nifti_name, '.nii']));
   
    tmp = spm_read_vols(V);
    tmp = permute(tmp, [2 1 3 4]);
    tmp = flip(tmp, 1);
    tmp= flip(tmp, 2);
    for angles=1:NbAngles
        data_to_process(:,:,:,angles) = tmp(:,:,:,angles);
    end
end


%% calculate T10map for all slices
data_in_vector = reshape(data_to_process, [size(data_to_process,1)*size(data_to_process,2)*size(data_to_process,3), size(data_to_process,4)]);
maxim=max(data_in_vector) * seuil_du_fit / 100;
fit_result = zeros([size(data_in_vector,1), 2]);
erreur = zeros([size(data_in_vector,1), 1]);

parfor i=1:size(data_in_vector,1)
    vector_to_process = data_in_vector(i,:);
    if  max(vector_to_process(:))>= maxim
        x = vector_to_process ./ tan(flip_angles'/180*pi)/ mean(vector_to_process);
        y = vector_to_process ./ sin(flip_angles'/180*pi)/ mean(vector_to_process);
        %p = polyfit(x,y,1);
        X = [ones(length(x'),1) x'];
        p = X\y';
        % Calcul de R² pour estimer la qualite du fit
        yfit = p(2)*x + p(1);
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        rsq = 1 - SSresid/SStotal;
        
        fit_result(i,1) = p(2); % Pente de la courbe
        erreur(i) = rsq*100;    % R²*100
    end
end
toc
warning('on', 'all');
fit_result(:,2) = erreur;
fit_result=reshape(fit_result,[size(data_to_process,1),size(data_to_process,2),size(data_to_process,3), 2]);
fit_result=permute(fit_result, [1 2 4 3]);
fit_result(fit_result(:,:,1,:,1)<0.) = 0.;
T1 = -FA(1).acq.tr ./ log(fit_result(:,:,1,:,1)); % en ms

% T10map.reco.data = zeros([size(data_to_process,1),size(data_to_process,2),2,size(data_to_process,3)]);
T10map.reco.data = zeros([size(data_to_process,1),size(data_to_process,2),1,size(data_to_process,3)]);

% 1er echo = carte T1
T10map.reco.data(:,:,1,:,1) = T1;
T10map.reco.unit(1,1) = {'ms'};
T10map.reco.echo_label(1,1) = {'carte T10'};
% % 2eme echo = valeur de R²
% T10map.reco.data(:,:,2,:,1) = fit_result(:,:,2,:,1);
% T10map.reco.unit(1,2) = {''};
% T10map.reco.echo_label(1,2) = {'Coefficient de determination'};

% T10map.reco.no_echoes = 2;
T10map.reco.no_echoes = 1;
T10map.reco.no_expts  = 1;
T10map.reco.globalmin=0;
T10map.reco.globalmax=max(T10map.reco.data(:));

T10map.reco.texte = 'T1 Fit Multi Angle';
T10map.reco.date = date;
files = [];
for angles = 1:NbAngles
    files = horzcat(files, sprintf('##$FA%d° file%d=%s\n',flip_angles(angles), angles,FA_filenames{angles}{:}));
end
ParamConfig=sprintf('##$QuantifMethod=T1map from %d angle data\n##$''A*[1-exp(-t/T1)]''\n%s\n',...
    NbAngles,files);
T10map.reco.paramQuantif = ParamConfig;
T10map.reco=orderfields(T10map.reco);
T10map.reco = orderfields(T10map.reco);

T10map.clip=[0 5000 1];

T10map.reco.fov_offsets = T10map.reco.fov_offsets(:,1,:);
T10map.reco.fov_orientation = T10map.reco.fov_orientation(:,1,:);
T10map.reco.fov_phase_orientation = T10map.reco.fov_phase_orientation(1,:);




