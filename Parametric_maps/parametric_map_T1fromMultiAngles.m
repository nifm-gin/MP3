function T10map = parametric_map_T1fromMultiAngles(FA5_filename, FA15_filename, FA35_filename, add_parameters)
% generate a T1 map a multi angle scans (3 angles FA 5, 15 and 35)
% this code come from the the T1_GE_MULTIPLE_ANGLES_UVASC_MODULE coded by
% I. Tropres and adapted by B. Lemasson
tic
% no additional parameter for now
seuil_du_fit=str2double(add_parameters{:}(1));
% load the FA5 file
fid=fopen(FA5_filename ,'r');
if fid>0
    fclose(fid);
    FA5 = load(FA5_filename);
    FA5 = FA5.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the T10 map because there is\n##$ Somthing wrong with the data\n##$FA15=%s\n##$',...
        FA5_filename);
    msgbox(warning_text, 'T10 map warning') ;
    T10map = [];
    return
end
% load the FA15 file
fid=fopen(FA15_filename ,'r');
if fid>0
    fclose(fid);
    FA15 = load(FA15_filename);
    FA15 = FA15.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the T10 map because there is\n##$ Somthing wrong with the data\n##$FA15=%s\n##$',...
        FA15_filename);
    msgbox(warning_text, 'T10 map warning') ;
    T10map = [];
    return
end
% load the FA35 file
fid=fopen(FA35_filename ,'r');
if fid>0
    fclose(fid);
    FA35 = load(FA35_filename);
    FA35 = FA35.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ Somthing wrong with the data\n##$FA35=%s\n##$',...
        FA35_filename);
    msgbox(warning_text, 'T10 map warning') ;
    T10map = [];
    return
end
temp_thickness = [FA5.reco.thickness FA15.reco.thickness  FA35.reco.thickness];
temp_slice_nbr = [FA5.reco.no_slices FA15.reco.no_slices  FA35.reco.no_slices];
temp_resolution = [FA5.reco.no_samples FA15.reco.no_samples  FA35.reco.no_samples];

% check data compatibility (slice thickness and slice number)
if  length(find(temp_thickness == FA5.reco.thickness)) ~= numel(temp_thickness)
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ a slice thickness missmatch between\n##$FA5=%s\n##$ and \n##$FA15=%s\n##$ and\n##$FA35=%s\n',...
        FA5_filename,FA15_filename,FA35_filename);
    msgbox(warning_text, 'T10map warning') ;
    T10map = [];
    return
end
if length(find(temp_resolution == FA5.reco.no_samples)) ~= numel(temp_resolution)
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ a resolution missmatch between\n##$FA5=%s\n##$ and \n##$FA15=%s\n##$ and\n##$FA35=%s\n',...
        FA5_filename,FA15_filename,FA35_filename);
    msgbox(warning_text, 'T10map warning') ;
    T10map = [];
    return
end

if length(find(temp_slice_nbr == FA5.reco.no_slices)) ~= numel(temp_slice_nbr)
    warning_text = sprintf('##$ Can not calculate the T10map because there is\n##$ a resolution missmatch between\n##$FA5=%s\n##$ and \n##$FA15=%s\n##$ and\n##$FA35=%s\n',...
        FA5_filename,FA15_filename,FA35_filename);
    msgbox(warning_text, 'T10map warning') ;
    T10map = [];
    return
end

% if Realign Raw scan first? option is set to Yes :
% need to code : realign estimate then reslice (spm_reslice(V))


% --> Convert data to nii
% --> Run realign and reslice pipeline in SPM_12
% --> Convert nii resliced to uvascim struct
if strcmp(add_parameters{:}(4), 'Yes')
    %convert FA5 scan
    [PATHSTR, NAME_FA5, ~] = fileparts(FA5_filename);
    NAME_FA5 = 'data1';
    tmpFA5.uvascim.image = FA5;
    %     NHdr=CreateNiiHdr_from_uvascim(FA5_filename, tmp1, NAME_FA5);
    dataFA5_4d = permute(tmpFA5.uvascim.image.reco.data, [1 2 4 3]);
    %     WriteNii(fullfile(PATHSTR, [NAME_FA5 '.nii']), NHdr, data1_4d)
    %     clear tmp data_4d
    %convert FA15 scan
    %     [PATHSTR, NAME_FA15, ~] = fileparts(FA15_filename);
    %     NAME_FA15 = 'data2';
    tmpFA15.uvascim.image = FA15;
    %     NHdr=CreateNiiHdr_from_uvascim(FA15_filename, tmp2, NAME_FA15);
    dataFA15_4d = permute(tmpFA15.uvascim.image.reco.data, [1 2 4 3]);
    %     WriteNii(fullfile(PATHSTR, [NAME_FA15 '.nii']), NHdr, data2_4d)
    %     clear tmp data_4d
    %convert FA35 scan
    %     [PATHSTR, NAME_FA35, ~] = fileparts(FA35_filename);
    NAME_FA35 = 'data3';
    tmpFA35.uvascim.image = FA35;
    tmpFA35.uvascim.image.reco.data(:,:,2,:) =tmpFA15.uvascim.image.reco.data;
    tmpFA35.uvascim.image.reco.data(:,:,3,:) =tmpFA5.uvascim.image.reco.data;
    
    NHdr=CreateNiiHdr_from_uvascim(FA35_filename, tmpFA35, NAME_FA35);
    dataFA35_4d = permute(tmpFA35.uvascim.image.reco.data, [1 2 4 3]);
    dataFA35_4d(:,:,:,2) = dataFA15_4d;
    dataFA35_4d(:,:,:,3) = dataFA5_4d;
    
    
    WriteNii(fullfile(PATHSTR, [NAME_FA35 '.nii']), NHdr, dataFA35_4d)
    clear tmp data_4d
    
    
    % Now run the realign (estimate/reslice) form spm12
    %     P = spm_select('ExtList', pwd, 'data1.nii', 1:32);
    %     P =   char({ [NAME_FA5, '.nii,1']
    %         [NAME_FA15, '.nii,1']
    %         [NAME_FA35, '.nii,1']});
    %     spm_realign(P)
    %     spm_reslice(P)
    %
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {
        { fullfile(PATHSTR, [NAME_FA35, '.nii,1'])
        fullfile(PATHSTR, [NAME_FA35, '.nii,2'])
        fullfile(PATHSTR, [NAME_FA35, '.nii,3'])}
        }';
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 0;
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
    
    % Load resliced data
    % FA 5
    %     V =spm_vol(fullfile(PATHSTR, [NAME_FA5, '.nii']));
    %     tmp = spm_read_vols(V);
    %     tmp = permute(tmp, [1 2 4 3]);
    %     tmp = permute(tmp, [2 1 3 4]);
    %     tmp = flip(tmp, 1);
    %     tmp= flip(tmp, 2);
    %     FA5.reco.data = tmp;
    % FA 15
    %      V =spm_vol(fullfile(PATHSTR, ['r',NAME_FA15, '.nii']));
    %      tmp = spm_read_vols(V);
    %     tmp = permute(tmp, [1 2 4 3]);
    %     tmp = permute(tmp, [2 1 3 4]);
    %     tmp = flip(tmp, 1);
    %     tmp= flip(tmp, 2);
    %     FA15.reco.data = tmp;
    % FA 35
%     P = spm_select('ExtList', pwd, [NAME_FA35, '.nii'], 1:3);
%     spm_reslice(P)
    V =spm_vol(fullfile(PATHSTR, ['r', NAME_FA35, '.nii']));
    tmp = spm_read_vols(V);
    tmp = permute(tmp, [1 2 4 3]);
    tmp = permute(tmp, [2 1 3 4]);
    tmp = flip(tmp, 1);
    tmp= flip(tmp, 2);
    FA5.reco.data(:,:,1,:) = tmp(:,:,3,:);
    FA15.reco.data(:,:,1,:) = tmp(:,:,2,:);
    FA35.reco.data(:,:,1,:) = tmp(:,:,1,:);
    
    close(SPMinter)
    close(SPMgraph)
end

T10map = FA15;
flip_angles = zeros(3,1,'single');
data_to_process = zeros(FA5.reco.no_samples,FA5.reco.no_views,FA5.reco.no_slices,3);

%On poole les donnees, on recupere les angles, et on verifie les dimensions
if strcmp(FA5.texte(1:27),'# === DATA DESCRIPTION FILE')
    % Philips data
    liste = textscan(FA5.texte, '%s','delimiter', char(13));
    lgn = size(liste{1,1},1);
    parinfo = liste{1,1}{lgn-3,1}; % donne la derniere ligne des parametres
    champs = textscan(parinfo,'%s');
    flip_angles(1) = str2num(champs{1,1}{36,1}); %#ok<*ST2NM>
    %     else % Bruker data
    %         flip_angles(m)=scan_acqp('##$PVM_ExcPulseAngle=',uvascim.image(images_choisies(m)).texte,1);
    if numel(size((FA5.reco.data))) == 5
        data_to_process(:,:,:,1) = squeeze(FA5.reco.data(:,:,1,:,1));
    else
        data_to_process(:,:,:,1)  = squeeze(FA5.reco.data(:,:,1,:));
    end
end
if strcmp(FA15.texte(1:27),'# === DATA DESCRIPTION FILE')
    % Philips data
    liste = textscan(FA15.texte, '%s','delimiter', char(13));
    lgn = size(liste{1,1},1);
    parinfo = liste{1,1}{lgn-3,1}; % donne la derniere ligne des parametres
    champs = textscan(parinfo,'%s');
    flip_angles(2) = str2num(champs{1,1}{36,1}); %#ok<*ST2NM>
    %     else % Bruker data
    %         flip_angles(m)=scan_acqp('##$PVM_ExcPulseAngle=',uvascim.image(images_choisies(m)).texte,1);
    if numel(size((FA15.reco.data))) == 5
        data_to_process(:,:,:,2) = squeeze(FA15.reco.data(:,:,1,:,1));
    else
        data_to_process(:,:,:,2)  = squeeze(FA15.reco.data(:,:,1,:));
    end
end
if strcmp(FA35.texte(1:27),'# === DATA DESCRIPTION FILE')
    % Philips data
    liste = textscan(FA35.texte, '%s','delimiter', char(13));
    lgn = size(liste{1,1},1);
    parinfo = liste{1,1}{lgn-3,1}; % donne la derniere ligne des parametres
    champs = textscan(parinfo,'%s');
    flip_angles(3) = str2num(champs{1,1}{36,1}); %#ok<*ST2NM>
    %     else % Bruker data
    %         flip_angles(m)=scan_acqp('##$PVM_ExcPulseAngle=',uvascim.image(images_choisies(m)).texte,1);
    if numel(size((FA35.reco.data))) == 5
        data_to_process(:,:,:,3) = squeeze(FA35.reco.data(:,:,1,:,1));
    else
        data_to_process(:,:,:,3)  = squeeze(FA35.reco.data(:,:,1,:));
    end
end

% calculate T10map for all slices

warning('off', 'all');
data_in_vector = reshape(data_to_process, [size(data_to_process,1)*size(data_to_process,2)*size(data_to_process,3), size(data_to_process,4)]);
maxim=max(data_in_vector) * seuil_du_fit / 100;
% data_in_vector(data_in_vector(:,1) + data_in_vector(:,2) + data_in_vector(:,3) <  (mean(data_in_vector(:,1) + data_in_vector(:,2) + data_in_vector(:,3)) *seuil_du_fit),:)= nan;
fit_result = zeros([size(data_in_vector,1), 2]);

%parfor i=1:size(data_in_vector,1)
for i=1:size(data_in_vector,1)
    vector_to_process = data_in_vector(i,:);
    if  max(vector_to_process(:))>= maxim
        x = vector_to_process ./ tan(flip_angles'/180*pi) / 1.0e10;
        y = vector_to_process ./ sin(flip_angles'/180*pi) / 1.0e10;
        fit_result(i,:) = polyfit(x,y,1);
        if  fit_result(i,1) <0
            fit_result(i,1) = 0;
        end
    else
        fit_result(i,:) = [0 0];
    end
end
toc
warning('on', 'all');
fit_result=reshape(fit_result,[size(data_to_process,1),size(data_to_process,2),size(data_to_process,3), 2]);
fit_result=permute(fit_result, [1 2 4 3]);
fit_result(fit_result(:,:,1,:,1)<0.) = 0.;
T1 = -FA35.acq.tr ./ log(fit_result(:,:,1,:,1)); % en ms
T1nul = isnan(T1) | (T1<0.);
T1(T1nul) = 0.;

% T10map.reco.data = zeros([size(data_to_process,1),size(data_to_process,2),2,size(data_to_process,3)]);
T10map.reco.data = zeros([size(data_to_process,1),size(data_to_process,2),1,size(data_to_process,3)]);
% 1er echo = carte T1
T10map.reco.data(:,:,1,:,1) = T1;
T10map.reco.unit(1,1) = {'ms'};
T10map.reco.echo_label(1,1) = {'carte T10'};
% % 2eme echo = carte des pixels exclus
% T10map.reco.data(:,:,2,:,1) = T1nul;
% T10map.reco.unit(1,2) = {''};
% T10map.reco.echo_label(1,2) = {'pixels exclus'};

% T10map.reco.no_echoes = 2;
T10map.reco.no_echoes = 1;
T10map.reco.no_expts  = 1;
T10map.reco.globalmin=0;
T10map.reco.globalmax=max(T10map.reco.data(:));

T10map.reco.texte = 'T1 Fit Multi Angle';
T10map.reco.date = date;

ParamConfig=sprintf('##$QuantifMethod=T1map from Multi angle data\n##$''A*[1-exp(-t/T1)]'')\n##$##$FA5 file=%s\n##$FA15* file=%s\n##$FA35 file=%s\n\n',...
    FA5_filename,FA15_filename, FA35_filename);
T10map.reco.paramQuantif = ParamConfig;
T10map.reco=orderfields(T10map.reco);
T10map.reco = orderfields(T10map.reco);

T10map.clip=[0 5000 1];

T10map.reco.fov_offsets = T10map.reco.fov_offsets(:,1,:);
T10map.reco.fov_orientation = T10map.reco.fov_orientation(:,1,:);
T10map.reco.fov_phase_orientation = T10map.reco.fov_phase_orientation(1,:);




