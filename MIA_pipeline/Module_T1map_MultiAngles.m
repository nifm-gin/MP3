function [files_in,files_out,opt] = Module_T1map_MultiAngles(files_in,files_out,opt)
%
% SYNTAX:
% [IN,OUT,OPT] = PSOM_TEMPLATE_BRICK(IN,OUT,OPT)
%
% _________________________________________________________________________
% INPUTS:
%
% IN        
%   (string) a file name of a 3D+t fMRI dataset .
%
% OUT
%   (structure) with the following fields:
%       flag_test
%   CORRECTED_DATA
%       (string, default <BASE NAME FMRI>_c.<EXT>) File name for processed 
%       data.
%       If OUT is an empty string, the name of the outputs will be 
%       the same as the inputs, with a '_c' suffix added at the end.
%
%   MASK
%       (string, default <BASE NAME FMRI>_mask.<EXT>) File name for a mask 
%       of the data. If OUT is an empty string, the name of the 
%       outputs will be the same as the inputs, with a '_mask' suffix added 
%       at the end.
%
% OPT           
%   (structure) with the following fields.  
%
%   TYPE_CORRECTION       
%      (string, default 'mean_var') possible values :
%      'none' : no correction at all                       
%      'mean' : correction to zero mean.
%      'mean_var' : correction to zero mean and unit variance
%      'mean_var2' : same as 'mean_var' but slower, yet does not use as 
%      much memory).
%
%   FOLDER_OUT 
%      (string, default: path of IN) If present, all default outputs 
%      will be created in the folder FOLDER_OUT. The folder needs to be 
%      created beforehand.
%
%   FLAG_VERBOSE 
%      (boolean, default 1) if the flag is 1, then the function prints 
%      some infos during the processing.
%
%   FLAG_TEST 
%      (boolean, default 0) if FLAG_TEST equals 1, the brick does not do 
%      anything but update the default values in IN, OUT and OPT.
%           
% _________________________________________________________________________
% OUTPUTS:
%
% IN, OUT, OPT: same as inputs but updated with default values.
%              
% _________________________________________________________________________
% SEE ALSO:
% NIAK_CORRECT_MEAN_VAR
%
% _________________________________________________________________________
% COMMENTS:
%
% _________________________________________________________________________
% Copyright (c) <NAME>, <INSTITUTION>, <START DATE>-<END DATE>.
% Maintainer : <EMAIL ADDRESS>
% See licensing information in the code.
% Keywords : PSOM, documentation, template, brick

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
    % define every option needed to run this module
    fields   = {'RefInput', 'InputToReshape','NbInput', 'NbOutput', 'threshold'  , 'flag_test' , 'folder_out', 'output_filename_ext', 'OutputSequenceName'};
    defaults = {1,1, 1, 1, 5, true, '', 'T1map', 'AllName'};
    opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
    
    % list of everything displayed to the user associated to their 'type'
    user_parameter_list = {'Select at least 3 scans with different flip angles'; 'Parameters'; '   .Output filename extension' ; '   .Threshold'};
    user_parameter_type = {'XScan'; ''; 'char'; 'numeric'};
    parameter_default = {''; ''; 'T1map'; 5};
    psom_parameter_list = {''; ''; 'output_filename_ext'; 'threshold'};
    scans_input_DOF = {{'SequenceName'}; ''; ''; ''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF'};
    opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, scans_input_DOF, 'VariableNames', VariableNames);
    % So for no input file is selected and therefore no output

%% Benjamin Modifications
%      % define every option needed to run this module
%       %%   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
%     module_option(:,1)   = {'folder_out',''};
%     module_option(:,2)   = {'flag_test',true};
%     module_option(:,3)   = {'threshold',5};
%     module_option(:,4)   = {'output_filename_ext','T2map'};
%     module_option(:,5)   = {'OutputSequenceName','AllName'};
%     opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
% 
% %     % define every option needed to run this module
% %     fields   = {'threshold'  , 'flag_test' , 'folder_out', 'output_filename_ext', 'OutputSequenceName'};
% %     defaults = {5, true, '', 'T2map', 'AllName'};
% %     opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
%     
% %% list of everything displayed to the user associated to their 'type'
%      % --> user_parameter(1,:) = user_parameter_list
%      % --> user_parameter(2,:) = user_parameter_type
%      % --> user_parameter(3,:) = parameter_default
%      % --> user_parameter(4,:) = psom_parameter_list
%      % --> user_parameter(5,:) = Help : text data which describe the parameter (it
%      % will be display to help the user)
%     user_parameter(:,1)   = {'Select a Multi Spin Echo scan as input','1Scan','','',''};
%     user_parameter(:,2)   = {'Parameters','','','',''};
%     user_parameter(:,3)   = {'   .Output filename extension','char','T2map','output_filename_ext',''};
%     user_parameter(:,4)   = {'   .Threshold','numeric', 5,'threshold',''};
%     VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Help'};
%     opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', 'VariableNames', VariableNames);
% 
% %     % list of everything displayed to the user associated to their 'type'
% %     user_parameter_list = {'Select a Multi Spin Echo scan as input'; 'Parameters'; '   .Output filename extension' ; '   .Threshold'};
% %     user_parameter_type = {'1Scan'; ''; 'char'; 'numeric'};
% %     parameter_default = {''; ''; 'T2map'; 5};
% %     psom_parameter_list = {''; ''; 'output_filename_ext'; 'threshold'};
% %     VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields'};
% %     opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, 'VariableNames', VariableNames);
%    
% % So for no input file is selected and therefore no output
%%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'T1map'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
%%%%%%%%
opt.NameOutFiles = {opt.output_filename_ext};
%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('T1map:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
if ~ischar(files_in.In1{1}) 
    error('files in should be a char');
end

[path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
if ~strcmp(ext_nii, '.nii')
     error('First file need to be a .nii');  
end


if isfield(opt,'threshold') && (~isnumeric(opt.threshold))
    opt.threshold = str2double(opt.threshold);
    if isnan(opt.threshold)
        disp('The threshold used was not a number')
        return
    end
end


%% Building default output names
if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
    opt.folder_out = path_nii;    
end

if isempty(files_out)
   files_out.In1 = {cat(2,opt.folder_out,filesep,name_nii,'_',opt.output_filename_ext,ext_nii)};
end

%% If the test flag is true, stop here !
if opt.flag_test == 1
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate a T1 map a multi angle scans 
% this code come from the the T1_GE_MULTIPLE_ANGLES_UVASC_MODULE coded by
% I. Tropres and adapted by B. Lemasson
tic
% no additional parameter for now

% load the FA files
NbAngles = size(files_in.In1,2);
for angles=1:NbAngles
    fid=fopen(files_in.In1{angles} ,'r');
    if fid>0
        fclose(fid);
        input.image(angles)  = niftiread(files_in.In1{angles});
        input.nifti_header(angles) = niftiinfo(files_in.In1{angles});
        [path, name, ~] = fileparts(files_in.In1{angles});
        jsonfile = [path, '/', name, '.json'];
        fid = fopen(jsonfile, 'r');
        raw = fread(fid, inf, 'uint8=>char');
        fclose(fid);
        %raw = reshape(raw, 1,length(raw));
        input.json(angles) = jsondecode(raw);

%         FAuvasc = load(FA_filenames{angles}{:});
%         FA(angles) = FAuvasc.uvascim.image;
    else
        warning_text = sprintf('##$ Can not calculate the T1 map because there is\n##$ something wrong with the data \n##$ Could not open the files');
        msgbox(warning_text, 'T10 map warning') ;
        return
    end
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
maxim=max(data_in_vector) * opt.threshold / 100;
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




