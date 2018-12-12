function [files_in,files_out,opt] = Module_Realign_Est(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
%  
%     %%   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','Prefix'};
    module_option(:,4)   = {'Function','nmi'};
    module_option(:,5)   = {'Separation','4 2'};
    module_option(:,6)   = {'Tolerence','0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001'};
    module_option(:,7)   = {'Hist_Smooth','7 7'};
    module_option(:,8)   = {'Interpolation','4th Degree B-Spline'};
    module_option(:,9)   = {'Wrapping','No wrap'};
    module_option(:,10)   = {'Masking','Dont mask images'};
    module_option(:,11)   = {'output_filename_ext','RealignEst'};
    module_option(:,12)   = {'RefInput',1};
    module_option(:,13)   = {'InputToReshape',1};
    module_option(:,14)   = {'Table_in', table()};
    module_option(:,15)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
%   
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
    'Within-subject registration using a rigid-body model and image reslicing.'
    ''
    'The registration method used here is based on work by Collignon et al.'
    'The original interpolation method described in this paper has been changed in order to give a smoother cost function.'
    'The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
    ''
    'A copy of the ''Source'' and '' Other Images'' is performed before the registration process'
    'This module DO NOT overwrite the orignal files' 
    ''
    'Registration parameters are stored in the headers of the "source" and the "other" images.'
    }'};
    user_parameter(:,2)   = {'Reference Image','1Scan','','', {'SequenceName'},'Mandatory',...
         'This is the image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
    user_parameter(:,3)   = {'Parameters','','','','','',''};
    user_parameter(:,4)   = {'    Estimation Options','', '','','','',...
        'Various registration options, which are passed to the Powell optimisation algorithm.'};
    user_parameter(:,5)   = {'       .Objective Function','cell',{'mi','ncc', 'nmi', 'ecc'},'Function','','',...
        'Registration involves finding parameters that either maximise or minimise some objective function. For inter-modal registration, use Mutual Information, Normalised Mutual Information, or Entropy Correlation Coefficient. For within modality, you could also use Normalised Cross Correlation.'};
    user_parameter(:,6)   = {'       .Tolerances','char','','Separation','','',...
        'The average distance between sampled points (in mm).  Can be a vector to allow a coarse registration followed by increasingly fine ones.'};
    user_parameter(:,7)  = {'       .Separation','numeric','','Tolerence','','',...
        'The accuracy for each parameter.  Iterations stop when differences between successive estimates are less than the required tolerance.'};
    user_parameter(:,8)  = {'       .Histogram Smoothing','numeric','','Hist_Smooth','','',...
        'Gaussian smoothing to apply to the 256x256 joint histogram. Other information theoretic coregistration methods use fewer bins, but Gaussian smoothing seems to be more elegant.'};
    user_parameter(:,9)  = {'       .Filename prefix','char','','output_filename_ext','','',...
        'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''RealignEst''.'};


    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional', 'Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)', 'VariableNames', VariableNames);
%%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
  
end
%%%%%%%%


if strcmp(files_out, '')
    [Path_In1, Name_In1, ~] = fileparts(files_in.In1{1});
    tags1 = opt.Table_in(opt.Table_in.Path == [Path_In1, filesep],:);
    tags1 = tags1(tags1.Filename == Name_In1,:);
    assert(size(tags1, 1) == 1);
    tags_out_In1 = tags1;
    tags_out_In1.IsRaw = categorical(0);
    tags_out_In1.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_In1.SequenceName = categorical(cellstr([opt.output_filename_ext, char(tags_out_In1.SequenceName)]));
    tags_out_In1.Filename = categorical(cellstr([char(tags_out_In1.Patient), '_', char(tags_out_In1.Tp), '_', char(tags_out_In1.SequenceName)]));
    f_out = [char(tags_out_In1.Path), char(tags_out_In1.Patient), '_', char(tags_out_In1.Tp), '_', char(tags_out_In1.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    opt.Table_out = tags_out_In1;
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Coreg_Est:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end


%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end

[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FixedImInfo = niftiinfo(files_in.In1{1});
% [path, name, ~] = fileparts(files_in.In1{1});
% FixedImJsonfile = [path, filesep, name, '.json'];
% fid = fopen(FixedImJsonfile, 'r');
% raw = fread(fid, inf, 'uint8=>char');
% fclose(fid);
% %raw = reshape(raw, 1,length(raw));
% FixedImJSON = jsondecode(raw);

info = niftiinfo(files_in.In1{1});

if length(info.ImageSize) ~= 4
    error(['The file ', files_in.In1{1}, ' is not a 4D volume.'])
end

[path, name, ~] = fileparts(files_in.In1{1});
%     TmpFile = [opt.folder_out, '/TMP_FILE_Module_Susceptibility', file,'.nii'];
%     OutTmpFile = [opt.folder_out, '/rTMP_FILE_Module_Susceptibility', file, '.nii'];
%     copyfile(files_in.In1{1},  TmpFile);

%     Scan_to_realign_nii_filename = cell(length(keepscans), 1);
%     for xx = 1:numel(keepscans)
%         Scan_to_realign_nii_filename{xx} =  fullfile([TmpFile ',' num2str(keepscans(xx))]);
%     end
 %matlabbatch{1}.spm.spatial.realign.estwrite.data  = {Scan_to_realign_nii_filename};
%%
matlabbatch{1}.spm.spatial.realign.estwrite.data  = {files_in.In1};
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
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = opt.output_filename_ext;
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch 



delete(strrep(files_out.In1{1}, '.nii', '.mat'))
delete()
%delete(strrep(OutTmpFile, ['rTMP_FILE_Module_Susceptibility', file, '.nii'], ['meanTMP_FILE_Module_Susceptibility', file, '.nii']))
%delete(strrep(OutTmpFile, ['rTMP_FILE_Module_Susceptibility', file, '.nii'], ['rp_TMP_FILE_Module_Susceptibility', file, '.txt']))

jsonfile = [path, '/', name, '.json'];  
J = ReadJson(jsonfile);
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)




    









