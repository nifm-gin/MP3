function [files_in,files_out,opt] = Module_FSL_BET(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

%     %   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'output_filename_ext','Brain_BET'};
    module_option(:,4)   = {'OutputSequenceName','AllName'};
    module_option(:,5)   = {'Fractional_intensity_threshold',0.35};
    module_option(:,6)   = {'Fill_Holes','No'};
    module_option(:,7)   = {'RefInput',1};
    module_option(:,8)   = {'InputToReshape',1};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)   = {'Table_out', table()};
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
        'Brain Extraction Tool from FSL software'
        ''
        'This module needs a brain anatomical image as input'
        'Output --> a brain mask saved as ROI'
        ''
        'WARNING'
        'FSL 6.0.1 needs to be install by the user first in /usr/local/fsl'
        }'};
    
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename','char','Brain_BET','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Fractional Intensity Threshold','numeric',0.35,'Fractional_intensity_threshold','', '',''};
     user_parameter(:,6)   = {'   .Fill holes in the brain?','cell',{'Yes', 'No'},'Fill_Holes','', '',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    %%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
    
end
%%%%%%%%


if isempty(files_out)
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    opt.Table_out.Type = categorical(cellstr('ROI'));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end








%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_FSL_BET','Bad syntax, type ''help %s'' for more info.',mfilename)
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


% set FSL environment

setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

% set FSL-BET command
cmd = strcat('/usr/local/fsl/bin/bet', {' '}, files_in.In1{:}, {' '}, files_out.In1{:},  ' -f', {' '}, num2str(opt.Fractional_intensity_threshold), {' '}, '-g 0 -n -m -t');
% execute the command
system(cmd{:});

% if ~exist strrep(files_out.In1{:}, '.nii', '_mask.nii.gz') 'file'
%     error('Module_FSL_BET:brick','Module failed')
% end
%unzip the mask
gunzip(strrep(files_out.In1{:}, '.nii', '_mask.nii.gz'))

% the Mask needs to be open in order to check the size of the nii. Indeed
% ROI can not be > 3D
N = niftiread(strrep(files_out.In1{:}, '.nii', '_mask.nii'));
info = niftiinfo(strrep(files_out.In1{:}, '.nii', '_mask.nii'));
if length(size(N))>3
    N = squeeze(N(:,:,:,1,1,1,1));
    info.ImageSize = size(N);
    info.PixelDimensions =  info.PixelDimensions(1:length(size(N)));
end
% if the option Fill_Holes is check exectute this loop
if strcmp(opt.Fill_Holes, 'Yes')
    for i=1:size(N,3)
        N(:,:,i) = imfill(N(:,:,i), 'holes');
    end    
end
% update the nii header
info.Datatype = class(N);
% write the nii 
niftiwrite(N, strrep(files_out.In1{:}, '.nii', '_mask.nii'), info);
% rename the output file (BET add by default '_mask' the the file name)
movefile(strrep(files_out.In1{:}, '.nii', '_mask.nii'), strrep(files_out.In1{:}, '.nii', '.nii'), 'f')
