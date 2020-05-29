function [files_in,files_out,opt] = Module_FSL_HD_BET(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename_ext','_hd-BET'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'tta','No'};
    module_option(:,6)   = {'CPU_GPU','cpu'};
    module_option(:,7)   = {'Mode','fast'};
    module_option(:,8)   = {'RefInput',1};
    module_option(:,9)   = {'InputToReshape',1};
    module_option(:,10)   = {'Table_in', table()};
    module_option(:,11)   = {'Table_out', table()};
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
        'Automated deep-learning based brain extraction from multi-sequence MRI. If you are using HD-BET, please cite the following publication: '
        'Isensee F, Schell M, Tursunova I, Brugnara G, Bonekamp D, Neuberger U, Wick A, Schlemmer HP, Heiland S, Wick W, Bendszus M, Maier-Hein KH, Kickingereder P.'
        'Automated brain extraction of multi-sequence MRI using artificial neural networks. Hum Brain Mapp. 2019; 1?13. https://doi.org/10.1002/hbm.24750'
        ''
        'This module needs a brain anatomical image as input (FLAIR/T1c/T1 or T2w)'
        'Output :' 
        '   --> a brain mask saved as ROI'
        '   --> input image masked'
        ''
        'WARNING'
        'First you need to install all Python code needed. Please follow instruction details here'
        'https://github.com/NeuroAI-HD/HD-BET'
        }'};
    
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',...
        {'file: can only be 3D. No support for 4d images'
         'can be either a pre- or postcontrast T1-w, T2-w or FLAIR MRI sequence.' 
         'Other modalities might work as well. Input images must match the orientation of standard MNI152 template! Use fslreorient2std 2 upfront to ensure that this is the case'}'
                        };
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extention (same as the ROI name)','char','_hd-BET','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Select the processor type you are using?','char','cpu','CPU_GPU','', '',...
        {'used to set on which device the prediction will run. Must be either int or str.'
        'Use int for GPU id or ''cpu'' to run on CPU. When using CPU you should consider disabling tta.'}'
        };
   user_parameter(:,6)   = {'   .Fast or accurate?','cell',{'fast', 'accurate'},'Mode','', '',...
       {'can be either ''fast'' or ''accurate''. Fast will use only one set of parameters whereas accurate will use the'
        'five sets of parameters that resulted from our cross-validation as an ensemble'}'
        };

    user_parameter(:,7)   = {'   .Test time data augmentation (mirroring)','cell',{'Yes', 'No'},'tta','', '',...
        {'whether to use test time data augmentation (mirroring).'
        'Disable this if you are using CPU to speed things up!'}'
        };
 
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
    % add the input image masked in the table_out
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
   %11047-20180321-FLAIR_20200520-103707002_hd-BET.nii
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    
    % add the mask in the table_out
    opt.Table_out(2,:) = opt.Table_out(1,:);
    opt.Table_out.Type(2) = categorical(cellstr('ROI'));
    opt.Table_out.SequenceName(2) = categorical(cellstr([char(opt.Table_out.SequenceName(1)), '_mask']));
    opt.Table_out.Filename(2) = categorical(cellstr([char(opt.Table_out.Patient(2)), '_', char(opt.Table_out.Tp(2)), '_', char(opt.Table_out.SequenceName(2))]));
    f_out2 = [char(opt.Table_out.Path(2)), char(opt.Table_out.Patient(2)), '_', char(opt.Table_out.Tp(2)), '_', char(opt.Table_out.SequenceName(2)), '.nii'];
    files_out.In2{1} =f_out2;
end


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_FSL_HD_BET','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
%setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

% command to excute the algorithm
if ~strcmp(opt.CPU_GPU, 'cpu') && isnan(str2double(opt.CPU_GPU))
     warndlg('Wrong processor type selected in the HD_BET module')
end
if strcmp(opt.CPU_GPU, 'cpu')
    device = '-device cpu';
else
    device = ['-device '  opt.CPU_GPU]; 
end
if strcmp(opt.Mode, 'fast')
    mode = '-mode fast';
else
    mode = '-mode accurate';
end
if strcmp(opt.tta, 'Yes')
   tta = '-tta 1';
else
   tta = '-tta 0' ;
end

cmd = strcat('hd-bet', {' -i '}, files_in.In1{:}, {' -o '}, files_out.In1{:},{' '}, device, {' '}, mode,{' '}, tta);
% execute the command
system(cmd{:});



%11047_20180321_FLAIR_hd-BET.nii_mask.nii
%% something is wrong with the transformation matrice between the input and
%% outputs. Therefore we have to correct that here in matlab
% first, we get the input nii information 
info = niftiinfo(files_in.In1{:});
%unzip the image masked
gunzip(strrep(files_out.In1{:}, '.nii', '.nii.gz'))
% load the nii
output_image = niftiread(files_out.In1{:});
% write the nii using the input header
% update the input header
info.Filename =files_out.In1{:};
% crop file before saving it (remove empty slices)
[output_image,FinalMat] = CropNifti(output_image,info.Transform.T');
info.Transform = affine3d(FinalMat');
info.ImageSize = size(output_image); 
%save the nifti updated (croped and correct matric)
niftiwrite(output_image,files_out.In1{:}, info);

% same for the mask
% unzip the mask % 
gunzip(strrep(files_out.In1{:}, '.nii', '.nii_mask.nii.gz'))
% load the nii
output_mask = niftiread(strrep(files_out.In1{:}, '.nii', '.nii_mask.nii'));
% reload and update the input header
clear info
info = niftiinfo(files_in.In1{:});
% update the dataype in the input header (because the datatype of the mask
% is not the same as the input data)
info.Filename =files_out.In2{:};
info.Datatype = 'int32';
% % crop file before saving it (remove empty slices)
[output_mask,FinalMat] = CropNifti(output_mask,info.Transform.T');
info.Transform = affine3d(FinalMat');
info.ImageSize = size(output_mask); 
% save the nifti updated (croped and correct matric)
niftiwrite(output_mask,files_out.In2{:}, info);

%% Json processing for the ouput 1 (the scan masked)
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
if isfile(jsonfile)
    J = ReadJson(jsonfile);
    
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    
    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end

