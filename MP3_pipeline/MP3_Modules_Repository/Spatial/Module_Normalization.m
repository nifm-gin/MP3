function [files_in,files_out,opt] = Module_Normalization(files_in,files_out,opt)

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
    module_option(:,3)   = {'OutputSequenceName','Extension'};
    module_option(:,4)   = {'RefInput',1};
    module_option(:,5)   = {'InputToReshape',1};
    module_option(:,6)   = {'Table_in', table()};
    module_option(:,7)   = {'Table_out', table()};
   % module_option(:,8)   = {'Operation', 'Addition'};
    module_option(:,8)   = {'Output_filename_ext','_Norm'};
    module_option(:,9)   = {'Value_of_normalization',1};
    module_option(:,10)  = {'norm_type', 'Mean'};

    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','', '','Description of the module'}  ;
    user_parameter(:,2)   = {'Select the scan to normalize','1Scan','','',{'SequenceName'}, 'Mandatory',''};
 %   user_parameter(:,3)   = {'Select the operation you would like to apply','cell', {'Addition', 'Subtraction', 'Multiplication', 'Division', 'Percentage'},'Operation','', '',''};
    user_parameter(:,3)   = {'Select one ROI','1ROI','','',{'SequenceName'}, 'Mandatory',''};
    norm_options          = {'Mean','Sqrt of square sum'};
    user_parameter(:,4)   = {'   .Normalize by','cell', norm_options, 'norm_type', '', '', {'Normalization to apply'}};  
    user_parameter(:,5)   = {'   .Output filename extension','char','_Norm','Output_filename_ext','','',...
        {'Specify the string to be added to the first filename.'
        'Default filename extension is ''_Norm''.'}'};
    user_parameter(:,6)   = {'   .Value of normalization','numeric',1,'Value_of_normalization','','',...
        {'Specify the value used to normalize the scan'}'};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    %%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
    
end
%%%%%%%%

if isempty(files_out)
    opt.Table_out = opt.Table_in(1,:);
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.Output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Normalize:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%% load input Nii file
input(1).nifti_header = spm_vol(files_in.In1{1});
%input(1).nifti_header.pinfo(1 ) = 1;
input(2).nifti_header = spm_vol(files_in.In2{1});


scan = read_volume(input(1).nifti_header, input(1).nifti_header, 0, 'Axial');
ROI = read_volume(input(2).nifti_header, input(1).nifti_header, 0, 'Axial');
% check the compatibity bewteen scan and ROI
if ~strcmp(class(scan), class(ROI))
    % if mismatch, ROI is converted into the class of the scan
    ROI = cast(ROI, class(scan));
end

%% code to perform the normalization of the scan by the ROI and using the value gave by the user
ROI(ROI == 0) = NaN;
tmp = scan.*ROI;
if strcmp(opt.norm_type, 'Mean')
    scan = scan/nanmean(tmp(:));
    scan =  scan .* opt.Value_of_normalization;
elseif strcmp(opt.norm_type, 'Sqrt of square sum')
    scan = scan./sqrt(sum(scan.^2,4));    
else
    error('Unknown mean normalization type :%s\n', opt.norm_type)
end

%%

% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(scan, input(1).nifti_header, 'Axial');
end
% save the new files (.nii & .json)
% update the header before saving the new .nii

info2 = niftiinfo(files_in.In1{1});

info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(scan);
info2.MultiplicativeScaling = 1;


% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% % so far copy the .json file of the first input
% copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))

%% Json Processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

% 


