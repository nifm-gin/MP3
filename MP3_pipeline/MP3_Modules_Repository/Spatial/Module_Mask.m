function [files_in,files_out,opt] = Module_Mask(files_in,files_out,opt)

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
    module_option(:,4)   = {'OutOfTheMask','NaN'};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    module_option(:,9)   = {'output_filename_ext','_Masked'};
    module_option(:,10)  = {'Output_orientation','First input'};

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
    user_parameter(:,2)   = {'Select the scan to mask','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the mask (ROI)','1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename extension','char','','output_filename_ext','','',...
        {'Specify the string to be added to the first filename.'
        'Default filename extension is ''_Masked''.'}'};
    user_parameter(:,5)   = {'   .Output orientation','cell',{'Scan', 'ROI'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = Scan'
        '--> Output orientation = ROI'
        }'};
    user_parameter(:,6)   = {'   .Value out of the mask','char','','OutOfTheMask','','',''};
    
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
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Arithmetic:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
if isfield(files_in, 'In2')
    input(2).nifti_header = spm_vol(files_in.In2{1});
end


if strcmp(opt.Output_orientation, 'Scan') || ~isfield(files_in, 'ROI')
    ref_scan = 1;
else
    ref_scan = 2;
end

Scan = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
Types{1} = class(Scan);




if isfield(files_in, 'In2')
    ROI = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
    Types{2} = class(ROI);
    % if there is a missmatch in the data type we have to convert data type
    % before doing any operation except the multiplication to a scalar  
    CommonType = FindCommonDatatype(Types);       
    Scan = cast(Scan, CommonType);
    ROI = cast(ROI, CommonType);
    
end




if length(size(Scan)) > length(size(ROI))
    OutputImages = Scan .* repmat(ROI, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]);    
else
    OutputImages =  repmat(Scan, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]) .* ROI;    
end

OutputImages(~ROI) = str2double(opt.OutOfTheMask);




% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header, 'Axial');
end

info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});

nifti_header_output = info;
nifti_header_output.Filename = files_out.In1{1};
nifti_header_output.Filemoddate = char(datetime('now'));
[OutputImages_reoriented, FinalMat] = CropNifti(OutputImages_reoriented, nifti_header_output.Transform.T');
nifti_header_output.Datatype = class(OutputImages_reoriented);
nifti_header_output.Transform = affine3d(FinalMat');
nifti_header_output.ImageSize = size(OutputImages_reoriented); 
nifti_header_output.PixelDimensions = info.PixelDimensions(1:length(nifti_header_output.ImageSize));
nifti_header_output.MultiplicativeScaling = 1;

% % save the new .nii file
 niftiwrite(OutputImages_reoriented, files_out.In1{1}, nifti_header_output);

%% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
if isfile(jsonfile)
    J = ReadJson(jsonfile);
    
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    
    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end



