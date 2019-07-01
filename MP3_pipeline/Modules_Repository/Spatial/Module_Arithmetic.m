function [files_in,files_out,opt] = Module_Arithmetic(files_in,files_out,opt)

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
    module_option(:,4)   = {'Constant',1};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    module_option(:,9)   = {'Operation', 'Addition'};
    module_option(:,10)   = {'output_filename_ext','_Arith'};
    module_option(:,11)  = {'Output_orientation','First input'};

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
    user_parameter(:,2)   = {'Select the first scan','1ScanOr1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the operation you would like to apply','cell', {'Addition', 'Subtraction', 'Multiplication (Between Scans)', 'Division', 'Percentage',...
        'Union', 'Intersection', 'Moyenne temporelle',...
        'Addition (Between a Scan and a Scalar)', 'Soubtraction (Between a Scan and a Scalar)', 'Multiplication (Between a Scan and a Scalar)', 'Division (Between a Scan and a Scalar)'},'Operation','', '',''};
    user_parameter(:,4)   = {'Select the second scan','1ScanOr1ROI','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,5)   = {'   .Output filename extension','char','_Smooth','output_filename_ext','','',...
        {'Specify the string to be added to the first filename.'
        'Default filename extension is ''_Arith''.'}'};
    user_parameter(:,6)   = {'   .Select the constant to multiply your scan with','numeric', '', 'Constant', '', '', ''};
    user_parameter(:,7)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        }'};
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


if strcmp(opt.Output_orientation, 'First input') || ~isfield(files_in, 'In2')
    ref_scan = 1;
else
    ref_scan = 2;
end
input1 = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
Types{1} = class(input1);

if isfield(files_in, 'In2')
    input2 = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
    Types{2} = class(input2);
    % if there is a missmatch in the data type we have to convert data type
    % before doing any operation except the multiplication to a scalar
    if ~strcmp(opt.Operation, 'Multiplication (Between Scan 1 and Scalar)')   
        CommonType = FindCommonDatatype(Types);       
        input1 = cast(input1, CommonType);
        input2 = cast(input2, CommonType);
    end
    
end

switch opt.Operation
    case 'Addition'
        OutputImages = input1 + input2;
    case 'Subtraction'
        OutputImages = input1 - input2;
        if strcmp(char(unique(opt.Table_in.Type)), 'ROI')
            OutputImages(OutputImages<0.5) = 0;
            OutputImages(OutputImages>0.5) = 1;
        end
    case 'Multiplication (Between Scans)'
       
        if  sum(strcmp(cellstr(opt.Table_in.Type), 'ROI')) > 0 %sum(opt.Table_in.Type == 'ROI') > 0
            if length(size(input1)) > length(size(input2))
                OutputImages = input1 .* repmat(input2, [1 1 1 size(input1,4) size(input1,5) size(input1,6) size(input1,7)]);    
            else
                OutputImages =  repmat(input1, [1 1 1 size(input1,4) size(input1,5) size(input1,6) size(input1,7)]) .* input2;    
            end
           % OutputImages = OutputImages + abs(min(OutputImages(:)));
            OutputImages(OutputImages == 0) = min(OutputImages(:));
           OutputImages(isnan(OutputImages)) = min(OutputImages(:));
        else
            OutputImages = input1 .* input2;
        end
    case 'Division'
        OutputImages = input1 ./ input2;
        OutputImages(isinf(OutputImages)) = nan;
    case 'Percentage'
        OutputImages = ((input2 - input1) ./ input1 .* 100);
        OutputImages(isinf(OutputImages)) = nan;
%     case 'Concentration'
%         question = strcat('Please enter the relaxivity (mM-1.sec-1)');
%         relaxivity  = inputdlg(question,'Relaxivity',1,{'6'});
%         if isempty(relaxivity)
%             return
%         end
%         result_map.reco.data = (1./(image1.uvascim.image.reco.data/1000) - 1/(image2.uvascim.image.reco.data/1000)) / str2double(relaxivity{:}) ;
%         Methode_info = ['relaxivity = ' relaxivity{:} ' mM-1.sec-1'];
    case 'Union'
        OutputImages = double(input1 | input2);
    case 'Intersection'
        OutputImages = double(input1 & input2);
    case 'Multiplication (Between a Scan and a Scalar)'
        OutputImages = input1 .* opt.Constant;
    case 'Soubtraction (Between a Scan and a Scalar)'
        OutputImages = input1 - opt.Constant;
    case 'Addition (Between a Scan and a Scalar)'
        OutputImages = input1 + opt.Constant;
    case 'Division (Between a Scan and a Scalar)'
        OutputImages = input1 ./ opt.Constant;
    case 'Moyenne temporelle'
        OutputImages = mean(input1,4);
end

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



