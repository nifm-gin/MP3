function [files_in,files_out,opt] = Module_Export_data4DL(files_in,files_out,opt)

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
    module_option(:,3)   = {'RefInput',1};
    module_option(:,4)   = {'InputToReshape',1};
    module_option(:,5)   = {'Table_in', table()};
    module_option(:,6)   = {'Table_out', table()};
    module_option(:,7)   = {'Output_filename','Data_4DL'};
    module_option(:,8)   = {'OutputSequenceName','AllName'};
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
         {'This module aims to reshape images in order to generate nifti in a specific nifti format'
         '    - inputs: several scans and 1 ROI'
         '    - output: one 4D nifti (image volume in 3D + each input concatenated in the 4th dimenssion per patient and per time point'}
        };
    user_parameter(:,2)   = {'   .Scan of reference','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the scan that will be use as reference. Every other scan selected below will be reoriented to this reference space'};
    user_parameter(:,3)   = {'   .Input Scans','XScan','','', {'SequenceName'},'Mandatory',...
         'Please select evey scan (other than the reference scan) that you would like to concatenate to the reference scan'};
    user_parameter(:,4)   = {'   .ROI','1ROI','','',{'SequenceName'},'Mandatory',...
         'Please select one ROI. This ROI will be used to extract evey voxel in this ROI'};
    user_parameter(:,5)   = {'   .Output filename','char','Data_4DL','Output_filename','', '',''};

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
    [Path_In, Name_In, ~] = fileparts(files_in.In1{1});
    tags = opt.Table_in(opt.Table_in.Path == [Path_In, filesep],:);
    tags = tags(tags.Filename == Name_In,:);
    assert(size(tags, 1) == 1);
    tags_out_In = tags;
    tags_out_In.IsRaw = categorical(0);
    tags_out_In.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_In.SequenceName = categorical(cellstr(opt.Output_filename));
    tags_out_In.Filename = categorical(cellstr([char(tags_out_In.Patient), '_', char(tags_out_In.Tp), '_', char(tags_out_In.SequenceName)]));
    f_out = [char(tags_out_In.Path), char(tags_out_In.Patient), '_', char(tags_out_In.Tp), '_', char(tags_out_In.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    opt.Table_out = tags_out_In;
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Export_data4DL:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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




%% load the scan of reference
scan_of_reference.header = spm_vol(files_in.In1{1});
scan_of_reference.data=  read_volume(scan_of_reference.header, scan_of_reference.header, 0, 'Axial'); %si tout marche, on peut virer cette ligne

%% load the ROI
ROI.header = spm_vol(files_in.In3{1});
ROI.data=  read_volume(ROI.header, scan_of_reference.header, 0, 'Axial'); %si tout marche, on peut virer cette ligne



 

%% load all the input scan (other than the scan of reference)
for i=1:length(files_in.In2)
    other_scan(i).header =  spm_vol(files_in.In2{i});
    other_scan(i).data = read_volume(other_scan(i).header, scan_of_reference.header, 0, 'Axial'); % si tout marche, on peut virer cette ligne
    % read json file
    other_scan(i).json = spm_jsonread(strrep(files_in.In2{i}, '.nii', '.json'));
    
end
% all data need to be in the same class
% we use the function FindCommonDatatype to define the common class
Types = cell(1,length(other_scan)+1);
Types{1} = class(scan_of_reference.data);
for i=1:length(other_scan)
    Types{i+1} = class(other_scan(i).data);
end
CommonType = FindCommonDatatype(Types);

% the ROI needs to have the same class as the data
output_data = cast(scan_of_reference.data, CommonType) .* cast(ROI.data, CommonType);


  % read json file
scan_of_reference.json = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));

output_json = scan_of_reference.json;

% complet the json (for each echoes, repetition, spin-echoes, diffusion
% values..
if ~isfield(output_json, 'SpinEchoTime')
    output_json.SpinEchoTime = output_json.EchoTime;
    output_json.SpinEchoTime.format = [];
    output_json.SpinEchoTime.description = 'SpinEcho time';
    output_json.SpinEchoTime.units = 'ms';
    output_json.SpinEchoTime.type = 'float';
    output_json.SpinEchoTime.value= 0;
end

if numel(size(scan_of_reference.data)) <= 4
    output_json.SequenceName.value = repmat(cellstr(opt.Table_in.SequenceName(1)), [size(scan_of_reference.data,4), 1]);
    if isfield(output_json, 'EchoTime') && numel(output_json.EchoTime.value) ~= size(scan_of_reference.data, 4)
        output_json.EchoTime.value = repmat(scan_of_reference.json.EchoTime.value, [size(scan_of_reference.data,4), 1]);
    end
    if isfield(output_json, 'RepetitionTime') && numel(output_json.RepetitionTime.value) ~= size(scan_of_reference.data, 4)
        output_json.RepetitionTime.value =repmat(scan_of_reference.json.RepetitionTime.value, [size(scan_of_reference.data,4), 1]);
    end
    if isfield(output_json, 'InversionTime') && numel(output_json.InversionTime.value) ~= size(scan_of_reference.data, 4)
        output_json.InversionTime.value =repmat(scan_of_reference.json.InversionTime.value, [size(scan_of_reference.data,4), 1]);
    end
     if isfield(output_json, 'SpinEchoTime') && numel(output_json.SpinEchoTime.value) ~= size(scan_of_reference.data, 4)
        output_json.SpinEchoTime.value =repmat(output_json.SpinEchoTime.value, [size(scan_of_reference.data,4), 1]);
    end
elseif numel(size(scan_of_reference.data)) == 5
    disp('not coded yet')
    return
end


for i=1:length(files_in.In2)
    if size(other_scan(i).data, 5) > 1 % for 5d scan
        
        %         [a, b, c, d, e ] = size(other_scan(i).data);
        %         output_data(:,:,:,size(output_data,4)+1:size(output_data,4)+length(other_scan(i).header)) =  reshape(other_scan(i).data,[a, b, c, d * e ]);
        disp('not coded for 5D data')
        output_data = [];
        return
        
    else
        % concatanate data but first the ROI and the scan to concatanate
        % need to have the same class as the scan of reference
        output_data(:,:,:,size(output_data,4)+1:size(output_data,4)+length(other_scan(i).header)) =  cast(other_scan(i).data, CommonType)  .* cast(ROI.data, CommonType);
        % update json information
        output_json.SequenceName.value = [output_json.SequenceName.value' repmat(cellstr(opt.Table_in.SequenceName(i+1)), [size(other_scan(i).data,4), 1])']';
        if isfield(other_scan(i).json, 'EchoTime.')
            if length(other_scan(i).json.EchoTime.value)  ~= size(other_scan(i).data,4)
                output_json.EchoTime.value = [output_json.EchoTime.value' repmat(other_scan(i).json.EchoTime.value, [size(other_scan(i).data,4), 1])']';
            else
                output_json.EchoTime.value = [output_json.EchoTime.value' other_scan(i).json.EchoTime.value']';
            end
        end
        if isfield(other_scan(i).json, 'RepetitionTime')
            if  length(other_scan(i).json.RepetitionTime.value)  ~= size(other_scan(i).data,4)
                output_json.RepetitionTime.value = [output_json.RepetitionTime.value' repmat(other_scan(i).json.RepetitionTime.value, [size(other_scan(i).data,4), 1])']';
            else
                output_json.RepetitionTime.value = [output_json.RepetitionTime.value' other_scan(i).json.RepetitionTime.value']';
            end
        end
        if isfield(other_scan(i).json, 'InversionTime')
            if  length(other_scan(i).json.InversionTime.value)  ~= size(other_scan(i).data,4)
                output_json.InversionTime.value = [output_json.InversionTime.value' repmat(other_scan(i).json.InversionTime.value, [size(other_scan(i).data,4), 1])']';
            else
                output_json.InversionTime.value = [output_json.InversionTime.value' other_scan(i).json.InversionTime.value']';
            end
        end
        if isfield(other_scan(i).json, 'SpinEchoTime')
            if   length(other_scan(i).json.SpinEchoTime.value)  ~= size(other_scan(i).data,4)
                output_json.SpinEchoTime.value = [output_json.SpinEchoTime.value' repmat(other_scan(i).json.SpinEchoTime.value, [size(other_scan(i).data,4), 1])']';
            else
                output_json.SpinEchoTime.value = [output_json.SpinEchoTime.value' other_scan(i).json.SpinEchoTime.value']';
            end
        end
    end
    
    
end




% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)
if ~exist('OutputImages_reoriented', 'var')
    output_data = write_volume(output_data, scan_of_reference.header, 'Axial');
end


% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.In1{1});

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(output_data);
info2.ImageSize = size(output_data);
info2.PixelDimensions = ones([1,length(size(output_data))]);
info2.PixelDimensions(1:length(info.PixelDimensions)) = info.PixelDimensions;

% save the new .nii file
niftiwrite(output_data, files_out.In1{1}, info2);


%% Json processing
% [path, name, ~] = fileparts(files_in.In1{1});
% jsonfile = [path, '/', name, '.json'];
% J = ReadJson(jsonfile);

output_json = KeepModuleHistory(output_json, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

WriteJson(output_json, strrep(files_out.In1{1}, '.nii', '.json'))





%fill output_data data matrix with the data


