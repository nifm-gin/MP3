function [files_in,files_out,opt] = Module_T1map_MultiAngles(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
     % define every option needed to run this module
         % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'threshold',5};
    module_option(:,4)   = {'output_filename_ext','T1map'};
    module_option(:,5)   = {'OutputSequenceName','AllName'};
    module_option(:,6)   = {'RefInput',1};
    module_option(:,7)   = {'InputToReshape',1};
    module_option(:,8)   = {'Table_in', table()};
    module_option(:,9)   = {'Table_out', table()};
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
    user_parameter(:,1)   = {'Select at least 3 scans with different flip angles','XScan','','',{'SequenceName'},'Mandatory',''};
    user_parameter(:,2)   = {'Parameters','','','','','',''};
    user_parameter(:,3)   = {'   .Output filename extension','char','T1map','output_filename_ext','','',''};
    user_parameter(:,4)   = {'   .Threshold','numeric', 5,'threshold','', '',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);

% So for no input file is selected and therefore no output
%
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
    error('T1map:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs


% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii');  
% end
% 
% 
% 
% %% Building default output names
% if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
%     opt.folder_out = path_nii;    
% end


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

% generate a T1 map a multi angle scans 
% this code come from the the T1_GE_MULTIPLE_ANGLES_UVASC_MODULE coded by
% I. Tropres and adapted by B. Lemasson

% no additional parameter for now

% load the FA files
NbAngles = size(files_in.In1,2);

for angles=1:NbAngles
    fid=fopen(files_in.In1{angles} ,'r');
    if fid>0
        fclose(fid);
        %tmp =  niftiread(files_in.In1{angles});
        input(angles).nifti_header = spm_vol(files_in.In1{angles});
        data_to_process(:,:,:,angles) = read_volume(input(angles).nifti_header, input(1).nifti_header, 0, 'Axial');
        input(angles).json = spm_jsonread(strrep(files_in.In1{angles}, '.nii', '.json'));
        flip_angles(angles) = input(angles).json.FlipAngle.value;
        TR(angles) = input(angles).json.RepetitionTime.value;
    else
        warning_text = sprintf('##$ Can not calculate the T1 map because there is\n##$ something wrong with the data \n##$ Could not open the files');
        msgbox(warning_text, 'T10 map warning') ;
        return
    end
end

% Order the data (from the smallest to hightest flip angle)
[~,flip_angles_index] = sort(flip_angles);
flip_angles = flip_angles(flip_angles_index);
TR = TR(flip_angles_index);
data_to_process =  data_to_process(:,:,:,flip_angles_index);

%% calculate T1map
data_in_vector = reshape(data_to_process, [size(data_to_process,1)*size(data_to_process,2)*size(data_to_process,3), size(data_to_process,4)]);
maxim=max(data_in_vector) * opt.threshold / 100;
fit_result = NaN([size(data_in_vector,1), 2]);
erreur = zeros([size(data_in_vector,1), 1]);

parfor i=1:size(data_in_vector,1)
    vector_to_process = data_in_vector(i,:);
    if  max(vector_to_process(:))>= maxim
        x = vector_to_process ./ tan(flip_angles/180*pi)/ mean(vector_to_process);
        y = vector_to_process ./ sin(flip_angles/180*pi)/ mean(vector_to_process);
        X = [ones(length(x'),1) x'];
        p = X\y';
        % Calculate R2 --> in order to estimate the fit quality
        yfit = p(2)*x + p(1);
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        rsq = 1 - SSresid/SStotal;
        
        fit_result(i,1) = p(2); % slope of the curve
        erreur(i) = rsq*100;    % R2*100
    end
end

fit_result(:,2) = erreur;
fit_result=reshape(fit_result,[size(data_to_process,1),size(data_to_process,2),size(data_to_process,3), 2]);
fit_result(fit_result(:,:,:,1)<0.) = nan;
T1map = -mean(TR)./ log(fit_result(:,:,:,1)); % en ms

% transform the T1map matrix in order to match to the nii hearder of the
% first input (rotation/translation)
T1map = write_volume(T1map, input(1).nifti_header, 'Axial');


%% save the new files (.nii & .json)
% first load the header of the first input
nii_header = niftiinfo(files_in.In1{1});
% update the header before saving the new .nii
nii_header.Filename = files_out.In1{1};
nii_header.Filemoddate = char(datetime('now'));
nii_header.Datatype = class(T1map);
nii_header.ImageSize = size(T1map);
%nii_header.Description = [nii_header.Description, ' Modified by the Module T1map_MultiAngles'];
% save the new .nii file
niftiwrite(T1map, files_out.In1{1}, nii_header)

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



