function [files_in,files_out,opt] = Module_Apply_DL_model(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(opt)
    %% Here define every parameter needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    
    
    module_option(:,1)   = {'model_folder_filename',  'MP3/data/DL_models.model.tar'};
    module_option(:,2)   = {'prefix',           'DL_'};
    module_parameters(:,1)  = {'OutputSequenceName','AllName'};
    
    
%     module_parameters(:,1)   = {'output_filename_ext','_Smooth'};
%     module_parameters(:,2)   = {'OutputSequenceName','Extension'};
%     module_parameters(:,3)   = {'Type','2-D Gaussian'};
%     module_parameters(:,4)   = {'Sigma',1};
%     
    
    %% System parameters : Do not modify without understanding the behaviour of the software.
    
    system_parameters(:,1)   = {'RefInput',1};
    system_parameters(:,2)   = {'InputToReshape',1};
    
    
    %% Initialisation parameters : Do not modify without understanding the behaviour of the software.
    
    initialisation_parameters(:,1)   = {'folder_out',''};
    initialisation_parameters(:,2)   = {'flag_test',true};
    initialisation_parameters(:,3)   = {'Table_in', table()};
    initialisation_parameters(:,4)   = {'Table_out', table()};
    
    Parameters = [module_parameters, system_parameters, initialisation_parameters];
    
    opt.Module_settings  = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    
    %% Each line displayed to the user :
    
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'Models must be .tar files stored under MP3/data/DL_models'
        }'};
    
    user_parameter(:,2)   = {'Select the scan to use','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    model_files	= dir(fullfile(s{1}, 'data/DL_models/'));
    model_files    = model_files(~[model_files.isdir]);
    [~, model_names, ~] = cellfun(@fileparts, {model_files.name}, 'UniformOutput', false);
    opt.Module_settings.folder = fullfile(s{1}, 'data/DL_models/');
    if isempty(model_files), model_files(1).name = ' '; end
    user_parameter(:,3)   = {'   .Model to use','cell', model_names, 'model_folder_filename','','Mandatory',...
        {'Select the model to use for data evaluation (.tar)'}};
%     user_parameter(:,3)   = {'Parameters','','','','', '', ''};
%     user_parameter(:,4)   = {'   .Output filename extension','char','','output_filename_ext','', '',''};
%     user_parameter(:,5)   = {'   .Type','cell', {'2-D Gaussian', '3-D Gaussian'},'Type','', '',''};
%     user_parameter(:,6)   = {'   .Sigma','numeric','','Sigma','', '',''};
    
    % Concatenate these user_parameters, and store them in opt.table
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', ...
        user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    %%
    
    % Initialize to an empty string the names of the input and output
    % files.
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
    
end
%%%%%%%%

% Here we generate the names of the files_out, thanks to the names of the
% files_in and some of the parameters selected by the user.
% This paragraph is peculiar to the number and the kind of the files_out.
% Please modify it to automatically generate the name of each output file.
if isempty(files_out)
    opt.Table_out = opt.Table_in;
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
    error('Smoothing:module','Bad syntax, type ''help %s'' for more info.',mfilename)
end


%% If the test flag is true, stop here !
% It's the end of the initialization part.
% At the execution, opt.flag_test will be set to 0 and most of the
% initialization part will not be computed.

if opt.flag_test == 1
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now the core of the module (the operations on the files) starts %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if the files_in exist
[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end


% Load the data, header, and associated json of the inputs
N = niftiread(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);


%% Process your data

Informations = whos('N');
FilteredImages = zeros(size(N), Informations.class);
Size = size(N);
NbDim = length(Size);

if NbDim >4
    Dim_To_Merge = Size(4:end);
    NewDim = prod(Dim_To_Merge);
    NewN = reshape(N, Size(1), Size(2), Size(3), NewDim);
else
    NewN = N;
end

Sigma = opt.Sigma;
if strcmp(opt.Type, '2-D Gaussian')
    for i=1:size(NewN,3)
        for j=1:size(NewN,4)
            FilteredImages(:,:,i,j) = imgaussfilt(NewN(:,:,i,j), Sigma);
        end
    end
elseif strcmp(opt.Type, '3-D Gaussian')
    for j=1:size(NewN,4)
        FilteredImages(:,:,:,j) = imgaussfilt3(NewN(:,:,:,j), Sigma);
    end
    
end

% Get the output data
NewFilteredImages = reshape(FilteredImages, Size);


% Get the output header
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));



% Save informations about the module in the json that will be associated to
% the file out.
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

%% Write the output files

% Write the nifti, thanks to the data, the header, and the file_out name.
niftiwrite(NewFilteredImages, files_out.In1{1}, info2)

% Write the associated json.
[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

%% It's already over !

