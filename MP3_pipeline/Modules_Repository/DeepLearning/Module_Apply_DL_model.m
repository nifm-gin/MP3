function [files_in,files_out,opt] = Module_Apply_DL_model(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(opt)
    %% Here define every parameter needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    
    
    module_parameters(:,1)   = {'model_folder_filename',  'MP3/data/DL_models.model.tar'};
%     module_parameters(:,2)   = {'model_mask', 'Mask'};
    module_parameters(:,2)   = {'OutputSequenceName','AllName'};
    module_parameters(:,3)   = {'output_filename_BVf','BVf'};
    module_parameters(:,4)   = {'output_filename_VSI','VSI'};
    module_parameters(:,5)   = {'output_filename_SO2','SO2'};   
    
    %% System parameters : Do not modify without understanding the behaviour of the software.
    
    system_parameters(:,1)   = {'RefInput',1};
    system_parameters(:,2)   = {'InputToReshape',1};
    
    
    %% Initialisation parameters : Do not modify without understanding the behaviour of the software.
    
    initialisation_parameters(:,1)   = {'folder_out',''};
    initialisation_parameters(:,2)   = {'flag_test',true};
    initialisation_parameters(:,3)   = {'Table_in', table()};
    initialisation_parameters(:,4)   = {'Table_out', table()};
    
    Parameters = [module_parameters, system_parameters, initialisation_parameters];
    
    opt.Module_settings  = psom_struct_defaults(struct(), Parameters(1,:), Parameters(2,:));
    opt.NameOutFiles = {opt.Module_settings.output_filename_BVf, opt.Module_settings.output_filename_VSI, opt.Module_settings.output_filename_SO2};
    
    %% Each line displayed to the user :
    
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'Models must be .tar files stored under MP3/data/DL_models'
        }'};
    
    user_parameter(:,2)   = {'Select the scan to use','1Scan','','',{'SequenceName'}, 'Mandatory',''};
%     user_parameter(:,3)   = {'   .Mask to apply', '1ROI', '', 'model_mask', {'SequenceName'}, 'Mandatory', 'Mask to restrain the computation to the ROI'};
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    model_files	= dir(fullfile(s{1}, 'data/DL_models/'));
    model_files    = model_files(~[model_files.isdir]);
    [~, model_names, ~] = cellfun(@fileparts, {model_files.name}, 'UniformOutput', false);
    opt.Module_settings.folder = fullfile(s{1}, 'data/DL_models/');
    if isempty(model_files), model_files(1).name = ' '; end
    user_parameter(:,3)   = {'   .Model to use','cell', model_names, 'model_folder_filename','','Mandatory',...
        {'Select the model to use for data evaluation (.tar)'}};
    
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
%%%%%%%
% Here we generate the names of the files_out, thanks to the names of the
% files_in and some of the parameters selected by the user.
% This paragraph is peculiar to the number and the kind of the files_out.
% Please modify it to automatically generate the name of each output file.

opt.NameOutFiles = {opt.output_filename_BVf, opt.output_filename_VSI, opt.output_filename_SO2};

if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in(1,:);
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(strcat(cellstr(opt.NameOutFiles{i}), '_', opt.model_folder_filename));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), '_', opt.NameOutFiles{i}]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
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
% set python command
cmd = sprintf('env -i /home/delphina/data_ssd/anaconda3/bin/python /home/delphina/data_ssd/Git/Drone/EvaluationScript.py -i %s --outFileBVf %s --outFileVSI %s --outFileSO2 %s -n %s%s.tar',...
    files_in.In1{1}, files_out.In1{1}, files_out.In1{2}, files_out.In1{3}, opt.folder, opt.model_folder_filename);
% cmd = strcat('/usr/local/fsl/bin/bet', {' '}, files_in.In1{:}, {' '}, files_out.In1{:},  ' -f', {' '}, num2str(opt.Fractional_intensity_threshold), {' '}, '-g 0 -n -m -t');
% execute the command
system(cmd);

% Get the output header
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));

% Save informations about the module in the json that will be associated to
% the file out.
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
J.ProtocolName.value = opt.model_folder_filename;
J.SequenceName.value = opt.model_folder_filename;

%% Write the output files

% Write the nifti, thanks to the data, the header, and the file_out name.
niftiwrite(NewFilteredImages, files_out.In1{1}, info2)

% Write the associated json.
[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

%% It's already over !

