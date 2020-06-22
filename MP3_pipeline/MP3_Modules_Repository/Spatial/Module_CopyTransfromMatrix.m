function [files_in,files_out,opt] = Module_CopyTransfromMatrix(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename_ext','_MatCopied'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
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
        'Copy the transform matrix of a scan to another scan.'
        }'};
    user_parameter(:,2)   = {'Select the scan to copy the matrix from','1ScanOr1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the scan to copy the matrix to','1ScanOr1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'Parameters','','','','', '', ''};
    user_parameter(:,5)   = {'   .Output filename extension','char','_Smooth','output_filename_ext','', '',''};
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
    opt.Table_out = opt.Table_in(2,:);
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
    error('Smoothing:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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



%N_from = niftiread(files_in.In1{1});
info_from = niftiinfo(files_in.In1{1});
%[path, name, ~] = fileparts(files_in.In1{1});
%jsonfile = [path, '/', name, '.json'];
%J_from = ReadJson(jsonfile);

N_to = niftiread(files_in.In2{1});
info_to = niftiinfo(files_in.In2{1});
[path, name, ~] = fileparts(files_in.In2{1});
jsonfile = [path, '/', name, '.json'];
if exist(jsonfile, 'file')
    J_to = ReadJson(jsonfile);
    J_to = KeepModuleHistory(J_to, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J_to, jsonfile)
end
info2 = info_to;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
%info2.Description = [info.Description, 'Modified by Smoothing Module'];
info2.Transform.T = info_from.Transform.T;
niftiwrite(N_to, files_out.In1{1}, info2)



