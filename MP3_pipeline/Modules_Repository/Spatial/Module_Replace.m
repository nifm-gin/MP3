function [files_in,files_out,opt] = Module_Replace(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename_ext','_Replaced'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'InitialValue',0};
    module_option(:,6)   = {'FinalValue',0};
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
        'This module is usefull if you would like to replace a value by another in a scan.'
        }'};
    user_parameter(:,2)   = {'Select one scan as input','1ScanOr1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Replaced','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Initial Value','numeric',0,'InitialValue','', '','The value that need to be replaced by the final value.'};
    user_parameter(:,6)   = {'   .Final Value','numeric',0,'FinalValue','', '','The Value by which replace the the initial value'};
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



N = niftiread(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
N = N*info.MultiplicativeScaling;
info.MultiplicativeScaling = 1;

NewImages = N;

NewImages(NewImages == opt.InitialValue) = opt.FinalValue;



info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
%info2.Description = [info.Description, 'Modified by Smoothing Module'];

niftiwrite(NewImages, files_out.In1{1}, info2)


    
if strcmp(char(opt.Table_out.Type), 'Scan')
    [path, name, ~] = fileparts(files_in.In1{1});
    jsonfile = [path, filesep, name, '.json'];
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, filesep, name, '.json'];
    WriteJson(J, jsonfile)
elseif strcmp(char(opt.Table_out.Type), 'Cluster')
    [path, name, ~] = fileparts(files_in.In1{1});
    matfilein = [path, filesep, name, '.mat'];
    [path, name, ~] = fileparts(files_out.In1{1});
    matfilefin = [path, filesep, name, '.mat'];
    if exist(matfilein, 'file')
        copyfile(matfilein, matfilefin);
    else
        Sorry = 'Sorry this cluster doesn''t have an associated .mat file';
        save(matfilefin, 'Sorry');
    end
end

