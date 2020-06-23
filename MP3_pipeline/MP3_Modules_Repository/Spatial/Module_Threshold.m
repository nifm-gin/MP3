function [files_in,files_out,opt] = Module_Mask(files_in,files_out,opt)
%Module_Template Template of a MP3 module, with comments and explanations.
%   This template can be modified and pasted into the Module_Repository
%   folder of the MP3 source code in order to be used in the MP3 GUI.

%   MP3 is an open source software aimed to support medical researchers all
%   along their studies. Indeed, MP³ offers a graphical interface to
%   convert, visualize, compute and analyze medical images. It allows to
%   use well known processes or to develop yours and to apply them in
%   complex pipelines. All the post processing you used to apply on each
%   subset of your data can now be stored, reproduced and parallelized.
%   Although it was initially developed for MRI data, this software can be
%   used for all medical images that matches one of the required formats.
%   It has for example been used on scanner data. MP³ has been developed at
%   the Grenoble Institute of Neurosciences, in France, and is downoadable
%   at: https://github.com/nifm-gin/MP3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(opt)
    %% Here define every parameter needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    
    module_parameters(:,1)   = {'output_filename_ext','_Masked'};
    module_parameters(:,2)   = {'OutputSequenceName','Extension'};
    module_parameters(:,3)   = {'Threshold','1'};
    
    
    %% System parameters : Do not modify without understanding the behaviour of the software.
    
    system_parameters(:,1)   = {'RefInput',1};
    system_parameters(:,2)   = {'InputToReshape',1};
    
    
    %% Initialisation parameters : Do not modify without understanding the behaviour of the software.
    
    initialisation_parameters(:,1)   = {'folder_out',''};
    initialisation_parameters(:,2)   = {'flag_test',true};
    initialisation_parameters(:,3)   = {'Table_in', table()};
    initialisation_parameters(:,4)   = {'Table_out', table()};
    
    Parameters = [module_parameters, system_parameters, initialisation_parameters];
    
    opt.Module_settings = psom_struct_defaults(struct(),Parameters(1,:),Parameters(2,:));
    
    
    %% Each line displayed to the user :
    
    % --> user_parameter(1,:) = user_parameter_list: String displayed to
    % the user
    
    % --> user_parameter(2,:) = user_parameter_type: Type of the needed
    % parameter. Types available so far:
    % Classic types:
    %   - Text: In order to display text in the Description window.
    %   - char: In order to ask the user to type a string.
    %   - cell: In order to ask the user to choose one value between
    %   several.
    %   - numeric: In order to ask the user to type a number.
    % Scans:
    %   - 1Scan: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   Scan.
    %   - XScan: In order to ask the user to select one or several values
    %   of the 'SequenceName' tag among the ones of its database whose type
    %   is Scan.
    %   - 1ROI: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   ROI.
    %   - XROI: In order to ask the user to select one or several values of
    %   the 'SequenceName' tag among the ones of its database whose type is
    %   ROI.
    %   - 1ScanOr1ROI: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   Scan, ROI, or Cluster.
    %   - XScanOrXROI: In order to ask the user to select one or several
    %   values of the 'SequenceName' tag among the ones of its database
    %   whose type is Scan, ROI, or Cluster.
    %   - 1Cluster: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   Cluster.
    %   - XCluster: In order to ask the user to select one or several
    %   values of the 'SequenceName' tag among the ones of its database
    %   whose type is Cluster.
    %   - 1Scan1TPXP: In order to ask the user to select one value of the
    %   'SequenceName' tag, one value of the 'Tp' tag, and one or several
    %   values of the 'Patient' tag among the ones of its database whose
    %   type is Scan. (Only use in the Coregistrations modules so far, in
    %   order to coregister scans between timepoints or patients.
    
    % --> user_parameter(3,:) = parameter_default: for parameters of type
    % 'cell', here lies the different values that will be displayed to the
    % user.
    
    % --> user_parameter(4,:) = psom_parameter_list: The name of the
    % parameter in the PSOM structure: 'module_parameters'
    
    % --> user_parameter(5,:) = Scans_input_DOF: Degrees of Freedom for
    % the user to choose the scan: for scans/ROI/Clusters inputs, please
    % type here the names of the tags selectables by the user. In fact, if
    % you have never heard of the tags or tried to modify them, dont touch
    % them and let the defaults values :
    %   - 1Scan, XScan, 1ROI, XROI, 1ScanOr1ROI, XScanOrXROI, 1Cluster,
    %   XCluster: {'SequenceName'}
    %   - 1Scan1TPXP: {'SequenceName', 'Tp', 'Patient'}
    
    % --> user_parameter(6,:) = IsInputMandatoryOrOptional: For
    % scans/ROI/Clusters, please type here if the input is Mandatory or
    % Optional. If none, the input is set as Optional.
    
    % --> user_parameter(7,:) = Help : text data which describe the
    % parameter (it will be display to help the user)
    
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'This module create a ROI with all the voxels whose value is higher than Threshold.'
        }'};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Threshold','numeric','','Threshold','', '',''};
    
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
    opt.Table_out.Type = categorical({'ROI'});
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

Output = cast(N>str2double(opt.Threshold), 'int16');
Output = Output(:,:,:,1,1,1,1);
% Get the output header
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = 'int16';


% Save informations about the module in the json that will be associated to
% the file out.
%J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

%% Write the output files

% Write the nifti, thanks to the data, the header, and the file_out name.
niftiwrite(Output, files_out.In1{1}, info2)

% Write the associated json.
%[path, name, ~] = fileparts(files_out.In1{1});
%jsonfile = [path, '/', name, '.json'];
%WriteJson(J, jsonfile)

%% It's already over !

