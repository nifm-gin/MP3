function [files_in,files_out,opt] = Module_Prov(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

	%   % define every option needed to run this module
	% --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values    
    module_option(:,1)   = {'RefInput',         1};
    module_option(:,2)   = {'InputToReshape',   1};
    module_option(:,3)   = {'Table_in',         table()};
    module_option(:,4)   = {'Table_out',        table()};
    module_option(:,5)   = {'OutputSequenceName','Extension'};
    
    opt.Module_settings  = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    
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
        ''
        'o'
        ''
        }'};
    
    user_parameter(:,2)   = {'Select anat map','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select all other maps','XScan','','',{'SequenceName'}, 'Optionnal',''};
    
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


% Load label map
info_ref = niftiinfo([char(opt.Table_in(1,:).Path) char(opt.Table_in(1,:).Filename) '.nii']);

for m = 2:size(opt.Table_in,1)
    info    = niftiinfo([char(opt.Table_in(m,:).Path) char(opt.Table_in(m,:).Filename) '.nii']);
    info.Transform.T(4,:) = info_ref.Transform.T(4,:);
    
    niftiwrite(niftiread([char(opt.Table_in(m,:).Path) char(opt.Table_in(m,:).Filename) '.nii']), [char(opt.Table_in(m,:).Path) char(opt.Table_in(m,:).Filename) '.nii'], info)
end




