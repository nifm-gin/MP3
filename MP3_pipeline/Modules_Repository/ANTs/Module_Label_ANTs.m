function [files_in,files_out,opt] = Module_Label_ANTs(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

	%   % define every option needed to run this module
	% --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'rois', ''};
    module_option(:,2)   = {'output_extension', '_report'};
    module_option(:,3)   = {'rotate_angle',     0};
    
    module_option(:,4)   = {'RefInput',         1};
    module_option(:,5)   = {'InputToReshape',   1};
    module_option(:,6)   = {'Table_in',         table()};
    module_option(:,7)   = {'Table_out',        table()};
    module_option(:,8)   = {'OutputSequenceName','Extension'};
    
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
        'ANTs Atlas Matching:'
        '    ANTS was created by:'
        '    Brian B. Avants, Nick Tustison and Gang Song'
        '    Penn Image Computing And Science Laboratory'
        '    University of Pennsylvania'
        '    https://github.com/ANTsX/ANTs'
        }'};
    
    user_parameter(:,2)   = {'Select label map','1ScanOr1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select maps as input','XScan','','',{'SequenceName'}, 'Optionnal',''};
    user_parameter(:,4)   = {'Parameters','','','','', '', ''};
    user_parameter(:,5)   = {'   .ROIs','char','','rois','', 'Mandatory',...
        {'Select all ROIs to extract separeted by '','' character (no space)'}};
    user_parameter(:,6)   = {'   .Ouput extension','char','','output_extension','', '',...
        {'Give the output extension'}};
    user_parameter(:,7)   = {'   .Rotation angle','char','','rotate_angle','', '',...
        {'Give the rotation angle if needed'}};
    
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
label_map   = niftiread(files_in.In1{1});
label_list	= unique(label_map(:));
label_list  = label_list(label_list ~= 0);

% Is needed, rotate label map
label_map = imrotate3(label_map,opt.rotate_angle,[0 0 1]);

% Extract parameters
nb_maps     = length(files_in.In2);
nb_rois     = length(label_list);

% Create and init struct
Report    	= struct();
%metadata
Report.name	= string(opt.Table_in.Patient(1));
Report.time	= string(opt.Table_in.Tp(1));
Report.model = string(opt.Table_in.Group(1));
Report.group = string();
%data
Report.data.labels{1} = string(opt.Table_in.SequenceName(2:end));
Report.data.labels{2} = string(label_list);
Report.data.labels{3} = {'Mean','Std','Skewness','Kurtosis',...
                       'Median','25thPercentile','75thPercentile',...
                       'Density'};
Report.data.values = nan(length(Report.data.labels{1}),...%Maps
                         length(Report.data.labels{2}),...%ROIs
                         length(Report.data.labels{3}));  %Summ
                   
                   
% Complete data values
for map = 1:nb_maps
    
    if exist(files_in.In2{map},'file')
        
        % Load map
        complete_map = niftiread(files_in.In2{map});
        %TODO: improve this rescale
        complete_map = imresize3(complete_map, size(label_map)); 
                
        for lab = 1:nb_rois
            
            % Hist
            h   = complete_map(label_map == label_list(lab));
            
            % Resume
            Report.data.values(map,lab,1)   = nanmean(h);
            Report.data.values(map,lab,2)   = nanstd(h);
            Report.data.values(map,lab,3)   = skewness(h);
            Report.data.values(map,lab,4)   = kurtosis(h);
            Report.data.values(map,lab,5)   = nanmedian(h);
            Report.data.values(map,lab,6)   = prctile(h,25);
            Report.data.values(map,lab,7)   = prctile(h,75);
            Report.data.values(map,lab,8)   = sum(~isnan(h));
        end
    end
end

folder = [fileparts(opt.folder_out) filesep 'Reports'];
if ~exist(folder,'dir'), mkdir(folder); end
    
save([folder filesep char(opt.Table_in.Patient(1)) '_' ...
      char(opt.Table_in.Tp(1)) opt.output_extension '.mat'], 'Report')


