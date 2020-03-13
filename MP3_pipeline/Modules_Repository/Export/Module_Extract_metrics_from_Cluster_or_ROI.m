function [files_in,files_out,opt] = Module_Extract_metrics_from_Cluster_or_ROI(files_in,files_out,opt)

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
    module_option(:,3)   = {'RefInput',         1};
    module_option(:,4)   = {'InputToReshape',   1};
    module_option(:,5)   = {'Table_in',         table()};
    module_option(:,6)   = {'Table_out',        table()};
    module_option(:,7)   = {'OutputSequenceName','Extension'};
    
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
        'The goal of this modul is to extract several metrics from a Labels apply to several MR parameters'
        '   The user have to select a Label map and 1 or several MR parameters'
        '   Then both a .csv and a .mat files will be created in the Reports directory of the project'
        '   So far the metrics extacted are :'
        '       Mean, Std, Skewness , Kurtosis' 
        '       Median, 25th Percentil,  75th Percentile and the Density' 
        }'};
    
    user_parameter(:,2)   = {'Select Cluster or ROI map','1ROIOr1Cluster','','',{'SequenceName'}, 'Mandatory','Please select ONE Label map '};
    user_parameter(:,3)   = {'Select MR maps','XScan','','',{'SequenceName'}, 'Mandatory','Please select ONE or SEVERAL MR parameters of interest'};
    user_parameter(:,4)   = {'Parameters','','','','', '', ''};
    user_parameter(:,5)   = {'   .Select labels','char','','rois','', '',...
        {
        'Select every Label you would like to extract separeted by '','' character (no space)'
        '    For example : 3,40,22'
        'Please leave this field empty if you would like to extract EVERY label'
        }'};
    user_parameter(:,6)   = {'   .Ouput extension','char','','output_extension','', '',...
        {'please enter the name of the report generated'}};
    user_parameter(:,7)   = {'   .Mask','1ROI','','',{'SequenceName'}, '',...
        {'Optional : user can add and additional mask to remove a specific areas'}};
    
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
label_vol   = spm_vol(files_in.In1{1});
label_map   = read_volume(label_vol,label_vol,0);
label_list	= unique(label_map(:));
label_list  = label_list(label_list ~= 0);

% Extract parameters
nb_maps     = length(files_in.In2);
nb_rois     = length(label_list);

% Load mask map
if isfield(files_in,'In3')
    mask_vol    = spm_vol(files_in.In3{1});
end

% Create and init struct
Report    	= struct();

%metadata
Report.name	= string(opt.Table_in.Patient(1));
Report.time	= string(opt.Table_in.Tp(1));
Report.model = string(opt.Table_in.Group(1));
Report.group = string();

%data
if isfield(files_in,'In3')
    Report.data.labels{1} = string(opt.Table_in.SequenceName(2:end-1));
else
    Report.data.labels{1} = string(opt.Table_in.SequenceName(2:end));
end
Report.data.labels{2} = string(label_list);
Report.data.labels{3} = {'Mean','Std','Skewness','Kurtosis',...
                       'Median','TwentyFifthPercentile','SeventyFifthPercentile',...
                       'Density', 'Vol_mmxmmxmm'};
Report.data.values = nan(length(Report.data.labels{1}),...%Maps
                         length(Report.data.labels{2}),...%ROIs
                         length(Report.data.labels{3}));  %Summ
                   
                   
% Complete data values
for map = 1:nb_maps
    
    if exist(files_in.In2{map},'file')
        
        % Load MR parameter
        map_vol_header = niftiinfo(files_in.In2{map});
        map_vol = spm_vol(files_in.In2{map});
        complete_map = read_volume(map_vol, map_vol, 0); % 1 for trilinear 
        
        if isfield(files_in,'In3')
             mask_map    = read_volume(mask_vol,map_vol,0); % 0 for nearest neighbour
            complete_map(logical(mask_map)) = nan;
        end
        
        % Load label map
        label_map   = read_volume(label_vol,map_vol,0);
        
        for lab = 1:nb_rois
            
            % Hist
            h   = complete_map(label_map == label_list(lab));
            % remove NaN
            v = ~isnan(h);
            h = h(v);
            if isempty(h)
                Report.data.values(map,lab,1)   = NaN;
                Report.data.values(map,lab,2)   = NaN;
                Report.data.values(map,lab,3)   = NaN;
                Report.data.values(map,lab,4)   = NaN;
                Report.data.values(map,lab,5)   = NaN;
                Report.data.values(map,lab,6)   = NaN;
                Report.data.values(map,lab,7)   = NaN;
                Report.data.values(map,lab,8)   = NaN;
                Report.data.values(map,lab,9)   = NaN;
            else
                
                % Resume
                Report.data.values(map,lab,1)   = nanmean(h);
                Report.data.values(map,lab,2)   = nanstd(h);
                Report.data.values(map,lab,3)   = skewness(h);
                Report.data.values(map,lab,4)   = kurtosis(h);
                Report.data.values(map,lab,5)   = nanmedian(h);
                Report.data.values(map,lab,6)   = prctile(h,25);
                Report.data.values(map,lab,7)   = prctile(h,75);
                Report.data.values(map,lab,8)   = sum(v);
                Report.data.values(map,lab,9)   = sum(v)*prod(map_vol_header.PixelDimensions(1:3));
            end
        end
    end
end

folder = [fileparts(opt.folder_out) filesep 'Reports'];
if ~exist(folder,'dir'), mkdir(folder); end

% save the Report as .mat file
save([folder filesep char(opt.Table_in.Patient(1)) '_' ...
    char(opt.Table_in.Tp(1)) opt.output_extension '.mat'], 'Report')

% save the Report as .csv file
Report_csv = table();
warning('off')
for i = 1:length(Report)
    % find the MR_parameters
    for j = 1 : length(Report.data.labels{1})   
        % find the Labels name
        for z = 1: length( Report.data.labels{2})
            % find the Metric value
            [rowsn , ~] = size(Report_csv);
            for y = 1: length( Report.data.labels{3})
                % store all the patient information 
                Report_csv.name(rowsn+1,:)	= categorical(cellstr(Report.name{i}));
                Report_csv.time(rowsn+1,:)	= categorical(cellstr(Report.time{i}));
                Report_csv.model(rowsn+1,:) = categorical(cellstr(Report.model{i}));
                Report_csv.ClusterOrROI(rowsn+1,:) = categorical(cellstr(opt.Table_in.SequenceName(1)));
                % store MR_parameters
                Report_csv.MR_parameter(rowsn+1,:) = categorical(cellstr(Report.data.labels{1}{j}));
                % store Labels name
                Report_csv.Label_name(rowsn+1,:) = categorical(cellstr(Report.data.labels{2}{z}));
                % store Metric name and value
                Report_csv.(Report.data.labels{3}{y})(rowsn+1,:) = Report.data.values(j, z, y);             
            end
        end
    end
end
warning('on')

writetable(Report_csv, [folder filesep char(opt.Table_in.Patient(1)) '_' ...
    char(opt.Table_in.Tp(1)) opt.output_extension '.csv'])

