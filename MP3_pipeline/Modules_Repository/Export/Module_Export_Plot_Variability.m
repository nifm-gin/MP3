function [files_in,files_out,opt] = Module_Export_Plot_Variability(files_in,files_out,opt)

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
    module_option(:,7)   = {'FieldsToStudy','Manufacturer, Manufacturer_sModel'};
    module_option(:,8)   = {'AutomaticJobsCreation', 'No'};
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
         {'This module plot the mean and std of the selected scan in the selected ROI for each session file. Each different values of the selected JSON Field are displayed with a different color while the different Group values differs by the marker used. '}
        };
    user_parameter(:,2)   = {'   .Scan','XScan','','', {'SequenceName'},'Mandatory',...
         'One figure for each scan'};
    user_parameter(:,3)   = {'   .ROI','1ROI','','',{'SequenceName'},'Optional',...
         'You can also select a ROI on which the mean and std will be computed'};
    user_parameter(:,4)   = {'   .JSON field to study','char','','FieldsToStudy','', '','The names of the fields on which base the variability, separated by commas. Be careful: if there is too much fields, each session scan will displayed in a different color ....'};
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


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Export_Values_VoxelByVoxel:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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





database = opt.Table_in;
FieldsToStudy = opt.FieldsToStudy;
FieldsToStudy = FieldsToStudy(~isspace(FieldsToStudy));
FieldsToStudy = strsplit(FieldsToStudy, ',');



%% Add Tag about scanner model
for i=1:size(database,1)
    if strcmp(char(database.Type(i)), 'Scan')
        json_fname = [char(database.Path(i)), filesep, char(database.Filename(i)), '.json'];
        J = ReadJson(json_fname);
        AllFields = {};
        for k=1:length(FieldsToStudy)
           if isfield(J, FieldsToStudy{k}) 
               Val = J.(FieldsToStudy{k}).value;
               if iscell(Val)
                   Val = Val{1};
               end
               if ~ischar(Val) && isnumeric(Val)
                   Val = num2str(Val);
               elseif ~ischar(Val)
                   Val = 'Non_Supported_Type';
               end
               AllFields = [AllFields, {Val}];
           else 
               AllFields = [AllFields, {'Unknown'}];
           end
        end
        AllFields = strjoin(AllFields, ' - ');
        database.VariabilityField(i) = categorical({AllFields});
    else
        database.VariabilityField(i) = categorical({'Undefined'});
    end
end





%% Compute the mean + std of each feature_map on each ROI and plot those 2 values in function of the scanner model and the group of the patient.
features_maps = cellstr(unique(opt.Table_in.SequenceName(opt.Table_in.Type == 'Scan',:)));% For each different feature, a new graph.
%features_maps = {'CoregEstCT_Dn_Resampled05_Clipped_Smooth_Clipped_Masked', 'CoregEstCT_Dn_Resampled05_Clipped_Smooth_Clipped_Masked_LocalMean'};
ro = unique(opt.Table_in.SequenceName(opt.Table_in.Type == 'ROI',:));
if ~isempty(ro)
    ROIs = cellstr(ro);
end
%ROIs = {'Brain_BET'};
Markers = {'o', 'x', '*', '+', '.', 's', 'd', '^', 'v', '>', '<', 'p', 'h'}; % If several Groups, please add more markers


for i=1:length(features_maps)
    Datab = database(database.SequenceName == features_maps{i},:);
    Means = [];
    STDs = [];
    Models = [];
    Group = [];
    ID = [];
    for k=1:height(Datab)
        f_fm = [char(Datab.Path(k)), filesep, char(Datab.Filename(k)), '.nii'];
        if ~exist(f_fm) && exist(strrep(f_fm, '.nii', '.nii.gz'))
            gunzip(strrep(f_fm, '.nii', '.nii.gz'))
        end
        h_fm = spm_vol(f_fm);
        if ~isempty(ro)
            db = database(database.Patient == Datab.Patient(k),:);
            db = db(db.Tp == Datab.Tp(k), :);
            db = db(db.SequenceName == ROIs{1}, :);
            if height(db)~=1
                warning(['error for ', char(Datab.Patient(k)), ' - ', char(Datab.Tp(k))]);
                continue
            end
            f_ROI = [char(db.Path), filesep, char(db.Filename), '.nii'];
            if ~exist(f_ROI) && exist(strrep(f_ROI, '.nii', '.nii.gz'))
                gunzip(strrep(f_ROI, '.nii', '.nii.gz'))
            end
            h_ROI = spm_vol(f_ROI);
            Vol = read_volume(h_fm, h_ROI, 0, 'Axial');
            ROI = read_volume(h_ROI, h_ROI, 0, 'Axial');
            Vol = cast(Vol, 'single');
            M = nanmean(Vol(logical(ROI)));
            S = nanstd(Vol(logical(ROI)));
        else
            Vol = read_volume(h_fm, h_fm, 0, 'Axial');
            Vol = cast(Vol, 'single');
            M = nanmean(Vol(:));
            S = nanstd(Vol(:));
        end
        Means = [Means, M];
        STDs = [STDs, S];
        Models = [Models, {char(Datab.VariabilityField(k))}];
        Group = [Group, {char(Datab.Group(k))}];
        ID = [ID, {[char(Datab.Patient(k)), ' - ', char(Datab.Tp(k))]}];

    end
    
    Types = unique(Models);
    Groups = unique(Group);
    Colormap = jet(length(Types));
    f = figure;
    for l = 1:length(Means)
        [~, ind] = max(contains(Types, Models(l)));
        Color = Colormap(ind,:);
        [~, ind2] = max(contains(Groups, Group(l)));
        Marker = Markers{ind2}; 
        plot(Means(l), STDs(l), Marker, 'Color', Color, 'DisplayName', [Models{l}, ' - ', Group{l}], 'MarkerSize', 10, 'linewidth', 2);
        f.Children(1).Children(1).UserData = ID{l};
        hold on
    end
    if ~isempty(ro)
        title([features_maps{i}, ' - ', ROIs{1} ], 'Interpreter', 'none')
    else
        title(features_maps{i}, 'Interpreter', 'none')
    end
    LH = [];
    for ii=1:length(Types)
        LH(ii) = plot(nan, nan, 's', 'Color', Colormap(ii,:), 'MarkerFaceColor', Colormap(ii,:));
    end
    for iii=1:length(Groups)
        LH(length(Types)+ iii) = plot(nan, nan, Markers{iii}, 'Color', [0 0 0]);
    end
    
    
    lgd = legend(LH, [Types, Groups], 'Interpreter', 'none');
    title(lgd, strjoin(FieldsToStudy, ' - '));
    xlabel('Mean')
    ylabel('STD')
end

