function [files_in,files_out,opt] = Module_Export_Values_ROIByROI(files_in,files_out,opt)

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
    module_option(:,7)   = {'Output_filename','Data_exported_ByROI'};
    module_option(:,8)   = {'OutputSequenceName','AllName'};
    module_option(:,9)   = {'AutomaticJobsCreation', 'No'};
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
         {'This module aims to export data ROI-by-ROI.'
         '    - inputs: one or several scans and 1 ROI'
         '    - output: cvs file which contains values for each ROI (means/SD). One line = the mean +/- SD for on ROI'}
        };
    user_parameter(:,2)   = {'   .First scan (used as the reference for the resolution)','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the scan that will be use as reference. Every other scan selected below will be reoriented to this reference space'};
    user_parameter(:,3)   = {'   .Other scans (if needed)','XScan','','', {'SequenceName'},'Optional',...
         'Please select evey scan (other than the reference scan) that you would like to export using the resolution of the scan of reference'};
    user_parameter(:,4)   = {'   .ROI','XROI','','',{'SequenceName'},'Mandatory',...
         'Please select one ROI. This ROI will be used to extract values of every voxel in this ROI'};
    user_parameter(:,5)   = {'   .Output filename','char','Data_exported_ByROI','Output_filename','', '',''};

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
databScans = opt.Table_in(opt.Table_in.Type == categorical(cellstr('Scan')),:);
parameters_to_export = unique(databScans.SequenceName);
databROIs = opt.Table_in(opt.Table_in.Type == categorical(cellstr('ROI')),:);
Patient_listing = unique(opt.Table_in.Patient);
Tp_listing = unique(opt.Table_in.Tp);

for i=1:size(parameters_to_export)
    if find(contains(files_in.In1(1),char(parameters_to_export(i)))) ==1
        scan_of_ref = parameters_to_export((i))  ;
    end
end
        
for x = 1:numel(Patient_listing)
    for y = 1:numel(Tp_listing)
        output_data = [];
        sub_database = opt.Table_in(opt.Table_in.Patient == Patient_listing(x) & opt.Table_in.Tp == Tp_listing(y), :);
        sub_databaseScan = sub_database(sub_database.Type == categorical(cellstr('Scan')),:);
        % figure out which scan is the scan of reference and which scan(s)
        % corresponds to the other scan(s)
        sub_databaseScanofRef = sub_databaseScan(sub_databaseScan.SequenceName == scan_of_ref,:);
        sub_databaseOtherScan = sub_databaseScan(sub_databaseScan.SequenceName ~= scan_of_ref,:);
        
        sub_databaseROI = sub_database(sub_database.Type == categorical(cellstr('ROI')),:);
        if ~isempty(sub_databaseScanofRef)
            %% load the scan of reference
            scan_of_reference.header =   spm_vol([char(sub_databaseScanofRef.Path) char(sub_databaseScanofRef.Filename) '.nii']);
            scan_of_reference.data=  read_volume(scan_of_reference.header, scan_of_reference.header, 0, 'Axial');
            % save all the data in output_data
            output_data = scan_of_reference.data;
            
            % start to create the CSV file
            PatientName = {};
            GroupName = {};
            TimePoint = {};
            ScanOfRef = {};
            Parameters_name= {};
            if length(size(output_data)) == 3
                GroupName{numel(GroupName)+1}  = char(sub_databaseScanofRef.Group);
                PatientName{numel(PatientName)+1}  = char(sub_databaseScanofRef.Patient);
                TimePoint{numel(TimePoint)+1}  = char(sub_databaseScanofRef.Tp);
                ScanOfRef{numel(Parameters_name)+1}  = char(sub_databaseScanofRef.SequenceName);
                Parameters_name{numel(Parameters_name)+1}  = strcat(char(sub_databaseScanofRef.SequenceName), '_mean');
                Parameters_name{numel(Parameters_name)+1}  = strcat(char(sub_databaseScanofRef.SequenceName), '_SD');
            else
                for j =  1:size(output_data, 5)
                    for i = 1:size(output_data, 4)   
                        GroupName{numel(GroupName)+1}  = char(sub_databaseScanofRef.Group);
                        PatientName{numel(PatientName)+1}  = char(sub_databaseScanofRef.Patient);
                        TimePoint{numel(TimePoint)+1}  = char(sub_databaseScanofRef.Tp);
                        ScanOfRef{numel(Scan_name)+1} = char(char(sub_databaseScanofRef.SequenceName));
                        Parameters_name{numel(Parameters_name)+1} = strcat(char(sub_databaseScanofRef.SequenceName), '_4thDim', num2str(i), '_5thDim', num2str(j), '_mean');   
                        Parameters_name{numel(Parameters_name)+1} = strcat(char(sub_databaseScanofRef.SequenceName), '_4thDim', num2str(i), '_5thDim', num2str(j), '_SD');   
                    end
                end
            end
            
            
            
            %% load other scan if needed
            if ~isempty(sub_databaseOtherScan)
                for i=1:size(sub_databaseOtherScan,1)
                    other_scan.header =  spm_vol([char(sub_databaseOtherScan.Path(i)) char(sub_databaseOtherScan.Filename(i)) '.nii']);
                    other_scan.data = read_volume(other_scan.header, scan_of_reference.header, 0, 'Axial');
                    output_data(:,:,:,size(output_data,4)+1:size(output_data,4)+length(other_scan.header)) =  reshape(other_scan.data, [size(output_data,1), size(output_data,2), size(output_data,3), size(other_scan.data,4)*size(other_scan.data,5)]);

                    % save the name of the voxel (head of the column)
                    if length(size(other_scan.data)) == 3
                        GroupName{numel(GroupName)+1}  = char(sub_databaseScanofRef.Group);
                        PatientName{numel(PatientName)+1}  = char(sub_databaseScanofRef.Patient);
                        TimePoint{numel(TimePoint)+1}  = char(sub_databaseScanofRef.Tp);
                        ScanOfRef{numel(ScanOfRef)+1}  = char(sub_databaseScanofRef.SequenceName);
                        Parameters_name{numel(Parameters_name)+1}  = strcat(char(sub_databaseOtherScan.SequenceName(i)), '_mean');
                        Parameters_name{numel(Parameters_name)+1}  = strcat(char(sub_databaseOtherScan.SequenceName(i)), '_SD');

                    else
                        for j =  1:size(other_scan.data, 5)
                            for ii = 1:size(other_scan.data, 4)  
                                GroupName{numel(GroupName)+1}  = char(sub_databaseScanofRef.Group);
                                PatientName{numel(PatientName)+1}  = char(sub_databaseScanofRef.Patient);
                                TimePoint{numel(TimePoint)+1}  = char(sub_databaseScanofRef.Tp);
                                ScanOfRef{numel(ScanOfRef)+1}  = char(sub_databaseScanofRef.SequenceName);
                                Parameters_name{numel(Parameters_name)+1}  = strcat(char(sub_databaseOtherScan.SequenceName(i)), '_4thDim', num2str(ii), '_5thDim', num2str(j), '_mean');
                                Parameters_name{numel(Parameters_name)+1}  = strcat(char(sub_databaseOtherScan.SequenceName(i)), '_4thDim', num2str(ii), '_5thDim', num2str(j), '_SD');

                            end
                        end
                    end
                end
            end
            
            
            

            %% load ROI one-by-one
            for i=1:size(sub_databaseROI,1)
                current_cvs_table = table;
                % first load the ROI
                ROI.header = spm_vol([char(sub_databaseROI.Path(i)) char(sub_databaseROI.Filename(i)) '.nii']);
                ROI.data=  read_volume(ROI.header, scan_of_reference.header, 0, 'Axial');
                % the ROI needs to have the same class as the data
                ROI.data = cast( ROI.data, class(output_data));
                % mask the output_data using the current ROI
                current_output_data = output_data.* ROI.data;
                
                % then create a cvs_table_template
                current_cvs_table(1,1:size([{'GroupName'}, {'PatientName'}, {'TimePoint'}, {'ScanOfRef'},  {'ROI_name'},  {'ROI_volume_mm3'},  Parameters_name(:)'  ],2)) = ...
                    [GroupName(1), PatientName(1), TimePoint(1), ScanOfRef(1),  cellstr(sub_databaseROI.SequenceName(i)), num2cell(NaN(1, size(Parameters_name,2)+1))]; % +1 for the ROI volume
                current_cvs_table.Properties.VariableNames = [{'GroupName'}, {'PatientName'}, {'TimePoint'}, {'ScanOfRef'},  {'ROI_name'},  {'ROI_volume_mm3'},  Parameters_name(:)'  ];
                
                
                % update the current_cvs_table with the measures
                parameter_in_vect = reshape(current_output_data, [size(current_output_data,1)*size(current_output_data,2)*size(current_output_data,3) size(current_output_data,4)]);
                parameter_in_vect(parameter_in_vect ==0) = NaN;
                computed_values = [];
                for ii = 1:size(parameter_in_vect,2)
                    computed_values(length(computed_values)+1) = nanmean(parameter_in_vect(:,ii));
                    computed_values(length(computed_values)+1) = nanstd(parameter_in_vect(:,ii));
                end
                %calcuate the ROI volume
                ROI_info = niftiinfo([char(sub_databaseROI.Path(i)) char(sub_databaseROI.Filename(i)) '.nii']);
                voxel_volume = abs(prod(diag(scan_of_reference.header.mat(1:3,1:3))));
                %voxel_volume = prod(ROI_info.raw.pixdim(2:4));
                current_cvs_table.ROI_volume_mm3 = sum(ROI.data(:))*voxel_volume;
                current_cvs_table(:,7:end) = num2cell(computed_values);
                % concatanate the current_cvs_table into one big table
                if ~exist('cvs_table', 'var')
                    cvs_table = current_cvs_table;
                else
                    cvs_table = outerjoin(cvs_table,current_cvs_table,'MergeKeys', true);
                end
            end
        end
    end
end

% generate the output filename
[Path_In, Name_In, ~] = fileparts(files_in.In1{1});
tags = opt.Table_in(strcmp(cellstr(opt.Table_in.Path),[Path_In, filesep]),:);
tags = tags(strcmp(cellstr(tags.Filename), Name_In),:);

tags_out_In = tags;
tags_out_In.SequenceName = categorical(cellstr(opt.Output_filename));
tags_out_In.Filename = categorical(cellstr([char(tags_out_In.Patient), '_', char(tags_out_In.Tp), '_', char(tags_out_In.SequenceName)]));
folder_out_silesep = strsplit(opt.folder_out, filesep);
folder_out_silesep(end) = {'Data_to_export'};
folder_out_silesep(end+1) = {''};

if exist(strjoin(folder_out_silesep, filesep),'dir') ~= 7
    [status, ~, ~] = mkdir(strjoin(folder_out_silesep, filesep));
    if status == false
        error('Cannot create the Data_to_export folder to save the csv files.')
    end
end
f_out = [strjoin(folder_out_silesep, filesep), opt.Output_filename, '.csv'];

% save as cvs
writetable(cvs_table, f_out)
