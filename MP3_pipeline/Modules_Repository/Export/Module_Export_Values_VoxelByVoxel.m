function [files_in,files_out,opt] = Module_Export_Values_VoxelByVoxel(files_in,files_out,opt)

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
    module_option(:,7)   = {'Output_filename','Voxel_values'};
    module_option(:,8)   = {'OutputSequenceName','AllName'};
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
         {'This module aims to export data voxel-by-voxel.'
         '    - inputs: one or several scans and 1 ROI'
         '    - output: cvs file which contains values for each voxel. One line = values for one voxel'}
        };
    user_parameter(:,2)   = {'   .First scan (used as the reference for the resolution)','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the scan that will be use as reference. Every other scan selected below will be reoriented to this reference space'};
    user_parameter(:,3)   = {'   .Other scans (if needed)','XScan','','', {'SequenceName'},'Optional',...
         'Please select evey scan (other than the reference scan) that you would like to export using the resolution of the scan of reference'};
    user_parameter(:,4)   = {'   .ROI','XROI','','',{'SequenceName'},'Mandatory',...
         'Please select one ROI. This ROI will be used to extract values of every voxel in this ROI'};
    user_parameter(:,5)   = {'   .Output filename','char','Voxel_values','Output_filename','', '',''};

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




%% load the scan of reference
scan_of_reference.header = spm_vol(files_in.In1{1});
scan_of_reference.data=  read_volume(scan_of_reference.header, scan_of_reference.header, 0, 'Axial'); %si tout marche, on peut virer cette ligne


output_data = scan_of_reference.data;% .* ROI.data;

voxel_name = {};
Scan_name = {};
TimePoint = {};
PatientName= {};
if length(size(output_data)) == 3
     voxel_name{numel(voxel_name)+1}  = char(opt.Table_in.SequenceName(1));
     Scan_name{numel(Scan_name)+1}  = char(opt.Table_in.SequenceName(1));
     TimePoint{numel(TimePoint)+1}  = char(opt.Table_in.Tp(1));
     PatientName{numel(PatientName)+1}  = char(opt.Table_in.Patient(1));

else
    for j =  1:size(output_data, 5)
        for i = 1:size(output_data, 4)
            voxel_name{numel(voxel_name)+1} = strcat(char(opt.Table_in.SequenceName(1)), '_4thDim', num2str(i), '_5thDim', num2str(j));
            Scan_name{numel(Scan_name)+1} = strcat(char(opt.Table_in.SequenceName(1)));
            TimePoint{numel(TimePoint)+1}  = char(opt.Table_in.Tp(1));
            PatientName{numel(PatientName)+1}  = char(opt.Table_in.Patient(1));
        end
    end
end


%% load all the input scan (other than the scan of reference)
if isfield(files_in, 'In2')
for i=1:length(files_in.In2)
    other_scan(i).header =  spm_vol(files_in.In2{i});
    other_scan(i).data = read_volume(other_scan(i).header, scan_of_reference.header, 0, 'Axial'); 
    
    % save the name of the voxel (head of the column)
    if length(size(other_scan(i).data)) == 3
        voxel_name{numel(voxel_name)+1} = char(opt.Table_in.SequenceName(i+1));
        Scan_name{numel(Scan_name)+1}  = char(opt.Table_in.SequenceName(i+1));
        TimePoint{numel(TimePoint)+1}  = char(opt.Table_in.Tp(i+1));
        PatientName{numel(PatientName)+1}  = char(opt.Table_in.Patient(i+1));
    else
        for j =  1:size(other_scan(i).data, 5)
            for ii = 1:size(other_scan(i).data, 4)
                voxel_name{numel(voxel_name)+1} = strcat(char(opt.Table_in.SequenceName(i+1)), '_4thDim', num2str(ii), '_5thDim', num2str(j));
                Scan_name{numel(Scan_name)+1} = strcat(char(opt.Table_in.SequenceName(i+1)));
                TimePoint{numel(TimePoint)+1}  = char(opt.Table_in.Tp(i+1));
                PatientName{numel(PatientName)+1}  = char(opt.Table_in.Patient(i+1));
            end
        end
    end
    data_to_add = other_scan(i).data; %;.* ROI.data;
     output_data(:,:,:,size(output_data,4)+1:size(output_data,4)+length(other_scan(i).header)) =  reshape(data_to_add, [size(output_data,1), size(output_data,2), size(output_data,3), size(data_to_add,4)*size(data_to_add,5)]); 
   
end
end
% Save data as csv file
% Patient_name / TimePoint / Scan_name / ROI_name /
% Coordonate_x / Coordonate_y / Coordonate_z / Name_of_the_voxel
cvs_table_template = table;


%% load the ROIs one-by-one and create the corresponding table 
% generate the ROI_table
ROI_table = opt.Table_in(strcmp(cellstr(opt.Table_in.Type), 'ROI'),:);
for i=1:numel(files_in.In3)
    % first load the ROI
    ROI.header = spm_vol(files_in.In3{i});
    ROI.data=  read_volume(ROI.header, scan_of_reference.header, 0, 'Axial');
    % mask the output_data using the current ROI
    % but first the ROI needs to have the same class as the data
    ROI.data = cast(ROI.data, class(output_data));
    current_output_data = output_data.* ROI.data;
    % then create a cvs_table_template
    cvs_table_template(1,1:size([{'Ref_Patient_name'}, {'Ref_TimePoint'}, {'Ref_Scan_name'},  {'ROI_name'}, {'Coordonate_x'}, {'Coordonate_y'}, {'Coordonate_z'},  voxel_name(:)'  ],2)) = ...
        [PatientName(1), TimePoint(1), Scan_name(1),  cellstr(ROI_table.SequenceName(i)), NaN, NaN, NaN,  num2cell(NaN(1, size(voxel_name,2)))];
    cvs_table_template.Properties.VariableNames = [{'Ref_Patient_name'}, {'Ref_TimePoint'}, {'Ref_Scan_name'},  {'ROI_name'}, {'Coordonate_x'}, {'Coordonate_y'}, {'Coordonate_z'},  voxel_name(:)'  ];
    % create 3 matrices for the voxel's coordonates
    [coord_x, coord_y, coord_z] = ndgrid(1:size(current_output_data,1), 1:size(current_output_data,2), 1:size(current_output_data,3));
    % reshape the coordonate and the data
    coord_x = reshape(coord_x, [size(current_output_data,1)*size(current_output_data,2)*size(current_output_data,3), 1]);
    coord_y = reshape(coord_y, [size(current_output_data,1)*size(current_output_data,2)*size(current_output_data,3), 1]);
    coord_z = reshape(coord_z, [size(current_output_data,1)*size(current_output_data,2)*size(current_output_data,3), 1]);
    current_output_data = reshape(current_output_data, [size(current_output_data,1)*size(current_output_data,2)*size(current_output_data,3), size(current_output_data,4)]);
    % create a binary variable containing which contains the information of any
    % voxel with non-zeros values
    index_voxel_to_save = any(current_output_data,2);
    % generate a variable with any voxel to save (coordonates and values of
    % each parameters)
    vox_in_table = [coord_x(index_voxel_to_save), coord_y(index_voxel_to_save), coord_z(index_voxel_to_save), current_output_data(index_voxel_to_save,:)];
    
    % design the final table (scan info + voxels values)
    current_cvs_table = repmat(cvs_table_template, [size(vox_in_table,1),1]);
    current_cvs_table(:,5:end) = num2cell(vox_in_table);
    
    % concatanate the current_cvs_table into one big table
    if ~exist('cvs_table', 'var')
        cvs_table = current_cvs_table;
    else
        cvs_table = vertcat(cvs_table,current_cvs_table);
    end
    %cvs_filename = fullfile(strrep(opt.folder_out, '/Tmp', 'Data_to_export'), strcat(cellstr(opt.Output_filename), '.csv'));  
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
f_out = [strjoin(folder_out_silesep, filesep), char(tags_out_In.Patient), '_', char(tags_out_In.Tp), '_', char(tags_out_In.SequenceName), '.csv'];

% save as cvs
writetable(cvs_table, f_out)
