function [files_in,files_out,opt] = Module_COMBAT_Harmonization(files_in,files_out,opt)

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
    module_option(:,3)   = {'OutputSequenceName','Extension'};
    module_option(:,4)   = {'variability_parameter','Manufacturer_sModel'};
    module_option(:,5)   = {'biological_info','Manufacturer, PatientWeight'};
    module_option(:,6)   = {'output_extension','_harmonized'};
    module_option(:,7)   = {'Parametric_Adjustments','No'};
    module_option(:,8)   = {'AutomaticJobsCreation', 'No'};
    module_option(:,9)   = {'RefInput',2};
    module_option(:,10)   = {'InputToReshape',2};
    module_option(:,11)   = {'Table_in', table()};
    module_option(:,12)   = {'Table_out', table()};
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
    'Combat harmonization : https://github.com/Jfortin1/ComBatHarmonization'
    }'};

    user_parameter(:,2)   = {'Input Scan','1Scan','','', {'SequenceName'},'Mandatory',...
         'The scan on which the others will be opened (same resolution, same referential)'};
    user_parameter(:,3)   = {'Input ROI','1ROI','','', {'SequenceName'},'Optional',...
         'The ROI that will contains the voxels to harmonize.'};
    user_parameter(:,4)   = {'Parameters','','','','','',''};
    user_parameter(:,5)   = {'   .Parametric Adjustments ?','cell',{'No', 'Yes'},'Parametric_Adjustments','','',...
        ''};
    user_parameter(:,6)   = {'   .Output extension','char','','output_extension','','',...
        ''};
    user_parameter(:,7)   = {'   .Variability json field','char','','variability_parameter','','',...
        'For each session, This field will represent the variability to reduce.'};
    user_parameter(:,8)   = {'   .Biological json fields','char','','biological_info','','',...
        'For each session, This field will represent the biological variability to keep.'};


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

if isempty(files_out)
    db = opt.Table_in(opt.Table_in.Type == 'Scan',:);
    for i=1:size(db,1)
        tag1 = db(i,:);
        tag1.IsRaw = categorical(0);   
        tag1.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            tag1.SequenceName = categorical(cellstr(opt.output_extension));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            tag1.SequenceName = categorical(cellstr([char(tag1.SequenceName), opt.output_extension]));
        end
        tag1.Filename = categorical(cellstr([char(tag1.Patient), '_', char(tag1.Tp), '_', char(tag1.SequenceName)]));
        f_out = [char(tag1.Path), char(tag1.Patient), '_', char(tag1.Tp), '_', char(tag1.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; tag1];
    end
end



%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_COMBAT_Harmonization:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

% Variability parameters
Variability_param = opt.variability_parameter;
Variability_param = Variability_param(~isspace(Variability_param));
batch = cell(1, length(files_in.In1));
batch_tmp = cell(1, length(files_in.In1));

% Biological informations
Biological_param = opt.biological_info;
Biological_param = Biological_param(~isspace(Biological_param));
Biological_param = strsplit(Biological_param, ',');
Mod = nan(length(files_in.In1), length(Biological_param));
Mod_tmp = cell(length(files_in.In1), length(Biological_param));

% The first file is the volume of reference, on which we open all the
% others
h_ref = spm_vol(files_in.In1{1});



%The ROI containing the voxels to harmonize.
h_first_roi = spm_vol(files_in.In2{1});
ActualROI = read_volume(h_first_roi, h_ref, 0, 'Axial');
for i=2:length(files_in.In2)
    h_r = spm_vol(files_in.In2{i});
    new_ROI = read_volume(h_r, h_ref, 0, 'Axial');
    ActualROI = ActualROI | new_ROI;
end


DATA = nan(sum(ActualROI(:)), length(files_in.In1));
Size_ref = h_ref.dim;
J = cell(1,length(files_in.In1));
info = cell(1,length(files_in.In1));


for i=1:length(files_in.In1)
    h = spm_vol(files_in.In1{i});
    info{i} = niftiinfo(files_in.In1{i});
    Vol = read_volume(h, h_ref, 0, 'Axial');
    DATA(:,i) = Vol(logical(ActualROI));
    jsonfile = strrep(files_in.In1{i}, '.nii', '.json');
    J{i} = ReadJson(jsonfile);
    J{i} = KeepModuleHistory(J{i}, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
    for j=1:length(Biological_param) %% Pour le moment les valeurs sont considérées comme categorical, pas continues
        Mod_tmp(i,j) = J{i}.(Biological_param{j}).value(1);
    end
    batch_tmp(1,i) = J{i}.(Variability_param).value(1);
end


btch = dummyvar(categorical(batch_tmp)');
Vec = 1:length(unique(batch_tmp));
batch = sum(Vec.*btch,2);

for j=1:length(Biological_param) %% Pour le moment les valeurs sont considérées comme categorical, pas continues
    md = dummyvar(categorical(Mod_tmp(:,j)));
    V = 1:length(unique(Mod_tmp(:,j)));
    Mod(:,j) = sum(V.*md,2);
end

if strcmp(opt.Parametric_Adjustments, 'Yes')
    parametric = 1;
else
    parametric = 0;
end

batch([2,4]) = 2;
Mod([1, end]) = 2;
bayesdata = combat(DATA, batch', Mod, parametric);

for i=1:length(files_out.In1)
    N = nan(Size_ref);
    N(logical(ActualROI)) = bayesdata(:,i);
    D = write_volume(N,h_ref, 'Axial');
    niftiwrite(D, files_out.In1{i}, info{i});
    [path, name, ~] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J{i}, jsonfile)
end
