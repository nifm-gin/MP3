function [files_in,files_out,opt] = Module_T2map(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
  
     % define every option needed to run this module
      %%   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'threshold',5};
    module_option(:,4)   = {'output_filename_ext','T2map'};
    module_option(:,5)   = {'OutputSequenceName','AllName'};
    module_option(:,6)   = {'RefInput',1};
    module_option(:,7)   = {'InputToReshape',1};
    module_option(:,8)   = {'Table_in', table()};
    module_option(:,9)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Select a Multi Spin Echo scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,2)   = {'Parameters','','','','', '',''};
    user_parameter(:,3)   = {'   .Output filename extension','char','T2map','output_filename_ext','', '',''};
    user_parameter(:,4)   = {'   .Threshold','numeric', 5,'threshold','', '',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'T2map'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename_ext};


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
    error('T2map:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end


% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii');  
% end
% 
% %% Building default output names
% if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
%     opt.folder_out = path_nii;    
% end


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


%% load input Nii file
N = niftiread(files_in.In1{1});
% load nifti info
info = niftiinfo(files_in.In1{1});

%% load input JSON file
J = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));

% Get information from the JSON data
EchoTime = J.EchoTime.value;

% reshape the data to a vector matric (speed the fitting process)
%data_to_fit = reshape(double(data.img), [size(data.img,1)*size(data.img, 2)*size(data.img,3) numel(EchoTime)]);
data_to_fit = reshape(double(N), [size(N,1)*size(N, 2)*size(N,3) numel(EchoTime)]);

%% create empty structures
T2map_tmp = NaN(size(data_to_fit,1),1);


% define the threshold and variables
maxim=max(data_to_fit(:)) * opt.threshold/100;
t2init_Cte = EchoTime(1) - EchoTime(end-1);


%init matlabpool
% schd = parcluster();
% poolobj = parpool('local', schd.NumWorkers);

parfor voxel_nbr = 1:size(data_to_fit,1)
    tmp_voxel_data=data_to_fit(voxel_nbr,:);
    if max(tmp_voxel_data(:))>= maxim
        %% fit data
        t2init=(t2init_Cte)/log(tmp_voxel_data(end-1)/tmp_voxel_data(1));
        if t2init<=0 || isnan(t2init)
            t2init=30;
        end
        [aaa, ~,  convergence]=levenbergmarquardt('AB_t2s',EchoTime', abs(tmp_voxel_data),[t2init max(abs(tmp_voxel_data))*1.5]);
        % the algorithm converged
        if convergence == -1
            % to remove when good input data
            if isreal(aaa(1))
                T2map_tmp(voxel_nbr)=aaa(1);
%                 M0map_tmp(voxel_nbr)=aaa(2);
%                 T2_Error_map_tmp(voxel_nbr)=bbb(1);
%                 M0_Error_map_tmp(voxel_nbr)=bbb(2);
            end
        end
    end
end
% delete(poolobj);

% [~,filename,~] = fileparts(MSE_map_filename);

% reshape matrix
OutputImages=reshape(T2map_tmp,[size(N,1) size(N, 2) size(N,3)]);
OutputImages(OutputImages < 0) = -1;
OutputImages(OutputImages > 5000) = -1;
OutputImages(isnan(OutputImages)) = -1;


%% need to update the json structure here before saving it with the T2map
%spm_jsonwrite(strrep(files_out.filename, '.nii', '.json'), data.json);

% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);
info2.Description = [info.Description, 'Modified by T2map Module'];

% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);

% % so far copy the .json file of the first input
% copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))

%% Json Processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)


