function [files_in,files_out,opt] = Module_DCE_phenomeno(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','Extension'};
    module_option(:,4)   = {'output_filename_ext_AUC','AUC'};
    module_option(:,5)   = {'output_filename_ext_Max','Max'};
    module_option(:,6)   = {'output_filename_ext_pc_rehau','pc_rehau'};
    module_option(:,7)   = {'output_filename_ext_TTP','TTP'};
    module_option(:,8)   = {'Las_dyn_before_bolus', 'Auto'};
    module_option(:,9)   = {'End_analysis',180};
    module_option(:,10)   = {'RefInput', 1};
    module_option(:,11)   = {'InputToReshape',1};
    module_option(:,12)   = {'Table_in', table()};
    module_option(:,13)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    opt.NameOutFiles = {opt.Module_settings.output_filename_ext_AUC, opt.Module_settings.output_filename_ext_Max, opt.Module_settings.output_filename_ext_pc_rehau, opt.Module_settings.output_filename_ext_TTP};
 
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Select one PERF scan as input','1Scan','','',{'SequenceName'},'Mandatory',''};
    user_parameter(:,2)   = {'Parameters','','','','','',''};
    user_parameter(:,3)   = {'   .Output filename extension AUC','char','AUC','output_filename_ext_AUC','','',''};
    user_parameter(:,4)   = {'   .Output filename extension Max','char', 'Max','output_filename_ext_Max','','',''};
    user_parameter(:,5)   = {'   .Output filename extension pc_rehau','char','pc_rehau','output_filename_ext_pc_rehau','','',''};
    user_parameter(:,6)   = {'   .Output filename extension TTP','char','TTP','output_filename_ext_TTP','','',''};
    user_parameter(:,7)   = {'   .Last dyn before bolus','char','Auto','Las_dyn_before_bolus','','',''};
    user_parameter(:,8)   = {'   .End of the analysis (in sec)','numeric','','End_analysis','','',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
%     
%%
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%

opt.NameOutFiles = {opt.output_filename_ext_AUC, opt.output_filename_ext_Max, opt.output_filename_ext_pc_rehau, opt.output_filename_ext_TTP};




if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in;
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), '_', opt.NameOutFiles{i}]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
end







%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_DCE_phenomeno:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, ext_nii, files_in);  
% end
% 
% 
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



N = niftiread(files_in.In1{1});
if size(N,4) == 1
   warndlg([files_in.In1{1} ' is not a 4d image'], 'Warning');
    return
end

N = double(N);
info = niftiinfo(files_in.In1{1});
[path, name, ext] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
J = jsondecode(raw');

%Informations = whos('N');





%% Reshape to process a 4 dimensions scan
% Size = size(N);
% NbDim = length(Size);
% if NbDim >4
%     Dim_To_Merge = Size(4:end);
%     NewDim = prod(Dim_To_Merge);
%     N = reshape(N, Size(1), Size(2), Size(3), NewDim);
%     FilteredImages = reshape(FilteredImages, size(N));
% end
%% Processing



% generate a 4 phenomenological maps from a DCE scan: 
%  i) max enhance 
%  ii)time-to-peak
%  iii) precentage of enhancement 
%  iv) exclued pixels 
%
% this code come from the permeability module 
% additional_parameters correspond to the size of the windows used for the loating mean and the beggining of the bolus 

DCE = N;
%TR = J.RepetitionTime.value/1000; %(From ms to s)
%TR=TR*19;%%Test Clement





% timing = sscanf(scan_acqp('##$PVM_ScanTimeStr=',data.uvascim.image.texte,2),'%dh%dm%ds%dms');
% duree_tot=timing(1)*3600+timing(2)*60+timing(3)+timing(4)/1000; %en seconde
% TR = duree_tot / data.uvascim.image.reco.no_expts;


duree_tot = J.ScanTime.value/1000; %% Conversion from ms to s
TR = duree_tot / size(N, 4);




repetition_nbr = size(DCE,4);
if ~strcmp(opt.Las_dyn_before_bolus, 'Auto')
    debut = str2double(opt.Las_dyn_before_bolus);
    fin = opt.End_analysis;
else
    mean_signal =max(reshape(DCE, [size(DCE,1)*size(DCE,2)*size(DCE,3) size(DCE,4)]),[], 1);
    mean_baseline = nanmean(mean_signal(1:3));
    sd_baseline = nanstd(mean_signal(1:3));
    debut = find(mean_signal>(mean_baseline+2*sd_baseline), 1)-1;
    repetion_time = 1:repetition_nbr;
    repetion_time = repetion_time*TR;
    [~, fin] = min(abs(repetion_time -(debut*TR + opt.End_analysis)));
end
% tmp_data = permute(reco.data, [1 2 4 3]);
data_in_vector = reshape(DCE, [size(DCE,1)*size(DCE,2)*size(DCE,3), size(DCE,4)]);


maxi = NaN(size(data_in_vector,1),1);
ttp  = NaN(size(data_in_vector,1),1);
rehaus  = NaN(size(data_in_vector,1),1);
AUC  = NaN(size(data_in_vector,1),1);


for voxel_nbr=1:size(data_in_vector,1)% 
%     if data_in_vector(voxel_nbr,1) ~= 0 && data_in_vector(voxel_nbr,2) ~= 0
        tmp = squeeze(data_in_vector(voxel_nbr,:)); %
        moy_ss_gd = mean(tmp(1:debut));
        %std_ss_gd = std(tmp(1:debut));
        
        %if moy_ss_gd && sum(tmp>(moy_ss_gd+2*std_ss_gd))>0
        % moving average -------------------------------
        for m=2:repetition_nbr-1
            tmp(m) = (data_in_vector(voxel_nbr,m-1)+data_in_vector(voxel_nbr,m)+data_in_vector(voxel_nbr,m+1))/3;
        end
        % maxi et ttp -------------------------------------
        [maxi(voxel_nbr),ttp(voxel_nbr)] = max(tmp(debut:fin));
        ttp(voxel_nbr) = ttp(voxel_nbr) * TR;
        % rehaussement ------------------------------------
        rehaus(voxel_nbr) = ((max(tmp)-moy_ss_gd)/moy_ss_gd)*100;
        % maxi(voxel_nbr) = data(i,j,k,l,ttp(i,j,k,l)); % Max intensity with
        % the filter
        % AUC Area under the curve ----------------------------
        AUC(voxel_nbr) = 0;
        for m=debut:fin
            AUC(voxel_nbr) = AUC(voxel_nbr) + (data_in_vector(voxel_nbr,m)-moy_ss_gd);
        end
%     else
%         exclus(voxel_nbr) = 1;
%     end
end


% maxi
Max = reshape(maxi,[size(DCE,1),size(DCE,2),size(DCE,3)]);
Max(Max < 0) = -1;
Max(Max > 5000) = -1;
Max(isnan(Max)) = -1;
% ttp (time-to-peak(
TTP =  reshape(ttp,[size(DCE,1),size(DCE,2),size(DCE,3)]);
TTP(TTP < 0) = -1;
TTP(TTP > 5000) = -1;
TTP(isnan(TTP)) = -1;
% rehaus = (intensite-MoySansGd)/MoySansGd
pc_rehau = reshape(rehaus,[size(DCE,1),size(DCE,2),size(DCE,3)]);
pc_rehau(pc_rehau < 0) = -1;
pc_rehau(pc_rehau > 5000) = -1;
pc_rehau(isnan(pc_rehau)) = -1;
% Area under the curve (AUC)
AUC = reshape(AUC,[size(DCE,1),size(DCE,2),size(DCE,3)]);
AUC(AUC < 0) = -1;
AUC(AUC > 5000) = -1;
AUC(isnan(AUC)) = -1;

%%
mapsVar = {AUC, Max, pc_rehau, TTP};
for i=1:length(mapsVar)
    info2 = info;
    info2.Filename = files_out.In1{i};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(mapsVar{i});
    info2.PixelDimensions = info.PixelDimensions(1:length(size(mapsVar{i})));
    info2.ImageSize = size(mapsVar{i});

    %info2.Description = [info.Description, 'Modified by Susceptibility Module'];
    niftiwrite(mapsVar{i}, files_out.In1{i}, info2)

    %% Json Processing
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
    [path, name, ~] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)

end