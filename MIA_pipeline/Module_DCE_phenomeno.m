function [files_in,files_out,opt] = Module_DCE_phenomeno(files_in,files_out,opt)
% This is a template file for "brick" functions in NIAK.
%
% SYNTAX:
% [IN,OUT,OPT] = PSOM_TEMPLATE_BRICK(IN,OUT,OPT)
%
% _________________________________________________________________________
% INPUTS:
%
% IN        
%   (string) a file name of a 3D+t fMRI dataset .
%
% OUT
%   (structure) with the following fields:
%       flag_test
%   CORRECTED_DATA
%       (string, default <BASE NAME FMRI>_c.<EXT>) File name for processed 
%       data.
%       If OUT is an empty string, the name of the outputs will be 
%       the same as the inputs, with a '_c' suffix added at the end.
%
%   MASK
%       (string, default <BASE NAME FMRI>_mask.<EXT>) File name for a mask 
%       of the data. If OUT is an empty string, the name of the 
%       outputs will be the same as the inputs, with a '_mask' suffix added 
%       at the end.
%
% OPT           
%   (structure) with the following fields.  
%
%   TYPE_CORRECTION       
%      (string, default 'mean_var') possible values :
%      'none' : no correction at all                       
%      'mean' : correction to zero mean.
%      'mean_var' : correction to zero mean and unit variance
%      'mean_var2' : same as 'mean_var' but slower, yet does not use as 
%      much memory).
%
%   FOLDER_OUT 
%      (string, default: path of IN) If present, all default outputs 
%      will be created in the folder FOLDER_OUT. The folder needs to be 
%      created beforehand.
%
%   FLAG_VERBOSE 
%      (boolean, default 1) if the flag is 1, then the function prints 
%      some infos during the processing.
%
%   FLAG_TEST 
%      (boolean, default 0) if FLAG_TEST equals 1, the brick does not do 
%      anything but update the default values in IN, OUT and OPT.
%           
% _________________________________________________________________________
% OUTPUTS:
%
% IN, OUT, OPT: same as inputs but updated with default values.
%              
% _________________________________________________________________________
% SEE ALSO:
% NIAK_CORRECT_MEAN_VAR
%
% _________________________________________________________________________
% COMMENTS:
%
% _________________________________________________________________________
% Copyright (c) <NAME>, <INSTITUTION>, <START DATE>-<END DATE>.
% Maintainer : <EMAIL ADDRESS>
% See licensing information in the code.
% Keywords : PSOM, documentation, template, brick

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
    % define every option needed to run this module
    %fields   = {'Type', 'HSize', 'Sigma', 'flag_test' , 'folder_out', 'output_filename_ext'};
%     fields   = {'RefInput', 'InputToReshape', 'NbInput', 'NbOutput', 'folder_out', 'flag_test', 'OutputSequenceName', 'output_filename_ext_CBV', 'output_filename_ext_CBF', 'output_filename_ext_MTT', 'output_filename_ext_TMAX', 'output_filename_ext_TTP', 'output_filename_ext_T0'};
%     defaults = {1,1,1,6,'', true, 'AllName','CBV', 'CBF', 'MTT', 'TMAX', 'TTP', 'T0'};
%     opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
%     
%     %opt.NameOutFiles = {opt.output_filename_ext_CBV, output_filename_ext_CBF, output_filename_ext_MTT, output_filename_ext_TMAX, output_filename_ext_TTP, output_filename_ext_T0};
%     
%     % list of everything displayed to the user associated to their 'type'
%     user_parameter_list = {'Select one PERF scan as input'; 'Parameters'; '   .Output filename extension CBV'; '   .Output filename extension CBF'; '   .Output filename extension MTT'; '   .Output filename extension TMAX'; '   .Output filename extension TTP'; '   .Output filename extension T0'};%; ''; ''};
%     user_parameter_type = {'1Scan'; ''; 'char'; 'char'; 'char'; 'char'; 'char'; 'char'};%; 'logical'; 'char'};
%     parameter_default = {''; ''; 'CBV'; 'CBF'; 'MTT'; 'TMAX'; 'TTP'; 'T0'};%; '1'; ''};
%     psom_parameter_list = {''; ''; 'output_filename_ext_CBV'; 'output_filename_ext_CBF'; 'output_filename_ext_MTT'; 'output_filename_ext_TMAX'; 'output_filename_ext_TTP'; 'output_filename_ext_T0'};%; 'flag_test'; 'folder_out'};
%     scans_input_DOF = {{'SequenceName'}; ''; ''; '';''; '';'';''};
%     VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF'};
%     %opt.table = table(categorical(user_parameter_list), categorical(user_parameter_type), categorical(parameter_default), categorical(psom_parameter_list), 'VariableNames', VariableNames);
%     opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, scans_input_DOF,'VariableNames', VariableNames);
%     %opt.parameter_link_psom = {'output_filename_ext', '   .Output filename extension'; 'Type', '   .Type'; 'HSize','   .HSize'; 'Sigma', '   .Sigma'};
%     opt.NameOutFiles = {'CBV', 'CBF', 'MTT', 'TMAX', 'TTP', 'T0'};
%% Benjamin Modifications
     % define every option needed to run this module
%       %%   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
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
% 
%  %% list of everything displayed to the user associated to their 'type'
%      % --> user_parameter(1,:) = user_parameter_list
%      % --> user_parameter(2,:) = user_parameter_type
%      % --> user_parameter(3,:) = parameter_default
%      % --> user_parameter(4,:) = psom_parameter_list
%      % --> user_parameter(5,:) = Help : text data which describe the parameter (it
%      % will be display to help the user)
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
    error('Smoothing:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
if ~ischar(files_in.In1{1}) 
    error('files in should be a char');
end

[path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
if ~strcmp(ext_nii, '.nii')
     error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, ext_nii, files_in);  
end

if isfield(opt,'threshold') && (~isnumeric(opt.threshold))
    opt.threshold = str2double(opt.threshold);
    if isnan(opt.threshold)
        disp('The threshold used was not a number')
        return
    end
end

%% Options
% %fields   = {'Type', 'HSize', 'Sigma', 'flag_test' , 'folder_out', 'output_filename_ext'};
% %defaults = {5, 3, 1, true, '', '_Smooth'};
% fields   = {'folder_out', 'flag_test', 'output_filename_ext', 'Type', 'HSize', 'Sigma'};
% defaults = {'', false, '_Smooth', {'gaussian'}, 3, 1};
% if nargin < 5
%     opt = psom_struct_defaults(struct(),fields,defaults);
% else
%     opt = psom_struct_defaults(opt,fields,defaults);
% end



%% Check the output files structure
%fields    = {'filename'};
%defaults  = {'output defaults' };
%files_out = psom_struct_defaults(files_out,fields,defaults);

%% Building default output names
if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
    opt.folder_out = path_nii;    
end

%if isempty(files_out.filename)
%    files_out.filename = cat(2,opt.folder_out,filesep,name_nii,opt.output_filename_ext,ext_nii);
%end
% 
% if strcmp(files_out, '')
%     files_out.In1{1} = cat(2,opt.folder_out,filesep,name_nii, '_',opt.output_filename_ext_CBV,ext_nii);
%     files_out.In1{2} = cat(2,opt.folder_out,filesep,name_nii, '_',opt.output_filename_ext_CBF,ext_nii);
%     files_out.In1{3} = cat(2,opt.folder_out,filesep,name_nii, '_',opt.output_filename_ext_MTT,ext_nii);
%     files_out.In1{4} = cat(2,opt.folder_out,filesep,name_nii, '_',opt.output_filename_ext_TMAX,ext_nii);
%     files_out.In1{5} = cat(2,opt.folder_out,filesep,name_nii, '_',opt.output_filename_ext_TTP,ext_nii);
%     files_out.In1{6} = cat(2,opt.folder_out,filesep,name_nii, '_',opt.output_filename_ext_T0,ext_nii);
% end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
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
J = jsondecode(raw);

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
TR = J.RepetitionTime.value/1000; %(From ms to s)
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
    %info2.raw.datatype = 16;
    %info2.BitsPerPixel = 32;
    info2.ImageSize = size(mapsVar{i});
    %info2.raw.dim(1) = 3;
    %info2.raw.dim(5) = 1;
    %info2.raw.bitpix = 32;
    %info2.raw = struct();
    info2.Description = [info.Description, 'Modified by Susceptibility Module'];
    niftiwrite(mapsVar{i}, files_out.In1{i}, info2)
    
    %info2 = info;
    %info2.Filename = files_out.In1{i};
    %info2.Filemoddate = char(datetime('now'));
    %info2.Description = [info.Description, 'Modified by Susceptibility Module'];

    %eval(['niftiwrite(',maps{i},', files_out.In1{i})'])
    JMod = jsonencode(J);
    [path, name, ext] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    fidmod = fopen(jsonfile, 'w');
    fwrite(fidmod, JMod, 'uint8');
    fclose(fidmod);
end



%% Initial function
% 
% % generate a 4 phenomenological maps from a DCE scan: 
% %  i) max enhance 
% %  ii)time-to-peak
% %  iii) precentage of enhancement 
% %  iv) exclued pixels 
% %
% % this code come from the permeability module 
% % additional_parameters correspond to the size of the windows used for the loating mean and the beggining of the bolus 
% 
% % load the DCE file
% fid=fopen(filename_DCE ,'r');
% if fid>0
%     fclose(fid);
%     data = load(filename_DCE);
%     reco = data.uvascim.image.reco;
% else
%      warning_text = sprintf('##$ Can not calculate the DCE-phenomeno map because there is\n##$ Somthing wrong with the data\n##$DCE=%s\n##$',...
%         filename_DCE);
%     msgbox(warning_text, 'DCE-phenomeno map warning') ;
%     maxi_struc = []; ttp_struc = []; rehaus_struc = []; AUC_struc = []; exclus_struc = [];
%     return
% end
% 
% 
% % check if human data (par/rec data)
% if strcmp(data.uvascim.image.filename(end-3:end), '.PAR')
%     TR = data.uvascim.image.acq.tr;
% else
%     timing = sscanf(scan_acqp('##$PVM_ScanTimeStr=',data.uvascim.image.texte,2),'%dh%dm%ds%dms');
%     duree_tot=timing(1)*3600+timing(2)*60+timing(3)+timing(4)/1000; %en seconde
%     TR = duree_tot / data.uvascim.image.reco.no_expts;
% end
% repetition_nbr = size(reco.data,5);
% if ~strcmp(add_parameters{:}(1), 'Auto')
%     debut = str2double(add_parameters{:}(1));
%     fin = str2double(add_parameters{:}(2));
% else
%     mean_signal =max(reshape(data.uvascim.image.reco.data, [size(data.uvascim.image.reco.data,1)*size(data.uvascim.image.reco.data,2)*size(data.uvascim.image.reco.data,4) size(data.uvascim.image.reco.data,5)]),[], 1);
%     mean_baseline = nanmean(mean_signal(1:3));
%     sd_baseline = nanstd(mean_signal(1:3));
%     debut = find(mean_signal>(mean_baseline+2*sd_baseline), 1)-1;
%     repetion_time = 1:repetition_nbr;
%     repetion_time = repetion_time*TR;
%     [~, fin] = min(abs(repetion_time -(debut*TR + str2double(add_parameters{:}(2)))));
% end
% % tmp_data = permute(reco.data, [1 2 4 3]);
% data_in_vector = reshape(reco.data, [size(reco.data,1)*size(reco.data,2)*size(reco.data,4), size(reco.data,5)]);
% 
% 
% maxi = NaN(size(data_in_vector,1),1);
% ttp  = NaN(size(data_in_vector,1),1);
% rehaus  = NaN(size(data_in_vector,1),1);
% AUC  = NaN(size(data_in_vector,1),1);
% exclus  = NaN(size(data_in_vector,1),1);
% 
% 
% for voxel_nbr=1:size(data_in_vector,1)% 
% %     if data_in_vector(voxel_nbr,1) ~= 0 && data_in_vector(voxel_nbr,2) ~= 0
%         tmp = squeeze(data_in_vector(voxel_nbr,:)); %
%         moy_ss_gd = mean(tmp(1:debut));
%         std_ss_gd = std(tmp(1:debut));
%         
%         %if moy_ss_gd && sum(tmp>(moy_ss_gd+2*std_ss_gd))>0
%         % moving average -------------------------------
%         for m=2:repetition_nbr-1
%             tmp(m) = (data_in_vector(voxel_nbr,m-1)+data_in_vector(voxel_nbr,m)+data_in_vector(voxel_nbr,m+1))/3;
%         end
%         % maxi et ttp -------------------------------------
%         [maxi(voxel_nbr),ttp(voxel_nbr)] = max(tmp(debut:fin));
%         ttp(voxel_nbr) = ttp(voxel_nbr) * TR;
%         % rehaussement ------------------------------------
%         rehaus(voxel_nbr) = ((max(tmp)-moy_ss_gd)/moy_ss_gd)*100;
%         % maxi(voxel_nbr) = data(i,j,k,l,ttp(i,j,k,l)); % Max intensity with
%         % the filter
%         
%         % AUC Area under the curve ----------------------------
%         AUC(voxel_nbr) = 0;
%         for m=debut:fin
%             AUC(voxel_nbr) = AUC(voxel_nbr) + (data_in_vector(voxel_nbr,m)-moy_ss_gd);
%         end
% %     else
% %         exclus(voxel_nbr) = 1;
% %     end
% end
% 
% % maxi = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% % ttp = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% % rehaus = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% % exclus = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% % AUC = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% % 
% % h = waitbar(0, 'DCE-phenomeno maps in progress...');
% % steps_tot= reco.no_slices*reco.no_samples*reco.no_views;
% % current_step = 1;
% % 
% % 
% % for i=1:size(reco.data,1) % lignes
% %     for j=1:size(reco.data,2) % colonnes
% %         for k=1:size(reco.data,3) % echo
% %             for l=1:size(reco.data,4) % coupes
% %                 waitbar(current_step/steps_tot, h)
% %                 tmp = squeeze(reco.data(i,j,k,l,:)); % donn�es temporelles
% %                 if tmp(1) ~= 0 && tmp(2) ~= 0
% %                     moy_ss_gd = mean(tmp(1:debut));
% %                     std_ss_gd = std(tmp(1:debut));
% %                     %if moy_ss_gd && sum(tmp>(moy_ss_gd+2*std_ss_gd))>0
% %                     % moyenne flottante -------------------------------
% %                     for m=2:size(tmp)-1
% %                         tmp(m) = (reco.data(i,j,k,l,m-1)+reco.data(i,j,k,l,m)+reco.data(i,j,k,l,m+1))/3;
% %                     end
% %                     % maxi et ttp -------------------------------------
% %                     [maxi(i,j,k,l),ttp(i,j,k,l)] = max(tmp(debut:fin));
% %                     ttp(i,j,k,l) = ttp(i,j,k,l) * TR;
% %                     % rehaussement ------------------------------------
% %                     rehaus(i,j,k,l) = ((max(tmp)-moy_ss_gd)/moy_ss_gd)*100;
% %                     % maxi(i,j,k,l) = data(i,j,k,l,ttp(i,j,k,l)); % pour r�cup�rer l'intensit� max et pas celle filtr�e
% %                     % AUC sous la courbe ----------------------------
% %                     AUC(i,j,k,l) = 0;
% %                     for m=debut:fin
% %                         AUC(i,j,k,l) = AUC(i,j,k,l) + (reco.data(i,j,k,l,m)-moy_ss_gd);
% %                     end
% %                     %else
% %                     %   exclus(i,j,k,l) = 1;
% %                     %end
% %                 else
% %                     exclus(i,j,k,l) = 1;
% %                 end
% %                 current_step = current_step+1;
% %             end
% %         end
% %     end
% % end
% 
% % save imformation
% 
% 
% tempstruct=data.uvascim.image;
% tempstruct.reco.data = zeros(size(reco.data,1),...
%     size(reco.data,2),1,size(reco.data,4),1,'single');
% tempstruct.reco.no_echoes = 1;
% tempstruct.reco.displayedecho = 1;
% tempstruct.reco.displayedslice = 1;
% tempstruct.reco.displayedexpt = 1;
% tempstruct.reco.no_expts = 1;
% tempstruct.reco.echotime = NaN('single');
% 
% %Adapt the fov offsets and orientation infos
% tmpfov=reco.fov_offsets(:,1,:,1);
% tmpori=reco.fov_orientation(:,1,:,1);
% tmpphaselabel=reco.phaselabel(1,:,1);
% tempstruct.reco.fov_offsets = [];
% tempstruct.reco.fov_orientation = [];
% tempstruct.reco.scaling_factor = [];
% tempstruct.reco.scaling_offset = [];
% tempstruct.reco.label = {};
% tempstruct.reco.phaselabel = {};
% 
% for m_slice=1:reco.no_slices,
%     tempstruct.reco.fov_offsets(:,:,m_slice)=squeeze(tmpfov(:,m_slice));
%     tempstruct.reco.fov_orientation(:,:,m_slice)=squeeze(tmpori(:,m_slice));
%     tempstruct.reco.label(1,m_slice)=reco.label(1,m_slice,1);
%     tempstruct.reco.phaselabel(1,m_slice)=reco.phaselabel(1,m_slice,1);
%     tempstruct.reco.scaling_factor(1,m_slice) = reco.scaling_factor(1,m_slice,1);
%     tempstruct.reco.scaling_offset(1,m_slice) = reco.scaling_offset(1,m_slice,1);
% end
% tempstruct.reco.fov_phase_orientation=reco.fov_phase_orientation(:,:,1);
% 
% 
% %delete useless field
% tempstruct.reco = rmfield(tempstruct.reco, 'angAP');
% tempstruct.reco = rmfield(tempstruct.reco, 'angFH');
% tempstruct.reco = rmfield(tempstruct.reco, 'angRL');
% tempstruct.reco = rmfield(tempstruct.reco, 'angulation');
% tempstruct.reco = rmfield(tempstruct.reco, 'echotime');
% tempstruct.reco = rmfield(tempstruct.reco, 'iminfos');
% tempstruct.reco = rmfield(tempstruct.reco, 'reco_meth');
% 
% ParamConfig=sprintf('##$QuantifMethod=phenomenological maps from a DCE scan\n##$Beggining of the injection=%s\n##$End of the analysis=%s\n##$Raw scan used=%s',...
%         num2str(debut),...
%         num2str(fin),...
%         filename_DCE);
% tempstruct.reco.paramQuantif = ParamConfig;
% 
% tempstruct.reco = rmfield(tempstruct.reco, 'data');
% maxi_struc = tempstruct;
% ttp_struc = tempstruct;
% rehaus_struc = tempstruct;
% AUC_struc = tempstruct;
% exclus_struc = tempstruct;
% 
% % maxi
% maxi_struc.reco.data = reshape(maxi,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
% maxi_struc.reco.unit = {'a.u.'};
% maxi_struc.reco.echo_label = {'DCE-Max'};
% % ttp (time-to-peak(
% ttp_struc.reco.data =  reshape(ttp,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
% ttp_struc.reco.unit = {'s'};
% ttp_struc.reco.echo_label = {'DCE-ttp'};
% % rehaus = (intensite-MoySansGd)/MoySansGd
% rehaus_struc.reco.data = reshape(rehaus,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
% rehaus_struc.reco.unit = {'%'};
% rehaus_struc.reco.echo_label = {'DCE-rehaus'};
% % Area under the curve (AUC)
% AUC_struc.reco.data = reshape(AUC,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]); 
% AUC_struc.reco.unit = {'a.u.'};
% AUC_struc.reco.echo_label = {'DCE-AUC'};
% % Exclued voxels
% exclus_struc.reco.data = reshape(exclus,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
% exclus_struc.reco.unit = {''};
% exclus_struc.reco.echo_label = {'DCE-exclus'};
% 
% % close(h);
