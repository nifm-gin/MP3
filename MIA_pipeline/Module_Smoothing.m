function [files_in,files_out,opt] = Module_Smoothing(files_in,files_out,opt)
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
%     %fields   = {'Type', 'HSize', 'Sigma', 'flag_test' , 'folder_out', 'output_filename_ext'};
%     fields   = {'RefInput', 'InputToReshape', 'NbInput', 'NbOutput', 'folder_out', 'flag_test', 'output_filename_ext', 'OutputSequenceName','Type', 'HSize', 'Sigma'};
%     defaults = {1,1,1,1,'', true, '_Smooth', 'Extension', 'gaussian', '3', '1'};
%     opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
%     
%     % list of everything displayed to the user associated to their 'type'
%     user_parameter_list = {'Select one scan or more as input'; 'Parameters'; '   .Output filename extension'; '   .Type';  '   .HSize';  '   .Sigma'};%; ''; ''};
%     user_parameter_type = {'XScan'; ''; 'char'; 'cell'; 'numeric'; 'numeric'};%; 'logical'; 'char'};
%     parameter_default = {''; ''; '_Smooth'; {'gaussian'}; '3'; '1'};%; '1'; ''};
%     psom_parameter_list = {''; ''; 'output_filename_ext'; 'Type'; 'HSize'; 'Sigma'};%; 'flag_test'; 'folder_out'};
%     scans_input_DOF = {{'SequenceName'}; ''; ''; ''; ''; ''};
%     VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF'};
%     %opt.table = table(categorical(user_parameter_list), categorical(user_parameter_type), categorical(parameter_default), categorical(psom_parameter_list), 'VariableNames', VariableNames);
%     opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, scans_input_DOF, 'VariableNames', VariableNames);
%     %opt.parameter_link_psom = {'output_filename_ext', '   .Output filename extension'; 'Type', '   .Type'; 'HSize','   .HSize'; 'Sigma', '   .Sigma'};

%% Benjamin Modifications
%       %%   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'output_filename_ext','Smooth'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'Type','gaussian'};
    module_option(:,6)   = {'HSize',3};
    module_option(:,7)   = {'Sigma',1};
    module_option(:,8)   = {'RefInput',1};
    module_option(:,9)   = {'InputToReshape',1};
    module_option(:,10)   = {'NbInput',1};
    module_option(:,11)   = {'NbOutput',1};
    module_option(:,12)   = {'Table_in', table()};
    module_option(:,13)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
% 
%     %% list of everything displayed to the user associated to their 'type'
%      % --> user_parameter(1,:) = user_parameter_list
%      % --> user_parameter(2,:) = user_parameter_type
%      % --> user_parameter(3,:) = parameter_default
%      % --> user_parameter(4,:) = psom_parameter_list
%      % --> user_parameter(5,:) = Help : text data which describe the parameter (it
%      % will be display to help the user)
    user_parameter(:,1)   = {'Select one scan or more as input','XScan','','',{'SequenceName'}, ''};
    user_parameter(:,2)   = {'Parameters','','','Execution_Mode','', ''};
    user_parameter(:,3)   = {'   .Output filename extension','char','Smooth','output_filename_ext','', ''};
    user_parameter(:,4)   = {'   .Type','cell', {'gaussian'},'Type','', ''};
    user_parameter(:,5)   = {'   .HSize','numeric',3,'HSize','', ''};
    user_parameter(:,6)   = {'   .Sigma','numeric',1,'Sigma','', ''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', 'VariableNames', VariableNames);
%%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%


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
    opt.threshold = opt.threshold;
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

if strcmp(files_out, '')
    files_out.In1 = {cat(2,opt.folder_out,filesep,name_nii,'_',opt.output_filename_ext,ext_nii)};
end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N = niftiread(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
[path, name, ext] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
J = jsondecode(raw);

Informations = whos('N');
FilteredImages = zeros(size(N), Informations.class);
Size = size(N);
NbDim = length(Size);


if NbDim >4
    Dim_To_Merge = Size(4:end);
    NewDim = prod(Dim_To_Merge);
    NewN = reshape(N, Size(1), Size(2), Size(3), NewDim);
else
    NewN = N;
end



h = fspecial(opt.Type, opt.HSize, opt.Sigma);
for i=1:size(NewN,3)
    for j=1:size(NewN,4)
        FilteredImages(:,:,i,j) = imfilter(NewN(:,:,i,j), h, 'replicate');
    end
end



NewFilteredImages = reshape(FilteredImages, Size);



info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Description = [info.Description, 'Modified by Smoothing Module'];

niftiwrite(NewFilteredImages, files_out.In1{1}, info2)
JMod = jsonencode(J);
[path, name, ext] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
fidmod = fopen(jsonfile, 'w');
fwrite(fidmod, JMod, 'uint8');
fclose(fidmod);
%opt.Success = 1;



% %% load input Nii file
% data.hdr = spm_vol(files_in{1});
% data.img = spm_read_vols(data.hdr);
% 
% %% load input JSON fil
% data.json = spm_jsonread(files_in{2});
% 
% % Get information from the JSON data
% EchoTime = data.json.EchoTime;
% 
% % reshape the data to a vector matric (speed the fitting process)
% data_to_fit = reshape(double(data.img), [size(data.img,1)*size(data.img, 2)*size(data.img,3) numel(EchoTime)]);
% 
% %% create empty structures
% Smooth_tmp = NaN(size(data_to_fit,1),1);
% % M0map_tmp = NaN(size(data_to_fit,1),1);
% % T2_Error_map_tmp = NaN(size(data_to_fit,1),1);
% % M0_Error_map_tmp = NaN(size(data_to_fit,1),1);
% 
% 
% 
% %% Old Code
% 
% if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%     % find the parameter position for this patient (if exist)
%     parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
%     map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
%     if ~isempty(map_filename) && sum(parameter_nbr) == 1
%         % load the file
%         fid=fopen(map_filename{:} ,'r');
%         if fid>0
%             fclose(fid);
%             load(map_filename{:});
%         else
%             logbook{numel(logbook)+1} = sprintf('Somthing wrong with the data\n##$%s\n##$',...
%                 map_filename);
%             %continue
%         end
%         h = fspecial(add_parameters{1},str2double(add_parameters{2}),str2double(add_parameters{3}));
%         [~, ~, E, Z] = size(uvascim.image.reco.data);
%         for s=1:Z
%             for e=1:E
%                 uvascim.image.reco.data(:,:,e,s) = imfilter(uvascim.image.reco.data(:,:,e,s),h,'replicate');
%             end
%         end
%         % save info
%         if isfield(uvascim.image.reco, 'paramQuantif')
%             ParamConfig=sprintf('##$Original file name=%s\n##$Original file info=%s\n##$MathOperation=Smoothing type %s\n##$Hsize=%s\n##$Sigma=%s\n',...
%                 map_filename{:},...
%                 uvascim.image.reco.paramQuantif,...
%                 add_parameters{1},...
%                 add_parameters{2},...
%                 add_parameters{3});
%         else
%             ParamConfig=sprintf('##$Original file name=%s\n##$MathOperation=Smoothing type %s\n##$Hsize=%s\n##$Sigma=%s\n',...
%                 map_filename{:},...
%                 add_parameters{1},...
%                 add_parameters{2},...
%                 add_parameters{3});
%         end
%         uvascim.image.reco.paramQuantif = ParamConfig;
%         
%         % save file
%         [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
%         scan_name = [name file_name_extension ext];
%         save(fullfile([path, scan_name]),'uvascim');
%         
%         % add to the database and update handles.parameters and handles.clips
%         handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
%         guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
%         MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
%         clear uvascim;
%         logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
%         
%         
%     else
%         logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped'];
%     end
% else
%     logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
% end
% %%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % define the threshold and variables
% maxim=max(data_to_fit(:)) * opt.threshold/100;
% 
% t2init_Cte = EchoTime(1) - EchoTime(end-1);
% 
% 
% %init matlabpool
% % schd = parcluster();
% % poolobj = parpool('local', schd.NumWorkers);
% 
% parfor voxel_nbr = 1:size(data_to_fit,1)
%     tmp_voxel_data=data_to_fit(voxel_nbr,:);
%     if max(tmp_voxel_data(:))>= maxim
%         %% fit data
%         t2init=(t2init_Cte)/log(tmp_voxel_data(end-1)/tmp_voxel_data(1));
%         if t2init<=0 || isnan(t2init)
%             t2init=30;
%         end
%         [aaa, ~,  convergence]=levenbergmarquardt('AB_t2s',EchoTime', abs(tmp_voxel_data),[t2init max(abs(tmp_voxel_data))*1.5]);
%         % the algorithm converged
%         if convergence == -1
%             % to remove when good input data
%             if isreal(aaa(1))
%                 T2map_tmp(voxel_nbr)=aaa(1);
% %                 M0map_tmp(voxel_nbr)=aaa(2);
% %                 T2_Error_map_tmp(voxel_nbr)=bbb(1);
% %                 M0_Error_map_tmp(voxel_nbr)=bbb(2);
%             end
%         end
%     end
% end
% % delete(poolobj);
% 
% % [~,filename,~] = fileparts(MSE_map_filename);
% 
% % reshape matrix
% T2map.img=reshape(T2map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% T2map.img(T2map.img < 0) = -1;
% T2map.img(T2map.img > 5000) = -1;
% T2map.img(isnan(T2map.img)) = -1;
% % M0map.img=reshape(M0map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% % T2_Error_map.img=reshape(T2_Error_map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% % M0_Error_map.img=reshape(M0_Error_map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% % save the T2 map
% %% if spm_function is used
% % T2map.hdr = spm_vol([files_in{1}, ', 1']);
% % T2map.hdr.fname = files_out.filename;
% % spm_write_vol(T2map.hdr, T2map.img);
% 
% %% if matlab function is used
% T2map.hdr = niftiinfo(files_in{1});
% T2map.hdr = update_nifti_hdr(T2map.hdr, T2map.img, files_out.filename);
% niftiwrite(single(T2map.img), files_out.filename, T2map.hdr)
% 
% %% need to update the json structure here before saving it with the T2map
% spm_jsonwrite(strrep(files_out.filename, '.nii', '.json'), data.json);
% 
% % save the M0map map
% M0map.hdr = spm_vol([MSE_map_filename, ', 1']);
% M0map.hdr.fname = char(strcat(filename, '-M0map.nii'));
% M0map.img(isnan(M0map.img)) = -1;
% spm_write_vol(M0map.hdr, M0map.img);
% 
% % save the T2_Error_map 
% T2_Error_map.hdr = spm_vol([MSE_map_filename, ', 1']);
% T2_Error_map.hdr.fname = char(strcat(filename, '-T2_Error.nii'));
% T2_Error_map.img(T2_Error_map.img < 0) = -1;
% T2_Error_map.img(T2_Error_map.img > 50) = -1;
% T2_Error_map.img(isnan(T2_Error_map.img)) = -1;
% spm_write_vol(T2_Error_map.hdr, T2_Error_map.img);
% 
% % save the M0_Error_map map
% M0_Error_map.hdr = spm_vol([MSE_map_filename, ', 1']);
% M0_Error_map.hdr.fname = char(strcat(filename, '-M0_Error.nii'));
% M0_Error_map.img(isnan(M0_Error_map.img)) = -1;
% spm_write_vol(M0_Error_map.hdr, M0_Error_map.img);
