function [files_in,files_out,opt] = Module_Susceptibility(files_in,files_out,opt)
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
    module_option(:,3)   = {'OutputSequenceName','AllName'};
    module_option(:,4)   = {'output_filename_ext_CBV','CBV'};
    module_option(:,5)   = {'output_filename_ext_CBF','CBF'};
    module_option(:,6)   = {'output_filename_ext_MTT','MTT'};
    module_option(:,7)   = {'output_filename_ext_TMAX','TMAX'};
    module_option(:,8)   = {'output_filename_ext_TTP','TTP'};
    module_option(:,9)   = {'output_filename_ext_T0','T0'}; 
    module_option(:,10)   = {'RefInput', 1};
    module_option(:,11)   = {'InputToReshape',1};
    module_option(:,12)   = {'Table_in', table()};
    module_option(:,13)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    opt.NameOutFiles = {opt.Module_settings.output_filename_ext_CBV, opt.Module_settings.output_filename_ext_CBF, opt.Module_settings.output_filename_ext_MTT, opt.Module_settings.output_filename_ext_TMAX, opt.Module_settings.output_filename_ext_TTP, opt.Module_settings.output_filename_ext_T0};
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
    user_parameter(:,3)   = {'   .Output filename extension CBV','char','CBV','output_filename_ext_CBV','','',''};
    user_parameter(:,4)   = {'   .Output filename extension CBF','char', 'CBF','output_filename_ext_CBF','','',''};
    user_parameter(:,5)   = {'   .Output filename extension MTT','char','MTT','output_filename_ext_MTT','','',''};
    user_parameter(:,6)   = {'   .Output filename extension TMAX','char','TMAX','output_filename_ext_TMAX','','',''};
    user_parameter(:,7)   = {'   .Output filename extension TTP','char','TTP','output_filename_ext_TTP','','',''};
    user_parameter(:,8)   = {'   .Output filename extension T0','char','T0','output_filename_ext_T0','','',''};
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

opt.NameOutFiles = {opt.output_filename_ext_CBV, opt.output_filename_ext_CBF, opt.output_filename_ext_MTT, opt.output_filename_ext_TMAX, opt.output_filename_ext_TTP, opt.output_filename_ext_T0};




if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in;
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), opt.NameOutFiles{i}]));
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
meansignal = mean(squeeze(N),4);
volume_mask = meansignal>max(N(:))*0.01;
[aif,scores] = extraction_aif_volume(squeeze(N),volume_mask);

if sum(cell2num(scores(:,5))) > 10
    error('No computation because the AIF is not good enough');
else
    [~,~,~,TMAX,TTP,T0, CBV,CBF,MTT,~,~,~,~] = deconvolution_perfusion_gui(aif,squeeze(N),J.RepetitionTime.value(1)*10^(-3),J.EchoTime.value*10^(-3));
    %maps = {'CBV','CBF','MTT','TMAX','TTP','T0'};
    mapsVar = {CBV, CBF, MTT, TMAX, TTP, T0};

end


%% Reshape to the input scan size
%FilteredImages = reshape(FilteredImages, Size);


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



% % check if DSC analysis scan exist already
% if  sum(strcmp(['CBV' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%     match = strcmp('DSC', [handles.ready_to_go.couple.map_wanted]') == 1;
%     map_name_selected = handles.ready_to_go.couple(match).map_name;
%     add_parameters = handles.ready_to_go.couple(match).additional_information;
%     DSC_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
%     filename_DSC= strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DSC_nbr));
%     
%     if ~isempty(filename_DSC) && sum(DSC_nbr) == 1
%         
%         fid=fopen(filename_DSC{:} ,'r');
%         if fid>0
%             fclose(fid);
%             load(filename_DSC{:});
%             data = uvascim.image;
%             data.reco.data= double(data.reco.data);
%         end
%         if ~strcmp(add_parameters{1}, 'None')
%             if ~isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
%                 Mask_nbr = strcmp(add_parameters{:}(1), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
%                 Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
%                 if ~isempty(Mask_filename)
%                     fid=fopen(Mask_filename{:} ,'r');
%                     if fid>0
%                         fclose(fid);
%                         Mask = load(Mask_filename{:});
%                         Mask = Mask.uvascroi;
%                         if ~isempty(Mask)
%                             [~, ~, E, Z,repet] = size(data.reco.data);
%                             
%                             tmp_data = data.reco.data;
%                             for slice = 1:Z
%                                 for repetition = 1:repet
%                                     for roi_number = 1:numel(Mask)
%                                         if abs(data.reco.fov_offsets(3,1,slice,repetition) - Mask(roi_number).fov_offsets(3)) < 1e-5
%                                             data.reco.data(:,:,:,slice,repetition) = tmp_data(:,:,:,slice,repetition).*repmat(Mask(roi_number).value,[1 1 E 1]);
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         else
%             
%         end
%         
%         % algo to detect the AIF
%         meansignal = mean(squeeze(data.reco.data),4);
%         volume_mask = meansignal>max(data.reco.data(:))*0.01;
%         [aif,scores] = extraction_aif_volume(squeeze(data.reco.data),volume_mask);
%         figure; plot(aif); legend(handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DSC_nbr), 'Location','NorthEast');
%         if sum(cell2num(scores(:,5))) > 10
%             %skip the process because the AIF is not
%             %good enough
%             logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC skip (No correct AIF)'];
%         else
%             
%             drawnow
%             % compute the parameteric maps using the AIF
%             % detected automatically
%             [~,~,~,TMAX,TTP,T0, CBV,CBF,MTT,~,~,~,~] = deconvolution_perfusion_gui(aif,squeeze(data.reco.data(:,:,1,:,:)),data.acq.tr*10^(-3),data.acq.echotime*10^(-3),handles); %#ok<ASGLU>
%             maps = {'CBV','CBF','MTT','TMAX','TTP','T0'};
%             if ~isempty(CBV)
%                 for xx = 1:numel(maps)
%                     uvascim.image = data;
%                     uvascim.image.reco.texte=maps{xx};
%                     uvascim.image.reco.no_expts = 1;
%                     uvascim.image.reco = rmfield(uvascim.image.reco, 'data');
%                     uvascim.image.reco.data(:,:,1,:) = eval(maps{xx});
%                     uvascim.image.reco.globalmax = max(uvascim.image.reco.data(:));
%                     uvascim.image.reco.globalmin = min(uvascim.image.reco.data(:));
%                     uvascim.image.clip=[0 250 1];
%                     uvascim.image.reco.unit = '%';
%                     ParamConfig=sprintf('##$QuantifMethod=Automatique deconvolution\n##$DSC=%s\n##$AIF=\n%s\n##$AIFScores=\n%s\n%s\n%s\n%s\n##$Warning=\n%d %d %d %d% d\n##$Warning Message=\n%s\n%s\n%s\n%s',...
%                         filename_DSC{:},num2str(aif), num2str([scores{:,1}]), num2str([scores{:,2}]), num2str([scores{:,3}]), num2str([scores{:,4}]),...
%                         [scores{:,5}], [scores{:,6}], [scores{:,7}], [scores{:,8}], [scores{:,9}]);
%                     uvascim.image.reco.paramQuantif = ParamConfig;
%                     uvascim.image.reco=orderfields(uvascim.image.reco);
%                     
%                     %Adapt the fov offsets and orientation infos
%                     uvascim.image.reco.fov_offsets=uvascim.image.reco.fov_offsets(:,:,:,1);
%                     uvascim.image.reco.fov_orientation=uvascim.image.reco.fov_orientation(:,:,:,1);
%                     uvascim.image.reco.label=uvascim.image.reco.label(:,:,1);
%                     uvascim.image.reco.phaselabel=uvascim.image.reco.phaselabel(:,:,1);
%                     uvascim.image.reco.fov_phase_orientation=uvascim.image.reco.fov_phase_orientation(:,:,1);
%                     uvascim.image.reco.scaling_factor=uvascim.image.reco.scaling_factor(:,:,1);
%                     uvascim.image.reco.scaling_offset=uvascim.image.reco.scaling_offset(:,:,1);
%                     
%                     % save the new maps andupdate MIA database
%                     % save file
%                     scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' maps{xx}  file_name_extension];
%                     save(fullfile(path, [scan_name '.mat']),'uvascim');
%                     
%                     % add to the database and update handles.parameters and handles.clips
%                     handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [maps{xx} file_name_extension], [scan_name '.mat']);
%                     guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
%                     
%                     clear uvascim;
%                 end
%                 clear CBV CBF MTT TMAX TTP T0;
%                 
%                 logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC Done'];
%             else
%                 logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC Skiped'];
%             end
%         end
%     else
%         logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC Skiped'];
%     end
% else
%     logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC already computed'];
% end