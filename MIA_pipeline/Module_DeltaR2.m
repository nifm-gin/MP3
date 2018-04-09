function [files_in,files_out,opt] = Module_DeltaR2(files_in,files_out,opt)
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
    %
    %     %%   % define every option needed to run this module
    %     % --> module_option(1,:) = field names
    %     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','AllName'};
    module_option(:,4)   = {'RefInput',1};
    module_option(:,5)   = {'InputToReshape',1};
    module_option(:,6)   = {'Table_in', table()};
    module_option(:,7)   = {'Table_out', table()};
    module_option(:,8)   = {'output_filename_ext','DeltaR2'};
    module_option(:,9)   = {'Output_orientation','First input'};
    
    
    
    
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    %
    %      additional_info_name = {'Operation', 'Day of the 2nd scan', 'Name of the 2nd scan'};
    %         additional_info_data = {'Addition', 'first', handles.Math_parameters_list{1}};
    %         additional_info_format = {{'Addition', 'Subtraction', 'Multiplication', 'Division', 'Percentage', 'Concentration'},...
    %             ['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list};
    %         additional_info_editable = [1 1 1];
    
    
    %     %% list of everything displayed to the user associated to their 'type'
    %      % --> user_parameter(1,:) = user_parameter_list
    %      % --> user_parameter(2,:) = user_parameter_type
    %      % --> user_parameter(3,:) = parameter_default
    %      % --> user_parameter(4,:) = psom_parameter_list
    %      % --> user_parameter(5,:) = Help : text data which describe the parameter (it
    %      % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','', '','Description of the module'}  ;
    user_parameter(:,2)   = {'Select the Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','DeltaR2','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''DeltaR2''.'}'};
    user_parameter(:,5)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        }'};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
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
    opt.Table_out = opt.Table_in(opt.RefInput,:);
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
    error('Module_Coreg_Est:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
%if ~ischar(files_in)
%    error('files in should be a char');
%end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load input Nii file
input(1).nifti_header = spm_vol(files_in.In1{1});
input(2).nifti_header = spm_vol(files_in.In2{1});

J1 = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
J2 = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));


if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
input_pre = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
input_post = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);

%% WTF
% %%% ? What's the difference between se_echotime and RAW_data_pre.echotime ?
% %%% I don't think those informations are available in the nifti files.
% se_echotime = scan_acqp('##$SpinEchoTime=',data_pre.uvascim.image.texte,1);
% if isnan(se_echotime)
%     se_echotime = scan_acqp('##$PVM_EchoTime=',data_pre.uvascim.image.texte,1);
% end
% se_echo_pos = abs(RAW_data_pre.echotime - se_echotime) <= 2;
% se_echo_pos = find(se_echo_pos == 1);
%%
%nb_echo_used=opt.nb_echos;
se_echotime = J1.EchoTime.value;
se_echo_pos = 1:length(J1.EchoTime.value);


if size(input_pre,3) ~= size(input_post,3)
    deltaR2_map=zeros(size(input_pre));
    deltaR2_slice_nbr = 0;
    for i = 1:size(input_pre, 3)
        for j = 1:size(input_post, 3)
            deltaR2_slice_nbr = deltaR2_slice_nbr+1;
            % Compute the CMRO2 map each slice with the same offset
            temp_avant=squeeze(input_pre(:,:,i,se_echo_pos));
            temp_apres=squeeze(input_post(:,:,j,se_echo_pos));
            index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
            index_apres=find(temp_apres<(1e-3*mean(temp_apres(:))));
            
            %To avoid division by small numbers
            temp_apres(index_apres)=1;
            warning off %#ok<WNOFF>
            deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
            warning on %#ok<WNON>
            deltaR2(index_avant)=0;
            deltaR2(index_apres)=0;
            deltaR2_map(:,:,deltaR2_slice_nbr,1)=deltaR2;
            
        end
    end
    
    if deltaR2_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ no slice offset match between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
            filename_SO2map,filename_CBF);
        msgbox(warning_text, 'deltaR2 map warning') ;
        return
    end
else
    % Compute the deltaR2 map for all slices
    deltaR2_map.reco='';
    deltaR2_map=zeros([size(input_pre,1), size(input_pre, 2), size(input_pre, 3)]);     % DeltaR2 map
    for m_slice=1:size(input_pre, 3)       
        temp_avant=squeeze(input_pre(:,:,m_slice,se_echo_pos));
        temp_apres=squeeze(input_post(:,:,m_slice,se_echo_pos));
        index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
        index_apres=find(temp_apres<(1e-3*mean(temp_apres(:)))); 
        %To avoid division by small numbers
        temp_apres(index_apres)=1;
        warning off %#ok<WNOFF>
        %deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
        A = -(1./se_echotime)';
        B = log(temp_apres ./ temp_avant);
        B = permute(B, [3,2,1]);
        for i=1:size(B,3)
            deltaR2(i,:) = A*squeeze(B(:,:,i));
        end
        %deltaR2= A.*B ;%ms-1
        %deltaR2 = permute(deltaR2, [3,2,1]);
        warning on %#ok<WNON>
        %deltaR2(index_avant)=0;
        %deltaR2(index_apres)=0;
        deltaR2_map(:,:,m_slice,1)=deltaR2;
        %Compute the err map
        %proc.err(:,:,1,m_slice)=tmpim;
    end
    %Adapt the fov offsets and orientations infos
end

%Set the dimensions of the proc structure


%Complete the reco structure




% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)

OutputImages = deltaR2_map;
OutputImages(OutputImages < 0) = -1;
OutputImages(OutputImages > 5000) = -1;
OutputImages(isnan(OutputImages)) = -1;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
end





% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});


info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);
%info2.Description = [info.Description, 'Modified by T2map Module'];


% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 


%% Initial function 
% fid=fopen(filename_pre ,'r');
% if fid>0
%     fclose(fid);
%     data_pre = load(filename_pre);
%     RAW_data_pre = data_pre.uvascim.image.reco;
% else
%     warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ Somthing wrong with the data\n##$MGESEpre=%s\n##$',...
%         filename_pre);
%     msgbox(warning_text, 'deltaR2 map warning') ;
%     deltaR2_map = [];
%     return
% end
% fid=fopen(filename_post ,'r');
% if fid>0
%     fclose(fid);
%     data_post = load(filename_post);
%     RAW_data_post = data_post.uvascim.image.reco;
% else
%     warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ Somthing wrong with the data\n##$MGESEpost=%s\n##$',...
%         filename_post);
%     msgbox(warning_text, 'deltaR2 map warning') ;
%     deltaR2_map = [];
%     return
% end
% % %% smooth data
% % for i=1:size(RAW_data_pre,1)
% %     if sum(data_in_vector(i,:)) ~= 0
% %         data_in_vector(i,:) = smooth(data_in_vector(i,:), 'lowess')';
% %     end
% % end
% % %% smooth data
% % for i=1:size(data_in_vector,1)
% %     if sum(data_in_vector(i,:)) ~= 0
% %         data_in_vector(i,:) = smooth(data_in_vector(i,:), 'lowess')';
% %     end
% % end
% 
% 
% se_echotime = scan_acqp('##$SpinEchoTime=',data_pre.uvascim.image.texte,1);
% if isnan(se_echotime)
%     se_echotime = scan_acqp('##$PVM_EchoTime=',data_pre.uvascim.image.texte,1);
% end
% se_echo_pos = abs(RAW_data_pre.echotime - se_echotime) <= 2;
% se_echo_pos = find(se_echo_pos == 1);
% 
% nb_echo_used=str2double(add_parameters{:}(1));
% final_res = add_parameters{:}(2);
% if strcmp(final_res, 'Original')
%    rescale = 0;
% else
%     rescale = 1;
%     final_res = str2double(add_parameters{:}(2));
% end
% 
% 
% % check data compatibility (slice thickness and slice number)
% if RAW_data_pre.thickness ~= RAW_data_post.thickness
%     warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ a slice thickness missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
%         filename_pre,filename_post);
%     msgbox(warning_text, 'deltaR2 map warning') ;
%     deltaR2_map = [];
%     return
% end
% if RAW_data_pre.no_samples ~= RAW_data_post.no_samples && rescale == 0
%      warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ a resolution missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
%         filename_pre,filename_post);
%     msgbox(warning_text, 'deltaR2 map warning') ;
%     deltaR2_map = [];
%     return
% end
% if RAW_data_pre.no_slices ~= RAW_data_post.no_slices
%     deltaR2_slice_nbr = 0;
%     deltaR2_map.reco.fov_offsets = [];
%     deltaR2_map.reco.fov_orientation = [];
%     deltaR2_map.reco.label = {};
%     deltaR2_map.reco.phaselabel = {};
%     deltaR2_map.reco.fov_phase_orientation = [];
%     for i = 1:size(RAW_data_pre.data, 4)
%         for j = 1:size(RAW_data_post.data, 4)
%             if abs(RAW_data_pre.fov_offsets(3,1,i) - RAW_data_post.fov_offsets(3,1,j)) < 1e-5
%                 deltaR2_slice_nbr = deltaR2_slice_nbr+1;
%                 % Compute the CMRO2 map each slice with the same offset
%                 if nb_echo_used == 1
%                     if rescale == 0
%                         temp_avant=squeeze(RAW_data_pre.data(:,:,se_echo_pos,i));
%                         temp_apres=squeeze(RAW_data_post.data(:,:,se_echo_pos,j));
%                     else
%                         temp_avant=imresize(squeeze(RAW_data_pre.data(:,:,se_echo_pos,i)),[final_res final_res],'bilinear');
%                         temp_apres=imresize(squeeze(RAW_data_post.data(:,:,se_echo_pos,j)),[final_res final_res],'bilinear');
%                     end
%                 else
%                     if rescale == 0
%                         temp_avant=squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,i),3));
%                         temp_apres=squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,j),3));
%                     else
%                         temp_avant=imresize(squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,i),3)),[final_res final_res],'bilinear');
%                         temp_apres=imresize(squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,j),3)),[final_res final_res],'bilinear');
%                     end
%                 end
%                 index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
%                 index_apres=find(temp_apres<(1e-3*mean(temp_apres(:))));
%                 
%                 %To avoid division by small numbers
%                 temp_apres(index_apres)=1;
%                 warning off %#ok<WNOFF>
%                 deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
%                 warning on %#ok<WNON>
%                 deltaR2(index_avant)=0;
%                 deltaR2(index_apres)=0;
%                 deltaR2_map.reco.data(:,:,1,deltaR2_slice_nbr)=deltaR2;
%                 
%                 %Compute the wrongpix map : wrongpix=1 if temp_apres or temp_avant < 0.001*mean(image)
%                 if rescale == 0
%                     tmpim=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views);
%                 else
%                     tmpim=zeros(final_res,final_res);
%                 end
%                 tmpim(index_avant)=1;
%                 tmpim(index_apres)=1;
%                 tmpim(deltaR2<0)=1;
%                 deltaR2_map.reco.wrongpix(:,:,1,deltaR2_slice_nbr)=tmpim;
%                 % Update the SO2map structure
%                 deltaR2_map.reco.fov_offsets(:,1,deltaR2_slice_nbr,1) = RAW_data_pre.fov_offsets(:,1,i,1);
%                 deltaR2_map.reco.fov_orientation(:,1,deltaR2_slice_nbr,1) = RAW_data_pre.fov_orientation(:,1,i,1);
%                 deltaR2_map.reco.label(1,deltaR2_slice_nbr,1) = RAW_data_pre.label(1,i,1);
%                 deltaR2_map.reco.phaselabel(1,deltaR2_slice_nbr,1) = RAW_data_pre.phaselabel(1,i,1);
%                 deltaR2_map.reco.fov_phase_orientation(1,deltaR2_slice_nbr,1) = RAW_data_pre.fov_phase_orientation(1,i,1);    
%             end
%         end
%     end
%     if deltaR2_slice_nbr == 0
%         warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ no slice offset match between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
%             filename_SO2map,filename_CBF);
%         msgbox(warning_text, 'deltaR2 map warning') ;
%         return
%     end
%     deltaR2_map.reco.no_slices=deltaR2_slice_nbr;
% else
%     % Compute the deltaR2 map for all slices
%     deltaR2_map.reco='';
%     if rescale == 0
%         deltaR2_map.reco.data=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views,1,RAW_data_pre.no_slices);     % DeltaR2 map
%         deltaR2_map.reco.wrongpix=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views,1,RAW_data_pre.no_slices);  % excluded voxels
%     else
%         deltaR2_map.reco.data=zeros(final_res,final_res,1,RAW_data_pre.no_slices);     
%         deltaR2_map.reco.wrongpix=zeros(final_res,final_res,1,RAW_data_pre.no_slices);  
%     end
%     
%     for m_slice=1:RAW_data_pre.no_slices
%         if nb_echo_used == 1
%             if rescale == 0
%                 
%                 temp_avant=squeeze(RAW_data_pre.data(:,:,se_echo_pos,m_slice));
%                 temp_apres=squeeze(RAW_data_post.data(:,:,se_echo_pos,m_slice));
%             else
%                 
%                 temp_avant=imresize(squeeze(RAW_data_pre.data(:,:,se_echo_pos,m_slice)),[final_res final_res],'bilinear');
%                 temp_apres=imresize(squeeze(RAW_data_post.data(:,:,se_echo_pos,m_slice)),[final_res final_res],'bilinear');
%             end
%         else
%             if rescale == 0
%                 temp_avant=squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3));
%                 temp_apres=squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3));
%             else
%                 temp_avant=imresize(squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3)),[final_res final_res],'bilinear');
%                 temp_apres=imresize(squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3)),[final_res final_res],'bilinear');
%             end
%         end
%         index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
%         index_apres=find(temp_apres<(1e-3*mean(temp_apres(:))));
%         
%         %To avoid division by small numbers
%         temp_apres(index_apres)=1;
%         warning off %#ok<WNOFF>
%         deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
%         warning on %#ok<WNON>
%         deltaR2(index_avant)=0;
%         deltaR2(index_apres)=0;
%         deltaR2_map.reco.data(:,:,1,m_slice)=deltaR2;
%         
%         %Compute the err map
%         %proc.err(:,:,1,m_slice)=tmpim;
%         
%         %Compute the wrongpix map : wrongpix=1 if temp_apres or temp_avant < 0.001*mean(image)
%         if rescale == 0
%             tmpim=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views);
%         else
%             tmpim=zeros(final_res,final_res);
%         end
%         tmpim(index_avant)=1;
%         tmpim(index_apres)=1;
%         tmpim(deltaR2<0)=1;
%         deltaR2_map.reco.wrongpix(:,:,1,m_slice)=tmpim;
%     end
%     %Adapt the fov offsets and orientations infos
%     deltaR2_map.reco.fov_offsets=RAW_data_pre.fov_offsets(:,1,:,1);
%     deltaR2_map.reco.fov_orientation=RAW_data_pre.fov_orientation(:,1,:,1);
%     deltaR2_map.reco.label=RAW_data_pre.label(1,:,1);
%     deltaR2_map.reco.phaselabel=RAW_data_pre.phaselabel(1,:,1);
%     deltaR2_map.reco.fov_phase_orientation=RAW_data_pre.fov_phase_orientation(1,:,1);
%     deltaR2_map.reco.no_slices=RAW_data_pre.no_slices;
% end
% 
% %Set the dimensions of the proc structure
% if rescale == 0
%     deltaR2_map.reco.no_samples=RAW_data_pre.no_samples;
%     deltaR2_map.reco.no_views=RAW_data_pre.no_views;
% else
%     deltaR2_map.reco.no_samples = final_res;
%     deltaR2_map.no_views = final_res;
% end
% deltaR2_map.reco.no_echoes=1; %Number of parameters stored here
% deltaR2_map.reco.no_expts=1;
% 
% %Complete the reco structure
% deltaR2_map.reco.texte='DeltaR2';
% deltaR2_map.reco.unit='ms-1';
% deltaR2_map.reco.date=date;
% deltaR2_map.reco.displayedecho=1;
% deltaR2_map.reco.displayedslice=1;
% deltaR2_map.reco.displayedexpt=1;
% deltaR2_map.reco_number = 1;
% deltaR2_map.scan_number = 104;
% deltaR2_map.clip=[0 0.4 1];
% deltaR2_map.reco.globalmax = max(deltaR2_map.reco.data(:));
% deltaR2_map.reco.globalmin = min(deltaR2_map.reco.data(:));
% deltaR2_map.reco.thickness = data_pre.uvascim.image.reco.thickness;
% 
% ParamConfig=sprintf('##$QuantifMethod=''-(1/se_echotime) * log(temp_apres ./ temp_avant)''\n##$Spin Echo time=%d\n##$Spin Echo position=%d\n##$Number of echo used=%s\n##$T2pre=%s\n##$T2post=%s\n\n##$Sequence pre info\n%s\n\n##$Sequence post info\n%s',...
%     se_echotime, se_echo_pos, add_parameters{:}{1},filename_pre,filename_post, [RAW_data_pre.iminfos{:}], [RAW_data_post.iminfos{:}]);
% deltaR2_map.reco.paramQuantif = ParamConfig;
% deltaR2_map.reco=orderfields(deltaR2_map.reco);
% 
% %complete the structure
% deltaR2_map.acq = data_pre.uvascim.image.acq;
% deltaR2_map.filename = data_pre.uvascim.image.filename;