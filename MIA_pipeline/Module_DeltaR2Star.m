function [files_in,files_out,opt] = Module_DeltaR2Star(files_in,files_out,opt)
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
    module_option(:,8)   = {'output_filename_ext','DeltaR2Star'};
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
    user_parameter(:,2)   = {'Select the T2_star_Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the T2_star_Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','DeltaR2Star','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''DeltaR2Star''.'}'};
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


if size(input_pre, 3) ~= size(input_post, 3)
    deltaR2star_slice_nbr = 0;
    
    for i = 1:size(input_pre, 3)
        for j = 1:size(input_post, 3)
            deltaR2star_slice_nbr = deltaR2star_slice_nbr+1;
            r2savant=input_pre(:,:,i,1);
            r2sapres=input_post(:,:,j,1);
            index_avant=find(r2savant<10e-2);% looking for the T2*pre < 100�s
            warning off					%To avoid division by zero warning
            r2savant=1./r2savant;       % R2*pre (ms-1)
            warning on
            r2savant(index_avant)=0;
            index_apres=find(r2sapres<10e-2);% olooking for the T2*post  < 100�s
            warning off					%To avoid division by zero warning
            r2sapres=1./r2sapres;	% R2*post (ms-1)
            warning on %#ok<*WNON>
            r2sapres(index_apres)=0;
            DeltaR2StarMap(:,:,deltaR2star_slice_nbr,1)=r2sapres-r2savant;
            index_nan=find(input_pre(:,:,j,1)==NaN); %#ok<FNAN>
            DeltaR2StarMap(index_nan)=0;
        end
    end
else
    %Initialisation des variables
    r2savant=input_pre(:,:,:,1);
    r2sapres=input_post(:,:,:,1);
    DeltaR2StarMap=zeros(size(input_pre));     %param�tres (DeltaR2*)
    index_avant=find(r2savant<10e-2);% on cherche les T2* < 100�s
    warning off					%To avoid division by zero warning
    r2savant=1./r2savant;	% R2*avant en ms-1
    warning on
    r2savant(index_avant)=0;
    index_apres=find(r2sapres<10e-2);% on cherche les T2* < 100�s
    warning off					%To avoid division by zero warning
    r2sapres=1./r2sapres;	% R2*apres en ms-1
    warning on %#ok<*WNON>
    r2sapres(index_apres)=0;
    DeltaR2StarMap(:,:,:,1)=r2sapres-r2savant;
    index_nan=find(input_pre==NaN); %#ok<FNAN>
    DeltaR2StarMap(index_nan)=0;
end

OutputImages = DeltaR2StarMap;
% OutputImages(OutputImages < 0) = -1;
% OutputImages(OutputImages > 5000) = -1;
% OutputImages(isnan(OutputImages)) = -1;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
end

%% Initial function 
% fid=fopen(filename_pre ,'r');
% if fid>0
%     fclose(fid);
%     data_pre = load(filename_pre);
%     RAW_T2star_pre = data_pre.uvascim.image.reco;
% else
%     warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ Somthing wrong with the data\n##$T2*pre=%s\n##$',...
%         filename_pre);
%     msgbox(warning_text, 'deltaR2star map warning') ;
%     deltaR2star_map = [];
%     return
% end
% fid=fopen(filename_post ,'r');
% if fid>0
%     fclose(fid);
%     data_post = load(filename_post);
%     RAW_T2star_post = data_post.uvascim.image.reco;
% else
%     warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ Somthing wrong with the data\n##$T2*post=%s\n##$',...
%         filename_pre);
%     msgbox(warning_text, 'deltaR2star map warning') ;
%     deltaR2star_map = [];
%     return
% end
% final_res = add_parameters{:}(1);
% if strcmp(final_res, 'Original')
%    rescale = 0;
% else
%     rescale = 1;
%     final_res = str2double(add_parameters{:}(1));
% end
% 
% check data compatibility (slice thickness and slice number)
% if RAW_T2star_pre.thickness ~= RAW_T2star_post.thickness
%     warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ a slice thickness missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
%         filename_pre,filename_post);
%     msgbox(warning_text, 'deltaR2star map warning') ;
%      deltaR2star_map = [];
%     return
% end
% if RAW_T2star_pre.no_samples ~= RAW_T2star_post.no_samples && rescale == 0
%      warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ a resolution missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
%         filename_pre,filename_post);
%     msgbox(warning_text, 'deltaR2star map warning') ;
%      deltaR2star_map = [];
%     return
% end
% if RAW_T2star_pre.no_slices ~= RAW_T2star_post.no_slices 
%     deltaR2star_slice_nbr = 0;
%     deltaR2star_map.reco.fov_offsets = [];
%     deltaR2star_map.reco.fov_orientation = [];
%     deltaR2star_map.reco.label = {};
%     deltaR2star_map.reco.phaselabel = {};
%     deltaR2star_map.reco.fov_phase_orientation = [];
%    
%     for i = 1:size(RAW_T2star_pre.data, 4)
%         for j = 1:size(RAW_T2star_post.data, 4)
%             if abs(RAW_T2star_pre.fov_offsets(3,1,i) - RAW_T2star_post.fov_offsets(3,1,j)) < 1e-5
%                 deltaR2star_slice_nbr = deltaR2star_slice_nbr+1; 
%                 if rescale == 0
%                      r2savant=RAW_T2star_pre.data(:,:,1,i);
%                      r2sapres=RAW_T2star_post.data(:,:,1,j);
%                      wrongpix_pre=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,1);  
%                     wrongpix_post=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,1);  
%                 else
%                     r2savant=imresize(RAW_T2star_pre.data(:,:,:,i),[final_res final_res],'bilinear');
%                     r2sapres=imresize(RAW_T2star_post.data(:,:,:,i),[final_res final_res],'bilinear');
%                     wrongpix_pre=zeros(final_res,final_res,1,1);  
%                     wrongpix_post=zeros(final_res,final_res,1,1); 
%                 end
% 
%                 index_avant=find(r2savant<10e-2);% looking for the T2*pre < 100�s
%                 warning off					%To avoid division by zero warning
%                 r2savant=1./r2savant;       % R2*pre (ms-1)
%                 warning on
%                 r2savant(index_avant)=0;
% 
%                 index_apres=find(r2sapres<10e-2);% olooking for the T2*post  < 100�s
%                 warning off					%To avoid division by zero warning
%                 r2sapres=1./r2sapres;	% R2*post (ms-1)
%                 warning on %#ok<*WNON>
%                 r2sapres(index_apres)=0;
%                 
%                 deltaR2star_map.reco.data(:,:,1,deltaR2star_slice_nbr)=r2sapres-r2savant;
%                 index_nan=find(RAW_T2star_pre.data(:,:,1,j)==NaN); %#ok<FNAN>
%                 deltaR2star_map.reco.data(index_nan)=0;
%                 Compute the error map
%                 proc.err=;
%                 
%                 Compute the wrongpix map
%                 if rescale == 0
%                     wrongpix_pre(find(RAW_T2star_pre.wrongpix(:,:,j)>0)) = 1;%#ok<FNDSB>
%                     wrongpix_post(find(RAW_T2star_post.wrongpix(:,:,j)>0)) = 1;%#ok<FNDSB>
%                 else
%                     wrongpix_pre(find(imresize(RAW_T2star_pre.wrongpix(:,:,j),[final_res final_res],'bilinear')>0)) = 1;%#ok<FNDSB>
%                     wrongpix_post(find(imresize(RAW_T2star_post.wrongpix(:,:,j),[final_res final_res],'bilinear')>0)) = 1;%#ok<FNDSB>    
%                 end
%                 deltaR2star_map.reco.wrongpix(:,:,1,deltaR2star_slice_nbr)=or(wrongpix_pre,wrongpix_post) ; %Propagation of existing wronpix %#ok<FNDSB>
%                 deltaR2star_map.reco.wrongpix(index_avant)=1;						%Pixels where r2savant has been set to zero
%                 deltaR2star_map.reco.wrongpix(index_apres)=1;						%Pixels where r2sapres has been set to zero
%                 deltaR2star_map.reco.wrongpix(index_nan)=1;						%Pixels where proc.data=NaN
%                 
%                 Update the deltaR2star_map structure
%                 deltaR2star_map.reco.fov_offsets(:,1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.fov_offsets(:,1,i,1);
%                 deltaR2star_map.reco.fov_orientation(:,1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.fov_orientation(:,1,i,1);
%                 deltaR2star_map.reco.label(1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.label(1,i,1);
%                 deltaR2star_map.reco.phaselabel(1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.phaselabel(1,i,1);
%                 deltaR2star_map.reco.fov_phase_orientation(1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.fov_phase_orientation(1,i,1);
%                 
%             end
%         end
%     end
%     if deltaR2star_slice_nbr == 0
%         warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ no slice offset match between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
%            filename_pre,filename_post);
%         msgbox(warning_text, 'deltaR2star map warning') ;
%         return
%     end
%     deltaR2star_map.reco.no_slices=deltaR2star_slice_nbr;
% else
%     Initialisation des variables
%     deltaR2star_map.reco='';
%     if rescale == 0
%         r2savant=RAW_T2star_pre.data(:,:,1,:);
%         r2sapres=RAW_T2star_post.data(:,:,1,:);
%         deltaR2star_map.reco.data=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,RAW_T2star_pre.no_slices);     %param�tres (DeltaR2*)
%         deltaR2star_map.reco.wrongpix=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,RAW_T2star_pre.no_slices);  %carte de pixel � exclure
%     else
%         r2savant=imresize(RAW_T2star_pre.data(:,:,:,:),[final_res final_res],'bilinear');
%         r2sapres=imresize(RAW_T2star_post.data(:,:,:,:),[final_res final_res],'bilinear');
%         deltaR2star_map.reco.data=zeros(final_res,final_res,1,RAW_T2star_pre.no_slices);     %param�tres (DeltaR2*)
%         deltaR2star_map.reco.wrongpix=zeros(final_res,final_res,1,RAW_T2star_pre.no_slices);  %carte de pixel � exclure
%     end
%     
%     index_avant=find(r2savant<10e-2);% on cherche les T2* < 100�s
%     warning off					%To avoid division by zero warning
%     r2savant=1./r2savant;	% R2*avant en ms-1
%     warning on
%     r2savant(index_avant)=0;
%     
%     index_apres=find(r2sapres<10e-2);% on cherche les T2* < 100�s
%     warning off					%To avoid division by zero warning
%     r2sapres=1./r2sapres;	% R2*apres en ms-1
%     warning on %#ok<*WNON>
%     r2sapres(index_apres)=0;
%     
%     deltaR2star_map.reco.data(:,:,1,:)=r2sapres-r2savant;
%     index_nan=find(RAW_T2star_pre.data==NaN); %#ok<FNAN>
%     deltaR2star_map.reco.data(index_nan)=0;
%     Compute the error map
%     proc.err=;
%     
%     Compute the wrongpix map
%     if rescale == 0
%         deltaR2star_map.reco.wrongpix(find(RAW_T2star_pre.wrongpix>0)) = 1;%#ok<FNDSB>
%         deltaR2star_map.reco.wrongpix(find(RAW_T2star_post.wrongpix>0)) = 1;%#ok<FNDSB>
%     else
%         for i = 1:RAW_T2star_pre.no_slices
%             wrongpix_pre = imresize(RAW_T2star_pre.wrongpix(:,:,i),[final_res final_res],'bilinear');
%             wrongpix_pre(find(wrongpix_pre>0)) = 1;%#ok<FNDSB>
%             index_nan_pre=isnan(wrongpix_pre);
%             wrongpix_pre(index_nan_pre)=0;
%              wrongpix_post = imresize(RAW_T2star_post.wrongpix(:,:,i),[final_res final_res],'bilinear');
%             wrongpix_post(find(wrongpix_post>0)) = 1;%#ok<FNDSB>
%             index_nan_post=isnan(wrongpix_post); 
%             wrongpix_post(index_nan_post)=0;
%             deltaR2star_map.reco.wrongpix(:,:,1,i)=or(wrongpix_pre,wrongpix_post) ;
%         end
%     end
% 
%     deltaR2star_map.reco.wrongpix(index_avant)=1;						%Pixels where r2savant has been set to zero
%     deltaR2star_map.reco.wrongpix(index_apres)=1;						%Pixels where r2sapres has been set to zero
%     deltaR2star_map.reco.wrongpix(index_nan)=1;						%Pixels where proc.data=NaN
%     
%     deltaR2star_map.reco.no_slices=RAW_T2star_pre.no_slices;
%     Adapt the fov offsets and orientations infos
%     deltaR2star_map.reco.fov_offsets=RAW_T2star_pre.fov_offsets(:,1,:,1);
%     deltaR2star_map.reco.fov_orientation=RAW_T2star_pre.fov_orientation(:,1,:,1);
%     deltaR2star_map.reco.label=RAW_T2star_pre.label(1,:,1);
%     deltaR2star_map.reco.phaselabel=RAW_T2star_pre.phaselabel(1,:,1);
%     deltaR2star_map.reco.fov_phase_orientation=RAW_T2star_pre.fov_phase_orientation(1,:,1);
% end
% 
% Set the dimensions of the proc structure
% if rescale == 0
%     deltaR2star_map.reco.no_samples=RAW_T2star_pre.no_samples;
%     deltaR2star_map.reco.no_views=RAW_T2star_pre.no_views;
% else
%     deltaR2star_map.reco.no_samples=final_res;
%     deltaR2star_map.reco.no_views=final_res;
% end
% deltaR2star_map.reco.no_echoes=1; %Number of parameters stored here
% deltaR2star_map.reco.no_expts=1;
% 
% Complete the reco structure
% deltaR2star_map.reco.texte='DeltaR2*';
% deltaR2star_map.reco.unit='ms-1';
% deltaR2star_map.reco.date=date;
% deltaR2star_map.reco.displayedecho=1;
% deltaR2star_map.reco.displayedslice=1;
% deltaR2star_map.reco.displayedexpt=1;
% deltaR2star_map.reco_number = 1;
% deltaR2star_map.scan_number = 104;
% deltaR2star_map.clip=[0 0.2 1];
% deltaR2star_map.reco.globalmax = max(deltaR2star_map.reco.data(:));
% deltaR2star_map.reco.globalmin = min(deltaR2star_map.reco.data(:));
% deltaR2star_map.reco.thickness = data_pre.uvascim.image.reco.thickness;
% deltaR2star_map.reco=orderfields(deltaR2star_map.reco);
% if iscell(RAW_T2star_pre.paramQuantif)
%     RAW_T2star_pre.paramQuantif = RAW_T2star_pre.paramQuantif{1};
% end
% if iscell(RAW_T2star_post.paramQuantif)
%     RAW_T2star_post.paramQuantif = RAW_T2star_post.paramQuantif{1};
% end
% ParamConfig=sprintf('##$QuantifMethod=''(1/T2*post)-(1/T2*pre)''\n##$final_res=%s\n##$T2*pre=%s\n##$T2*post=%s\n\n##$T2*pre scan info\n%s\n\n##$T2*post scan info\n%s',...
%     add_parameters{:}{1},filename_pre,filename_post, RAW_T2star_pre.paramQuantif, RAW_T2star_post.paramQuantif);
% deltaR2star_map.reco.paramQuantif = ParamConfig;
% deltaR2star_map.reco=orderfields(deltaR2star_map.reco);
% 
% 
% complete the structure
% deltaR2star_map.acq = data_pre.uvascim.image.acq;
% deltaR2star_map.filename = data_pre.uvascim.image.filename;

%% 


% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.Description = [info.Description, 'Modified by the DeltaR2 Module'];

% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 


