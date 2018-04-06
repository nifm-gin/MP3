function [files_in,files_out,opt] = Module_R2Prim(files_in,files_out,opt)
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
    module_option(:,8)   = {'output_filename_ext','R2Prim'};
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
    user_parameter(:,2)   = {'Select the T2_map scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the T2*corr3D scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','R2Prim','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''R2Prim''.'}'};
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
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), '_', opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_R2Prim:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

JT2Map = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
JT2StarCorr3D = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));


if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
%T2Map = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
%T2StarCorr3D = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);
T2Map = niftiread(files_in.In1{1});
T2StarCorr3D = niftiread(files_in.In2{1});

% check data compatibility (slice thickness and slice number)
% if JT2StarCorr3D.SliceThickness.value ~= JT2Map.SliceThickness.value
%     warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ a slice thickness missmatch between\n##$T2*map=%s\n##$ and \n##$Tmap=%s',...
%         T2star_filename,T2map_filename);
%     msgbox(warning_text, 'R2prim map warning') ;
%     R2prim_map = [];
%     return
% end


R2prim_map = T2StarCorr3D;
if size(T2Map,3) ~= size(T2StarCorr3D,3)
    for i = 1:size(T2StarCorr3D, 3)
        for j = 1:size(T2Map, 3)
            R2prim_slice_nbr = R2prim_slice_nbr+1;
            % Compute the CMRO2 map each slice with the same offset
            R2prim_map(:,:,R2prim_slice_nbr,1)=1./(RAW_T2star.data(:,:,i,1)*10^-3)-1./(RAW_T2map.data(:,:,j,1)*10^-3);
        end
    end
    if R2prim_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ no slice offset match between\n##$T2*map=%s\n##$ and \n##$T2map=%s',...
            files_in.In2{1},files_in.In1{1});
        msgbox(warning_text, 'R2prim map warning') ;
        return
    end
else
    final_res = max([size(T2Map,1), size(T2Map, 2), size(T2StarCorr3D, 1), size(T2StarCorr3D, 2)]);
    if size(T2StarCorr3D,1) ~= final_res || size(T2StarCorr3D,2) ~= final_res
        for i = 1:size(T2StarCorr3D, 3)
            tmp(:,:,i,1)= imresize(T2StarCorr3D(:,:,i,1),[final_res final_res],'bilinear');
        end
        T2StarCorr3D = tmp;
    end
    if size(T2Map,1) ~= final_res || size(T2Map,2) ~= final_res
        for i = 1:size(T2Map,3)
            tmp(:,:,i,1)= imresize(T2Map(:,:,i,1),[final_res final_res],'bilinear');
        end
        T2Map = tmp;
    end
    clear tmp
    
    % Compute the R2prim map for all slices
    R2prim_map=1./(T2StarCorr3D*10^-3)-1./(T2Map*10^-3);
end

%replace negative value of R2prim by the interpolation of their neibourgh 
R2prim_map_m = convn(R2prim_map,ones(3,3,3)./26,'same');
R2prim_map(R2prim_map(:)<=0) =R2prim_map_m(R2prim_map(:)<=0);


OutputImages = R2prim_map;
% OutputImages(OutputImages < 0) = -1;
% OutputImages(OutputImages > 5000) = -1;
% OutputImages(isnan(OutputImages)) = -1;
% if ~exist('OutputImages_reoriented', 'var')
%     OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
% end


% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(1)]){1});

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);

% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 


