function [files_in,files_out,opt] = Module_Arithmetic(files_in,files_out,opt)
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
    module_option(:,3)   = {'OutputSequenceName','Extension'};
    module_option(:,4)   = {'RefInput',1};
    module_option(:,5)   = {'InputToReshape',1};
    module_option(:,6)   = {'Table_in', table()};
    module_option(:,7)   = {'Table_out', table()};
    module_option(:,8)   = {'Operation', 'Addition'};
    module_option(:,9)   = {'output_filename_ext','_Arith'};
    module_option(:,10)   = {'Output_orientation','First input'};
    
    
    
    
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
     user_parameter(:,1)   = {'Description','Text','','','', 'Description of the module'}  ;
    user_parameter(:,2)   = {'Select the first scan','1ScanOr1ROI','','',{'SequenceName'}, ''};
    user_parameter(:,3)   = {'Select the operation you would like to apply','cell', {'Addition', 'Subtraction', 'Multiplication', 'Division', 'Percentage'},'Operation','', ''};
    user_parameter(:,4)   = {'Select the second scan','1ScanOr1ROI','','',{'SequenceName'}, ''};
    user_parameter(:,5)   = {'   .Output filename extension','char','_Smooth','output_filename_ext','',...
        {'Specify the string to be added to the first filename.'
        'Default filename extension is ''_Arith''.'}'};
    user_parameter(:,6)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        }'};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', 'VariableNames', VariableNames);
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
    opt.Table_out = opt.Table_in(1,:);
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

if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
input1 = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
input2 = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);



switch opt.Operation
    case 'Addition'
        OutputImages = input1 + input2;
    case 'Subtraction'
        OutputImages = input2 - input1;
    case 'Multiplication'
       
        if sum(opt.Table_in.Type == 'ROI') > 0
            if length(size(input1)) > length(size(input2))
                OutputImages = input1 .* repmat(input2, [1 1 1 size(input1,4) size(input1,5) size(input1,6) size(input1,7)]);    
            else
                OutputImages =  repmat(input1, [1 1 1 size(input1,4) size(input1,5) size(input1,6) size(input1,7)]) .* input2;    
            end
        else
            OutputImages = input1 .* input2;
        end
    case'Division'
        OutputImages = input2 ./ input1;
    case 'Percentage'
        OutputImages = ((input2 - input1) ./ input2 .* 100);
%     case 'Concentration'
%         question = strcat('Please enter the relaxivity (mM-1.sec-1)');
%         relaxivity  = inputdlg(question,'Relaxivity',1,{'6'});
%         if isempty(relaxivity)
%             return
%         end
%         result_map.reco.data = (1./(image1.uvascim.image.reco.data/1000) - 1/(image2.uvascim.image.reco.data/1000)) / str2double(relaxivity{:}) ;
%         Methode_info = ['relaxivity = ' relaxivity{:} ' mM-1.sec-1'];
end

% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)
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
info2.Description = [info.Description, 'Modified by the Arithmetic Module'];

% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 

