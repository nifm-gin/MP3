function [files_in,files_out,opt] = Module_Coreg_Est(files_in,files_out,opt)
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
    module_option(:,3)   = {'Execution_Mode','Through all sessions of one Patient'};
    module_option(:,4)   = {'OutputSequenceName','Prefix'};
    module_option(:,5)   = {'Function','nmi'};
    module_option(:,6)   = {'Separation','4 2'};
    module_option(:,7)   = {'Tolerence','0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001'};
    module_option(:,8)   = {'Hist_Smooth','7 7'};
    module_option(:,9)   = {'Interpolation','4th Degree B-Spline'};
    module_option(:,10)   = {'Wrapping','No wrap'};
    module_option(:,11)   = {'Masking','Dont mask images'};
    module_option(:,12)   = {'output_filename_ext','CoregEst'};
    module_option(:,13)   = {'RefInput',2};
    module_option(:,14)   = {'InputToReshape',1};
    module_option(:,15)   = {'Table_in', table()};
    module_option(:,16)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
%   
%     %% list of everything displayed to the user associated to their 'type'
%      % --> user_parameter(1,:) = user_parameter_list
%      % --> user_parameter(2,:) = user_parameter_type
%      % --> user_parameter(3,:) = parameter_default
%      % --> user_parameter(4,:) = psom_parameter_list
%      % --> user_parameter(5,:) = Help : text data which describe the parameter (it
%      % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
    'Within-subject registration using a rigid-body model and image reslicing.'
    ''
    'The registration method used here is based on work by Collignon et al.'
    'The original interpolation method described in this paper has been changed in order to give a smoother cost function.'
    'The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
    ''
    'A copy of the ''Source'' and '' Other Images'' is performed before the registration process'
    'This module DO NOT overwrite the orignal files' 
    ''
    'Registration parameters are stored in the headers of the "source" and the "other" images.'
    }'};
    user_parameter(:,2)   = {'Execution Mode','cell',{'All Database','Through all sessions of one Patient','Through Each Session'},'Execution_Mode','','',...
        'Please select the excecution mode for the coreg module'};
    user_parameter(:,3)   = {'Reference Image','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
         'This is the image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
    user_parameter(:,4)   = {'Source Image','1Scan','','',{'SequenceName'},'Mandatory',...
         'This is the image that is jiggled about to best match the reference.'};
    user_parameter(:,5)   = {'Other Images','XScanOrXROI','','',{'SequenceName'},'Optional',...
         'These are any images that need to remain in alignment with the source image.'};
    user_parameter(:,6)   = {'Parameters','','','','','',''};
    user_parameter(:,7)   = {'    Estimation Options','', '','','','',...
        'Various registration options, which are passed to the Powell optimisation algorithm.'};
    user_parameter(:,8)   = {'       .Objective Function','cell',{'mi','ncc', 'nmi', 'ecc'},'Function','','',...
        'Registration involves finding parameters that either maximise or minimise some objective function. For inter-modal registration, use Mutual Information, Normalised Mutual Information, or Entropy Correlation Coefficient. For within modality, you could also use Normalised Cross Correlation.'};
    user_parameter(:,9)   = {'       .Tolerances','char','','Separation','','',...
        'The average distance between sampled points (in mm).  Can be a vector to allow a coarse registration followed by increasingly fine ones.'};
    user_parameter(:,10)  = {'       .Separation','numeric','','Tolerence','','',...
        'The accuracy for each parameter.  Iterations stop when differences between successive estimates are less than the required tolerance.'};
    user_parameter(:,11)  = {'       .Histogram Smoothing','numeric','','Hist_Smooth','','',...
        'Gaussian smoothing to apply to the 256x256 joint histogram. Other information theoretic coregistration methods use fewer bins, but Gaussian smoothing seems to be more elegant.'};
    user_parameter(:,12)  = {'       .Filename prefix','char','','output_filename_ext','','',...
        'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''CoregEst''.'};


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


if strcmp(files_out, '')
    [Path_In2, Name_In2, ~] = fileparts(files_in.In2{1});
    tags2 = opt.Table_in(opt.Table_in.Path == [Path_In2, filesep],:);
    tags2 = tags2(tags2.Filename == Name_In2,:);
    assert(size(tags2, 1) == 1);
    tags_out_In2 = tags2;
    tags_out_In2.IsRaw = categorical(0);
    tags_out_In2.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_In2.SequenceName = categorical(cellstr([opt.output_filename_ext, char(tags_out_In2.SequenceName)]));
    tags_out_In2.Filename = categorical(cellstr([char(tags_out_In2.Patient), '_', char(tags_out_In2.Tp), '_', char(tags_out_In2.SequenceName)]));
    f_out = [char(tags_out_In2.Path), char(tags_out_In2.Patient), '_', char(tags_out_In2.Tp), '_', char(tags_out_In2.SequenceName), '.nii'];
    files_out.In2{1} = f_out;
    opt.Table_out = tags_out_In2;
    for i=1:length(files_in.In3)
        if ~isempty(files_in.In3{i})
            [Path_In3, Name_In3, ~] = fileparts(files_in.In3{i});
            tags3 = opt.Table_in(opt.Table_in.Path == [Path_In3, filesep],:);
            tags3 = tags3(tags3.Filename == Name_In3,:);
            assert(size(tags3, 1) == 1);
            tags_out_In3 = tags3;
            tags_out_In3.IsRaw = categorical(0);
            tags_out_In3.SequenceName = categorical(cellstr([opt.output_filename_ext, char(tags_out_In3.SequenceName)]));
            if tags_out_In3.Type == 'Scan'
                tags_out_In3.Path = categorical(cellstr([opt.folder_out, filesep]));
                f_out = [char(tags_out_In3.Path), char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-', char(tags_out_In3.SequenceName), '.nii'];
                 tags_out_In3.Filename = categorical(cellstr([char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-', char(tags_out_In3.SequenceName)]));
            else
                f_out = [char(tags_out_In3.Path), char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-ROI-', char(tags_out_In3.SequenceName), '.nii'];
                tags_out_In3.Filename = categorical(cellstr([char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-ROI-', char(tags_out_In3.SequenceName)]));
            end
            files_out.In3{i} = f_out;
            opt.Table_out = [opt.Table_out ; tags_out_In3];
        end
    end
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Coreg_Est:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
%if ~ischar(files_in) 
%    error('files in should be a char');
%end



if isfield(opt,'threshold') && (~isnumeric(opt.threshold))
    opt.threshold = str2double(opt.threshold);
    if isnan(opt.threshold)
        disp('The threshold used was not a number')
        return
    end
end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FixedImInfo = niftiinfo(files_in.In1{1});
[path, name, ~] = fileparts(files_in.In1{1});
FixedImJsonfile = [path, filesep, name, '.json'];
fid = fopen(FixedImJsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
FixedImJSON = jsondecode(raw);

% matlabbatch{1}.spm.spatial.coreg.estimate.ref = '<UNDEFINED>';
% matlabbatch{1}.spm.spatial.coreg.estimate.source = '<UNDEFINED>';
% matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};

% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[files_in.In1{1}, ',1']};
% First duplicate the source scan using the prefix string (user-defined)
% Otherwise the CoregEstimate will overwrite the file!!
copyfile(files_in.In2{1},  files_out.In2{1})
copyfile(strrep(files_in.In2{1},'.nii','.json'),  strrep(files_out.In2{1},'.nii','.json'))
% Use the files_out as source (the CoregEstimate will overwrite the file)
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[files_out.In2{1}, ',1']};
if ~isempty(files_in.In3{1})
    other = {};
    for i=1:length(files_in.In3)
        header = niftiinfo(files_in.In3{i});
        % First duplicate all the 'Other' scans using the prefix string (user-defined)
        % Otherwise the CoregEstimate will overwrite the raw files!!
        copyfile(files_in.In3{i},  files_out.In3{i})
        if isfile(strrep(files_in.In3{i},'.nii','.json'))
            copyfile(strrep(files_in.In3{i},'.nii','.json'),  strrep(files_out.In3{i},'.nii','.json'))
        end
        if length(header.ImageSize) > 3
            for j=1:header.ImageSize(4)
                other= [other, [files_in.In3{i}, ',', num2str(j)]];
            end
        else
            other = [other, [files_in.In3{i}, ',1']];
        end
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.other = other';
end

header = niftiinfo(files_in.In2{1});
if length(header.ImageSize) > 3
    for j=2:header.ImageSize(4)
        other= [other, [files_in.In2{1}, ',', num2str(j)]];
    end
end

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = opt.Function;
if strcmp(opt.Separation, 'Auto= [slice thickness voxel_size voxel_size/2]') 
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [FixedImJSON.SliceThickness.value, FixedImInfo.PixelDimensions(2)*10, FixedImInfo.PixelDimensions(3)/2*10];
else
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = str2num(opt.Separation);
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = str2num(opt.Tolerence); %#ok<*ST2NM>
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = str2num(opt.Hist_Smooth);


[SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});


close(SPMinter)
close(SPMgraph)


