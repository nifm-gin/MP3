function [files_in,files_out,opt] = Module_Coreg_Est_Res(files_in,files_out,opt)
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
    module_option(:,12)   = {'output_filename_ext','Coreg'};
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
%      % --> user_parameter(5,:) = Scan_Input_DOF
%      % --> user_parameter(6,:) = InputMandatoryOrOptional
%      % --> user_parameter(7,:) = Help : text data which describe the parameter (it
%      % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
    'Within-subject registration using a rigid-body model and image reslicing.'
    ''
    'The registration method used here is based on work by Collignon et al.'
    'The original interpolation method described in this paper has been changed in order to give a smoother cost function.'
    'The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
    ''
    'At the end of coregistration, the voxel-to-voxel affine transformation matrix is displayed, along with the histograms for the images in the original orientations, and the final orientations.'
    'The registered images are displayed at the bottom.'
    ''
    'Registration parameters are stored in the headers of the "source" and the "other" images. These images are also resliced to match the source image voxel-for-voxel.'
    }'};
% 
% 

% if isempty(opt)
    % define every option needed to run this module
    %fields   = {'Type', 'HSize', 'Sigma', 'flag_test' , 'folder_out', 'output_filename_ext'};
%     fields   = {'RefInput', 'InputToReshape', 'NbInput', 'NbOutput', 'folder_out', 'flag_test', 'Execution_Mode','output_filename_ext', 'OutputSequenceName', 'Function', 'Separation', 'Tolerence', 'Hist_Smooth', 'Interpolation', 'Wrapping'};
%     defaults = {2,1, -1, -1, '', true, 'All Database','_Coreg', 'Extension', 'mi', 'Auto= [slice thickness voxel_size voxel_size/2]', '0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001', '7 7', '4th Degree B-Spline', 'No wrap'};
%     opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
%     
%     % list of everything displayed to the user associated to their 'type'
%     user_parameter_list = {'Description';'Execution Mode';'Select one scan or more as input reference image'; 'Select one scan as input to compute and apply the coreg on'; 'Select one scan or more as input image to apply the coreg on'; 'Parameters'; '   .Output filename prefix';  '    Estimation Options';'       .Objective Function';  '       .Separation'; '       .Tolerances'; '       .Histogram Smoothing'; '    Reslice options';'       .Interpolation'; '       .Wrapping'};%; ''; ''};
%     user_parameter_type = {'Text'; 'cell';'1Scan1TPXP'; '1Scan'; 'XScan'; ''; 'char'; '';'cell'; 'char';'numeric';'numeric';'';'cell';'cell'};%; 'logical'; 'char'};
%     parameter_default = {help;{'All Database','Through all sessions of one Patient','Through Each Session'};'';'';''; ''; '';''; {'mi','ncc', 'nmi', 'ecc'}; ''; ''; ''; '';{'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'}; {'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'}};%; 'Dont Show';'Dont Show'};
%     psom_parameter_list = {'';'Execution_Mode';'';'';''; ''; 'output_filename_prefix'; '';'Function'; 'Separation';'Tolerence' ; 'Hist_Smooth';'';'Interpolation'; 'Wrapping'};%;'flag_test'; 'folder_out' };
%     scans_input_DOF = {''; '';{'SequenceName', 'Tp', 'Patient'}; {'SequenceName'}; {'SequenceName'}; ''; ''; '';''; '';'';'';'';'';''};
%     VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF'};
%     %opt.table = table(categorical(user_parameter_list), categorical(user_parameter_type), categorical(parameter_default), categorical(psom_parameter_list), 'VariableNames', VariableNames);
%     opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, scans_input_DOF, 'VariableNames', VariableNames);
    %opt.parameter_link_psom = {'output_filename_ext', '   .Output filename extension'; 'Type', '   .Type'; 'HSize','   .HSize'; 'Sigma', '   .Sigma'};

    %% Benjamin Modifications
%     }'};

%     user_parameter(:,3)   = {'Reference Image','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},...
%          'This is the image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
%     user_parameter(:,4)   = {'Source Image','1Scan','','',{'SequenceName'},...
%          'This is the image that is jiggled about to best match the reference.'};
%     user_parameter(:,5)   = {'Other Images','XScan','','',{'SequenceName'},...
%          'These are any images that need to remain in alignment with the source image.'};

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
    user_parameter(:,12)  = {'    Reslice options','','','','','',...
        'These images are resliced to the same dimensions, voxel sizes, orientation etc as the space defining image.'};
    user_parameter(:,13)  = {'       .Interpolation','cell',{'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'},'Interpolation','','',...
        'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation is slower because it uses more neighbours.'};
    user_parameter(:,14)  = {'       .Masking','cell',{'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'},'Wrapping','','',...
        'These are typically: No wrapping - for PET or images that have already been spatially transformed.  Wrap in  Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'};
   user_parameter(:,15)  = {'       .Wrapping','cell',{'Mask images', 'Dont mask images'},'Masking','','',...
       'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};         
    user_parameter(:,16)  = {'       .Filename prefix','char','','output_filename_ext','','',...
        'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''Coreg''.'};


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
    error('Smoothing:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
%if ~ischar(files_in) 
%    error('files in should be a char');
%end

[path_nii,name_nii,ext_nii] = fileparts(files_in.In1{1});
%if ~strcmp(ext_nii, '.nii')
%     error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, name_nii, files_in);  
%end

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
%if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
%    opt.folder_out = path_nii;    
%end

%if isempty(files_out.filename)
%    files_out.filename = cat(2,opt.folder_out,filesep,name_nii,opt.output_filename_ext,ext_nii);
%end

% if strcmp(files_out, '')
%     files_out = files_in;
%     files_out = rmfield(files_out, 'In1');
%     [path_nii,name_nii,ext_nii] = fileparts(files_in.In2{1});
%     files_out.In2 = {cat(2,path_nii,filesep,opt.output_filename_prefix, '_', name_nii,ext_nii)};
%     for i=1:length(files_out.In3)
%         [path_nii,name_nii,ext_nii] = fileparts(files_in.In3{i});
%         files_out.In3{i}= cat(2,path_nii,filesep,opt.output_filename_prefix,'_', name_nii,ext_nii);
%     end
% end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FixedImInfo = niftiinfo(files_in.In1{1});
[path, name, ext] = fileparts(files_in.In1{1});
FixedImJsonfile = [path, filesep, name, '.json'];
fid = fopen(FixedImJsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
FixedImJSON = jsondecode(raw);




matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[files_in.In1{1}, ',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[files_in.In2{1}, ',1']};
if ~isempty(files_in.In3{1})
    for i=1:length(files_in.In3)
        files_in.In3{i}= [files_in.In3{i}, ',1'];
    end
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = files_in.In3;
end

matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = opt.Function;
if strcmp(opt.Separation, 'Auto= [slice thickness voxel_size voxel_size/2]') 
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [FixedImJSON.SliceThickness, FixedImInfo.PixelDimensions(2)*10, FixedImInfo.PixelDimensions(3)/2*10];
    %matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [fixed.uvascim.image.reco.thickness fixed.uvascim.image.reco.fov(1)/fixed.uvascim.image.reco.no_samples  fixed.uvascim.image.reco.fov(1)/fixed.uvascim.image.reco.no_samples/2];
else
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = str2num(opt.Separation);
end
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = str2num(opt.Tolerence); %#ok<*ST2NM>
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = str2num(opt.Hist_Smooth);
%define options
% Type of interpolation
switch opt.Interpolation
    case 'Nearest neighbour'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
    case 'Trilinear'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7;
end
%  Type of Warpping
switch opt.Wrapping
    case 'No wrap'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    case 'Warp X'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 1 1];
end

%% always set mask to 0
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask= 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = opt.output_filename_ext;

[SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});

%date_str = date;
% if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
%     movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
%         [PATHSTR, filesep, NAME_scan_to_coreg(1:end-8) '_SPMcoreg.ps'], 'f'); 
% end




%% JSON de l'input 2
[path, name, ext] = fileparts(files_in.In2{1});
SpmOutputFile  = [path, filesep, opt.output_filename_ext, name, ext];
if exist(SpmOutputFile, 'file') ~=2
    error('Cannot find the file %s', SpmOutputFile);
end
movefile(SpmOutputFile,files_out.In2{1});

jsonfile = [path, filesep, name, '.json'];
fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);

J = jsondecode(raw);

JMod = jsonencode(J);
[path, name, ext] = fileparts(files_out.In2{1});
jsonfile = [path, filesep, name, '.json'];
fidmod = fopen(jsonfile, 'w');
fwrite(fidmod, JMod, 'uint8');
fclose(fidmod);

%% JSON de l'input 3
if isfield(files_out, 'In3')
    for i=1:length(files_out.In3)
        [path, name, ext] = fileparts(files_in.In3{i});
        SpmOutputFile  = [path, filesep, opt.output_filename_ext, name, '.nii'];
        if exist(SpmOutputFile, 'file') ~=2
            error('Cannot find the file %s', SpmOutputFile);
        end
        movefile(SpmOutputFile,files_out.In3{i})
        jsonfile = [path, filesep, name, '.json'];
        fid = fopen(jsonfile, 'r');
        raw = fread(fid, inf, 'uint8=>char');
        fclose(fid);
        J = jsondecode(raw);
        
        JMod = jsonencode(J);
        [path, name, ext] = fileparts(files_out.In3{i});
        jsonfile = [path, filesep, name, '.json'];
        fidmod = fopen(jsonfile, 'w');
        fwrite(fidmod, JMod, 'uint8');
        fclose(fidmod);
    end
end



close(SPMinter)
close(SPMgraph)


% N = niftiread(files_in);
% info = niftiinfo(files_in);
% [path, name, ext] = fileparts(files_in);
% jsonfile = [path, '/', name, '.json'];
% fid = fopen(jsonfile, 'r');
% raw = fread(fid, inf, 'uint8=>char');
% fclose(fid);
% %raw = reshape(raw, 1,length(raw));
% J = jsondecode(raw);
% 
% Informations = whos('N');
% FilteredImages = zeros(size(N), Informations.class);
% NbDim = length(size(N));
% if NbDim>4
%     warning('Too much dimensions. This module deals with at most 4 dimensions.')
% end
% h = fspecial(opt.Type,str2double(opt.HSize), str2double(opt.Sigma));
% for i=1:size(N,3)
%     for j=1:size(N,4)
%         FilteredImages(:,:,i,j) = imfilter(N(:,:,i,j), h, 'replicate');
%     end
% end
% 
% info2 = info;
% info2.Filename = files_out;
% info2.Filemoddate = char(datetime('now'));
% info2.Description = [info.Description, 'Modified by Smoothing Module'];
% 
% niftiwrite(FilteredImages, files_out, info2)
% JMod = jsonencode(J);
% [path, name, ext] = fileparts(files_out);
% jsonfile = [path, '/', name, '.json'];
% fidmod = fopen(jsonfile, 'w');
% fwrite(fidmod, JMod, 'uint8');
% fclose(fidmod);