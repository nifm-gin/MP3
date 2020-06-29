function [files_in,files_out,opt] = Module_Reslice(files_in,files_out,opt)

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
    module_option(:,5)   = {'Smoothing','No'};
    module_option(:,6)   = {'Interpolation','4th Degree B-Spline'};
    module_option(:,7)   = {'Wrapping','No wrap'};
    module_option(:,8)   = {'Masking','Dont mask images'};
    module_option(:,9)   = {'output_filename_ext','r'};
    module_option(:,10)   = {'RefInput',2};
    module_option(:,11)   = {'InputToReshape',1};
    module_option(:,12)   = {'Table_in', table()};
    module_option(:,13)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    %
    %% list of everything displayed to the user associated to their 'type'
    % --> user_parameter(1,:) = user_parameter_list
    % --> user_parameter(2,:) = user_parameter_type
    % --> user_parameter(3,:) = parameter_default
    % --> user_parameter(4,:) = psom_parameter_list
    % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
    % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional.
    % --> user_parameter(7,:) = Help : text data which describe the parameter (it
    % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'Reslice images to match voxel-for-voxel with an image defining some space.'
        'The resliced images are named the same as the originals except that they are prefixed by ?r?.'
        ''
        'This branch contains 3 tiems'
        '     Image defining space'
        '     Image to reslice'
        '     Reslice option'
        }'};
    user_parameter(:,2)   = {'Reference Image','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
        'This is the image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
    user_parameter(:,3)   = {'Image to reslice','1Scan','','',{'SequenceName'},'Mandatory',...
        'This is the image is reslice to match to the reference.'};
    user_parameter(:,4)   = {'Parameters','','','','','',''};
    user_parameter(:,5)   = {'    Pre-processing: Do you want to smooth your data first?','cell',{'No', 'Yes'},'Smoothing','','',...
        'If the spatial resolution of the reference image is lower than the image to reslice we recommand to smooth the image to reslice (by selecting yes) before the reslicing'};
    user_parameter(:,6)   = {'    Reslice Options','', '','','','',...
        {'Various registration options'
        'This branch contains 4 items:'
        '     Interpolation'
        '     Warpping'
        '     Masking'
        '     Filename prefix'}};
    user_parameter(:,7)  = {'       .Interpolation','cell',{'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'},'Interpolation','','',...
        'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation is slower because it uses more neighbours.'};
    user_parameter(:,8)  = {'       .Wrapping','cell',{'Mask images', 'Dont mask images'},'Masking','','',...
        'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};
    user_parameter(:,9)  = {'       .Masking','cell',{'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'},'Wrapping','','',...
        'These are typically: No wrapping - for PET or images that have already been spatially transformed.  Wrap in  Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'};
    user_parameter(:,10)  = {'       .Filename prefix','char','','output_filename_ext','','',...
        'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''r''.'};
    
    
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    
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
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Reslice:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

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

% set the reference image
matlabbatch{1}.spm.spatial.coreg.write.ref = {[files_in.In1{1}, ',1']};

% First duplicate the source scan because the output will be automatically
% gererate in the same directory (/Tmp) and the fill the souce image to SPM
% stucture
source_filename = fullfile(char(opt.Table_out.Path),[char(opt.Table_in.Patient(2)), '_', char(opt.Table_in.Tp(2)),'_', char(opt.Table_in.SequenceName(2)),'.nii']);
if ~strcmp(files_in.In2{1}, source_filename)
    copyfile(files_in.In2{1},  source_filename)
    copyfile(strrep(files_in.In2{1},'.nii','.json'), strrep(source_filename, '.nii', '.json'))
end

data_to_reslice = niftiread(source_filename);
data_to_reslice_info = niftiinfo(source_filename);

% % first, reshape the data if needed (in order to only have 4D data)
if length(data_to_reslice_info.ImageSize) > 4
    data_to_reslice = reshape(data_to_reslice, [data_to_reslice_info.ImageSize(1:3) prod(data_to_reslice_info.ImageSize(4:end))]);
    data_to_reslice_info_reshaped = data_to_reslice_info;
    data_to_reslice_info_reshaped.ImageSize = [data_to_reslice_info.ImageSize(1:3) prod(data_to_reslice_info.ImageSize(4:end))];
    data_to_reslice_info_reshaped.PixelDimensions = data_to_reslice_info.PixelDimensions(1:4);
    
    % overwrite the tmp file by the smoothed one
    niftiwrite(data_to_reslice,source_filename, data_to_reslice_info_reshaped)
end

if strcmp(opt.Smoothing, 'Yes')
    reference_image_info =  niftiinfo(files_in.In1{1});
    % the sigma depends of both resolution (image to reslice and reference
    % image) multiply by a factor 0.4
    Sigma = (reference_image_info.PixelDimensions(1:3)./data_to_reslice_info.PixelDimensions(1:3) ) * 0.4;
    for j=1:size(data_to_reslice,4)
        data_to_reslice(:,:,:,j) = imgaussfilt3(data_to_reslice(:,:,:,j), Sigma);
    end
    % overwrite the tmp file by the smoothed one
    niftiwrite(data_to_reslice,source_filename, data_to_reslice_info)

end


header = spm_vol(source_filename);
if numel(header) > 1
     other = {};
    for j=1:numel(header)
        other= [other, [source_filename, ',', num2str(j)]];
    end
      matlabbatch{1}.spm.spatial.coreg.write.source = other';
else
    matlabbatch{1}.spm.spatial.coreg.write.source = {[source_filename, ',1']};
end


%% Json Processing
[path, name, ~] = fileparts(files_in.In2{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

[path, name, ~] = fileparts(files_out.In2{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

%define options
% Type of interpolation
switch opt.Interpolation
    case 'Nearest neighbour'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    case 'Trilinear'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 7;
end
%  Type of Warpping
switch opt.Wrapping
    case 'No wrap'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    case 'Warp X'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 1 1];
end

%% always set mask to 0
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = opt.output_filename_ext;

%[SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});

[path,name,ext] = fileparts(source_filename);

% reshape data if needed (in case the dim > 4)

data_smoothed = niftiread(fullfile(path, [opt.output_filename_ext, name,ext]));
data_smoothed_info = niftiinfo(fullfile(path, [opt.output_filename_ext, name,ext]));
%% reshap data if needed (if dim >4)
if length(data_to_reslice_info.ImageSize) > 4
    data_smoothed = reshape(data_smoothed, [data_smoothed_info.ImageSize(1:3) data_to_reslice_info.ImageSize(4:end)]);
    data_smoothed_info.ImageSize = [data_smoothed_info.ImageSize(1:3) data_to_reslice_info.ImageSize(4:end)];
    tmp = zeros([1 length(data_to_reslice_info.ImageSize)]);
    tmp(1:length(data_smoothed_info.PixelDimensions))=data_smoothed_info.PixelDimensions ;
    data_smoothed_info.PixelDimensions = tmp;
    
    % overwrite the tmp file by the smoothed one
    niftiwrite(data_smoothed,fullfile(path, [opt.output_filename_ext, name,ext]), data_smoothed_info)
end


movefile(fullfile(path, [opt.output_filename_ext, name,ext]), files_out.In2{1});
% close(SPMinter)
% close(SPMgraph)

