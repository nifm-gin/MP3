function [files_in,files_out,opt] = Module_Coreg_Est(files_in,files_out,opt)

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
    module_option(:,6)   = {'Separation','4 2 1 0.5 0.1'};
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
    user_parameter(:,9)   = {'       .Tolerances','char','','Tolerence','','',...
        'The average distance between sampled points (in mm).  Can be a vector to allow a coarse registration followed by increasingly fine ones.'};
    user_parameter(:,10)  = {'       .Separation','numeric','','Separation','','',...
        {'The accuracy for each parameter.  Iterations stop when differences between successive estimates are less than the required tolerance.'
        'Default Value :'
        '     - Human data : 4 2'
        '     - Rat data : 4 2 1 0.5 0.1'}};
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
    if isfield(files_in, 'In3')
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
                    f_out = [char(tags_out_In3.Path), char(tags_out_In3.Patient), '_', char(tags_out_In3.Tp), '_', char(tags_out_In3.SequenceName), '.nii'];
                     tags_out_In3.Filename = categorical(cellstr([char(tags_out_In3.Patient), '_', char(tags_out_In3.Tp), '_', char(tags_out_In3.SequenceName)]));
                else
                    tags_out_In3.Path = categorical(cellstr([opt.folder_out, filesep]));
                    f_out = [char(tags_out_In3.Path), char(tags_out_In3.Patient), '_', char(tags_out_In3.Tp), '_ROI_', char(tags_out_In3.SequenceName), '.nii'];
                    tags_out_In3.Filename = categorical(cellstr([char(tags_out_In3.Patient), '_', char(tags_out_In3.Tp), '_ROI_', char(tags_out_In3.SequenceName)]));
                end
                files_out.In3{i} = f_out;
                opt.Table_out = [opt.Table_out ; tags_out_In3];
            end
        end
    end
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Coreg_Est:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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


matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[files_in.In1{1}, ',1']};
% First duplicate the source scan using the prefix string (user-defined)
% Otherwise the CoregEstimate will overwrite the file!!
copyfile(files_in.In2{1},  files_out.In2{1})
copyfile(strrep(files_in.In2{1},'.nii','.json'),  strrep(files_out.In2{1},'.nii','.json'))
%% Json Processing
[path, name, ~] = fileparts(files_in.In2{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In2{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)
% if the image of reference is > to 3D need to split the data
% header = niftiinfo(files_out.In2{1});
header = spm_vol(files_out.In2{1});
if numel(header) > 1
     other = {};
    for j=2:numel(header)
        other= [other, [files_out.In2{1}, ',', num2str(j)]];
    end
end
% Use the files_out as source (the CoregEstimate will overwrite the file)
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[files_out.In2{1}, ',1']};
if isfield(files_in, 'In3')
    if ~exist('other', 'var')
        other = {};
    end
    for i=1:length(files_in.In3)
        if ~isempty(files_in.In3{i})
            %header = niftiinfo(files_in.In3{i});
            header = spm_vol(files_in.In3{i});
            % First duplicate all the 'Other' scans using the prefix string (user-defined)
            % Otherwise the CoregEstimate will overwrite the raw files!!
            copyfile(files_in.In3{i},  files_out.In3{i})
            if isfile(strrep(files_in.In3{i},'.nii','.json'))
                %% Json Processing
                [path, name, ~] = fileparts(files_in.In3{i});
                jsonfile = [path, '/', name, '.json'];
                J = ReadJson(jsonfile);

                J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

                [path, name, ~] = fileparts(files_out.In3{i});
                jsonfile = [path, '/', name, '.json'];
                WriteJson(J, jsonfile)
            end
            if numel(header) > 1
                for j=1:numel(header)
                    other= [other, [files_out.In3{i}, ',', num2str(j)]];
                end
            else
                other = [other, [files_out.In3{i}, ',1']];
            end
        end
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.other = other';
end



matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = str2num(opt.Separation);

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = str2num(opt.Tolerence); %#ok<*ST2NM>
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = str2num(opt.Hist_Smooth);


% [SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});


%close(SPMinter)
%close(SPMgraph)


