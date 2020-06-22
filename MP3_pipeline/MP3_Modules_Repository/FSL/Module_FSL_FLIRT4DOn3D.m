function [files_in,files_out,opt] = Module_FSL_FLIRT4DOn3D(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

%     %   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','Prefix'};
    module_option(:,4)   = {'Registration_type','Affine (12 parameter model)'};
    module_option(:,5)   = {'Cost_Function','Normalised Mutual information'};
    module_option(:,6)   = {'Search','Incorrect Orentation'};
    module_option(:,7)   = {'Interpolation','Trilinear'};
    module_option(:,8)   = {'output_filename_ext','Coreg'};
    module_option(:,9)   = {'RefInput',2};
    module_option(:,10)   = {'InputToReshape',1};
    module_option(:,11)   = {'Table_in', table()};
    module_option(:,12)   = {'Table_out', table()};
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
    'FLIRT - FMRIB''s Linear Image Registration Tool V6.0. for register 4D or 5D images on 3D'
    ''
    'This function come from the FSL toolbox'
    ''
    'A copy of the ''Source'' and '' Other Images'' is performed before the registration process'
    'This module DO NOT overwrite the orignal files' 
    ''
    'Registration parameters are stored in the headers of the "source" and the "other" images.'
    }'};
    user_parameter(:,2)   = {'Reference Image','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
         'This is the 3D image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
    user_parameter(:,3)   = {'Source Image','1Scan','','',{'SequenceName'},'Mandatory',...
         'This is the 4D/5D image that is jiggled about to best match the reference.'};
    user_parameter(:,4)   = {'Other Images','XScanOrXROI','','',{'SequenceName'},'Optional',...
         'These are any images that need to remain in alignment with the source image.'};
    user_parameter(:,5)   = {'Parameters','','','','','',...
        'Various registration options, which are passed to the Powell optimisation algorithm.'};
    user_parameter(:,6)   = {'       .Registration Type','cell',{'Affine (12 parameter model)','Traditional (9 parameter model)', 'Global Rescale (7 parameter model)', 'Rigid Body (6 parameter model)', 'Translation Only (3 parameter model)'},'Registration_type','','',...
        'Please select a type of registration (between translational only to affine)'};
    user_parameter(:,7)   = {'       .Cost Function','cell',{'Normalised Mutual information','Correlation Ratio', 'Mutual Information', 'Least Square (intra-modal)'},'Cost_Function','','',...
        'Registration involves finding parameters that either maximise or minimise some objective function. See FSL-FLIRT documentation for more information'};
    user_parameter(:,8)   = {'       .Search','cell',{'Incorrect Orentation','Already virtually aligned (no search)', 'Not aligned, but same orientation'},'Search','','',...
        'This option is used to help the algorithm to converge'};
     user_parameter(:,9)  = {'       .Interpolation','cell',{'Nearest neighbour', 'Trilinear', 'Spline'},'Interpolation','','',...
        'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation is slower because it uses more neighbours.'};
    user_parameter(:,10)  = {'       .Filename prefix','char','','output_filename_ext','','',...
        'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''Coreg''.'};
    

    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional', 'Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)', 'VariableNames', VariableNames);
 %%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
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
                    f_out = [char(tags_out_In3.Path), char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-', char(tags_out_In3.SequenceName), '.nii'];
                    tags_out_In3.Filename = categorical(cellstr([char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-', char(tags_out_In3.SequenceName)]));
                else
                    tags_out_In3.Path = categorical(cellstr([opt.folder_out, filesep]));
                    f_out = [char(tags_out_In3.Path), char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-ROI-', char(tags_out_In3.SequenceName), '.nii'];
                    tags_out_In3.Filename = categorical(cellstr([char(tags_out_In3.Patient), '-', char(tags_out_In3.Tp), '-ROI-', char(tags_out_In3.SequenceName)]));
                end
                files_out.In3{i} = f_out;
                opt.Table_out = [opt.Table_out ; tags_out_In3];
            end
        end
    end
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_FSL_FLIRT:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end

[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_FSL_FLIRT:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

switch opt.Registration_type
    case 'Affine (12 parameter model)'
        Registration_type = '12';
    case 'Traditional (9 parameter model)'
        Registration_type = '9';
    case 'Global Rescale (7 parameter model)'
        Registration_type = '7';
    case 'Rigid Body (6 parameter model)'
        Registration_type = '6';
    case 'Translation Only (3 parameter model)'
        Registration_type = '3';
end

switch opt.Search
    case 'Incorrect Orentation'
        search = '-searchrx -180 180 -searchry -180 180 -searchrz -180 180';
    case 'Already virtually aligned (no search)'
        search = '-searchrx 0 0 -searchry 0 0 -searchrz 0 0';
    case 'Not aligned, but same orientation'
        search = '-searchrx -90 90 -searchry -90 90 -searchrz -90 90';
end

switch opt.Cost_Function
    case 'Normalised Mutual information'
        Cost_Function = 'normmi';
    case 'Correlation Ratio'
        Cost_Function = 'corratio';
    case 'Mutual Information'
        Cost_Function = 'mutualinfo';
    case 'Least Square (intra-modal)'
         Cost_Function = 'leastsq';
end

switch opt.Interpolation
    case 'Trilinear'
        Interpolation = 'trilinear';
    case 'Nearest neighbour'
        Interpolation = 'nearestneighbour';
    case 'Spline'
        Interpolation = 'interp spline';
end

% set FSL environment
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

% Keep only the first volume of the 4D/5D image in files_in.In2
info = niftiinfo(files_in.In2{1});
N = niftiread(info);


% Prepare  FSL-FLIRT command
cmd = strcat('/usr/local/fsl/bin/flirt', {' -in '}, files_in.In2{:},...
    {' -ref '}, files_in.In1{:},...
    {' -out '}, files_out.In2{:},...
    {' -omat '}, strrep(files_out.In2{:}, '.nii', '.mat'),...
    {' -bins 256'},...
    {' -cost '}, Cost_Function,...
    {' '}, search,...
    {' -dof  '}, Registration_type, ...
    {' -interp '}, Interpolation);

% execute the command
system(cmd{:});

%unzip the ouput file
gunzip(strrep(files_out.In2{:}, '.nii', '.nii.gz'))

% copy the corresponding json file
%copyfile(strrep(files_in.In2{:},'.nii','.json'),  strrep(files_out.In2{:},'.nii','.json'))

%% Json Processing
[path, name, ~] = fileparts(files_in.In2{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In2{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)


%[Mov] = ExtractMovement('Patient1-J3-T2_EG_AXIAL_HEMO_20180504-163255801.nii','Patient1_J3_CoregADC.nii');
%ApplyMovement('ADC_moved.nii',mat')      

 % if the User wants to apply the same transformation to other scans
if isfield(files_in, 'In3')
    if ~exist('other', 'var')
        other = {};
    end
    % compute the transformation matrice which will be apply
%    [Mov] = ExtractMovement(files_out.In2{:}, files_in.In2{:});

    %
    for i=1:length(files_in.In3)
        if ~isempty(files_in.In3{i})
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
                                copyfile(strrep(files_in.In3{i},'.nii','.json'),  strrep(files_out.In3{i},'.nii','.json'))
            end
            % first apply reslice the 'other scan' to the image moved
%             info3 = niftiinfo(files_out.In3{i});
%             S = info3.ImageSize;
%             if length(S)>3
%                 for j=1:S(4)
%                     for k=1:S(5)
%                         other= [other, [files_out.In3{i}, ',', num2str(j), ',', num2str(k)]];
%                     end
%                 end
%             else
%                 other = [other, [files_out.In3{i}, ',1']];
%             end
            
            info3 = niftiinfo(files_out.In3{i});
            S = info3.ImageSize;
            if length(S)>4
               N = niftiread(info3);
               NewN = reshape(N, [S(1:3), S(4)*S(5)]);
               info3bis = info3;
               info3bis.ImageSize = size(NewN);
               info3bis.PixelDimensions = info3.PixelDimensions(1:4);
               niftiwrite(NewN, files_out.In3{i}, info3bis);                
            end

            header = spm_vol(files_out.In3{i});
            if numel(header) > 1
                for j=1:numel(header)
                    other= [other, [files_out.In3{i}, ',', num2str(j)]];
                end
            else
                other = [other, [files_out.In3{i}, ',1']];
            end
            
            
            P.which = 1; % spm_reslice's option
            spm_reslice([files_in.In2(1), other]', P);
            % delete unusefull file
            delete(files_out.In3{i})
            
            % Then apply the transformation calculated during the FLIRT process
            [filepath,name,ext] = fileparts(files_out.In3{i}); 

            
            cmd_other = strcat('/usr/local/fsl/bin/flirt',...
                {' -in '}, fullfile(filepath,['r', name, ext]),...   
             ' -applyxfm',...
             {' -init '}, strrep(files_out.In2{:}, '.nii', '.mat'),...
             {' -out '}, files_out.In3{i},...
             {' -paddingsize 0.0 -interp trilinear'},...
             {' -ref '}, files_in.In1{:});
         % execute the command
         system(cmd_other{:});
         
         % delete unusefull file
         delete(fullfile(filepath,['r', name, ext]))
         
         %unzip the ouput file
         gunzip(strrep(files_out.In3{i}, '.nii', '.nii.gz'))
         
         
        info4 = niftiinfo(files_out.In3{i});
        S2 = info4.ImageSize;
        if length(S)>4
           N = niftiread(info4);
           NewN = reshape(N, [S2(1:3), S(4:5)]);
           info4bis = info4;
           info4bis.ImageSize = size(NewN);
           info4bis.PixelDimensions = [info4.PixelDimensions(1:3), info3.PixelDimensions(4:5)];
           niftiwrite(NewN, files_out.In3{i}, info4bis);                
        end
         
         
         % image and the image moved
         %ApplyMovement(files_out.In3{i},Mov)
         
         
         
            
          %  ApplyMovement(files_out.In3{i},Mov)
%           {'/usr/local/fsl/bin/flirt' ...
%               '-in /Users/Lemasson/Data_non_synchro/2018_Collab_Nancy/MIA_data/Tmp/ADC_moved.nii'...
%               '-applyxfm -init /Users/Lemasson/Data_non_synchro/2018_Collab_Nancy/MIA_data/Tmp/Patient1_J3_CoregADC.mat'...
%               '-out /Users/Lemasson/Data_non_synchro/2018_Collab_Nancy/MIA_data/Tmp/toto.nii'...
%               '-paddingsize 0.0 -interp trilinear'...
%               '-ref /Users/Lemasson/Data_non_synchro/2018_Collab_Nancy/MIA_data/Tmp/Patient1-J3-T2_EG_AXIAL_HEMO_20180504-163255801.nii'}
          
        end
    end
end

