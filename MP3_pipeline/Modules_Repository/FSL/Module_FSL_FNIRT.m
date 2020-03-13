function [files_in,files_out,opt] = Module_FSL_FNIRT(files_in,files_out,opt)

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
    
    module_option(:,4)   = {'output_filename_prefix','Coreg_fnirt_'}; 
    module_option(:,5)   = {'output_filename_suffixCoef','_coef'};
    module_option(:,6)   = {'infwhm','2,1,0,0'};
    module_option(:,7)   = {'reffwhm','1,0.5,0,0'};
    module_option(:,8)   = {'warpres','0.25,0.25,0.25'};
    module_option(:,9)   = {'biasres','0.25,0.25,0.25'};
    module_option(:,10)   = {'regmod','bending_energy'};
    module_option(:,11)   = {'intmod','global_non_linear_with_bias'}; 
    
    module_option(:,12)   = {'RefInput',2};
    module_option(:,13)   = {'InputToReshape',1};
    module_option(:,14)   = {'Table_in', table()};
    module_option(:,15)   = {'Table_out', table()};
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
    'FNIRT - High resolution nonlinear registration with simultaneous modelling of intensities'
    ''
    'This function come from the FSL toolbox'
    ''
    'https://www.biorxiv.org/content/10.1101/646802v1.full'
    ''
    'A copy of the ''Source'' and '' Other Images'' is performed before the registration process'
    'This module DO NOT overwrite the orignal files' 
    ''
    'Registration parameters are stored in the headers of the "source" and the "other" images.'
    }'};
    user_parameter(:,2)   = {'Reference Image','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
         'This is the image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
    user_parameter(:,3)   = {'Source Image','1Scan','','',{'SequenceName'},'Mandatory',...
         'This is the image that is jiggled about to best match the reference.'};
    user_parameter(:,4)   = {'Other Images','XScanOrXROI','','',{'SequenceName'},'Optional',...
         'These are any images that need to remain in alignment with the source image.'};
    user_parameter(:,5)   = {'Parameters','','','','','',...
        'Various registration options, which are passed to the algorithm.'};
    
    user_parameter(:,6)  = {'       .Output filename prefix','char','','output_filename_prefix','','',...
        'Specify the string to be prepended to the filenames of the coregistered image file(s). Default prefix is ''Coreg_fnirt_''.'};
    
    user_parameter(:,7)  = {'       .Output filename suffix for the Coefficients ','char','Coeficients','output_filename_suffixCoef','','',...
        {'Specify the name of output file with field coefficients.'
        'Default filename is ''_coef''.'}'};
    user_parameter(:,8)  = {'       .FWHM (in mm) of gaussian smoothing kernel for input volume','numeric','','infwhm','','',...
        {'FWHM (in mm) of gaussian smoothing kernel for input volume'
        'Default Value :'
        '     - Human data : 6,4,2,2'
        '     - Mouse data : 2,1,0,0'}};
     user_parameter(:,9)  = {'       .FWHM (in mm) of gaussian smoothing kernel for ref volume','numeric','','reffwhm','','',...
        {'FWHM (in mm) of gaussian smoothing kernel for ref volume '
        'Default Value :'
        '     - Human data : 4,2,0,0'
        '     - Mouse data : 1,0.5,0,0'}};
    user_parameter(:,10)  = {'       .Resolution (in mm) of warp basis','numeric','','warpres','','',...
        {'(approximate) resolution (in mm) of warp basis in x-, y- and z-direction'
        'Default Value :'
        '     - Human data : 10,10,10'
        '     - Mouse data : 0.25,0.25,0.25'}};
     user_parameter(:,11)  = {'       .Resolution (in mm) of bias-field modelling local intensities','numeric','','biasres','','',...
        {'Resolution (in mm) of bias-field modelling local intensities'
        'Default Value :'
        '     - Human data : 50,50,50'
        '     - Mouse data : 0.25,0.25,0.25'}};
    
     user_parameter(:,12)  = {'       .Model for regularisation of warp-field','cell',{'membrane_energy', 'bending_energy'},'regmod','','',...
        {'Model for regularisation of warp-field [membrane_energy bending_energy]'
        'Default : bending_energy'}};
    
      user_parameter(:,13)  = {'       .Model for intensity-mapping ','cell',{'none', 'global_linear', 'global_non_linear', 'local_linear', 'global_non_linear_with_bias', 'local_non_linear'},'intmod','','',...
        {'Model for intensity-mapping [none global_linear global_non_linear local_linear global_non_linear_with_bias local_non_linear]'
        'Default : global_non_linear_with_bias'}};

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
    tags_out_In1 = tags2;
    tags_out_In1.IsRaw = categorical(0);
    tags_out_In1.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_In1.SequenceName = categorical(cellstr([opt.output_filename_prefix, char(tags_out_In1.SequenceName)]));
    tags_out_In1.Filename = categorical(cellstr([char(tags_out_In1.Patient), '_', char(tags_out_In1.Tp), '_', char(tags_out_In1.SequenceName)]));
    f_out = [char(tags_out_In1.Path), char(tags_out_In1.Patient), '_', char(tags_out_In1.Tp), '_', char(tags_out_In1.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    
    tags_out_In2 = tags_out_In1;
    tags_out_In2.SequenceName = categorical(cellstr([char(tags_out_In1.SequenceName) , opt.output_filename_suffixCoef])); %categorical(cellstr([char(opt.Table_in.SequenceName(i)), opt.output_filename_ext_Scan]));
    tags_out_In2.Filename = categorical(cellstr([char(tags_out_In2.Patient), '_', char(tags_out_In2.Tp), '_', char(tags_out_In2.SequenceName)]));
    f_out2 = [char(tags_out_In2.Path), char(tags_out_In2.Patient), '_', char(tags_out_In2.Tp), '_', char(tags_out_In2.SequenceName), '.nii'];
    files_out.In2{1} = f_out2;
    opt.Table_out = [tags_out_In1; tags_out_In2];
    
    
    if isfield(files_in, 'In3')
        for i=1:length(files_in.In3)
            if ~isempty(files_in.In3{i})
                [Path_In3, Name_In3, ~] = fileparts(files_in.In3{i});
                tags3 = opt.Table_in(opt.Table_in.Path == [Path_In3, filesep],:);
                tags3 = tags3(tags3.Filename == Name_In3,:);
                assert(size(tags3, 1) == 1);
                tags_out_In3 = tags3;
                tags_out_In3.IsRaw = categorical(0);
                tags_out_In3.SequenceName = categorical(cellstr([opt.output_filename_prefix, char(tags_out_In3.SequenceName)]));
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
% 
% switch opt.regmod
%     case 'Affine (12 parameter model)'
%         Registration_type = '12';
%     case 'Traditional (9 parameter model)'
%         Registration_type = '9';
%     case 'Global Rescale (7 parameter model)'
%         Registration_type = '7';
%     case 'Rigid Body (6 parameter model)'
%         Registration_type = '6';
%     case 'Translation Only (3 parameter model)'
%         Registration_type = '3';
% end
% 
% switch opt.intmod
%     case 'Incorrect Orentation'
%         search = '-searchrx -180 180 -searchry -180 180 -searchrz -180 180';
%     case 'Already virtually aligned (no search)'
%         search = '-searchrx 0 0 -searchry 0 0 -searchrz 0 0';
%     case 'Not aligned, but same orientation'
%         search = '-searchrx -90 90 -searchry -90 90 -searchrz -90 90';
% end
% set FSL environment
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

% Prepare  FSL-FLIRT command
cmd = strcat('/usr/local/fsl/bin/fnirt',...
    {' --ref='}, files_in.In1{:},...
    {' --in='}, files_in.In2{:},...
    {' --iout='}, files_out.In1{:},...
    {' --cout='}, files_out.In2{:},...
    {' --infwhm='}, opt.infwhm,...
    {' --reffwhm='}, opt.reffwhm,...
    {' --warpres='}, opt.warpres, ...
    {' --biasres='}, opt.biasres, ...
    {' --regmod='}, opt.regmod, ...
    {' --intmod='}, opt.intmod, ...
    {' --verbose'});
% execute the command
system(cmd{:});

%unzip the ouput file
gunzip(strrep(files_out.In1{:}, '.nii', '.nii.gz'))
gunzip(strrep(files_out.In2{:}, '.nii', '.nii.gz'))

% copy the corresponding json file
%copyfile(strrep(files_in.In2{:},'.nii','.json'),  strrep(files_out.In2{:},'.nii','.json'))

%% Json Processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

[path, name, ~] = fileparts(files_in.In2{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
[path, name, ~] = fileparts(files_out.In2{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)
% 
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
%             % Then apply reslice the 'other scan' to the image moved
            P.which = 1; % spm_reslice's option
            % Interpolation option 
            SPM_Inter = 0; % O = 'nearestneighbour';
            P.interp = SPM_Inter;
            spm_reslice({files_in.In2{:}, files_out.In3{i}}', P); %#ok<CCAT>
            % delete unusefull file
            delete(files_out.In3{i})
            
            % Then apply the transformation calculated during the FLIRT process
            [filepath,name,ext] = fileparts(files_out.In3{i});
            cmd_other= strcat('/usr/local/fsl/bin/applywarp',...
                {' --in='}, fullfile(filepath,['r', name, ext]),...
                {' --ref='}, files_in.In1{:},...
                {' --out='}, files_out.In3{i},...
                {' --warp='}, files_out.In2{:});
           % applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --in=my_structural --warp=my_nonlinear_transf --out=my_warped_structural

            % execute the command
            system(cmd_other{:});
            
            % delete unusefull file
            delete(fullfile(filepath,['r', name, ext]))
            
            %unzip the ouput file
            gunzip(strrep(files_out.In3{i}, '.nii', '.nii.gz'))
            
            
            % image and the image moved
            %ApplyMovement(files_out.In3{i},Mov)
 
        end
    end
end

