function [files_in,files_out,opt] = Module_DTI_1sens(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the module's parameters with default values
if isempty(opt)
    %% Define every option needed to run this module
    %     % --> module_option(1,:) = field names
    %     % --> module_option(2,:) = defaults values
    module_parameters(:,1)   = {'OutputSequenceName','AllName'};
    module_parameters(:,2)   = {'output_filename_ext_FA','dti_FA'};
    module_parameters(:,3)   = {'output_filename_ext_MD','dti_MD'};
    module_parameters(:,4)   = {'output_filename_ext_L1','dti_L1'};
    module_parameters(:,5)   = {'output_filename_ext_L2','dti_L2'};
    module_parameters(:,6)   = {'output_filename_ext_L3','dti_L3'};
    module_parameters(:,7)   = {'FSL_path',''};
    
     %% System parameters : Do not modify without understanding the behaviour of the software.
    
    system_parameters(:,1)   = {'RefInput',1};
    system_parameters(:,2)   = {'InputToReshape',1};
    
    
    %% Initialisation parameters : Do not modify without understanding the behaviour of the software.
    
    initialisation_parameters(:,1)   = {'folder_out',''};
    initialisation_parameters(:,2)   = {'flag_test', true};
    initialisation_parameters(:,3)   = {'Table_in', table()};
    initialisation_parameters(:,4)   = {'Table_out', table()};
    
    Parameters = [module_parameters, system_parameters, initialisation_parameters];
    
    opt.Module_settings = psom_struct_defaults(struct(),Parameters(1,:),Parameters(2,:));
    
     %% Each line displayed to the user :
    
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'Extraction of FA, MD and other diffusion maps from 2 opposite direction DTI scans'
        '*This module calls FSL'
        }'};
    user_parameter(:,2)   = {'Select the DTI  scan','1Scan','','',{'SequenceName'}, 'Mandatory','Please copy the .bvals and .bvecs files associated to this file into the Raw_data folder of this project'};
    user_parameter(:,3)   = {'Indicate path to FSL','char','','FSL_path','','Mandatory','Result of <which fsl> on any terminal (usually on Linux: /usr/local/fsl/bin/)'};
    user_parameter(:,4)   = {'Parameters','','','','', '', ''};
    user_parameter(:,5)   = {'   .Output filename extension FA','char', 'dti_FA','output_filename_ext_FA','','',''};
    user_parameter(:,6)   = {'   .Output filename extension MD','char','dti_MD','output_filename_ext_MD','','',''};
    user_parameter(:,7)   = {'   .Output filename extension L1','char','dti_L1','output_filename_ext_L1','','',''};
    user_parameter(:,8)   = {'   .Output filename extension L2','char','dti_L2','output_filename_ext_L2','','',''};
    user_parameter(:,9)   = {'   .Output filename extension L3','char','dti_L3','output_filename_ext_L3','','',''};
    
    % Concatenate these user_parameters, and store them in opt.table
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', ...
        user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    %%
    
    % Initialize to an empty string the names of the input and output
    % files.
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename_ext_FA, opt.output_filename_ext_MD, opt.output_filename_ext_L1, opt.output_filename_ext_L2, opt.output_filename_ext_L3};

if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in(1,:);
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), opt.NameOutFiles{i}]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
end


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_DTI_2sens:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%% load input Nii file
input(1).nifti_header = spm_vol(files_in.In1{1});

J1 = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));

inFile = read_volume(input(1).nifti_header, input(1).nifti_header, 0);

%nombval = char(table_out_tmp.Path)

setenv('FSLOUTPUTTYPE','NIFTI');

%% fslroi: extract b0_P and b0_A files
nombval=strrep(files_in.In1{1},'.nii', '.bval');
nombvec=strrep(files_in.In1{1},'.nii', '.bvec');


cd(opt.folder_out)
current_folder = fullfile(opt.folder_out, [char(opt.Table_in.Patient(1)), '_', char(opt.Table_in.Tp(1))]);
mkdir(current_folder)
cd(current_folder)

name_split = strsplit(files_in.In1{:},filesep);
copyfile(files_in.In1{1}, [current_folder, filesep, name_split{end}])
fidval = fopen([current_folder, filesep, strrep(name_split{end}, '.nii', '.bval')], 'w');
fprintf(fidval, '%d ', J1.Bval.value.');
fclose(fidval);




formatSpec = repmat('%f ', 1,size(J1.Bvec.value,2)+1);
formatSpec(end-1:end) = '';
formatSpec(end-1:end) = '\n';
%Test2 = compose(formatSpec, J1.Bvec.value.');


fidvec = fopen([current_folder, filesep, strrep(name_split{end}, '.nii', '.bvec')], 'w');
fprintf(fidvec, formatSpec, J1.Bvec.value.');
fclose(fidvec);

% 
% Space = repmat({' '}, size(J1.Bvec.value,1), size(J1.Bvec.value,2)-1);
% RetourChariot = repmat({'\n'}, 3,1);
% Test = [Space, RetourChariot].';
% Cell = num2str(J1.Bvec.value).';
% col_interleave = reshape([Cell(:) Test(:)]',2*size(Cell,1), [])';



if ~exist([current_folder '/acqparams.txt'],'file')
    system(['printf "0 1 0 0.05\n0 -1 0 0.05" > ' current_folder '/acqparams.txt']);
end

index = repmat('1 ', 1,size(J1.Bvec.value,2));
index(end) = '';

if ~exist([current_folder '/index.txt'],'file')
    system(['printf "' index '" > ' current_folder '/index.txt']);
end

%% Skull extraction
disp('---------Skull extraction----------');
    status = system([opt.FSL_path 'bet ' name_split{end} ' b0_ss.nii -f 0.1 -m -Z']);
    if status ==0
        disp('DONE');
    else
        warning('FAIL');
    end

%% Eddy 
%     system([FSLcommand 'eddy_correct ' strrep(indir,' ','\ ') Nomdti ' ' strrep(out_ec,' ','\ ') ' 1 - spline']);
%FSLcommand2= '/usr/share/fsl/5.0/bin/';
%system([FSLcommand2 'eddy --imain=' sinFile ' --mask=' strrep(outdir,' ','\ ') 'mask2_mask.nii --index=' strrep(outdir,' ','\ ') 'index.txt --acqp=' strrep(outdir,' ','\ ') 'acqparams.txt --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0  --slm=linear --topup=' strrep(outdir,' ','\ ') 'topup_results --flm=quadratic --out=' strrep(outdir,' ','\ ') 'DTI_esc.nii']);
%status=system([FSLcommand2 'eddy --imain=' sinFile ' --mask=' [outdir 'b0_ss_mask.nii'] ' --index=' [outdir 'index.txt'] ' --acqp=' [outdir 'acqparams.txt'] ' --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0 --resamp=jac --interp=spline --dont_peas --fep --repol --flm=quadratic --slm=linear --topup=' [outdir 'topup_results'] ' --out=' [outdir 'DTI_esc.nii']]);
disp('---------Eddy correction----------');
%status=system([opt.FSL_path 'eddy_openmp --imain=' files_in.In1{1} ' --mask=b0_ss_mask.nii --bvecs=' strrep(files_in.In1{1},'.nii', '.bvec') ' --bvals=' strrep(files_in.In1{1},'.nii', '.bval') ' --out=DTI_esc --flm=quadratic --slm=linear --fwhm=0 --fep --interp=spline --resamp=jac --dont_peas']);
status=system([opt.FSL_path 'eddy_openmp --imain=' name_split{end} ' --mask=b0_ss_mask.nii --index=index.txt --acqp=acqparams.txt --bvecs=' strrep(name_split{end},'.nii', '.bvec') ' --bvals=' strrep(name_split{end},'.nii', '.bval') ' --out=DTI_esc --flm=quadratic --slm=linear --fwhm=0 --fep --interp=spline --resamp=jac --dont_peas']);
if status ==0
    disp('DONE');
else
    warning('FAIL');
end

%% DTIFIT
disp('--------DTIFIT: Extraction de cartes MD et FA---------');
system([opt.FSL_path 'dtifit --data=DTI_esc.nii --out=dti --mask=b0_ss_mask.nii --bvecs='  strrep(name_split{end},'.nii', '.bvec') ' --bvals=' strrep(name_split{end},'.nii', '.bval') ' --save_tensor']);
if status == 0
    disp('DONE');
else
    warning('FAIL');
end

%% Save file and create Jsonfile
for i=1:length(opt.NameOutFiles)

    copyfile(fullfile(current_folder, [opt.NameOutFiles{i} '.nii']), files_out.In1{i});
    
    Ji = KeepModuleHistory(J1, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

    jsonfile = strrep(files_out.In1{i},'.nii','.json');
    WriteJson(Ji, jsonfile)
end

%rmdir(current_folder, 's');


