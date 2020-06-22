function [files_in,files_out,opt] = Module_DTI_2sens(files_in,files_out,opt)

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
    user_parameter(:,2)   = {'Select the DTI APP scan','1Scan','','',{'SequenceName'}, 'Mandatory','Please copy the .bvals and .bvecs files associated to this file into the Raw_data folder of this project'};
    user_parameter(:,3)   = {'Select the DTI APA scan','1Scan','','',{'SequenceName'}, 'Mandatory','Please copy the .bvals and .bvecs files associated to this file into the Raw_data folder of this project'};
    user_parameter(:,4)   = {'Indicate path to FSL','char','','FSL_path','','Mandatory','Result of <which fsl> on any terminal'};
    user_parameter(:,5)   = {'Parameters','','','','', '', ''};
    user_parameter(:,6)   = {'   .Output filename extension FA','char', 'dti_FA','output_filename_ext_FA','','',''};
    user_parameter(:,7)   = {'   .Output filename extension MD','char','dti_MD','output_filename_ext_MD','','',''};
    user_parameter(:,8)   = {'   .Output filename extension L1','char','dti_L1','output_filename_ext_L1','','',''};
    user_parameter(:,9)   = {'   .Output filename extension L2','char','dti_L2','output_filename_ext_L2','','',''};
    user_parameter(:,10)  = {'   .Output filename extension L3','char','dti_L3','output_filename_ext_L3','','',''};
    
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
input(2).nifti_header = spm_vol(files_in.In2{1});

J1 = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
J2 = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));

inFile_APP = read_volume(input(1).nifti_header, input(1).nifti_header, 0);
inFile_APA = read_volume(input(2).nifti_header, input(2).nifti_header, 0);

%nombval = char(table_out_tmp.Path)

setenv('FSLOUTPUTTYPE','NIFTI');

%% fslroi: extract b0_P and b0_A files
nombval=strrep(files_in.In1{1},'.nii', '.bval');
nombvec=strrep(files_in.In1{1},'.nii', '.bvec');

nombval2 = strrep(files_in.In2{1},'.nii', '.bval');
% bvaltest= importdata([indir nombval]);
% posb0= find(bvaltest==0);
% clear bvaltest
%system([FSLcommand 'fslroi ' strrep(indir,' ','\ ') inFile ' ' strrep(indir,' ','\ ') 'b0_P.nii '  num2str(posb0(1)-1) ' 1']);
%system([FSLcommand 'fslroi ' indir inFile ' ' outdir 'b0_P.nii 0 1']);

cd(opt.folder_out)
current_folder = fullfile(opt.folder_out, char(opt.Table_in.Patient(1)));
mkdir(current_folder)
cd(current_folder)

system([opt.FSL_path 'fslroi ' files_in.In1{1} ' ' current_folder '/b0_P.nii 0 1']);

% bvaltest= importdata([indir nombval2]);
% posb0A= find(bvaltest==0);
% clear bvaltest
%system([FSLcommand 'fslroi ' strrep(PathSujetS,' ','\ ') '/' NomDossierB0 '/' B0file ' ' strrep(PathSujetS,' ','\ ') '/' NomDossierB0 '/b0_A.nii ' num2str(posb0A(1)-1) ' 1']);
%system([FSLcommand 'fslroi ' indir APAfile ' ' outdir 'b0_A.nii 0 1']);
system([opt.FSL_path 'fslroi ' files_in.In2{1} ' ' current_folder '/b0_A.nii 0 1']);
 
clear matlabbatch
spm_jobman('initcfg');

%% Coregistration du b0_A sur le b0_P
% Pour assurer une coregistration parfaite pour l'estimation de la carte de champs par topup
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[current_folder '/b0_P.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[current_folder '/b0_A.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);
clear matlabbatch

spm_check_registration([current_folder '/b0_P.nii'],[current_folder '/rb0_A.nii'])
disp('Vérification de la coregistration, pressez une touche pour continuer');
        
%% fslmerge:merge b0_P and rb0_A into b0_PA
system([opt.FSL_path 'fslmerge -t ' current_folder '/b0_PA.nii ' current_folder '/b0_P.nii ' current_folder '/rb0_A.nii']); %merge les 2 b0 pour topup

%% Delete supplementary Philips volume
bvecf = importdata(nombvec,' ');
bvalf = importdata(nombval,' ');
name_split = strsplit(files_in.In1{:},filesep);

copyfile(nombval,fullfile(current_folder,strrep(name_split{end},'.nii', '.bval')));
copyfile(nombvec,fullfile(current_folder,strrep(name_split{end},'.nii', '.bvec')));

Vol = 1:size(bvecf,2);
idx = find(bvecf(1,:)==0 & bvecf(2,:)==0 & bvecf(3,:)==0 & bvalf(1,:)~=0);

if ~isempty(idx)
    Tempbval=fullfile(current_folder,strrep(files_in.In1{1},'.nii', '.bval'));
    Tempbvec=fullfile(current_folder,strrep(files_in.In1{1},'.nii', '.bvec'));
%     
    bvecf(:,idx(1))=[]; % élimine le volume supplémentaire
    Vol(idx(1))=[];
    
    bvalf(idx(1))=[];
    
    disp('Extraction du volume supplémentaire...');
    %     system([FSLcommand 'fslroi ' Indir NomFichierNii ' ' Outdir NomFichierNiiX ' 0 ' num2str(idx(2)-1)]); %Exclue le volume supplémentaire des Nii
    
    for j= 1 : size(Vol,2)
        matlabbatch{1}.spm.util.cat.vols{j,1} = [files_in.In1{1}  ',' num2str(Vol(j)) ];
    end
    
    matlabbatch{1}.spm.util.cat.dtype = 4;
    matlabbatch{1}.spm.util.cat.name = fullfile(opt.folder_out,name_split{end});
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
    
    
    disp('Extraction du volume supplémentaire OK');
    fmt= '%d'; % Définit le format pour l'écriture dans le fichier bvec et bval
    for j= 1 :(size(bvecf,2)-1)
        fmt=[fmt ' %6.4E'];
    end
    fmt=[fmt '\n'];
    
    fid=fopen(Tempbvec,'wt');
    fprintf(fid,fmt,bvecf');
    fclose(fid);
    
    fid2=fopen(Tempbval,'wt');
    fprintf(fid2,'%d ',bvalf');
    fclose(fid2);
    clear fid fid2
    
    disp('Modification des bvec et bval OK');
else
    copyfile(files_in.In1{:},fullfile(current_folder,name_split{end}));
end

%--------------------------------------------------
fid=fopen(fullfile(current_folder,'/index.txt'),'w');
ind=ones(1,size(bvecf,2));
fprintf(fid,'%d ',ind);
fclose(fid);
clear ind fid

%% Topup
if ~exist([current_folder '/acqparams.txt'],'file')
    system(['printf "0 -1 0 0.03327\n0 1 0 0.03327" > ' current_folder '/acqparams.txt']);
end

if exist([current_folder '/b0_PA.nii'], 'file')
    disp('Estimation de la carte de champs ...')
    v1 = spm_vol([current_folder '/b0_P.nii']);
    if mod(v1.dim(3),2)==0
        status=system([opt.FSL_path 'topup --imain=b0_PA.nii --datain=acqparams.txt --config=b02b0.cnf --out=topup_results']);
%     else
%         status=system([opt.FSL_path 'topup --imain=b0_PA.nii --datain=acqparams.txt --config=/home/veronica/Donnees/Patients/Park_DTI/b02b02.cnf --out=topup_results']);% b02b02.cnf sans subsampling dans le répertoire contenant tous les centres
    end
    if status==0
        disp('Estimation de la carte de champs  OK')
    else
        warning('L''estimation de la carte de champs par topup a échoué, Vérification requise (dbcont pour continuer)');
        dbstop in PreprocessingDTI_FSL.m
    end
    disp('Correction des distorsions sur les b0 des images de diffusions ...')
    status=system([opt.FSL_path 'applytopup --imain=b0_P.nii,rb0_A.nii --topup=topup_results --datain=acqparams.txt --inindex=1,2 --out=b0_corr.nii']); % Correction de la b0 pour servir de cible de coregistration pour les images anatomiques
    if status ==0
        disp('Correction des distorsions sur les b0 des images de diffusions  OK')
    else
        warning('La correction des distorsions sur la B0 a échoué, Vérification requise (dbcont pour continuer)');
    end
end
   
%% Skull extraction
disp('---------Skull extraction----------');
if exist([current_folder '/b0_corr.nii'],'file')
    status = system([opt.FSL_path 'bet b0_corr.nii b0_ss.nii -f 0.1 -m -Z']);
    if status ==0
        disp('DONE');
    else
        warning('FAIL');
    end
end

%% Eddy 
%     system([FSLcommand 'eddy_correct ' strrep(indir,' ','\ ') Nomdti ' ' strrep(out_ec,' ','\ ') ' 1 - spline']);
%FSLcommand2= '/usr/share/fsl/5.0/bin/';
%system([FSLcommand2 'eddy --imain=' sinFile ' --mask=' strrep(outdir,' ','\ ') 'mask2_mask.nii --index=' strrep(outdir,' ','\ ') 'index.txt --acqp=' strrep(outdir,' ','\ ') 'acqparams.txt --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0  --slm=linear --topup=' strrep(outdir,' ','\ ') 'topup_results --flm=quadratic --out=' strrep(outdir,' ','\ ') 'DTI_esc.nii']);
%status=system([FSLcommand2 'eddy --imain=' sinFile ' --mask=' [outdir 'b0_ss_mask.nii'] ' --index=' [outdir 'index.txt'] ' --acqp=' [outdir 'acqparams.txt'] ' --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0 --resamp=jac --interp=spline --dont_peas --fep --repol --flm=quadratic --slm=linear --topup=' [outdir 'topup_results'] ' --out=' [outdir 'DTI_esc.nii']]);
disp('---------Eddy correction----------');
status=system([opt.FSL_path 'eddy_openmp --imain=' name_split{end} ' --mask=b0_ss_mask.nii --index=index.txt --acqp=acqparams.txt --bvecs=' strrep(name_split{end},'.nii', '.bvec') ' --bvals=' strrep(name_split{end},'.nii', '.bval') ' --out=DTI_esc --topup=topup_results --flm=quadratic --slm=linear --fwhm=0 --fep --interp=spline --resamp=jac --dont_peas']);
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


