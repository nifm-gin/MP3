function [files_in,files_out,opt] = Module_Apply_clustering_GMM(files_in,files_out,opt)

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
    module_option(:,3)   = {'Clustering_model',''};
    module_option(:,4)   = {'Output_cluster_Name','Name_of_the_cluster_generated'};
    module_option(:,5)   = {'Anomaly_exclusion', 'No'};
    module_option(:,6)   = {'ROI_abnormal_voxels', 'Abnormal_voxels'};
    module_option(:,7)   = {'Classification_threshold_method', 'Normal/Abnormal'};
    module_option(:,8)   = {'NbClusters',5};
    module_option(:,9)   = {'Number_of_replicate',10};
    module_option(:,10)   = {'ROI_abnormal_voxels_Yes_No','Yes'};
    module_option(:,11)   = {'Threashold_value',Inf};
    module_option(:,12)   = {'RefInput',2};
    module_option(:,13)   = {'InputToReshape',2};
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
        'This module can apply a clustering model on seleced data'
        }'};
    
    user_parameter(:,2)   = {'Input Scans','XScan','','', {'SequenceName'},'Mandatory',...
        'The scans on which the clustering will be performed.'};
    user_parameter(:,3)   = {'ROI','1ROI','','',{'SequenceName'},'Mandatory',...
        'This ROI will select the pixels to apply the module to.'};
    user_parameter(:,4)   = {'Parameters','','','','','',''};
    handles_MP3 = guidata(findobj('Tag', 'MP3_GUI'));
    Cluster_listing    = dir(strcat(handles_MP3.database.Properties.UserData.MP3_ROI_path, '*.mat*'));
    if isempty(Cluster_listing)
        Cluster_listing(1).name = ' ';
    end
    user_parameter(:,5)   = {'   .Select a clustering model','cell', {Cluster_listing.name}, 'Clustering_model','','',...
        {'Please a clustering model amoung the ones already generated'}};
    user_parameter(:,6)   = {'   .Exclude abnormal voxels?','cell',{'No', 'Yes'},'Anomaly_exclusion','','',...
        ''};
    user_parameter(:,7)  = {'        .if Yes --> Would you like to save abnormal voxel as ROI?','cell',{'Yes', 'No'},'ROI_abnormal_voxels_Yes_No','','',...
        ''};
    user_parameter(:,8)  = {'             .if Yes --> Name of the ROI','char', 'Abnormal_voxels','ROI_abnormal_voxels','','',...
        ''};
    user_parameter(:,9)  = {'             .Exclusion method?','cell',{'Normal/Abnormal', 'SlopeHeuristic', '95%', 'Threashold value'},'Classification_threshold_method','','',...
        ''};
    user_parameter(:,10)  = {'                  .if SlopeHeuristic --> Number of replicate per iteration?','numeric',10,'Number_of_replicate','','',...
        'Please select the number of time you would like to replicate the clustering per iteration; For more inforamtion cf. fitgmdist function'};
    user_parameter(:,11)  = {'                  .if SlopeHeuristic --> Number of clusters?','numeric',5,'NbClusters','','',...
        'Number of clusters in which will be sorted the data. If the slope heuristic parameter is set to ''yes'', this number will represent the maximum number of clusters that will be tested by the algorithm.'};
    user_parameter(:,12)  = {'                  .if Threashold value --> which threashold?','numeric',Inf,'Threashold_value','','',...
        'Number of clusters in which will be sorted the data. If the slope heuristic parameter is set to ''yes'', this number will represent the maximum number of clusters that will be tested by the algorithm.'};
    
    user_parameter(:,13)   = {'   .Name of the resulting cluster','char','','Output_cluster_Name','','',...
        'This module will create one cluster type of file for each input scan. '};
    
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

if isempty(files_out)
    opt.Table_out = opt.Table_in(1,:);
    opt.Table_out.IsRaw = categorical(0);
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    opt.Table_out.SequenceName = categorical(cellstr(opt.Output_cluster_Name));
    opt.Table_out.Type = categorical(cellstr('Cluster'));
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
        opt.Table_out(2,:) = opt.Table_out(1,:);
        opt.Table_out.Type(2) = categorical(cellstr('ROI'));
        opt.Table_out.SequenceName(2) = categorical(cellstr(opt.ROI_abnormal_voxels));
        opt.Table_out.Filename(2) = categorical(cellstr([char(opt.Table_out.Patient(2)), '_', char(opt.Table_out.Tp(2)), '_', char(opt.Table_out.SequenceName(2))]));
        f_out = [char(opt.Table_out.Path(2)), char(opt.Table_out.Patient(2)), '_', char(opt.Table_out.Tp(2)), '_', char(opt.Table_out.SequenceName(2)), '.nii'];
        files_out.In2{1} = f_out;
        
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

Tag1 = 'Patient';
Tag2 = 'Tp';

databScans = opt.Table_in(opt.Table_in.Type == categorical(cellstr('Scan')),:);
databROIs = opt.Table_in(opt.Table_in.Type == categorical(cellstr('ROI')),:);
UTag1 = unique(databScans.(Tag1));
UTag2 = unique(databScans.(Tag2));
in_files = {};
roi_files = {};
for i=1:length(UTag1)
    for j=1:length(UTag2)
        DbRois = databROIs(strcmp(cellstr(databROIs.(Tag1)), cellstr(UTag1(i))),:);
       % DbRois = DbRois(DbRois.(Tag2) == UTag2(j),:);
        DbRois = DbRois(strcmp(cellstr(DbRois.(Tag2)), cellstr(UTag2(j))),:);

        if size(DbRois, 1) == 0
            continue
        end
        roi = [char(DbRois.Path(1)), char(DbRois.Filename(1)), '.nii'];
        datab = databScans(strcmp(cellstr(databScans.(Tag1)), cellstr(UTag1(i))),:);
        datab = datab(strcmp(cellstr(datab.(Tag2)), cellstr(UTag2(j))),:);
        fi = cell(size(datab,1),1);
        for k=1:size(datab,1)
            fi{k} = [char(datab.Path(k)), char(datab.Filename(k)), '.nii'];
        end
        in_files = [in_files ; {fi}];
        roi_files = [roi_files ; {roi}];
    end
end









% ROI_nifti_header = spm_vol(files_in.In2{1});
% ROI = read_volume(ROI_nifti_header, ROI_nifti_header, 0, 'axial');
% NbVox = int64(sum(sum(sum(ROI))));
% Data = zeros(NbVox, length(files_in.In1));
All_Data = {};
ROI_nifti_header = cell(length(roi_files),1);
ROI = cell(length(roi_files),1);
VecVoxToDelete = {};



for i=1:length(roi_files)
    roi_filename = roi_files{i};
    [~, fname] = fileparts(roi_filename);
    Files = in_files{i};
    Name_ROI = opt.Table_in.SequenceName(opt.Table_in.Filename == categorical(cellstr(fname)));
    ROI_nifti_header{i} = spm_vol(roi_filename);
    ROI{i} = read_volume(ROI_nifti_header{i}, ROI_nifti_header{i}, 0, 'axial');
    NbVox = int64(sum(sum(sum(ROI{i}))));
    %Data = zeros(NbVox, length(Files));
    %Data = zeros(NbVox, 1);
    Data = [];
    NameScans = {};
    for j=1:length(Files)
        [~, sname] = fileparts(Files{j});
        Name_Scan = opt.Table_in.SequenceName(opt.Table_in.Filename == categorical(cellstr(sname)));
        Name_Patient = opt.Table_in.Patient(opt.Table_in.Filename == categorical(cellstr(sname)));
        Name_TP = opt.Table_in.Tp(opt.Table_in.Filename == categorical(cellstr(sname)));
        Name_Group = opt.Table_in.Group(opt.Table_in.Filename == categorical(cellstr(sname)));
        nifti_header = spm_vol(Files{j});
        input{j} = read_volume(nifti_header, ROI_nifti_header{i}, 0, 'axial');
        % the ROI needs to have the same class as the data
        ROI_NaN = cast(ROI{i}, class(input{j}));
        ROI_NaN(ROI_NaN == 0) = NaN;
        input{j} = input{j}.*ROI_NaN;
        %% merge the 4th and the 5th dimension (mean data)
        %         Vec = mean(input{j},4);
        %         NameScans = [NameScans, {char(Name_Scan)}];
        %         Vec(isnan(Vec)) = [];
        %         Data(:,j) = Vec.';
        %% create one row for each 4th and 5th dimension
        for x = 1:size(input{j},4)
            for xx= 1:size(input{j},5)
                if x==1 && xx==1
                    suffix = '';
                elseif size(input{j},5) == 1
                    suffix = ['_Ech', num2str(x)];
                else
                    suffix = ['_Ech', num2str(x), '_Rep', num2str(xx)];
                end
                NameScans = [NameScans, {[char(Name_Scan), suffix]}];
                Vec = input{j}(:,:,:,x,xx);
                %Vec(isnan(Vec)) = [];
                Data(:,size(Data,2)+1) = Vec(:);
                
            end
        end
    end
    VecVoxToDelete{i} = sum(isnan(Data),2)~=0;
    Data(VecVoxToDelete{i},:) = []; % Exclude any voxel that have at least one NaN parameter.
    Names_Patients = repmat({char(Name_Patient)},size(Data,1),1);
    Names_TPs = repmat({char(Name_TP)},size(Data,1),1);
    Names_Groups = repmat({char(Name_Group)},size(Data,1),1);
    NameScans = [{'Group', 'Patient', 'Tp'}, NameScans];
    %Data(:,1) = [];
    if size(Data, 2) ~= length(Files)
        continue
    end
    Data = [Names_Groups, Names_Patients, Names_TPs, num2cell(Data)];
    
    
    Data = cell2table(Data, 'VariableNames', NameScans);
    
    All_Data = [All_Data, {Data}];
end

Size = [];
for i=1:length(All_Data)
    Size = [Size, size(All_Data{i},2)];
end
SizeMax = max(Size);
%SizeMaxVec = Size == SizeMax;


All_Data_Clean = {};
ROI_Clean = {};
VecVoxToDeleteClean = {};
Clust_Data_In = [];
ROI_nifti_header_Clean = {};
for i=1:length(All_Data)
    if size(All_Data{i},2) ~= SizeMax
        continue
    end
    All_Data_Clean = [All_Data_Clean, All_Data(i)];
    ROI_Clean = [ROI_Clean, ROI(i)];
    VecVoxToDeleteClean = [VecVoxToDeleteClean, VecVoxToDelete{i}];
    Clust_Data_In = [Clust_Data_In ; All_Data{i}];
    ROI_nifti_header_Clean = [ROI_nifti_header_Clean, ROI_nifti_header{i}];
end





%% apply trained model
trainedModel_loaded = load(strcat(strrep(opt.folder_out, 'Tmp', 'ROI_data'), filesep, opt.Clustering_model));

VariableNames = Clust_Data_In.Properties.VariableNames;
exit = 0;
if isfield(trainedModel_loaded, 'Informations')
    % case of GMM
    trainedModel = trainedModel_loaded.Informations.Modele;
    PredictorNames = trainedModel_loaded.Informations.Cartes;
    
    if length(PredictorNames) ~= length(VariableNames) - 3 || ...
            sum(strcmp(sort(PredictorNames), sort(VariableNames(4:end)))) ~= length(PredictorNames)
        exit = 1;
    else
        data = zeros([size(Clust_Data_In,1) length(PredictorNames)]);
        for i = 1:length(PredictorNames)
            data(:,i) = Clust_Data_In.(PredictorNames{i}); % 1='Group'; 2 = 'Patient' and 3 = 'Tp'
        end
        % Normalize data if needed. We used the mean +/- SD value used to generate the model
        if isfield(trainedModel_loaded.Informations, 'Normalization_mode') && ...
                strcmp(trainedModel_loaded.Informations.Normalization_mode,'All Database')
            data = (data-trainedModel_loaded.Informations.NanMean_VoxValues)./trainedModel_loaded.Informations.NanStd_VoxValues;
        end
        if isfield(trainedModel_loaded.Informations, 'Normalization_mode') && ...
                strcmp(trainedModel_loaded.Informations.Normalization_mode,'Patient-by-Patient')
            data = (data-nanmean(data))./nanstd(data);
        end
        
        %[IDX,NLOGL,POST,LOGPDF,MAHALAD] = CLUSTER(trainedModel,data)
        [ClusteredVox,~,~,LOGPDF,dist]= cluster(trainedModel, data);
        if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
            ClusteredVox_excluded = zeros(size(ClusteredVox));
        end
        
        %% test plot LOGPDF
        %         Img_LOGPDF = ROI_Clean{1};
        %         Img_LOGPDF(~VecVoxToDeleteClean{i}) = LOGPDF(1:1+size(All_Data_Clean{1},1)-1);
        %         figure;imshow3D(Img_LOGPDF);
        
        if strcmp(opt.Anomaly_exclusion, 'Yes')
            if strcmp(opt.Classification_threshold_method, 'Threashold value')
%                 LOGPDF_cutoff = opt.Threashold_value;
%                 ClusteredVox(LOGPDF<LOGPDF_cutoff) = 0;
%                 if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
%                     ClusteredVox_excluded(LOGPDF<LOGPDF_cutoff) = 1;
%                 end
%                 disp(['cutoff_abnormality : ', num2str(LOGPDF_cutoff)])
                 %% test with the distance
                for i=1:size(dist,1)
                    dist_from_best_cluster(i) = dist(i,ClusteredVox(i));
                end
                dist_cutoff = opt.Threashold_value;
                ClusteredVox(dist_from_best_cluster>dist_cutoff) = 0;
                if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
                    ClusteredVox_excluded(dist_from_best_cluster>dist_cutoff) = 1;
                end
                disp(['cutoff_distance : ', num2str(dist_cutoff)])
            elseif strcmp(opt.Classification_threshold_method, '95%')
                %                 LOGPDF_cutoff = prctile(LOGPDF,5);
                %                 ClusteredVox(LOGPDF<LOGPDF_cutoff) = 0;
                %                 if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
                %                     ClusteredVox_excluded(LOGPDF<LOGPDF_cutoff) = 1;
                %                 end
                %                 disp(['cutoff_abnormality : ', num2str(LOGPDF_cutoff)])
                
                %% test with the distance
                for i=1:size(dist,1)
                    dist_from_best_cluster(i) = dist(i,ClusteredVox(i));
                end
                dist_cutoff = prctile(dist_from_best_cluster,95);
                ClusteredVox(dist_from_best_cluster>dist_cutoff) = 0;
                if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
                    ClusteredVox_excluded(dist_from_best_cluster>dist_cutoff) = 1;
                end
                disp(['cutoff_distance : ', num2str(dist_cutoff)])
                
            elseif strcmp(opt.Classification_threshold_method, 'Normal/Abnormal')
                options = statset ( 'maxiter', 1000);
                
                %% use the LOGPDF to split the abnormality in 2 clusters
                % split abdnomaly in 2
                modele_LOGPDF_k2 = fitgmdist(LOGPDF, 2, 'Options', options, 'Regularize', 1e-5, 'Replicates',  opt.Number_of_replicate);
                abnormality_separation_k2 =  cluster(modele_LOGPDF_k2, LOGPDF);
                
                min_of_max_abnormality = max(LOGPDF(abnormality_separation_k2 == find(modele_LOGPDF_k2.mu == min(modele_LOGPDF_k2.mu))));
                max_of_min_abnormality = min(LOGPDF(abnormality_separation_k2 == find(modele_LOGPDF_k2.mu == max(modele_LOGPDF_k2.mu))));
                LOGPDF_cutoff = (max_of_min_abnormality + min_of_max_abnormality) /2;
                disp(['cutoff_abnormality : ', num2str(LOGPDF_cutoff)])
                %
                % find the cluster with the lowest log-score. This cluster
                % corresponds to the less adequacy with the reference model
                % --> abdormal voxels
                ClusteredVox(abnormality_separation_k2 == find(modele_LOGPDF_k2.mu == min(modele_LOGPDF_k2.mu))) = 0;
                if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
                    ClusteredVox_excluded(abnormality_separation_k2 == find(modele_LOGPDF_k2.mu == min(modele_LOGPDF_k2.mu))) = 1;
                end
                
                
            elseif  strcmp(opt.Classification_threshold_method, 'SlopeHeuristic')
                options = statset ( 'maxiter', 1000);
                
                %% use the LOGPDF to split the abnormality in k classes (using the Heuristic slop)
                % then use the LOGPDF to spit the abnormality in 2 classes
                % finally the seperatation obained form the 2 clases separations to split the k classes in 2
                % (cf Alexis Arnaud et al. IEEE TRANSACTIONS ON MEDICAL IMAGING, VOL. 37, NO. 7, JULY 2018
                %% code to find tune the number of abnormality classes
                
                ptsheurist = opt.NbClusters + 5;
                
                
                %Vecteur pour stocker la logvraisemblance
                loglike = zeros(1,ptsheurist);
                
                %On stocke les modeles calcules pour ne pas avoir a les recalculer une
                %fois le nombre de classes optimal trouve.
                %modeles = cell(1,ptsheurist);
                Number_of_replicate = opt.Number_of_replicate;
                % find the number of abnormality classes
                parfor kk=1:ptsheurist
                    %La ligne suivante permet uniquement de suivre l'avancement du
                    %calcul des modeles
                    disp(strcat('Modele_', num2str(kk), '_started'))
                    %L'option "Replicate,10" signifie que l'on va calculer 10 fois le
                    %modele en modifiant l'initialisation. Le modele renvoye est celui
                    %de plus grande vraisemblance.
                    modeles_LOGPDF_kn{kk} = fitgmdist(LOGPDF, kk, 'Options', options, 'Regularize', 1e-5, 'Replicates', Number_of_replicate);
                    
                    loglike(kk) = -modeles_LOGPDF_kn{kk}.NegativeLogLikelihood;
                    
                    %La ligne suivante permet uniquement de suivre l'avancement du
                    %calcul des modeles
                    disp(strcat('Modele_', num2str(kk), '_done'))
                end
                
                %
                %Le vecteur alpha contient les coefficients directeurs des regressions
                %lineaires de la logvraisemblance en fonction du nombre de classes du
                %modele. Sa ieme composante contient le coefficient directeur de la
                %regression lineaire de la log vraisemblance en fonction du nombre de
                %classes du modele en ne prenant pas en compte les i-1 premiers points.
                alpha = zeros(opt.NbClusters,2);
                
                %Le vecteur eqbic contient pour chaque valeur alpha l'equivalent BIC
                %applique a chaque valeur de la log vraisemblance. On obtient donc une
                %matrice ou chaque ligne correspond a l'equivalent BIC applique en
                %chaque valeur de la log vraisemblance pour une valeur de alpha. On
                %passe ainsi d'une ligne a l'autre en modifiant alpha. Dans l'optique
                %de tracer les courbes uniquement a partir du point i, la matrice est initialisee a la valeur NaN.
                eqbic = NaN(opt.NbClusters,length(loglike));
                
                %Le vecteur eqbic2 est similaire au vecteur eqbic mais avec un autre
                %critere.
                %eqbic2 = NaN(opt.NbClusters,length(loglike));
                
                for j = 1:opt.NbClusters
                    %La regression lineaire
                    alpha(j,:) = polyfit(j:ptsheurist,loglike(j:end),1);
                    for i=j:length(loglike)
                        %eqbic2(j,i) = 2*alpha(j,1)*(i-1+NbCartes*i+(1+NbCartes)*NbCartes/2*i)-loglike(i);
                        eqbic(j,i) = 2*alpha(j,1)*i-loglike(i);
                    end
                end
                
                figure
                plot(eqbic.')
                [~,I] = nanmin(eqbic,2);  %Pour chacune des courbes de l'eqbic, l'indice pour lequel le minimum est atteint est considere comme etant le nombre optimal de clusters
                figure
                plot(0:opt.NbClusters,0:opt.NbClusters,'r')
                hold on
                plot(I,'b')
                k = 0;
                % On vient de tracer le nombre de clusters optimal pour chaque courbe,
                % donc pour chaque coefficient directeur, donc pour chaque point i
                % debut de la regression lineaire. Le nombre de classes optimal global
                % est le point k minimum pour lequel f(k) = k, soit l'intersection de
                % la courbe tracee et de le bissectrice du plan.
                for i=1:length(I)
                    if I(i) == i && k == 0
                        k = i;
                    end
                end
                
                if k == 0
                    warndlg('Cannot find an optimal number of cluster, try again and test higher numbers of clusters','Cannot find a number of clusters; k has been set to the default number of clusters asked');
                    %error('Cannot find an optimal number of cluster, try again and test higher numbers of clusters');
                    k = opt.NbClusters;
                    gmfit_LOGPDF_kn = modeles_LOGPDF_kn{k};
                    %return
                else
                    gmfit_LOGPDF_kn = modeles_LOGPDF_kn{k};
                end
                
                % define the PDF cutoff using the clustering with 2 classes
                abnormality_separation_k2 =  cluster(modeles_LOGPDF_kn{2}, LOGPDF);
                min_of_max_abnormality = max(LOGPDF(abnormality_separation_k2 == find(modeles_LOGPDF_kn{2}.mu == min(modeles_LOGPDF_kn{2}.mu))));
                max_of_min_abnormality = min(LOGPDF(abnormality_separation_k2 == find(modeles_LOGPDF_kn{2}.mu == max(modeles_LOGPDF_kn{2}.mu))));
                cutoff_abnormality = (max_of_min_abnormality + min_of_max_abnormality) /2;
                
                ClusteredVox_LOGPDF_kn = cluster(gmfit_LOGPDF_kn, LOGPDF);
                % find the closed mu_PDF to the cutoff_abnormality
                [~, threashold_cluster] = min(abs(gmfit_LOGPDF_kn.mu-cutoff_abnormality));
                mu_LOGPDS_threashold = gmfit_LOGPDF_kn.mu(threashold_cluster);
                % every abnormal cluster below the cutoff_cluster is consider
                % abnormal --> set to 0
                abnormal_clusters = find(gmfit_LOGPDF_kn.mu <= mu_LOGPDS_threashold);
                for z = 1:numel(abnormal_clusters)
                    ClusteredVox(ClusteredVox_LOGPDF_kn == abnormal_clusters(z)) = 0;
                    if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
                        ClusteredVox_excluded(ClusteredVox_LOGPDF_kn == abnormal_clusters(z)) = 1;
                    end
                end
            end
        end
        Clust_Data_In.Cluster = ClusteredVox;
    end
    
else
    % case of homemade models
    trainedModel_name = fieldnames(trainedModel_loaded);
    trainedModel = trainedModel_loaded.(trainedModel_name{:});
    data = table;
    PredictorNames = trainedModel.ClassificationTree.PredictorNames;
    if length(PredictorNames) ~= length(VariableNames) - 3
        exit = 1;
    else
        for i = 1:length(PredictorNames)
            data.(PredictorNames{i}) = Clust_Data_In.(VariableNames{3+i}); % 1='Group'; 2 = 'Patient' and 3 = 'Tp'
        end
    end
    
    [ClusteredVox, score] = trainedModel.predictFcn(data);
    ClusteredVox(max(score, [], 2) < (opt.Classification_threshold/100)) =  max(ClusteredVox)+1; %0;
    
    Clust_Data_In.Cluster = ClusteredVox;
end

if exit == 1
    errordlg('The scans selected does not correspond to the scans used to built this model')
    return
end
[MoyCartesVolume, ProbVolume, Ecart_Type_Global, Sign, MoyGlobal] = AnalyseClusterGMM(Clust_Data_In);


%On cree 2 strutures que l'on va sauvegarder avec chaque uvascroi. Ces
%structures contiennent les informations et les statistiques du clustering.
%Informations = struct('Cartes', {NameScans(4:end)} , 'Modele', gmfit, 'Sign', Sign, 'ROI', char(Name_ROI));
Informations = struct('Cartes', {NameScans(4:end)} , 'Modele', trainedModel, 'Sign', Sign, 'ROI', char(Name_ROI));

Statistiques = struct('MoyCartesVolume', MoyCartesVolume , 'ProbVolume', ProbVolume, 'Ecart_Type_Global', Ecart_Type_Global,'MoyGlobal', MoyGlobal);
% IA_patient_list = {handles.MIA_data.database.name};
% NomDossier = [];
% for i = 1:length(Informations.Cartes)
%     NomDossier = [NomDossier '_' char(Informations.Cartes(i))];
% end
% NomDossier2 = [num2str(length(Informations.Cartes)) 'Cartes' filesep NomDossier];
% logbook = {};
% answer = inputdlg('Comment voulez vous nommer ce clustering ?', 'Choix du nom du clustering', 1,{strcat(num2str(k),'C_',NomDossier)});
% if ~exist(strcat(handles.MIA_data.database(1).path,NomDossier2), 'dir')
%     mkdir(strcat(handles.MIA_data.database(1).path,NomDossier2));
% end

for i=1:length(files_out.In1)
    save(strrep(files_out.In1{i}, '.nii', '.mat'),'Informations', 'Statistiques');
end





ind = 1;
for i=1:length(All_Data_Clean)
    Cluster = ROI_Clean{i};
    %ROI_Clust = logical(ROI{i});
    Cluster(~VecVoxToDeleteClean{i}) = ClusteredVox(ind:ind+size(All_Data_Clean{i},1)-1);
    ROI_cluster_header = ROI_nifti_header_Clean{i};
    % On a fait le même traitement sur les files_out (tout début du code) que sur l'ouverture des ROI. Il y a donc tout à penser que l'ordre des fichiers correspondra.
    ROI_cluster_header.fname = files_out.In1{i};
    ROI_cluster_header = rmfield(ROI_cluster_header, 'pinfo');
    ROI_cluster_header = rmfield(ROI_cluster_header, 'private');
    
    Cluster = write_volume(Cluster,ROI_nifti_header_Clean{i}, 'axial');
    spm_write_vol(ROI_cluster_header, Cluster);
    
    if strcmp(opt.ROI_abnormal_voxels_Yes_No, 'Yes')
        ROI_abnormal_voxels = ROI_Clean{i};
        ROI_abnormal_voxels(~VecVoxToDeleteClean{i}) = ClusteredVox_excluded(ind:ind+size(All_Data_Clean{i},1)-1);
        % save the mask
        ROI_abnormal_voxels = single(ROI_abnormal_voxels);
        % copy nifti_header from the Cluster
        ROI_abnormal_voxels_header = ROI_cluster_header;
        % update the header with the new info
        ROI_abnormal_voxels_header.fname = files_out.In2{1};
        ROI_abnormal_voxels_header.Filemoddate = char(datetime('now'));
        ROI_abnormal_voxels_header.Datatype = class(ROI_abnormal_voxels);
        %reorient the volume
        ROI_abnormal_voxels = write_volume(ROI_abnormal_voxels,ROI_nifti_header_Clean{i}, 'axial');
        
        spm_write_vol(ROI_abnormal_voxels_header, ROI_abnormal_voxels);
        
        %% test saving LOGPDF
        %         LOGPDF_map = ROI_Clean{i};
        %         LOGPDF_map(~VecVoxToDeleteClean{i}) = LOGPDF(ind:ind+size(All_Data_Clean{i},1)-1);
        %         % save the mask
        %         LOGPDF_map = single(LOGPDF_map);
        %         % copy nifti_header from the Cluster
        %         ROI_abnormal_voxels_header = ROI_cluster_header;
        %         % update the header with the new info
        %         ROI_abnormal_voxels_header.fname = files_out.In2{1};
        %         ROI_abnormal_voxels_header.Filemoddate = char(datetime('now'));
        %         ROI_abnormal_voxels_header.Datatype = class(LOGPDF_map);
        %         %reorient the volume
        %         LOGPDF_map = write_volume(LOGPDF_map,ROI_nifti_header_Clean{i}, 'axial');
        %
        %         spm_write_vol(ROI_abnormal_voxels_header, LOGPDF_map);
    end
    
    ind = ind+size(All_Data_Clean{i},1);
    
end

