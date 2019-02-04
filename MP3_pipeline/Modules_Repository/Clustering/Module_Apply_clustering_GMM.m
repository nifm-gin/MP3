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
    module_option(:,3)   = {'Clustering_model_path','Please enter the path of the clustering model to apply'};
    module_option(:,4)   = {'output_cluster_Name','Name of the cluster generated'};
    module_option(:,5)   = {'Classification_threshold', 0};
    module_option(:,6)   = {'AutomaticJobsCreation', 'No'};
    module_option(:,7)   = {'RefInput',2};
    module_option(:,8)   = {'InputToReshape',2};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)   = {'Table_out', table()};
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
    user_parameter(:,5)   = {'   .Path of the clustering model','char','Please enter the path of the clustering model to apply','Clustering_model_path','','',...
        'Please enter the path of the clustering model to apply'};
    user_parameter(:,6)   = {'   .Threshold to apply to the classification','numeric',0,'Classification_threshold','','',...
        {'Threshold of the classification above witch the classifecation is validated'
        'For instance only voxel above 90% of correct classification is classify'}};
    user_parameter(:,7)   = {'   .Name of the resulting cluster','char','','output_cluster_Name','','',...
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

Tag1 = 'Patient';
Tag2 = 'Tp';
Table_out = table();
if strcmp(files_out, '')
    databScans = opt.Table_in(opt.Table_in.Type == categorical(cellstr('Scan')),:);
    databROIs = opt.Table_in(opt.Table_in.Type == categorical(cellstr('ROI')),:);
    UTag1 = unique(databScans.(Tag1));
    UTag2 = unique(databScans.(Tag2));
    out_file = {};
    in_files = {};
    Tailles = [];
    for i=1:length(UTag1)
        for j=1:length(UTag2)
            datab = databScans(databScans.(Tag1) == UTag1(i),:);
            datab = datab(datab.(Tag2) == UTag2(j),:);
            Tailles = [Tailles size(datab,1)];
        end
    end
    MaxTaille = max(Tailles);
    TailleBin = Tailles == MaxTaille;
    ind = 0;
    
    for i=1:length(UTag1)
        for j=1:length(UTag2)
            ind = ind+1;
            if TailleBin(ind)
                DbRois = databROIs(databROIs.(Tag1) == UTag1(i),:);
                DbRois = DbRois(DbRois.(Tag2) == UTag2(j),:);
                if size(DbRois, 1) == 0
                    continue
                end
                datab = databScans(databScans.(Tag1) == UTag1(i),:);
                datab = datab(datab.(Tag2) == UTag2(j),:);
                %                 fi = cell(size(datab,1),1);
                %                 for k=1:size(datab,1)
                %                     fi{k} = [char(datab.Path(k)), char(datab.Patient(k)), '_', char(datab.Tp(k)), '_', char(datab.SequenceName(k)), '.nii'];
                %                 end
                %                 in_files = [in_files ; fi];
                tags = databScans(1,:);
                tags.Patient = UTag1(i);
                tags.Tp = UTag2(j);
                tags.Type = categorical(cellstr('Cluster'));
                tags.IsRaw = categorical(1);
                Cluster_path = opt.folder_out; % strrep(opt.folder_out, 'Derived_data', 'ROI_data');
                tags.Path = categorical(cellstr([Cluster_path, filesep]));
                tags.SequenceName = categorical(cellstr([opt.output_cluster_Name]));
                tags.Filename = categorical(cellstr([char(tags.Patient), '_', char(tags.Tp), '_', char(tags.SequenceName)]));
                f_out = [char(tags.Path), char(tags.Patient), '_', char(tags.Tp), '_', char(tags.SequenceName), '.nii'];
                out_file = [out_file ; {f_out}];
                Table_out = [Table_out ; tags];
            end
        end
    end
    files_out.In1 = out_file;
    opt.Table_out = Table_out;
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



databScans = opt.Table_in(opt.Table_in.Type == categorical(cellstr('Scan')),:);
databROIs = opt.Table_in(opt.Table_in.Type == categorical(cellstr('ROI')),:);
UTag1 = unique(databScans.(Tag1));
UTag2 = unique(databScans.(Tag2));
in_files = {};
roi_files = {};
for i=1:length(UTag1)
    for j=1:length(UTag2)
        DbRois = databROIs(databROIs.(Tag1) == UTag1(i),:);
        DbRois = DbRois(DbRois.(Tag2) == UTag2(j),:);
        if size(DbRois, 1) == 0
            continue
        end
        roi = [char(DbRois.Path(1)), char(DbRois.Filename(1)), '.nii'];
        datab = databScans(databScans.(Tag1) == UTag1(i),:);
        datab = datab(datab.(Tag2) == UTag2(j),:);
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
        ROI_NaN = ROI{i};
        ROI_NaN(ROI_NaN == 0) = NaN;
        input{j} = read_volume(nifti_header, ROI_nifti_header{i}, 0, 'axial').*ROI_NaN;
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
trainedModel_loaded = load(opt.Clustering_model_path);
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
        [ClusteredVox, ~, score] = cluster(trainedModel, data);
        ClusteredVox(max(score, [], 2) < (opt.Classification_threshold/100)) = 0;
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
    ClusteredVox(max(score, [], 2) < (opt.Classification_threshold/100)) = 0;
    
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
    ind = ind+size(All_Data_Clean{i},1);
    
    ROI_cluster_header = ROI_nifti_header_Clean{i};
    % On a fait le même traitement sur les files_out (tout début du code) que sur l'ouverture des ROI. Il y a donc tout à penser que l'ordre des fichiers correspondra.
    ROI_cluster_header.fname = files_out.In1{i};
    ROI_cluster_header = rmfield(ROI_cluster_header, 'pinfo');
    ROI_cluster_header = rmfield(ROI_cluster_header, 'private');
    
    Cluster = write_volume(Cluster,ROI_nifti_header_Clean{i}, 'axial');
    spm_write_vol(ROI_cluster_header, Cluster);
end


% data_in_table.cluster(isnan(data_in_table.cluster)) = size(handles.color_RGB,1);
%
% Couleurs = handles.color_RGB;
% Cartes = nominal(data_in_table.Properties.VarNames(8:end-1));
%
% % A partir du classement des pixels, on calcule les statistiques du
% % clustering
% [MoyCartesTranches, ProbTranches, MoyCartesVolume, ProbVolume, Ecart_Type_Global, Sign, MoyGlobal] = AnalyseClusterGMM(data_in_table);
%
% ROI = char(unique(data_in_table.VOI));
%
% %On cree 2 strutures que l'on va sauvegarder avec chaque uvascroi. Ces
% %structures contiennent les informations et les statistiques du clustering.
% Informations = struct('Couleurs',Couleurs,'Cartes', Cartes , 'Modele', gmfit, 'Sign', Sign,'ROI',ROI);
% Statistiques = struct('MoyCartesTranches', MoyCartesTranches , 'ProbTranches', ProbTranches , 'MoyCartesVolume', MoyCartesVolume , 'ProbVolume', ProbVolume, 'Ecart_Type_Global', Ecart_Type_Global,'MoyGlobal', MoyGlobal);
% % IA_patient_list = {handles.MIA_data.database.name};
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
% save([strcat( handles.MIA_data.database(1).path,NomDossier2), filesep,answer{1} '-clusterGMM.mat'],'Informations', 'Statistiques');






