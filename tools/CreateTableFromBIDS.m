function [database] = CreateTableFromBIDS(folderBIDS, PathProject, SaveDbFlag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('folderBIDS')
    folderBIDS = '/home/cbrossard/Bureau/STEMRI/';
end
if ~exist('PathProject')
    PathProject = '/home/cbrossard/Documents/TestBIDS/';
end
if ~exist('SaveDbFlag')
    SaveDbFlag = 0;
end

if ~exist(PathProject)
    mkdir(PathProject);
end
if ~exist([PathProject, 'Raw_data'])
    mkdir([PathProject, 'Raw_data']);
end
if ~exist([PathProject, 'ROI_data'])
    mkdir([PathProject, 'ROI_data']);
end



database = table();

listPatients = dir([folderBIDS, 'sub-*']);
for i=1:length(listPatients)
    if listPatients(i).isdir
        listTP = dir([folderBIDS, listPatients(i).name, filesep, 'ses-*']);
        for j=1:length(listTP)
            if listTP(j).isdir
                % Unzip all .nii.gz files
                listScansToUnzip = dir([folderBIDS, listPatients(i).name, filesep, listTP(j).name, filesep, '**/*.nii.gz']);
                NameUncompressedFile = cell(length(listScansToUnzip),1);
                for ll=1:length(listScansToUnzip)
                    NameCompressedFile = [listScansToUnzip(ll).folder, filesep, listScansToUnzip(ll).name];
                    NameUncompressedFile{ll} = strrep(NameCompressedFile, '.nii.gz', '.nii');
                    if ~exist(NameUncompressedFile{ll})
                        gunzip(NameCompressedFile)
                    end
                end
                listScans = dir([folderBIDS, listPatients(i).name, filesep, listTP(j).name, filesep, '**/*.nii']);
                for k=1:length(listScans)
                    info = niftiinfo([listScans(k).folder, filesep, listScans(k).name]);
                    N = niftiread(info);
                    if isempty(N)
                        ProcessWrongNifti([listScans(k).folder, filesep, listScans(k).name]);
                        info = niftiinfo([listScans(k).folder, filesep, listScans(k).name]);
                        N = niftiread(info);
                    end
                    NbValeurs = unique(N);
                    if length(NbValeurs)<6
                        for l=1:length(NbValeurs)
                            if NbValeurs(l) == 0
                                continue
                            end
                            Tags = table();
                            Tags.Group = categorical(cellstr('Undefined'));
                            Tags.Patient = categorical(cellstr(listPatients(i).name(5:end)));
                            Tags.Tp = categorical(cellstr(listTP(j).name(5:end)));
                            Tags.Path = categorical(cellstr([PathProject, 'ROI_data', filesep]));
                            Tags.Filename = categorical(cellstr([listPatients(i).name(5:end), '_', listTP(j).name(5:end), '_', listScans(k).name(1:end-4), '_', num2str(NbValeurs(l))]));
                            Tags.Type = categorical(cellstr('ROI'));
                            Tags.IsRaw = categorical(1);
                            Tags.SequenceName = categorical(cellstr([listScans(k).name(1:end-4), '_', num2str(NbValeurs(l))]));
                            filename = [PathProject, 'ROI_data', filesep, listPatients(i).name(5:end), '_', listTP(j).name(5:end), '_', listScans(k).name(1:end-4), '_', num2str(NbValeurs(l)), '.nii'];
                            %copyfile([listScans(k).folder, filesep, listScans(k).name], filename)
                            if ~exist(filename)
                                N2=int16((N==NbValeurs(l)));
                                info2 = info;
                                info2.Datatype = 'int16';
                                [N3, Mat3] = CropROI(N2, info.Transform.T.');
                                info2.Transform.T = Mat3.';
                                info2.ImageSize = size(N3);
                                %N2 = permute(N2, [2,1,3]);
                                niftiwrite(N3, filename, info2)
                            end
                            database = [database; Tags];
                        end
                    else
                        Tags = table();
                        Tags.Group = categorical(cellstr('Undefined'));
                        Tags.Patient = categorical(cellstr(listPatients(i).name(5:end)));
                        Tags.Tp = categorical(cellstr(listTP(j).name(5:end)));
                        Tags.Path = categorical(cellstr([PathProject, 'Raw_data', filesep]));
                        Tags.Filename = categorical(cellstr([listPatients(i).name(5:end), '_', listTP(j).name(5:end), '_', listScans(k).name(1:end-4)]));
                        Tags.Type = categorical(cellstr('Scan'));
                        Tags.IsRaw = categorical(1);
                        Tags.SequenceName = categorical(cellstr(listScans(k).name(1:end-4)));
                        filename = [PathProject, 'Raw_data', filesep, listPatients(i).name(5:end), '_', listTP(j).name(5:end), '_', listScans(k).name(1:end-4), '.nii'];
                        infilename = [listScans(k).folder, filesep, listScans(k).name];
                        if ~exist(filename)
                            copyfile(infilename, filename)
                        end
                        jsonInfilename = strrep(infilename, '.nii', '.json');
%                         if exist(jsonInfilename)
%                             copyfile(jsonInfilename, strrep(filename, '.nii', '.json'))
%                         else
%                             CreateJsonFromNifti(filename, listScans(k).name(1:end-4), listTP(j).name(5:end))
%                         end
                        CreateJsonFromNifti(filename, listScans(k).name(1:end-4), listTP(j).name(5:end))
                        database = [database; Tags];
                    end
                    
                end
                
            end
            for ll=1:length(NameUncompressedFile)
                delete(NameUncompressedFile{ll})
            end
        end

        
    end
end


database.Properties.UserData.MIA_data_path = PathProject;
database.Properties.UserData.MIA_Raw_data_path = [PathProject, 'Raw_data/'];
database.Properties.UserData.MIA_Derived_data_path = [PathProject, 'Derived_data/'];
database.Properties.UserData.MIA_ROI_path = [PathProject, 'ROI_data/'];
database.Properties.UserData.db_filename = 'MIA_database.mat';


REF_BIDS_Folder = folderBIDS;
databasePwi = database(contains(cellstr(database.SequenceName),'pwi'),:);

for i=1:size(databasePwi,1)
    %% JSON automatically created by my script 'Import BIDS data'
    filename = [char(databasePwi.Path(i)), char(databasePwi.Filename(i)), '.json'];
    fid = fopen(filename, 'r');
    raw = fread(fid, inf, 'uint8=>char');
    fclose(fid);
    %raw = reshape(raw, 1,length(raw));
    J = jsondecode(raw);
    
    %% JSON where the valid data are stocked
    
    RecJSON = [REF_BIDS_Folder, 'sub-', char(databasePwi.Patient(i)), filesep, 'ses-', char(databasePwi.Tp(i)), filesep, 'pwi', filesep, '*.json'];
    listJSON = dir(RecJSON);
    for ll=1:length(listJSON)
        if contains(filename, listJSON(ll).name)
            listJSON = listJSON(ll);
            break
        end
    end
    assert(length(listJSON) == 1 )
    filename2 = [listJSON.folder, filesep, listJSON.name];
    fid2 = fopen(filename2, 'r');
    raw2 = fread(fid2, inf, 'uint8=>char');
    fclose(fid2);
    %raw = reshape(raw, 1,length(raw));
    J2 = jsondecode(raw2);
    
    %%
    J.RepetitionTime.value = J2.RepetitionTime*1000; %The J JSON specify the repetition time in ms, while the J2 JSON uses s.
    J.EchoTime.value = J2.EchoTime*1000; %The J JSON specify the echo time in ms, while the J2 JSON uses s.

    JMod = jsonencode(J);
    delete(filename)
    fidmod = fopen(filename, 'w');
    fwrite(fidmod, JMod, 'uint8');
    fclose(fidmod);
    
end















if SaveDbFlag
    save([PathProject, 'MIA_database.mat'], 'database')
end





end

