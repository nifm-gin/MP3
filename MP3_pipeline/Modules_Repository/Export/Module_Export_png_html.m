function [files_in,files_out,opt] = Module_Template(files_in,files_out,opt)
%Module_Template Template of a MP3 module, with comments and explanations.
%   This template can be modified and pasted into the Module_Repository
%   folder of the MP3 source code in order to be used in the MP3 GUI.

%   MP3 is an open source software aimed to support medical researchers all
%   along their studies. Indeed, MP³ offers a graphical interface to
%   convert, visualize, compute and analyze medical images. It allows to
%   use well known processes or to develop yours and to apply them in
%   complex pipelines. All the post processing you used to apply on each
%   subset of your data can now be stored, reproduced and parallelized.
%   Although it was initially developed for MRI data, this software can be
%   used for all medical images that matches one of the required formats.
%   It has for example been used on scanner data. MP³ has been developed at
%   the Grenoble Institute of Neurosciences, in France, and is downoadable
%   at: https://github.com/nifm-gin/MP3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(opt)
    %% Here define every parameter needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    
    module_parameters(:,1)   = {'output_filename','Export_Log.txt'};
    module_parameters(:,2)   = {'OutputSequenceName','AllName'};
    module_parameters(:,3)   = {'savingPath',  '/home/aurelien/Documents/ExportHTML'};
    module_parameters(:,4)   = {'dirName', 'Study0'};
    module_parameters(:,5)   = {'applyBbox', 'Yes'};
    module_parameters(:,6)   = {'tableName', 'Vf'};
    module_parameters(:,7)   = {'AutomaticJobsCreation', 'No'};
    module_parameters(:,8)   = {'makeMain', 'No'};
    module_parameters(:,9)   = {'dispSlice', 2};
    module_parameters(:,10)  = {'Cmap', 'gray'};
    module_parameters(:,11)  = {'Clipped', 'No'};   
    
    %% System parameters : Do not modify without understanding the behaviour of the software.
    
    system_parameters(:,1)   = {'RefInput',1};
    system_parameters(:,2)   = {'InputToReshape',1};
    
    
    %% Initialisation parameters : Do not modify without understanding the behaviour of the software.
    
    initialisation_parameters(:,1)   = {'folder_out',''};
    initialisation_parameters(:,2)   = {'flag_test',true};
    initialisation_parameters(:,3)   = {'Table_in', table()};
    initialisation_parameters(:,4)   = {'Table_out', table()};
    
    Parameters = [module_parameters, system_parameters, initialisation_parameters];
    
    opt.Module_settings = psom_struct_defaults(struct(),Parameters(1,:),Parameters(2,:));
    
    
    %% Each line displayed to the user :
    
    % --> user_parameter(1,:) = user_parameter_list: String displayed to
    % the user
    
    % --> user_parameter(2,:) = user_parameter_type: Type of the needed
    % parameter. Types available so far:
    % Classic types:
    %   - Text: In order to display text in the Description window.
    %   - char: In order to ask the user to type a string.
    %   - cell: In order to ask the user to choose one value between
    %   several.
    %   - numeric: In order to ask the user to type a number.
    % Scans:
    %   - 1Scan: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   Scan.
    %   - XScan: In order to ask the user to select one or several values
    %   of the 'SequenceName' tag among the ones of its database whose type
    %   is Scan.
    %   - 1ROI: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   ROI.
    %   - XROI: In order to ask the user to select one or several values of
    %   the 'SequenceName' tag among the ones of its database whose type is
    %   ROI.
    %   - 1ScanOr1ROI: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   Scan, ROI, or Cluster.
    %   - XScanOrXROI: In order to ask the user to select one or several
    %   values of the 'SequenceName' tag among the ones of its database
    %   whose type is Scan, ROI, or Cluster.
    %   - 1Cluster: In order to ask the user to select 1 value of the
    %   'SequenceName' tag among the ones of its database whose type is
    %   Cluster.
    %   - XCluster: In order to ask the user to select one or several
    %   values of the 'SequenceName' tag among the ones of its database
    %   whose type is Cluster.
    %   - 1Scan1TPXP: In order to ask the user to select one value of the
    %   'SequenceName' tag, one value of the 'Tp' tag, and one or several
    %   values of the 'Patient' tag among the ones of its database whose
    %   type is Scan. (Only use in the Coregistrations modules so far, in
    %   order to coregister scans between timepoints or patients.
    
    % --> user_parameter(3,:) = parameter_default: for parameters of type
    % 'cell', here lies the different values that will be displayed to the
    % user.
    
    % --> user_parameter(4,:) = psom_parameter_list: The name of the
    % parameter in the PSOM structure: 'module_parameters'
    
    % --> user_parameter(5,:) = Scans_input_DOF: Degrees of Freedom for
    % the user to choose the scan: for scans/ROI/Clusters inputs, please
    % type here the names of the tags selectables by the user. In fact, if
    % you have never heard of the tags or tried to modify them, dont touch
    % them and let the defaults values :
    %   - 1Scan, XScan, 1ROI, XROI, 1ScanOr1ROI, XScanOrXROI, 1Cluster,
    %   XCluster: {'SequenceName'}
    %   - 1Scan1TPXP: {'SequenceName', 'Tp', 'Patient'}
    
    % --> user_parameter(6,:) = IsInputMandatoryOrOptional: For
    % scans/ROI/Clusters, please type here if the input is Mandatory or
    % Optional. If none, the input is set as Optional.
    
    % --> user_parameter(7,:) = Help : text data which describe the
    % parameter (it will be display to help the user)
    
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'This module exports .png data designed to be used with the HTML viewer available on the NIFM-GIN git page'
        }'};
    user_parameter(:,2)   = {'Parameters','','','','', '', ''};
    user_parameter(:,3)   = {'    .Saving path','char','','savingPath','', '', {'Path to the directory where the results will be stored'}};
    user_parameter(:,4)   = {'    .Directory Name','char','','dirName','', '', {'Name of the directory that will be created'}};
    user_parameter(:,5)   = {'    .Apply bounding box','cell',{'No', 'Yes'},'applyBbox','','',{'Crop png to brain bounding box?'}};
    user_parameter(:,6)   = {'    .Table Name','char','','tableName','','Mandatory',{'Name of the table in HTML viewer'}};
    user_parameter(:,7)   = {'    .Create main HTML file?','cell',{'No', 'Yes'},'makeMain','','',{'Create main HTML page at the end of the module?'}};
    user_parameter(:,8)   = {'Scans','','','','', '', ''};
    user_parameter(:,9)   = {'    .Reference scan','1Scan','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,10)   = {'    .Scans to compare','XScan','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,11)  = {'    .ROI for BBox','1ROI','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,12)  = {'    .Slice to show','numeric','','dispSlice','', 'Optional',''};
    user_parameter(:,13)  = {'    .Colormap','cell',{'gray', 'hot', 'phase', 'jet', 'cool', 'spring'},'Cmap','', 'Optional',''};
    user_parameter(:,14)  = {'    .Is Ref clipped?','cell', {'No', 'Yes'}, 'Clipped','', 'Optional',''};
    
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

% Here we generate the names of the files_out, thanks to the names of the
% files_in and some of the parameters selected by the user.
% This paragraph is peculiar to the number and the kind of the files_out.
% % Please modify it to automatically generate the name of each output file.
% if isempty(files_out)
% %     files_out = files_in;
%     opt.Table_out = opt.Table_in;
%     opt.Table_out.IsRaw(:) = categorical(0);
%     opt.Table_out.Path(:) = categorical(cellstr([opt.folder_out, filesep]));
%     if strcmp(opt.OutputSequenceName, 'AllName')
%         opt.Table_out.SequenceName(:) = categorical(cellstr(opt.output_filename));
%     elseif strcmp(opt.OutputSequenceName, 'Extension')
%         opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
%     end
%
%     opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
%     f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
%     files_out.In1{1} = f_out;
% %
% end


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Smoothing:module','Bad syntax, type ''help %s'' for more info.',mfilename)
end


%% If the test flag is true, stop here !
% It's the end of the initialization part.
% At the execution, opt.flag_test will be set to 0 and most of the
% initialization part will not be computed.

if opt.flag_test == 1
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now the core of the module (the operations on the files) starts %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get required colormap
switch strcmp(opt.Cmap, 'phase')
    case 1
        C1 = [0:1/255:1];
        C2 = [0:0.5/128:0.5, 0.5-0.5/128:-0.5/127:0];
        C3 = [1:-1/255:0];
        Cmap = [C1', C2', C3'];
    otherwise
        Cmap = eval(sprintf([opt.Cmap '(256)']));
end
Q = [0, ones(1,255)];
% Check if the files_in exist
[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end

if isempty(opt.savingPath)
    folder_out_silesep = strsplit(opt.folder_out, filesep);
    folder_out_silesep(end) = {opt.dirName};
    folder_out_silesep(end+1) = {''};
else
    folder_out_silesep = strsplit(opt.savingPath, filesep);
    folder_out_silesep(end+1) = {opt.dirName};
    folder_out_silesep(end+1) = {''};
end

% Create the right repos
studyDir = fullfile(filesep, folder_out_silesep{:}, filesep);
htmlDir = fullfile(studyDir, 'HTML/');
resourcesDir = fullfile(studyDir, 'Resources/');
imgDir = fullfile(studyDir, 'img/');

if exist(studyDir, 'dir') ~= 7
    mkdir(studyDir);
end
if exist(resourcesDir, 'dir') ~= 7
    pathMP3 = split(mfilename('fullpath'),'MP3_pipeline',1);
    pathData = fullfile(pathMP3{1},'data/Resources/');
    copyfile(pathData, resourcesDir);
end
if exist(htmlDir, 'dir') ~= 7
    mkdir(htmlDir);
end
if exist(imgDir, 'dir') ~= 7
    mkdir(imgDir);
end

Patient_listing = unique(opt.Table_in.Patient);
% Allow different time points per patient
for a=1:numel(Patient_listing)
    Tp_listing{a} = unique(opt.Table_in.Tp(opt.Table_in.Patient == Patient_listing(a)));
end

Seq_listing = unique(opt.Table_in.SequenceName(opt.Table_in.Type == 'Scan'));

[~, refName, ~] = fileparts(files_in.In1{1});
refSeqName = opt.Table_in(opt.Table_in.Patient == Patient_listing(1) & opt.Table_in.Filename == refName, :).SequenceName;
%sortingTable = opt.Table_in(opt.Table_in.Patient == Patient_listing(1) & opt.Table_in.Filename == refName, :);
I = find(refSeqName == Seq_listing);

switch I>1
    case 1
        if I ~= numel(Seq_listing) %If I is not the last element
            sortingIdx = [I, 1:I-1, I+1:numel(Seq_listing)]';
        else
            sortingIdx = [I, 1:numel(Seq_listing)-1]';
        end
    otherwise
        sortingIdx = [1:numel(Seq_listing)]';
end
Seq_listing=Seq_listing(sortingIdx);

if strcmp(opt.applyBbox, 'Yes')
    BboxRoi = unique(opt.Table_in.SequenceName(opt.Table_in.Type == 'ROI'));
end
% Create a table that will contain png names
colNames = {};
varTypes = {};
for x = 1:numel(Patient_listing)
    for y = 1:numel(Tp_listing{x})
        colNames{end+1} = strcat(char(Patient_listing(x)), '_', char(Tp_listing{x}(y)));
        varTypes{end+1} = 'string';
    end
end

rowNames = {};
for i = 1:numel(Seq_listing)
    rowNames{end+1} = char(Seq_listing(i));
end

localTable = table('Size', [numel(rowNames), numel(colNames)],'VariableTypes', varTypes, 'VariableNames', colNames, 'RowNames',rowNames);
minTable = zeros(numel(rowNames), numel(colNames));
maxTable = zeros(numel(rowNames), numel(colNames));
colIndex = 0;
% fill the table and generate .png files
for x = 1:numel(Patient_listing)
    for y = 1:numel(Tp_listing{x})
        colIndex = colIndex+1;
        % process reference scan
        %sub_database = sortrows(opt.Table_in(opt.Table_in.Patient == Patient_listing(x) & opt.Table_in.Tp == Tp_listing{x}(y) & opt.Table_in.Type == 'Scan', :), {'SequenceName'}, {'ascend'});
        sub_database = opt.Table_in(opt.Table_in.Patient == Patient_listing(x) & opt.Table_in.Tp == Tp_listing{x}(y) & opt.Table_in.Type == 'Scan', :);
        sub_database=sub_database(sortingIdx,:);
        roi_database = opt.Table_in(opt.Table_in.Patient == Patient_listing(x) & opt.Table_in.Tp == Tp_listing{x}(y) & opt.Table_in.Type == 'ROI', :);
        % load the scan of reference
        scan_of_reference.header =   spm_vol([char(sub_database.Path(1)) char(sub_database.Filename(1)) '.nii']);
        scan_of_reference.data=  read_volume(scan_of_reference.header, scan_of_reference.header, 0, 'Axial');
        Slices = {scan_of_reference.header.dim};
        midSlice = opt.dispSlice;
        %midSlice = floor(Slices{1}(end)/2);
        
        % get bounding box if asked
        if strcmp(opt.applyBbox, 'Yes')
            roi.header = spm_vol([char(roi_database.Path(1)) char(roi_database.Filename(1)) '.nii']);
            roi.data = read_volume(roi.header, scan_of_reference.header, 0, 'Axial');
            switch scan_of_reference.header.dim(3) > 1
                case 0
                    roiSlice = roi.data;
                    Data = scan_of_reference.data;
                case 1
                    roiSlice = roi.data(:,:,midSlice);
                    Data = scan_of_reference.data(:,:,midSlice,1);
            end
            Data(roiSlice ==0 ) = nan;
            BBox = regionprops(roiSlice, 'BoundingBox');
            %pngData = imcrop(scan_of_reference.data(:,:,midSlice,1), BBox.BoundingBox);
            pngDataRaw = imcrop(Data, BBox.BoundingBox);
        else
            pngDataRaw = scan_of_reference.data(:,:,midSlice,1);
        end
        
        % Scale data on [0,1] before png conversion
        pngData = double(pngDataRaw - min(min(pngDataRaw)));
        pngData = double(pngData./max(max(pngData)));
        pngData = im2uint8(pngData);
        pngData(pngDataRaw==min(min(pngDataRaw(pngDataRaw~=0))))=1;
        if x == 1 && y == 1
            barWidth = floor(size(pngData,2)/20);
        end
        %         pngData = pngData - min(min(pngData));
        %         pngData = pngData/max(max(pngData));
        
        % Remove
        bar = linspace(255,0,size(pngData,1))';
        toWrite = [pngData, repmat(bar,[1,barWidth])];
        imwrite(toWrite, Cmap,[imgDir char(Patient_listing(x)) '_' char(Tp_listing{x}(y)) '_' char(sub_database.SequenceName(1)) '.png'], 'Transparency', Q);
        localTable(1, colIndex) = {[imgDir char(Patient_listing(x)) '_' char(Tp_listing{x}(y)) '_' char(sub_database.SequenceName(1)) '.png']};
        minTable(1, colIndex) = min(min(pngDataRaw(pngDataRaw~=0)));
        maxTable(1, colIndex) = max(max(pngDataRaw));
        % process following scans
        for i = 2:height(sub_database)
            scan.header =   spm_vol([char(sub_database.Path(i)) char(sub_database.Filename(i)) '.nii']);
            scan.data=  read_volume(scan.header, scan_of_reference.header, 0, 'Axial'); % opening in ref space
            %             Slices = {scan.header.dim};
            %             midSlice = floor(Slices{1}(end)/2);
            
            % get bounding box if asked
            if strcmp(opt.applyBbox, 'Yes')
                roi.header = spm_vol([char(roi_database.Path(1)) char(roi_database.Filename(1)) '.nii']);
                roi.data = read_volume(roi.header, scan_of_reference.header, 0, 'Axial');
                switch scan_of_reference.header.dim(3) > 1
                    case 0
                        roiSlice = roi.data;
                        Data = scan.data;
                    case 1
                        roiSlice = roi.data(:,:,midSlice);
                        Data = scan.data(:,:,midSlice,1);
                end
                Data(roiSlice ==0 ) = nan;
                BBox = regionprops(roiSlice, 'BoundingBox');
                %pngData = imcrop(scan_of_reference.data(:,:,midSlice,1), BBox.BoundingBox);
                pngDataRaw = imcrop(Data, BBox.BoundingBox);
            else
                pngDataRaw = scan_of_reference.data(:,:,midSlice,1);
            end
            % Scale data on [0,1] before png conversion
            pngData = double(pngDataRaw - min(min(pngDataRaw)));
            pngData = double(pngData./max(max(pngData)));
            pngData = im2uint8(pngData);
            bar = linspace(255,0,size(pngData,1))';
            pngData(pngDataRaw==min(min(pngDataRaw(pngDataRaw~=0))))=1;
            toWrite = [pngData, repmat(bar,[1,barWidth])];
            imwrite(toWrite, Cmap,[imgDir char(Patient_listing(x)) '_' char(Tp_listing{x}(y)) '_' char(sub_database.SequenceName(i)) '.png'], 'Transparency', Q)
            localTable(i, colIndex) = {[imgDir char(Patient_listing(x)) '_' char(Tp_listing{x}(y)) '_' char(sub_database.SequenceName(i)) '.png']};
            minTable(i, colIndex) = min(min(pngDataRaw(pngDataRaw~=0)));
            maxTable(i, colIndex) = max(max(pngDataRaw));
        end
    end
end

%% Make the HTML table
% Remove file if already exist
if exist([htmlDir 'Table' opt.tableName '.html'], 'file') == 2
    delete([htmlDir 'Table' opt.tableName '.html'])
end
htmlFid = fopen([htmlDir 'Table' opt.tableName '.html'], 'a');
if htmlFid < 0
    error('Unable to create html table')
end

fprintf(htmlFid, '<table>\n');
fprintf(htmlFid, '<tr>\n');
fprintf(htmlFid, '<th>  </th>\n');
for i = 1:size(localTable, 2)
    fprintf(htmlFid, ['<th> ' localTable.Properties.VariableNames{i} ' </th>\n']);
end
fprintf(htmlFid, '</tr>\n');

for l = 1:size(localTable, 1)
    fprintf(htmlFid, '<tr>\n');
    fprintf(htmlFid, ['<td class="myLabel"> <p class="rotate">' localTable.Properties.RowNames{l} '</p> </td>\n']);
    for i = 1:size(localTable, 2)
        if ~ismissing(localTable(l,i).Variables)
            switch (l==1 & strcmp(opt.Clipped, 'Yes'))
                case 1
                    fprintf(htmlFid, strcat('<td> ', strcat('>',num2str(maxTable(l,i))), '<br/> <img src="', localTable(l,i).Variables, '" id="', localTable.Properties.VariableNames{i},...
                    '_', localTable.Properties.RowNames{l}, '" onclick="singleClick(this);" /> <br/>', strcat('<',num2str(minTable(l,i))), ' </td>\n'));
                case 0
                    fprintf(htmlFid, strcat('<td> ', num2str(maxTable(l,i)), '<br/> <img src="', localTable(l,i).Variables, '" id="', localTable.Properties.VariableNames{i},...
                    '_', localTable.Properties.RowNames{l}, '" onclick="singleClick(this);" /> <br/>', num2str(minTable(l,i)), ' </td>\n'));
            end
        else
            fprintf(htmlFid, strcat('<td> </td>\n'));
        end
    end
    fprintf(htmlFid, '</tr>\n');
end
fprintf(htmlFid, '</table>\n');
fclose(htmlFid);

if strcmp(opt.makeMain, 'Yes')
    S = strsplit(studyDir,'/');
    mainName = sprintf('%s.html', S{end-1});
    mainFid = fopen([studyDir mainName], 'a');
    if mainFid < 0
        error('Unable to create html table')
    end
    fprintf(mainFid,...
        strcat('<html>\n',...
        '<head>\n',...
        '<title> Viewer </title>\n',...
        '<meta charset="utf-8"></meta>\n',...
        '<link rel="stylesheet" href="Resources/style.css" />\n', ...
        '<script src="Resources/jquery-3-4-1.js"></script>\n', ...
        '</head>'));
    fprintf(mainFid, ...
        strcat('<body>\n',...
        '<table class="selTable">\n',...
        '<td>Map to display:</td>\n',...
        '<td> <select id="mapToDisplay">'));
    % Retrieve the names of the displayable tables
    html_files	= dir(htmlDir);
    html_files    = html_files(~[html_files.isdir]);
    [~, html_names, ~] = cellfun(@fileparts, {html_files.name}, 'UniformOutput', false);
    for k=1:numel(html_names)
        [~, shortName, ~] = fileparts(html_names{k});
        fprintf(mainFid, strcat('<option value="', shortName, '">', shortName, '</option>'));
    end
    fprintf(mainFid, '</select> </td>\n');
    fprintf(mainFid, '<td><button type="button" id="show" onclick="showMap(this);">Show</button></td>\n');
    fprintf(mainFid, '<td><input type="range" min="-255" max="255" defaultValue="0" class="slider" id="brightRange" onmouseup="brightAdjust();" autocomplete="off"></td>');
    fprintf(mainFid, '<td><p id="demo"></p></td>');
    fprintf(mainFid, ...
        strcat('</table>\n', ...
        '</table\n>',...
        '<div id="includedContent"></div>\n',...
        '<script src="Resources/Interaction.js"></script>\n',...
        '<body style="background-color:black;">',...
        '</body>\n',...
        '</html>'));
    
end

%% It's already over !

