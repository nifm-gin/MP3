function varargout = MP3_pipeline(varargin)
% MP3_PIPELINE MATLAB code for MP3_pipeline.fig
%      MP3_PIPELINE, by itself, creates a new MP3_PIPELINE or raises the existing
%      singleton*.
%
%      H = MP3_PIPELINE returns the handle to a new MP3_PIPELINE or the handle to
%      the existing singleton*.
%
%      MP3_PIPELINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MP3_PIPELINE.M with the given input arguments.
%
%      MP3_PIPELINE('Property','Value',...) creates a new MP3_PIPELINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MP3_pipeline_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MP3_pipeline_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MP3_pipeline

% Last Modified by GUIDE v2.5 19-May-2020 15:05:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MP3_pipeline_OpeningFcn, ...
                   'gui_OutputFcn',  @MP3_pipeline_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT





% --- Executes just before MP3_pipeline is made visible.
function MP3_pipeline_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MP3_pipeline (see VARARGIN)

% Choose default command line output for MP3_pipeline
handles.output = hObject;

handles.MP3_data = varargin{3};
Spl = strsplit(mfilename('fullpath'), filesep);
%Spl{end} = ['Modules_Repository', filesep, '**', filesep, 'Module_*.m'];
Spl{end} = 'Modules';
path_modules = strjoin(Spl, filesep);
list_mod = dir(path_modules);

Mod_listing = {};
for i=3:length(list_mod)
    Mod_listing = [Mod_listing, list_mod(i).name];
    list_bis = dir([path_modules, filesep, list_mod(i).name, filesep]);

    for l=3:length(list_bis)
        [~, nam, ext] = fileparts(list_bis(l).name);
        if (~list_bis(l).isdir && ~strcmp(ext, '.m')) || (~list_bis(l).isdir && ~startsWith(nam, 'Module_')) || isempty(nam)% exist([list_mod(i).folder, filesep, list_mod(i).name])==
            continue
        end
        Mod_listing = [Mod_listing, ['   .', list_bis(l).name]];
        if list_bis(l).isdir
            %Mod_listing = [Mod_listing, list_mod(i).name];
            SubMod = dir([list_bis(l).folder, filesep, list_bis(l).name, filesep, 'Module_*.m']);
            for j=1:length(SubMod)
                Mod_listing = [Mod_listing, ['      ..', SubMod(j).name]];
            end
        end
    end
end

handles.Modules_listing = Mod_listing;


% handles.Modules_listing = {'Relaxometry', '   .T1map (Multi Inversion Time)', '   .T1map (Multi Angles)', '   .T2map', '   .Fit_T2_T2star',...
%                 '   .deltaR2', '   .deltaR2*', '   .MGE2Dfrom3D',...
%     'Perfusion', '   .Blood volume fraction (steady-state)', '   .Dynamic Susceptibility Contrast', '   .Vessel Size Imaging (steady-state)', ...
%                 '   .Vessel Density (steady-state)', '   .Cerebral blood flow (ASL)',  '   .Cerebral blood flow (ASL-Dynamic)',...
%                 '   .Inversion Efficiency (ASL_InvEff)',...
%      'Diffusion', '   .ADCmap',...
%      'Permeability', '   .Dynamic Contrast Enhancement (Phenomenology)', '   .Dynamic Contrast Enhancement (Quantitative)',...          
%      'Oxygenation', '   .R2prim', '   .SO2map', '   .CMRO2',...
%      'MRFingerprint', '   .Vascular MRFingerprint'...
%      'SPM', '   .SPM: Coreg (Est)', '   .SPM: Coreg (Est & Res)', '   .SPM: Reslice','   .SPM: Realign', ...
%      'Texture Analyses', '   .Texture Matlab',...
%      'Spatial', '   .Smoothing', '   .Shift images', '   .Arithmetic', '   .Normalization','   .Clip Image', '   .Brain Extraction (BET Function from FSL)',...
%      '   .Brain Mask (using PCNN3D function)', '   .FLIRT-FMRIB Linear Image Registration Tool (from FSL)', '   .Bias Estimation (MICO algorithm)', '   .Reshape (Extraction)',...
%      'Clustering', '   .Clustering GMM', ...
%      'Export', '   .Export data for deeplearing', ...
%      };
handles.Module_groups = {'Relaxometry','Perfusion', 'Diffusion', 'Permeability', 'Oxygenation', 'MRFingerprint', 'SPM', 'Spatial', 'Texture Analyses', 'Clustering', 'Export' };
 
set(handles.MP3_pipeline_module_listbox, 'String', handles.Modules_listing);
handles.Add_Tags_listing = handles.MP3_data.database.Properties.VariableNames;
handles.Add_list = handles.Add_Tags_listing;
set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_Tags_listing);
set(handles.MP3_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})
handles.Source_selected = handles.Add_Tags_listing{1};
%handles.MP3_pipeline_Unique_Values_Selection = 
handles.MP3_pipeline_TagsToPrint = {'Patient', 'Tp', 'SequenceName'};
handles.Remove_Tags_listing = {'NoMoreTags'};
handles.Select_Number_Workers.String = 'Max';
handles.Remove_selected = handles.Remove_Tags_listing{1};
%handles.MP3_data.database.IsRaw = categorical(handles.MP3_data.database.IsRaw);

handles.MP3_pipeline_TmpDatabase = handles.MP3_data.database;
handles.MP3_pipeline_Filtered_Table = handles.MP3_pipeline_TmpDatabase;
CellsColoured = DisplayColoredTable(handles.MP3_pipeline_TmpDatabase, handles.MP3_pipeline_TagsToPrint);
handles.MP3_pipeline_Filtering_Table.Data = CellsColoured;
% handles.MP3_pipeline_Filtering_Table.Data = cellstr(handles.MP3_pipeline_Filtered_Table{:,handles.MP3_pipeline_TagsToPrint});
% colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
% for i=1:size(handles.MP3_pipeline_Filtering_Table.Data, 1)
%     handles.MP3_pipeline_Filtering_Table.Data{i,1} = colergen(rgb2hex([1 0 0]), handles.MP3_pipeline_Filtering_Table.Data{i,1});
% end



handles.MP3_pipeline_Filtering_Table.ColumnName = handles.MP3_pipeline_TagsToPrint;

% handles.module_parameters_string = {};
% handles.module_parameters_fields = {};

MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
handles.FilterParameters = {};
%MP3_pipeline_Unique_Values_Tag_CellSelectionCallback(hObject, eventdata, handles)
% Update handles structure
handles.MP3_pipeline_Unique_Values_Selection = {};

% update the 'String' of MP3_pipeline_pushMP3Selection and MP3_pipeline_pushMP3TPSelection push button
%data_selected =  MP3('finddata_selected',handles.MP3_data);
data_selected = finddata_selected(handles.MP3_data);
%if data_selected == 0
set(handles.MP3_pipeline_pushMP3Selection, 'String', [char(handles.MP3_data.database.Patient(data_selected(1))) '-' char(handles.MP3_data.database.Tp(data_selected(1))) ' only'])
set(handles.MP3_pipeline_pushMP3TPSelection, 'String', ['All time point of :' char(handles.MP3_data.database.Patient(data_selected(1)))])

%set(findall(gcf,'-property','FontName'), 'FontName', 'Courier')
set(findall(handles.MP3_pipeline_manager_GUI.Children, 'Tag', 'MP3_pipeline_module_parameters'), 'FontName', 'Courier')
%set(findall(gcf,'-property','FontSize'), 'FontSize', 30)
guidata(hObject, handles);

function CellsColoured = DisplayColoredTable(InputTable, Tags)
colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
CellsColoured = cellstr(InputTable{:,Tags});
for i=1:size(CellsColoured, 1)
    Path =  InputTable.Path(i);
    Folders = strsplit(char(Path), filesep);
    Folder = Folders(end-1);
    switch Folder{1}
        case 'Raw_data'
            if ~strcmp(char(InputTable.Type(i)), 'Deleted')
                Color = [0 1 1];%Cyan
            else
                Color = [1 0 0];%Red
            end
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex(Color), CellsColoured{i,j});
            end
        case 'Derived_data'
            if ~strcmp(char(InputTable.Type(i)), 'Deleted')
                Color = [0 1 0];%Green
            else
                Color = [1 0 0];%Red
            end
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex(Color), CellsColoured{i,j});
            end
        case 'Tmp'
            if ~strcmp(char(InputTable.Type(i)), 'Deleted')
                Color = [1 1 0];%Yellow
            else
                Color = [1 0 0];%Red
            end
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex(Color), CellsColoured{i,j});
            end
        case 'ROI_data'
            if ~strcmp(char(InputTable.Type(i)), 'Deleted')
                Color = [1 0 1];%Magenta
            else
                Color = [1 0 0];%Red
            end
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex(Color), CellsColoured{i,j});
            end
            
    end
    
    %%Files to be deleted will be displayed in red!
    
    if strcmp(char(InputTable.Type(i)), 'Deleted')
        for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex([1 0 0]), CellsColoured{i,j});%Red
        end
    end
    
end


% UIWAIT makes MP3_pipeline wait for user response (see UIRESUME)
% uiwait(handles.MP3_pipeline_manager_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = MP3_pipeline_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function merged_struct = smart_pipeline_merge(old_pipeline, new_pipeline)
ModToRename = intersect(fieldnames(old_pipeline), fieldnames(new_pipeline));
for i=1:length(ModToRename)
    if psom_cmp_var(old_pipeline.(ModToRename{i}), new_pipeline.(ModToRename{i}))
        warndlg('The module %s already exist. This new occurence has been deleted.', ModToRename{i})
        new_pipeline = rmfield(new_pipeline, ModToRename{i});
    else
        j=1;
        New_Name = [ModToRename{i},'_' num2str(j)];
        while isfield(new_pipeline, New_Name) || isfield(old_pipeline, New_Name)
            j = j+1;
            New_Name = [ModToRename{i},'_' num2str(j)];
        end
        new_pipeline = renameStructField(new_pipeline, ModToRename{i}, New_Name);
    end
    
end
ModToRename = intersect(fieldnames(old_pipeline), fieldnames(new_pipeline));
assert(isempty(ModToRename));

merged_struct = psom_merge_pipeline(old_pipeline, new_pipeline);

function [hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles)
%% Update of Filtering/Filtered Tables
%[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
NewTable = handles.MP3_pipeline_TmpDatabase;
for i=1:length(handles.FilterParameters)
    Tag = handles.FilterParameters{1,i}{1};
    TagTable = table();
    for j=2:length(handles.FilterParameters{1,i})
        SelectedValue = handles.FilterParameters{1,i}{j};
        TagTable = [TagTable;NewTable(NewTable{:,Tag}==SelectedValue,:)];
        TagTable = unique(TagTable);
    end
    NewTable = TagTable;
end

handles.MP3_pipeline_Filtered_Table = NewTable;
CellsColoured = DisplayColoredTable(NewTable, handles.MP3_pipeline_TagsToPrint);
handles.MP3_pipeline_Filtering_Table.Data = CellsColoured;
if ~isempty(handles.MP3_pipeline_Filtering_Table.Data)
    % set ColumnWidth to auto
    column_name = handles.MP3_pipeline_TagsToPrint;
    column_data = cellstr(NewTable{:,handles.MP3_pipeline_TagsToPrint});
    merge_Data = [column_name;column_data];
    dataSize = size(merge_Data);
    % Create an array to store the max length of data for each column
    maxLen = zeros(1,dataSize(2));
    % Find out the max length of data for each column
    % Iterate over each column
    for i=1:dataSize(2)
        
        % 2 cases: the ColumnName is longer than any of the text contained
        % in this column or otherwise
        if max(cellfun('length',merge_Data(2:end,i))) >= length(merge_Data{1,i})
            maxLen(1,i) = max(cellfun('length',merge_Data(2:end,i)))*8;
        else
            maxLen(1,i) = length(merge_Data{1,i})*4;
        end
    end
    % Some calibration needed as ColumnWidth is in pixels
    cellMaxLen = num2cell(maxLen);
    % Set ColumnWidth of UITABLE
    set(handles.MP3_pipeline_Filtering_Table, 'ColumnWidth', cellMaxLen);
    
end




handles.MP3_pipeline_Filtering_Table.ColumnName = handles.MP3_pipeline_TagsToPrint;
%% Update of Unique Values
if strcmp(handles.Source_selected , 'NoMoreTags')
    TagValues = {};
else
    TagValues = unique(handles.MP3_pipeline_Filtered_Table{:,handles.Source_selected});
end
handles.MP3_pipeline_Unique_Values_Tag.Data = cellstr(TagValues);

%% Update of modules parameters
MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles)
%Coloredlistbox = DisplayColoredListbox(handles.MP3_pipeline_pipeline_listbox_Raw, handles);
%set(handles.MP3_pipeline_pipeline_listbox,'String', Coloredlistbox);
%[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
guidata(hObject, handles);



% --- Executes on button press in MP3_pipeline_add_module_button.
function MP3_pipeline_add_module_button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_add_module_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.MP3_pipeline_manager_GUI, 'pointer', 'watch');
drawnow;
if ~isfield(handles, 'new_module') 
    set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
    return
end
if ~isfield(handles.new_module, 'opt')
    set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
    return
end
[new_pipeline, output_database] = MP3_pipeline_generate_psom_modules(handles.new_module, handles.FilterParameters, handles.MP3_pipeline_TmpDatabase, handles.MP3_data.database.Properties.UserData.MP3_data_path, 1);

if isempty(fieldnames(new_pipeline)) && isempty(output_database)
    set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
    return
end


% if isfield(handles, 'tmp_database')
%     handles.tmp_database = [handles.tmp_database ; output_database];
% else
%     handles.tmp_database = output_database;
% end
    
if isfield(handles, 'MP3_pipeline_ParamsModules')
    old_modules = handles.MP3_pipeline_ParamsModules;
    %old_pipeline = handles.psom.pipeline;
else
    old_modules = struct();
    %old_pipeline = struct();
end

%merged_pipe = smart_pipeline_merge(old_pipeline, new_pipeline);


%handles.psom.pipeline = new_pipeline;
Name_New_Mod = handles.new_module.module_name;
j=2;
while ~isempty(intersect(fieldnames(old_modules), Name_New_Mod))
    Name_New_Mod = [handles.new_module.module_name, '_', num2str(j)];
    j=j+1;
end

SaveModule = struct();
SaveModule.Filters = handles.FilterParameters;
SaveModule.ModuleParams = handles.new_module;
SaveModule.OutputDatabase = output_database;
SaveModule.Jobs = new_pipeline;
handles.MP3_pipeline_ParamsModules.(Name_New_Mod) = SaveModule;

%handles.psom.Output_databases = setfield(old_databases, Name_New_Mod, output_database);
%handles.psom.Modules = setfield(old_modules, Name_New_Mod, new_pipeline);


%handles.MP3_pipeline_Filtered_Table = [handles.MP3_pipeline_Filtered_Table ; output_database];
%handles.MP3_pipeline_Filtering_Table.Data = cellstr(handles.MP3_pipeline_Filtered_Table{:,handles.MP3_pipeline_TagsToPrint});
%handles.psom.pipeline = merged_pipe;
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
%handles.MP3_pipeline_TmpDatabase = unique([handles.MP3_pipeline_TmpDatabase ; output_database]);


%handles.MP3_pipeline_pipeline_listbox_Raw = fieldnames(handles.MP3_pipeline_ParamsModules);
[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%guidata(hObject, handles);
%MP3_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)
%MP3_pipeline_Remove_Tag_Button_Callback(hObject, eventdata, handles)


module_listing = get(handles.MP3_pipeline_pipeline_listbox,'String');
set(handles.MP3_pipeline_pipeline_listbox,'String', [module_listing' {handles.new_module.module_name}]');

Coloredlistbox = DisplayColoredListbox(handles.MP3_pipeline_pipeline_listbox.String, handles);
set(handles.MP3_pipeline_pipeline_listbox,'String', Coloredlistbox);
[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
MP3_pipeline_JobsList_Callback(hObject, eventdata, handles)


set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
guidata(hObject, handles);

%MP3_pipeline_exectute_module_button_Callback(hObject, eventdata, handles)






% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.MP3_pipeline_manager_GUI)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.MP3_pipeline_manager_GUI,'Name') '?'],...
                     ['Close ' get(handles.MP3_pipeline_manager_GUI,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.MP3_pipeline_manager_GUI)



% --- Executes on selection change in MP3_pipeline_module_parameters.
function MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_module_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_module_parameters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_module_parameters


parameter_selected = get(handles.MP3_pipeline_module_parameters,'Value');
% display the help associated to the parameter selected
if isempty(parameter_selected) || ~isfield(handles, 'new_module')
    return
end
set(handles.MP3_pipeline_parameter_setup_text, 'String', handles.new_module.opt.table.Help{parameter_selected});

switch handles.new_module.opt.table.Type{parameter_selected}
    case 'Text'
        %Text = handles.new_module.opt.table.Default{1}';
        %disp('l MP3_pipeline ligne 225')

        table.ColumnFormat = {'char'};
        table.data = '';
        table.columnName = '';
        table.editable = false;
    case '1Scan1TPXP'
        SequenceType_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Scan')));
        TP_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.Tp(handles.MP3_pipeline_Filtered_Table.Type == 'Scan')));
        Patients_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.Patient(handles.MP3_pipeline_Filtered_Table.Type == 'Scan')));
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
           
            table.data(1:numel(TP_listing),3) = cellstr(TP_listing);
            table.data(1:numel(TP_listing),4) = {false};

            table.data(1:numel(Patients_listing),5) = cellstr(Patients_listing);
            table.data(1:numel(Patients_listing),6) = {false};
        else
            Def = handles.new_module.opt.table.Default{parameter_selected};
            Def2 = Def(:,2);
            Def = Def(~cellfun('isempty', Def2'),1:2);
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            NamesSelected = {Def{cell2mat(Def(:,2)),1}};
            IndSelected = find(ismember(SequenceType_listing, NamesSelected)) ;
            for i=1:length(IndSelected)
               table.data(IndSelected(i),2) = {true};
            end
            
            Def = handles.new_module.opt.table.Default{parameter_selected};
            if size(Def,2)<4
                table.data(1:numel(TP_listing),3) = cellstr(TP_listing);
                table.data(1:numel(TP_listing),4) = {false};
            else
                Def2 = Def(:,4);
                Def = Def(~cellfun('isempty', Def2'),3:4);
                table.data(1:numel(TP_listing),3) = cellstr(TP_listing);
                table.data(1:numel(TP_listing),4) = {false};
                TPSelected = {Def{cell2mat(Def(:,2)),1}};
                IndSelected = find(ismember(TP_listing, TPSelected)) ;
                for i=1:length(IndSelected)
                   table.data(IndSelected(i),4) = {true};
                end
            end
            
            Def = handles.new_module.opt.table.Default{parameter_selected};
            if size(Def,2)<6
                table.data(1:numel(Patients_listing),5) = cellstr(Patients_listing);
                table.data(1:numel(Patients_listing),6) = {false};
            else
                Def2 = Def(:,6);
                Def = Def(~cellfun('isempty', Def2'),5:6);
                table.data(1:numel(Patients_listing),5) = cellstr(Patients_listing);
                table.data(1:numel(Patients_listing),6) = {false};
                PatientSelected = {Def{cell2mat(Def(:,2)),1}};
                IndSelected = find(ismember(Patients_listing, PatientSelected)) ;
                for i=1:length(IndSelected)
                   table.data(IndSelected(i),6) = {true};
                end
            end
            
            %table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName', 'Tp', 'Patient'};
        table.columnName = {'SequenceName', 'Select ONE Input','Tp', 'Select ONE Input', 'Patient', 'Select Input'};
        handles.new_module.opt.ColumnNamesInput1Scan1TPXP = table.columnName;
        table.editable = [false true false true false true];
        table.ColumnFormat = {'char'}; 
        
%     case '1Scan'
%         SequenceType_listing = unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Scan'));
%         if isempty(handles.new_module.opt.table.Default{parameter_selected})
%             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
%             table.data(1:numel(SequenceType_listing),2) = {false};
%         else
%             Def = handles.new_module.opt.table.Default{parameter_selected};
%             NamesSelected = {Def{cell2mat(Def(:,2)),1}};
%             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
%             table.data(1:numel(SequenceType_listing),2) = {false};
%             for i=1:length(NamesSelected)
%                 IndSelected = find(strcmp(cellstr(SequenceType_listing), NamesSelected(i)));
%                 table.data(IndSelected,2) = {true};
%             end
%             %table.data = handles.new_module.opt.table.Default{parameter_selected};
%         end
%         %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
%         %% New Test
% %         SequenceType_listing = unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Scan'));
% %         if isempty(handles.new_module.opt.table.Default{parameter_selected})
% %             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
% %             table.data(1:numel(SequenceType_listing),2) = {false};
% %         else
% %             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
% %             table.data(1:numel(SequenceType_listing),2) = {false};
% %             NamesSelected = {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}};
% %             for i=1:length(NamesSelected)
% %                 IndSelected = find(strcmp(cellstr(SequenceType_listing), NamesSelected(i))) ;
% %                 table.data(IndSelected,2) = {true};
% %             end
% %             %table.data = handles.new_module.opt.table.Default{parameter_selected};
% %         end
%         %% Test
% %         SequenceType_listing = unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Scan'));
% %         table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
% %         table.data(1:numel(SequenceType_listing),2) = {false};
% %         if ~isempty(handles.new_module.opt.table.Default{parameter_selected}) && sum(strcmp(handles.new_module.opt.table.Default{parameter_selected}(2), cellstr(SequenceType_listing)))==1
% %             table.data(strcmp(handles.new_module.opt.table.Default{parameter_selected}(2), cellstr(SequenceType_listing)), 2) = {true};
% %         end
% %         
%         table.columnName = {'SequenceName', 'Select ONE Input'};
%         table.editable = [false true];
%         table.ColumnFormat = {'char'};
    case {'XScan', '1Scan'}
        SequenceType_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Scan')));
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            Def = handles.new_module.opt.table.Default{parameter_selected};
            NamesSelected = {Def{cell2mat(Def(:,2)),1}};
            if isempty(NamesSelected)
                NamesSelected = '';
            end
            %NamesSelected = {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}};
            IndSelected = find(ismember(SequenceType_listing, NamesSelected)) ;
            for i=1:length(IndSelected)
               table.data(IndSelected(i),2) = {true};
            end
            %table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        
        
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
        
    case {'1ScanOr1ROI', 'XScanOrXROI'}
        SequenceType_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type ~= 'Deleted')));
        handles.MP3_pipeline_Filtered_Table = handles.MP3_pipeline_Filtered_Table(handles.MP3_pipeline_Filtered_Table.Type ~= 'Deleted', :);
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            for i=1:numel(SequenceType_listing)
                table.data(i,3) = cellstr(unique(handles.MP3_pipeline_Filtered_Table.Type(handles.MP3_pipeline_Filtered_Table.SequenceName == SequenceType_listing(i))));
            end
        else
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            for i=1:numel(SequenceType_listing)
                table.data(i,3) = cellstr(unique(handles.MP3_pipeline_Filtered_Table.Type(handles.MP3_pipeline_Filtered_Table.SequenceName == SequenceType_listing(i))));
            end
            Def = handles.new_module.opt.table.Default{parameter_selected};
            NamesSelected = {Def{cell2mat(Def(:,2)),1}};
            %NamesSelected = {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}};
            if isempty(NamesSelected)
                NamesSelected = '';
            end
            IndSelected = find(ismember(SequenceType_listing, NamesSelected)) ;
            for i=1:length(IndSelected)
               table.data(IndSelected(i),2) = {true};
            end
            %table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        
        
        table.columnName = {'SequenceName', 'Select Input', 'Type'};
        table.editable = [false true false];
        table.ColumnFormat = {'char'};
%         if isempty(handles.new_module.opt.table.Default{parameter_selected})
%             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
%             table.data(1:numel(SequenceType_listing),2) = {false};
%             for i=1:numel(SequenceType_listing)
%                 table.data(i,3) = cellstr(unique(handles.MP3_pipeline_Filtered_Table.Type(handles.MP3_pipeline_Filtered_Table.SequenceName == SequenceType_listing(i))));
%             end
%         else
%             table.data = handles.new_module.opt.table.Default{parameter_selected};
%         end
%         %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
%         table.columnName = {'SequenceName', 'Select ONE Input', 'Type'};
%         table.editable = [false true false];
%         table.ColumnFormat = {'char'};  
%     case '1ROI'
%         SequenceType_listing = unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'ROI'));
%         if isempty(handles.new_module.opt.table.Default{parameter_selected})
%             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
%             table.data(1:numel(SequenceType_listing),2) = {false};
%         else
%             table.data = handles.new_module.opt.table.Default{parameter_selected};
%         end
%         %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
%         table.columnName = {'SequenceName', 'Select ONE Input'};
%         table.editable = [false true];
%         table.ColumnFormat = {'char'};
    case {'XROI', '1ROI'}
        SequenceType_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'ROI')));
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            Def = handles.new_module.opt.table.Default{parameter_selected};
            NamesSelected = {Def{cell2mat(Def(:,2)),1}};
            if isempty(NamesSelected)
                NamesSelected = '';
            end
            %NamesSelected = {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}};
            IndSelected = find(ismember(SequenceType_listing, NamesSelected)) ;
            for i=1:length(IndSelected)
               table.data(IndSelected(i),2) = {true};
            end
            %table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        
        
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
        
    case {'XCluster', '1Cluster'}
        SequenceType_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Cluster')));
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            Def = handles.new_module.opt.table.Default{parameter_selected};
            NamesSelected = {Def{cell2mat(Def(:,2)),1}};
            if isempty(NamesSelected)
                NamesSelected = '';
            end
            %NamesSelected = {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}};
            IndSelected = find(ismember(SequenceType_listing, NamesSelected)) ;
            for i=1:length(IndSelected)
               table.data(IndSelected(i),2) = {true};
            end
            %table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        
        
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
        
    case {'XROIOrXCluster', '1ROIOr1Cluster'}
        SequenceType_listing = cellstr(unique(handles.MP3_pipeline_Filtered_Table.SequenceName(handles.MP3_pipeline_Filtered_Table.Type == 'Cluster' | handles.MP3_pipeline_Filtered_Table.Type == 'ROI')));
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            Def = handles.new_module.opt.table.Default{parameter_selected};
            NamesSelected = {Def{cell2mat(Def(:,2)),1}};
            if isempty(NamesSelected)
                NamesSelected = '';
            end
            %NamesSelected = {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}};
            IndSelected = find(ismember(SequenceType_listing, NamesSelected)) ;
            for i=1:length(IndSelected)
               table.data(IndSelected(i),2) = {true};
            end
            %table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        
        
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
%      case 'XScanOrXROI'
%         SequenceType_listing = unique(handles.MP3_pipeline_Filtered_Table.SequenceName);
%         if isempty(handles.new_module.opt.table.Default{parameter_selected})
%             table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
%             table.data(1:numel(SequenceType_listing),2) = {false};
%         else
%             table.data = handles.new_module.opt.table.Default{parameter_selected};
%         end
%         %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
%         table.columnName = {'SequenceName', 'Select Input'};
%         table.editable = [false true];
%         table.ColumnFormat = {'char'};
    case 'cell'
        table.ColumnFormat = handles.new_module.opt.table.Default(parameter_selected);
        table.data = {getfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected})};
        table.columnName = handles.new_module.opt.table.PSOM_Fields{parameter_selected};
        table.editable = true;
        
    case {'char', 'numeric'}
        if strcmp(handles.new_module.opt.table.Default{parameter_selected}, 'Dont Show')
            table.data = '';
            table.ColumnFormat = {'char'};
            table.columnName = '';
            table.editable = false;
        else
            table.ColumnFormat = {'char'};
            table.data(1,1) = {getfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected})};
            table.columnName = handles.new_module.opt.table.PSOM_Fields{parameter_selected};
            table.editable = true;
        end
    case 'logical'
        if strcmp(handles.new_module.opt.table.Default{parameter_selected}, 'Dont Show')
            table.data = '';
            table.ColumnFormat = {'char'};
            table.columnName = '';
            table.editable = false;
        else
            table.ColumnFormat = {'logical'};
            table.data(1,1) = {getfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected})};
            %table.data(1,1)= handles.new_module.opt.table.Default(parameter_selected);
            table.columnName = handles.new_module.opt.table.PSOM_Fields{parameter_selected};
            table.editable = true;
        end
    case 'check'
        Temp =  {getfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected})};
        if ~iscell(Temp{1})
            Params = handles.new_module.opt.table.Default(parameter_selected);
            %table.data = Params;

            table.data(1:numel(Params{1}),1) = cellstr(Params{1});
            table.data(1:numel(Params{1}),2) = {false};
        else
            Params = {getfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected})};
            table.data = Params{1};
        end
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
    otherwise
        table.ColumnFormat = {'char'};
        table.data = '';
        table.columnName = '';
        table.editable = false;
end



%% update the setup table
set(handles.MP3_pipeline_parameter_setup_table, 'ColumnFormat', table.ColumnFormat);
% set names of the columns
set(handles.MP3_pipeline_parameter_setup_table, 'ColumnName', table.columnName);
% set data (default's parameters)
set(handles.MP3_pipeline_parameter_setup_table, 'Data', table.data);
% set each colomn editable
set(handles.MP3_pipeline_parameter_setup_table, 'columnEditable',  table.editable );
if ~isempty(table.data)
    % set ColumnWidth to auto

    merge_Data = [table.columnName; table.data];
    dataSize = size(merge_Data);
    % Create an array to store the max length of data for each column
    maxLen = zeros(1,dataSize(2));
    % Find out the max length of data for each column
    % Iterate over each column
    for i=1:dataSize(2)
        
        % 2 cases: the ColumnName is longer than any of the text contained
        % in this column or otherwise
        if max(cellfun('length',merge_Data(2:end,i))) >= length(merge_Data{1,i})
            maxLen(1,i) = max(cellfun('length',merge_Data(2:end,i)))*13;
        else
            maxLen(1,i) = length(merge_Data{1,i})*7;
        end
    end
    % Some calibration needed as ColumnWidth is in pixels
    cellMaxLen = num2cell(maxLen);
    % Set ColumnWidth of UITABLE
    set(handles.MP3_pipeline_parameter_setup_table, 'ColumnWidth', cellMaxLen);
    
end

guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function MP3_pipeline_module_parameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_module_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_setting_windows(hObject, eventdata, handles)




% --- Executes on selection change in MP3_pipeline_parameter_setup.
function MP3_pipeline_parameter_setup_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_parameter_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_parameter_setup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_parameter_setup



% --- Executes during object creation, after setting all properties.
function MP3_pipeline_parameter_setup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_parameter_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MP3_pipeline_clear_pipeline_button.
function [hObject, eventdata, handles] = MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_clear_pipeline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles, 'MP3_pipeline_ParamsModules')
    if ~strcmp(eventdata.Source.Tag, {'MP3_pipeline_Save_Module', 'MP3_pipeline_load_pipeline', 'MP3_pipeline_DeleteModule',...
            'MP3_pipeline_execute_button','MP3_pipeline_Module_UP', 'MP3_pipeline_Module_DOWN', 'MP3_pipeline_update_pipeline'})
        answer = questdlg('Are you sure you want to remove it?','Clean Pipeline', 'Yes', 'No', 'No');
        if strcmp(answer, 'No')
            return
        end
    end
    set(handles.MP3_pipeline_pipeline_listbox, 'String', '');
    if ~isempty(findobj('Tag', 'BioGraphTool'))
        close(findobj('Tag', 'BioGraphTool'));
    end
        handles = rmfield(handles, 'MP3_pipeline_ParamsModules');

end

%handles.tmp_database = table();
%handles.MP3_pipeline_TmpDatabase = handles.MP3_data.database;
%handles.MP3_pipeline_pipeline_listbox_Raw = {''};
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%set(handles.MP3_pipeline_JobsList, 'String', {''});
set(handles.MP3_pipeline_pipeline_listbox, 'String', {''});
set(handles.MP3_pipeline_pipeline_listbox, 'Value', 1);
[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
guidata(hObject, handles);



function edge_callbacks(hObject, eventdata, handles)
eventdata = [];
handles = guidata(findobj('Tag', 'MP3_pipeline_manager_GUI'));
module_list = get(handles.MP3_pipeline_module_popupmenu, 'String');
get(hObject, 'ID')
% sub_module = strfind(hObject.ID, '_');
% if ~isempty(sub_module)
%     module_name = hObject.ID(1:sub_module-1);
% else
%     module_name =hObject.ID;
% end
% idx = find(ismember(module_list, module_name));
% set(handles.MP3_pipeline_module_popupmenu, 'Value', idx)


function node_callbacks(hObject, ~, handles)
eventdata = [];
handles = guidata(findobj('Tag', 'MP3_pipeline_manager_GUI'));

pipeline_module_names = fieldnames(handles.psom.pipeline);
idx = strcmp(pipeline_module_names,hObject.ID);
module_selected = handles.psom.pipeline.(pipeline_module_names{idx});
module_selected.files_in
%% update handles.MP3_pipeline_module_parameters using the information the module_selected




%% test delete node
% handles.pipeline = rmfield(handles.pipeline, (char(hObject.ID)));
% update_setting_windows(hObject, eventdata, handles)
% if ~isempty(findobj('Tag', 'BioGraphTool'))
%     close(handles.biograph_fig);
% end
% handles = rmfield(handles, 'biograph_fig');
% handles = rmfield(handles, 'biograph_obj');
% handles.biograph_obj = psom_visu_dependencies(handles.pipeline);
% set(0, 'ShowHiddenHandles', 'on')
% handles.biograph_fig = gcf;
% set(handles.biograph_fig, 'Name', 'MP3 pipeline manager');
% guidata(findobj('Tag', 'MP3_pipeline_manager_GUI'), handles);

%%%%
%handles.new_module = handles.pipeline.(char(hObject.ID));
%update_setting_windows(hObject, eventdata, handles)


% switch node.ID
% case 'Node 1'
%     disp('Hello, I''m node 1');
% case 'Node 2'
%     disp('What''s up? This is node 2');
% case 'Node 3'
%     disp('Hi! You''ve clicked node 3');
% case 'Node 4'
%     disp('I''m node 4 and I don''t want to talk!');
% case 'Node 5'
%     disp('Who dares bother node 5??');
% end


% --- Executes when entered data in editable cell(s) in MP3_pipeline_parameter_setup_table.
function MP3_pipeline_parameter_setup_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_parameter_setup_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
parameter_selected = get(handles.MP3_pipeline_module_parameters,'Value');

%table_data = get(handles.MP3_pipeline_parameter_setup_table, 'Data');
if strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XScan') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XScanOrXROI') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XROI') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XCluster') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XROIOrXCluster') 
    %handles.new_module.opt.parameter_default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1};
    handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1Scan') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ScanOr1ROI') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ROI') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1Cluster') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ROIOr1Cluster') 
    if sum(cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2))) == 1 || sum(cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2))) == 0
        handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
        %handles.new_module.opt.table.Default{parameter_selected} = [handles.MP3_pipeline_parameter_setup_table.ColumnName(1); {handles.MP3_pipeline_parameter_setup_table.Data{cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,2)),1}}];
        
    end
% elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ScanOr1ROI')
%     if sum(cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,3))) == 1 || sum(cell2mat(handles.MP3_pipeline_parameter_setup_table.Data(:,3))) == 0
%         handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
%         
%     end
    
elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1Scan1TPXP')
%     if isempty(handles.new_module.opt.table.Default{parameter_selected})
%         handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
%     else
        A = handles.MP3_pipeline_parameter_setup_table.Data(:,2);
        A = A(~cellfun('isempty',A));
        %C = handles.new_module.opt.table.Default{parameter_selected};
        %D = C(:,2);
        %D = D(~cellfun('isempty',D));
        B = handles.MP3_pipeline_parameter_setup_table.Data(:,4);
        B = B(~cellfun('isempty',B));
        %E = C(:,4);
        %E = E(~cellfun('isempty',E));
        %if sum(cell2mat(A)) + sum(cell2mat(D)) == 1 && sum(cell2mat(B)) + sum(cell2mat(E)) == 0
        %    handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
        %elseif sum(cell2mat(A)) + sum(cell2mat(D)) == 0 && sum(cell2mat(B)) + sum(cell2mat(E)) == 1
        %    handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
        %end
        if (sum(cell2mat(A)) == 1 || sum(cell2mat(A)) == 0) && (sum(cell2mat(B)) == 1 || sum(cell2mat(B)) == 0)
            C = handles.MP3_pipeline_parameter_setup_table.Data(:,6);
            C = C(~cellfun('isempty',C));
            if any(cell2mat(C))
                DataToStore= handles.MP3_pipeline_parameter_setup_table.Data;
            elseif any(cell2mat(B))
                DataToStore= handles.MP3_pipeline_parameter_setup_table.Data(:,1:4);
            else
                DataToStore= handles.MP3_pipeline_parameter_setup_table.Data(:,1:2);
            end
            handles.new_module.opt.table.Default{parameter_selected} = DataToStore;
        end
%     end
    %handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data;
elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'check')
    handles.new_module.opt.Module_settings = setfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected},handles.MP3_pipeline_parameter_setup_table.Data); 
else
    handles.new_module.opt.Module_settings = setfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected},handles.MP3_pipeline_parameter_setup_table.Data{1,1}); 
    %handles.new_module.opt.table.Default{parameter_selected} = handles.MP3_pipeline_parameter_setup_table.Data{1,1};
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end
[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
    %patient_listing = table_data(:,1);
% patient_selected = patient_listing(find([table_data{:,2}]' == true));

% % case go back to all patient
% if table_data{1,2} == 1 && handles.new_module.files_in_filter_data{1,2} == 0
%     idex_patient = true(numel(handles.MP3_data.database.Patient),1);
%     table_data(1,2) = {true};
%     table_data(2:numel(unique(handles.MP3_data.database.Patient))+1,2) = {false};
%     % from all patient to a specific patient
% elseif sum(find([table_data{:,2}]' == true)) >1 && table_data{1,2} == 1
%     table_data{1,2} = false;
%     patient_selected = patient_listing(find([table_data{:,2}]' == true));
%     idex_patient = handles.MP3_data.database.Patient == patient_selected;
%     
%     % if all patient is selected
% elseif sum(find([table_data{:,2}]' == true))  == 1 && table_data{1,2} == 1
%     idex_patient = true(numel(handles.MP3_data.database.Patient),1);
%     % 1 or more patient (but not all)
% else
%     idex_patient = false(numel(handles.MP3_data.database.Patient),1);
%     for i=1:numel(patient_selected)
%         idex_patient = idex_patient | handles.MP3_data.database.Patient == patient_selected(i);
%     end
% end
% 
% tp_listing = table_data(:,3);
% tp_selected = tp_listing(find([table_data{:,4}]' == true));
% 
% % case go back to all time point
% if table_data{1,4} == 1 && handles.new_module.files_in_filter_data{1,4} == 0
%     idex_tp = true(numel(handles.MP3_data.database.Tp),1);
%     table_data(1,4) = {true};
%     table_data(2:numel(unique(handles.MP3_data.database.Tp))+1,4) = {false};
%     % from all time point to a specific time point
% elseif sum(find([table_data{:,4}]' == true)) >1 && table_data{1,4} == 1
%     table_data{1,4} = false;
%     tp_selected = tp_listing(find([table_data{:,4}]' == true));
%     idex_tp = handles.MP3_data.database.Tp == tp_selected;
%     
%     % if all time point is selected
% elseif sum(find([table_data{:,4}]' == true))  == 1 && table_data{1,4} == 1
%     idex_tp = true(numel(handles.MP3_data.database.Tp),1);
%     % if 1 or more time point (but not all)
% else
%     idex_tp = false(numel(handles.MP3_data.database.Tp),1);
%     for i=1:numel(tp_selected)
%         idex_tp = idex_tp | handles.MP3_data.database.Tp == tp_selected(i);
%     end
% end
% 
% SequenceName_listing =table_data(:,5);
% SequenceName_selected = SequenceName_listing(find([table_data{:,6}]' == true));
% if isempty(SequenceName_selected)
%     index_SequenceName =true(numel(handles.MP3_data.database.Tp,1));
% else
%      index_SequenceName = false(numel(handles.MP3_data.database.Tp),1);
%     for i=1:numel(SequenceName_selected)
%         index_SequenceName = index_SequenceName | handles.MP3_data.database.SequenceName == SequenceName_selected(i);
%     end
% end
% handles.new_module.files_in = char(handles.MP3_data.database.Filename(idex_patient & idex_tp & index_SequenceName));
% handles.new_module.files_in_index = find(idex_patient & idex_tp & index_SequenceName);
% handles.new_module.files_in_filter_data = table_data;
%guidata(hObject, handles);
MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)
%MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
  %             handles.new_module.files_in_filter_name = {'Patient Name', '', 'Time Point','', 'Sequence Name',''};
        %             Patient_listing = unique(handles.MP3_data.database.Patient);
        %             Tp_listing = unique(handles.MP3_data.database.Tp);
        %             SequenceName_listing = unique(handles.MP3_data.database.SequenceName);
        %             handles.new_module.files_in_filter_data = {'all', false, 'all', false, char(SequenceName_listing(1)),false};
        %             handles.new_module.files_in_filter_data(1:numel(Patient_listing)+1,1) = ['all' cellstr(Patient_listing)'];
        %             handles.new_module.files_in_filter_data(1,2) = {true};
        %             handles.new_module.files_in_filter_data(2:numel(Patient_listing),2) = {false};
        %
        %             handles.new_module.files_in_filter_data(1:numel(Tp_listing)+1,3) = ['all' cellstr(Tp_listing)'];
        %             handles.new_module.files_in_filter_data(1,4) = {true};
        %             handles.new_module.files_in_filter_data(2:numel(Patient_listing),4) = {false};
        %
        %             handles.new_module.files_in_filter_data(1:numel(SequenceName_listing),5) = cellstr(SequenceName_listing)';
        %             handles.new_module.files_in_filter_data(1:numel(SequenceName_listing),6) = {false};
        %
        %             handles.new_module.files_in_filter_format = {'char', 'logical','char', 'logical','char', 'logical' };
        %             handles.new_module.files_in_filter_editable = [0 1 0 1 0 1];
        %
        %             set(handles.MP3_pipeline_parameter_setup,  'String', cellstr(handles.MP3_data.database.nii));

function output_file_names = MP3_pipeline_generate_file_name(handles, database_indexes, output_extention)

output_file_names = [...
    char(handles.MP3_data.database.Patient(database_indexes)) ...
    repmat('-', [numel(database_indexes),1])...
    char(handles.MP3_data.database.Tp(database_indexes))...
    repmat('-', [numel(database_indexes),1])...
    repmat(output_extention, [numel(database_indexes),1])...
    repmat('_', [numel(database_indexes),1])...
    repmat(datestr(now,'yyyymmdd-HHMMSSFFF'), [numel(database_indexes),1])...
    repmat('.nii', [numel(database_indexes),1])
    ] ;
output_file_names  = strrep(cellstr(output_file_names), ' ', '');


% --- Executes on selection change in MP3_pipeline_add_tag_popupmenu.
function MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)

% hObject    handle to MP3_pipeline_add_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Source_selected = handles.MP3_pipeline_add_tag_popupmenu.String{handles.MP3_pipeline_add_tag_popupmenu.Value};
%[val, ind] = max(contains(handles.MP3_pipeline_Filtered_Table.Properties.VariableNames, handles.Source_selected));
if strcmp(handles.Source_selected, 'NoMoreTags')
    handles.MP3_pipeline_Unique_Values_Tag.Data = {};
    handles.MP3_pipeline_Unique_Values_Tag.ColumnName = {};
else

    TagValues = unique(handles.MP3_pipeline_Filtered_Table{:,handles.Source_selected});
%TagValues = [{'all'} ; unique(handles.MP3_pipeline_Filtering_Table.Data(:,ind))];
%TagValues = [{'all'} ; cellstr(char(unique(handles.MP3_data.database{:,handles.Source_selected})))];

    handles.MP3_pipeline_Unique_Values_Tag.Data = cellstr(TagValues);
    handles.MP3_pipeline_Unique_Values_Tag.ColumnName = {handles.Source_selected};
    % automatic column width
    if ~isempty(handles.MP3_pipeline_Filtering_Table.Data)
        % set ColumnWidth to auto
        column_name = get(handles.MP3_pipeline_Unique_Values_Tag,'columnName')';
        column_data = get(handles.MP3_pipeline_Unique_Values_Tag,'Data');
        %column_data = cellstr(NewTable{:,handles.MP3_pipeline_TagsToPrint});
        merge_Data = [column_name;column_data];
        dataSize = size(merge_Data);
        % Create an array to store the max length of data for each column
        maxLen = zeros(1,dataSize(2));
        % Find out the max length of data for each column
        % Iterate over each column
        for i=1:dataSize(2)
            
            % 2 cases: the ColumnName is longer than any of the text contained
            % in this column or otherwise
            if max(cellfun('length',merge_Data(2:end,i))) >= length(merge_Data{1,i})
                maxLen(1,i) = max(cellfun('length',merge_Data(2:end,i)))*5 + 35;
            else
                maxLen(1,i) = length(merge_Data{1,i})*7;
            end
        end
        % Some calibration needed as ColumnWidth is in pixels
        cellMaxLen = num2cell(maxLen);
        % Set ColumnWidth of UITABLE
        set(handles.MP3_pipeline_Unique_Values_Tag, 'ColumnWidth', cellMaxLen);
        
    end
    
end


guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_add_tag_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_add_tag_popupmenu


% --- Executes during object creation, after setting all properties.
function MP3_pipeline_add_tag_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_add_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');

end

function nii_json_fullfilename = fullfilename(handles, nii_index, ext)

nii_json_fullfilename = [char(handles.MP3_data.database.Path(nii_index)) char(handles.MP3_data.database.Filename(nii_index)) ext];


% --- Executes on selection change in MP3_pipeline_module_listbox.
function MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_module_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_module_listbox


% set(handles.MP3_pipeline_parameter_setup, 'Value', 1);
set(handles.MP3_pipeline_module_parameters, 'Value', 1);
module_selected = get(handles.MP3_pipeline_module_listbox, 'Value');
if isfield(handles, 'new_module')
    handles = rmfield(handles, 'new_module');
end

ismodule = 0;
% switch char(handles.Modules_listing(module_selected))
%     case handles.Module_groups	
%         module_parameters_string = [char(handles.Modules_listing(module_selected)) ' modules'];
%     case '   .SPM: Coreg (Est & Res)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Coreg_Est_Res('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Coreg_Est_Res(files_in,files_out,opt)';
%         handles.new_module.module_name = 'Module_Coreg_Est_Res';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%      case '   .SPM: Coreg (Est)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Coreg_Est('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Coreg_Est(files_in,files_out,opt)';
%         handles.new_module.module_name = 'Module_Coreg_Est';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;   
%     case '   .T2map'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T2map('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_T2map(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_T2map';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;   
% %     case '   .ADCmap'
% %         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_ADCmap('',  '', '');
% %         handles.new_module.command = '[files_in,files_out,opt] = Module_ADCmap(char(files_in),files_out,opt)';
% %         handles.new_module.module_name = 'Module_ADCmap';
% %         module_parameters_string = handles.new_module.opt.table.Names_Display;
% %         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
% %         ismodule = 1;
%     case '   .Fit_T2_T2star'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Fit_T2_T2star('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Fit_T2_T2star(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Fit_T2_T2star';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Smoothing'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Smoothing('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Smoothing(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Smoothing';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Shift images'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Shift_image('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Shift_image(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Shift_image';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;  
%         ismodule = 1;
%     case '   .Dynamic Susceptibility Contrast'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Susceptibility('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Susceptibility(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Susceptibility';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .T1map (Multi Angles)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T1map_MultiAngles('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_T1map_MultiAngles(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_T1map_MultiAngles';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Arithmetic'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Arithmetic('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Arithmetic(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Arithmetic';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Normalization'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Normalization('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Normalization(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Normalization';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%      case '   .Clip Image'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Clipping('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Clipping(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Clipping';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;    
%     case '   .FLIRT-FMRIB Linear Image Registration Tool (from FSL)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_FSL_FLIRT('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_FSL_FLIRT(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_FSL_FLIRT';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%         
%     case '   .Brain Extraction (BET Function from FSL)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_FSL_BET('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_FSL_BET(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_FSL_BET';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%      case '   .Brain Mask (using PCNN3D function)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Brain_Mask_PCNN3D('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Brain_Mask_PCNN3D(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Brain_Mask_PCNN3D';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;   
%         
%            
%     case '   .Inversion Efficiency (ASL_InvEff)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_ASL_InvEff('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_ASL_InvEff(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_ASL_InvEff';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .MGE2Dfrom3D'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_MGE2Dfrom3D('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_MGE2Dfrom3D(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_MGE2Dfrom3D';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .deltaR2'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_DeltaR2('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_DeltaR2(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_DeltaR2';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .deltaR2*'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_DeltaR2Star('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_DeltaR2Star(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_DeltaR2Star';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Dynamic Contrast Enhancement (Phenomenology)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_DCE_phenomeno('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_DCE_phenomeno(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_DCE_phenomeno';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .R2prim'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_R2Prim('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_R2Prim(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_R2Prim';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Blood volume fraction (steady-state)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_BVf('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_BVf(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_BVf';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Vessel Size Imaging (steady-state)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_VSI('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_VSI(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_VSI';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Vessel Density (steady-state)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Density('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Density(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Density';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .SO2map'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_SO2('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_SO2(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_SO2';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .CMRO2'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_CMRO2('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_CMRO2(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_CMRO2';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .T1map (Multi Inversion Time)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T1map_MIT('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_T1map_MIT(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_T1map_MIT';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Cerebral blood flow (ASL)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_CBF('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_CBF(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_CBF';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%      case '   .Texture Matlab'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Texture_matlab('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Texture_matlab(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Texture_matlab';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;     
%     case '   .Clustering GMM'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_ClusteringGMM('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_ClusteringGMM(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_ClusteringGMM';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Bias Estimation (MICO algorithm)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_MICO('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_MICO(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_MICO';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Reshape (Extraction)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Reshape('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Reshape(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Reshape';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .SPM: Realign'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Realign_Est('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Realign_Est(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Realign_Est';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case '   .Export data for deeplearing'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Export_data4DL('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_Export_data4DL(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_Export_data4DL';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     case 'Test'
%         Mod = 'Module_Test.m';
%         eval(['[handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = ', Mod, '('',  '', '')']);
%         handles.new_module.command = ['[files_in,files_out,opt] = ', Mod, '(char(files_in),files_out,opt)'];
%         handles.new_module.module_name = Mod;
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
%     otherwise
%         module_parameters_string = 'Not Implemented yet!!';    
%         set(handles.MP3_pipeline_parameter_setup_text, 'String', '');
% 
% end

Mod = char(handles.Modules_listing(module_selected));
    

if ~endsWith(Mod, '.m')
    module_parameters_string = '';
    set(handles.MP3_pipeline_parameter_setup_text, 'String', 'Not a .m file !');
elseif contains(Mod, '(') || contains(Mod, ')')
    module_parameters_string = '';
    set(handles.MP3_pipeline_parameter_setup_text, 'String', 'Please don''t use {''('', '')''} characters in the module file name.');
else
    if startsWith(Mod, '      ..')
        %Remove '   .'
        Mod = Mod(9:end);
    end
    %Remove '.m'
    Mod = Mod(1:end-2);
    
    eval(['[handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = ', Mod, '('''',  '''', '''');']);
    handles.new_module.command = ['[files_in,files_out,opt] = ', Mod, '(char(files_in),files_out,opt)'];
    handles.new_module.module_name = Mod;
    module_parameters_string = handles.new_module.opt.table.Names_Display;
    module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
    ismodule = 1;
end


set(handles.MP3_pipeline_module_parameters, 'String', char(module_parameters_string));
if ismodule
    MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)
    %TableSize = handles.MP3_pipeline_module_parameters.Position(3);
    %FigureSizeInPix = handles.MP3_pipeline_manager_GUI.Position(3);
    %TableSizeInPix = TableSize*FigureSizeInPix;
    %DotsPerInch = get(0,'ScreenPixelsPerInch');
    %DotsPerChar = handles.MP3_pipeline_module_parameters.FontSize;
    %TableSizeInChar = TableSizeInPix/DotsPerChar;
    handles.module_parameters_string = module_parameters_string;
    handles.module_parameters_fields = module_parameters_fields;
    
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
else
    table.data = '';
    table.columnName = '';
    table.editable = false;
    
    %% update the setup table
    % set names of the columns
    set(handles.MP3_pipeline_parameter_setup_table, 'ColumnName', table.columnName);
    % set data (default's parameters)
    set(handles.MP3_pipeline_parameter_setup_table, 'Data', table.data);
    % set each colomn editable
    set(handles.MP3_pipeline_parameter_setup_table, 'columnEditable',  table.editable );
end
   

%% save the data
%guidata(hObject, handles);
guidata(findobj('Tag', 'MP3_pipeline_manager_GUI'), handles);


function [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles)
ActualValues = cell(size(handles.module_parameters_string));
StrToDisplay = cell(size(handles.module_parameters_string));
LString = zeros(size(handles.module_parameters_string));
for i=1:length(handles.module_parameters_fields)
    if strcmp(handles.module_parameters_fields{i}, '') && any(any(strcmp(handles.new_module.opt.table.Default{i}, '')))
        ActualValues{i} = ' ';
    else
        if ~strcmp(handles.module_parameters_fields{i}, '') && isnumeric(handles.new_module.opt.Module_settings.(handles.module_parameters_fields{i}))
            ActualValues{i} = num2str(handles.new_module.opt.Module_settings.(handles.module_parameters_fields{i}));
        elseif any(contains(handles.new_module.opt.table.Type{i}, 'check'))
            ActualValues{i} = char('');
        elseif any(contains(handles.new_module.opt.table.Type{i}, 'Scan')) || any(contains(handles.new_module.opt.table.Type{i}, 'ROI')) || any(contains(handles.new_module.opt.table.Type{i}, 'Cluster'))
            if isempty(handles.new_module.opt.table.Default{i})
                Scan = [];
            else
                Val = handles.new_module.opt.table.Default{i}(:,2);
                Val = Val(~cellfun('isempty', Val'));
                Scan = handles.new_module.opt.table.Default{i}(cell2mat(Val),1);
                %Scan = handles.new_module.opt.table.Default{i}(cell2mat(handles.new_module.opt.table.Default{i}(:,2)),1);
            end
            if isempty(Scan)
                ActualValues{i} = ' ';
            elseif length(Scan)>1
                ActualValues{i} = '';
                for j=1:length(Scan)
                    ActualValues(i) = {[ActualValues{i},'  ', Scan{j}]};
                end
            else
                ActualValues(i) = Scan;
            end
        else
            ActualValues{i} = char(handles.new_module.opt.Module_settings.(handles.module_parameters_fields{i}));
        end
    end
    LString(i) = length(handles.module_parameters_string{i})+length(ActualValues{i});
end
%FinalLength = floor(TableSizeInChar);
FinalLength = 90;
for i=1:length(handles.module_parameters_fields)
    StrToDisplay{i} = [handles.module_parameters_string{i}, repmat(' ',1,FinalLength-LString(i)), ActualValues{i}];
end


%set(handles.MP3_pipeline_module_parameters, 'String', char(module_parameters_string));
set(handles.MP3_pipeline_module_parameters, 'String', char(StrToDisplay));



% --- Executes during object creation, after setting all properties.
function MP3_pipeline_module_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function [pipeline, output_database] = MP3_pipeline_generate_psom_modules(New_module, FilterParameters, TmpDatabase, MP3_path, warning)
% hObject    handle to MP3_pipeline_execute_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pipeline = struct();
Types = New_module.opt.table.Type;
ScanInputs = find(contains(Types, {'Scan', 'ROI', '1Cluster'}));
NbScanInput = length(ScanInputs);
%%% NOTE : Il serait peut-être plus judicieux de boucler sur tous les temps
%%% et patients de la databse de sortie du filtre grosses mailles, ce qui
%%% donnerait 3 matrices de tailles égales. Néanmoins, on perd un peu
%%% l'intérêt de ces histoires de matrices... Le soucis actuel c'est qu'on
%%% se retrouve avec des tailles de matrices incompatibles pour les 3
%%% inputs. Réfléchir là dessus. Edit : Trouvé ! :)

% for i=1:height(handles.new_module.opt.table)
%     if ~isempty(handles.new_module.opt.table.PSOM_Fields{i})
%         Value = handles.new_module.opt.Module_settings.(handles.new_module.opt.table.PSOM_Fields{i});
%         Type = handles.new_module.opt.table.Type{i};
%         if ~isa(Value, Type) && ~strcmp(Type,'cell')
%             text = [handles.new_module.opt.table.Names_Display{i}, ' value is not a ', Type, '. Please select a correct value.'];
%             warndlg(text, 'Wrong type of parameters');
%             pipeline = struct();
%             output_database = table();
%             return
%         end
%     end
% end



for i=1:length(FilterParameters)
    Tag = FilterParameters{1,i}{1};
    TagTable = table();
    for j=2:length(FilterParameters{1,i})
        SelectedValue = FilterParameters{1,i}{j};
        TagTable = [TagTable;TmpDatabase(TmpDatabase{:,Tag}==SelectedValue,:)];
        TagTable = unique(TagTable);
    end
    TmpDatabase = TagTable;
end




%% Build Tp*Patients Matrixes filled up with each inputs files.
DatabaseInput = cell(NbScanInput, 1);
MatricesInputs = cell(NbScanInput,1);
Tag1 = 'Tp';
Tag2 = 'Patient';
EmptyParams = cell(NbScanInput, 1);
for i=1:NbScanInput
    EmptyParams{i} = 0;
    Input = New_module.opt.table.Default{ScanInputs(i)};
    
    %% Other option
%     NbParameters = size(Input,2);
%     Datab = TmpDatabase;
%     Databtmp = table();
%     if NbParameters == 0
%         Datab = table();
%     end
%     for j=1:NbParameters
%         ParamsSelected = Input(2,j);
%         if ~isempty(ParamsSelected)
%             Databtmp = table();
%             for k = 1:length(ParamsSelected)
%                 if j==1
%                     SelectionScans{j,k,i} = ParamsSelected{k};
%                 end
%                 Selection{j,k,i} = ParamsSelected{k};
%                 Tag = New_module.opt.table.Scans_Input_DOF{ScanInputs(i)}{j};
%                 Databtmp = unique([Databtmp ; Datab(getfield(Datab,Tag)==ParamsSelected{k},:)]);
%             end
%             Datab = Databtmp;
%         else
%             EmptyParams{i} = EmptyParams{i} +1;
%         end
%     end
%% Initial Option
    NbParameters = size(Input,2)/2;
    Datab = TmpDatabase;
    Databtmp = table();
    if NbParameters == 0
        Datab = table();
    end
    for j=1:NbParameters
        A = Input(:,2*j);
        A = A(~cellfun('isempty',A));
        ParamsSelected = Input(cell2mat(A),2*j-1);
        if isempty(ParamsSelected) || NbParameters ~= length(New_module.opt.table.Scans_Input_DOF{ScanInputs(i)})
            EmptyParams{i} = 1;
        end
        if ~isempty(ParamsSelected)
            Databtmp = table();
            for k = 1:length(ParamsSelected)
                if j==1
                    SelectionScans{j,k,i} = ParamsSelected{k};
                end
                Selection{j,k,i} = ParamsSelected{k};
                Tag = New_module.opt.table.Scans_Input_DOF{ScanInputs(i)}{j};
                Databtmp = unique([Databtmp ; Datab(getfield(Datab,Tag)==ParamsSelected{k},:)]);
            end
            Datab = Databtmp;
        else
            %EmptyParams{i} = EmptyParams{i} +1;
            Datab = table();
        end
    end

    if ~exist('Databtmp', 'var')
        output_database = table();
        pipeline = struct();
        text = ['Please select at least one parameter (Input ',num2str(i),').'];
        warndlg(text)
        %set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
        return
    end
    DatabaseInput{i} = Databtmp;
    if ~isempty(Datab)
        UTag2 = unique(getfield(Datab, Tag2));
        UTag1 = unique(getfield(Datab, Tag1));
        Mat = cell(length(UTag2), length(UTag1));
        for m=1:length(UTag2)
            Datab2 = Datab(getfield(Datab, Tag2) == UTag2(m),:);
            for n=1:length(UTag1)
                Datab3 = Datab2(getfield(Datab2,Tag1) == UTag1(n),:);
                % Check if this database contains some entries which shares
                % all the fields but the Path one (One might be in the Tmp
                % folder). It means that the automatic generation doesn't
                % know which file take. The one already existing or the
                % virtual (temporary) one ? We currently don't deal with
                % this issue.
                if size(Datab3, 1) > 1
                    FNames = fieldnames(Datab3);
                    RestrictedDatab = Datab3(:,~strcmp(FNames(1:end-3), 'Path'));
                    if size(RestrictedDatab,1) ~= size(unique(RestrictedDatab), 1)
                        %UPDATE : We select the virtual one.
%                         output_database = table();
%                         pipeline = struct();
%                         text = 'It seems that you are trying to create a module which takes as input a scan that exists both virtually and concretely. It it not clear which one you want to use. Please delete the concretely one from your database, filter your database, or modify the module that creates the virtualy one.';
%                         warndlg(text)
%                         return
                        Path = Datab3.Path(1);
                        Path = strsplit(char(Path), filesep);
                        Path{end-1} = 'Tmp';
                        Path = strjoin(Path, filesep);
                        Datab3 = Datab3(Datab3.Path == Path,:);
                    end
                end
                for o=1:size(Datab3,1)
                    Mat{m,n,o} = [char(Datab3.Path(o)) char(Datab3.Filename(o)) '.nii'];
                end
            end
        end
    else
        Mat = {};
    end
    MatricesInputs{i} = Mat;
end


% On which input matrix size must we create the final matrixes ?
%RefInput = 2;


for i=1:length(ScanInputs)
    if strcmp(New_module.opt.table.IsInputMandatoryOrOptional{ScanInputs(i)}, 'Mandatory') && isempty(MatricesInputs{i})
        if warning
            warndlg('A mandatory scan is missing.', 'Missing Scan');
        end
        pipeline = struct();
        output_database = table();
        %set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
        return
    end
end

if ~isfield(New_module.opt.Module_settings, 'AutomaticJobsCreation')  || ...
        isfield(New_module.opt.Module_settings, 'AutomaticJobsCreation') && strcmp(New_module.opt.Module_settings.AutomaticJobsCreation, 'Yes')
    
    RefInput = New_module.opt.Module_settings.RefInput;
    RefDatab = DatabaseInput{RefInput};
    RefMat = MatricesInputs{RefInput};
    
    % Pour creer le tableau Tag1*Tag2 sur l'union des 2 databases plutot
    % que sur la database de reference :
%     In_to_reshape = New_module.opt.Module_settings.InputToReshape; % added
%     Datab_Reshap = DatabaseInput{In_to_reshape}; % added
%     UTag1 = unique(getfield([Datab_Reshap;RefDatab], Tag1)); % added
%     UTag2 = unique(getfield([Datab_Reshap;RefDatab], Tag2)); % added
    
    UTag1 = unique(getfield(RefDatab, Tag1)); 
    UTag2 = unique(getfield(RefDatab, Tag2));
    
    
    FinalMat = cell(NbScanInput,1);
    FinalMat{RefInput} = RefMat;
    
    for i=1:NbScanInput
        if i~=RefInput
            Mattmp = cell(length(UTag2), length(UTag1));
            for j=1:length(UTag2)
                for k=1:length(UTag1)
                    Databtmp = DatabaseInput{i};
                    if ~isempty(Databtmp)
                        Databtmp2 = Databtmp(getfield(Databtmp, Tag2) == UTag2(j), :);
                        Databtmp3 = Databtmp2(getfield(Databtmp2, Tag1) == UTag1(k), :);
                    if size(Databtmp3, 1) > 1
                        FNames = fieldnames(Databtmp3);
                        RestrictedDatab = Databtmp3(:,~strcmp(FNames(1:end-3), 'Path'));
                        if size(RestrictedDatab,1) ~= size(unique(RestrictedDatab), 1)
                            %UPDATE : We select the virtual one.
        %                         output_database = table();
        %                         pipeline = struct();
        %                         text = 'It seems that you are trying to create a module which takes as input a scan that exists both virtually and concretely. It it not clear which one you want to use. Please delete the concretely one from your database, filter your database, or modify the module that creates the virtualy one.';
        %                         warndlg(text)
        %                         return
                            dat = Databtmp3;
                            Path = Databtmp3.Path(1);
                            Path = strsplit(char(Path), filesep);
                            Path{end-1} = 'Tmp';
                            Path = strjoin(Path, filesep);
                            Databtmp3 = Databtmp3(Databtmp3.Path == Path,:);
                            %Now that we only selected the virtual scans,
                            %lets select the real scans that dont have a
                            %similar virtual scan.
                            [~, ind_diff] = setdiff(RestrictedDatab, Databtmp3(:,~strcmp(FNames(1:end-3), 'Path')));
                            Databtmp3 = [Databtmp3; dat(ind_diff,:)];
                        end
                    end
                        for l=1:size(Databtmp3,1)
                            Mattmp{j,k,l} = [char(Databtmp3.Path(l)) char(Databtmp3.Filename(l)) '.nii'];
                        end
                    else
                        Mattmp{j,k} = '';
                    end
                end
            end
            FinalMat{i} = Mattmp;
        end
    end
    
    %RefMat(strcmp('',RefMat)) = []
    %RefMat(:,~all(cellfun('isempty', RefMat(:, :, :)), 1))
    
    % A(:,find(all(cellfun(@isempty,A),1)) & find(all(cellfun(@isempty,A),3)),:) = [];
    % A(find(all(cellfun(@isempty,A),2)),:,:) = [];
    % A(:,:,find(all(cellfun(@isempty,A),2))) = [];
    
    
    
    
    % Which input must we adapt ?
    InputToReshape = New_module.opt.Module_settings.InputToReshape;
    if InputToReshape ~= RefInput
        InToReshape = FinalMat{InputToReshape};
        %InToReshape = InToReshape(~cellfun('isempty',InToReshape));
        InToReshape(:,find(all(cellfun(@isempty,InToReshape),1))) = [];
        
        % Ligne suivante commentee pour debeuger le pipeline de Ludovic . .. . .  YO
        % NO SE :/
        % Cette ligne commentee permet l'execution du pipeline de Ludovic mais
        % empeche du coup celle d'un module Coreg dont l'entrée 1 est
        % un fichier unique (1 Scan + 1 TP + 1 Patient selectionné)... A
        % investiguer. 
        
        %Commenter la ligne suivante permet d'executer des pipelines
        %lorsque l'input de reference manque dans l'un des TP alors que
        %l'input to reshape y est presente. Ex : aucun scan de ref à D0
        %pour tous les patients, alors que le scan a reshape en possede un.
        %InToReshape(find(all(cellfun(@isempty,InToReshape),2)),:) = [];
        
        %Pour régler ce probleme, on sépare les deux cas :
        InToReshape_bis = InToReshape;
        InToReshape_bis(find(all(cellfun(@isempty,InToReshape_bis),2)),:) = [];
        if size(InToReshape_bis) == [1,1]
            InToReshape = InToReshape_bis;
        end
        
        %if size(InToReshape,1) == 1 && size(InToReshape,2) == 1 && EmptyParams{InputToReshape} == 0
        if EmptyParams{InputToReshape} == 0
            NewIn = repmat(InToReshape, size(RefMat));
            assert(size(NewIn,1) == size(RefMat,1));
            assert(size(NewIn,2) == size(RefMat,2));
            assert(size(NewIn,3) == size(RefMat,3));
            FinalMat{InputToReshape} = NewIn;
        elseif size(InToReshape,1) == 1 && size(InToReshape,2) == 1 && EmptyParams{InputToReshape} == 1
%             string = ['Your ' Tag1 ' selection or your ' Tag2 ' selection is empty, but this multiple Tag is useless because it leads to a unique file. So, all the files from the input 2 will be coregistered on this unique file.'];
%             choice = questdlg(string, 'Unique file detected','Continue', 'Return', 'Return');
%             switch choice
%                 case 'Continue'
                    NewIn = repmat(InToReshape, size(RefMat));
                    assert(size(NewIn,1) == size(RefMat,1));
                    assert(size(NewIn,2) == size(RefMat,2));
                    assert(size(NewIn,3) == size(RefMat,3));
                    FinalMat{InputToReshape} = NewIn;
%                 case 'Return'
%                     output_database = table();
%                     pipeline = struct();
%                     %set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
%                     return
%             end
        elseif size(InToReshape,1) == 1 && EmptyParams{InputToReshape} == 1
            NewIn = repmat(InToReshape, size(RefMat,1), 1);
            assert(size(NewIn,1) == size(RefMat,1));
            assert(size(NewIn,2) == size(RefMat,2));
            assert(size(NewIn,3) == size(RefMat,3));
            FinalMat{InputToReshape} = NewIn;
        elseif size(InToReshape,2) == 1 && EmptyParams{InputToReshape} == 1
            NewIn = repmat(InToReshape, 1, size(RefMat,2));
            assert(size(NewIn,1) == size(RefMat,1));
            assert(size(NewIn,2) == size(RefMat,2));
            assert(size(NewIn,3) == size(RefMat,3));
            FinalMat{InputToReshape} = NewIn;
        end
    end
    
    NbModules = size(RefMat,1)*size(RefMat,2);
    for i=1:NbScanInput
        FinalMat{i} = reshape(FinalMat{i}, [NbModules,size(FinalMat{i},3)]);
    end
    
    
    
    output_database = table();
    % add module to a pipeline if feasible (all conditions checked)
    % and create output_database
    for i=1:NbModules
        table_in = table();
        Files_in = struct();
        % For the module i, check if all input contains the number of files as
        % the  number of scan selected by the user. if these numbers are
        % differents and the condered input is mantatory then we do not create
        % the module. If not, the module muste be created and added to the
        % pipeline
        B = zeros(1,NbScanInput);
        All_Selected_Scans = [];
        % here we check conditions input by input
        for l=1:NbScanInput
            AllScans = 1;
            if strcmp(New_module.opt.table.IsInputMandatoryOrOptional{ScanInputs(l)}, 'Mandatory')
                A = {FinalMat{l}{i,:}};
                if length({FinalMat{l}{i,:}}) ~= sum(~cellfun('isempty',SelectionScans(:,:,l))) %length(SelectionScans(:,:,l))
                    AllScans = 0;
                end
                B(l) = any(cellfun('isempty', A));
            end
            All_Selected_Scans = [All_Selected_Scans, AllScans];
        end
        if ~any(B) && all(All_Selected_Scans)
            %if ~isempty(FinalMat{InputToReshape}{i}) && ~isempty(FinalMat{RefInput}{i})
            for j=1:NbScanInput
                %A = {FinalMat{j}{i,:}};
                %B = cellfun('isempty', A);
                %if all(B)
                for k=1:size(FinalMat{j},2)
                    %                 if strcmp(handles.new_module.opt.table.IsInputMandatoryOrOptional{ScanInputs(j)}, 'Invariant') && isempty(FinalMat{j}{i,k})
                    %
                    %                 end
                    %                 if isempty(FinalMat{j}{i,k})
                    %                     FinalMat{j}{i,k} = '';
                    %                 end
                    %                  eval(['Files_in.In' num2str(j) '{' num2str(k) '} = FinalMat{' num2str(j) '}{' num2str(i) ',' num2str(k) '};']);
                    if ~isempty(FinalMat{j}{i,k})
                        eval(['Files_in.In' num2str(j) '{' num2str(k) '} = FinalMat{' num2str(j) '}{' num2str(i) ',' num2str(k) '};']);
                        [PATHSTR,NAME,~] = fileparts(FinalMat{j}{i,k});
                        databtmp = TmpDatabase(TmpDatabase.Filename == categorical(cellstr(NAME)),:);
                        databtmp = databtmp(databtmp.Path == categorical(cellstr([PATHSTR, filesep])),:);
                        % When dealing with ROIs, databtmp sometimes
                        % contains 2 times the ROI entry. I haven't fixed
                        % this issue so far, but it seems to have no
                        % importance as commenting this line don't reveal
                        % any unwanted issue.
                        %assert(size(databtmp, 1) == 1);
                        InTags = databtmp(1,:);
                        table_in = [table_in; InTags];
                    end
                end
            end
            New_module.opt.Module_settings.folder_out = [MP3_path, 'Tmp'];
            New_module.opt.Module_settings.Table_in = unique(table_in, 'stable');
            pipeline = psom_add_job(pipeline, [New_module.module_name, num2str(i)], New_module.module_name, Files_in, '', New_module.opt.Module_settings);
            Mod_Struct = getfield(pipeline, [New_module.module_name, num2str(i)]);
            output_database = [output_database; Mod_Struct.opt.Table_out];
        end
    end
    
else
    table_in = table();
    for i=1:length(DatabaseInput)
        Table = DatabaseInput{i};
        Files = cell(size(Table,1),1);
        for j=1:size(Table,1)
            Files{j} = [char(Table.Path(j)), char(Table.Filename(j)), '.nii'];
            table_in = [table_in ; DatabaseInput{i}];
        end
        eval(['Files_in.In' num2str(i) ' = Files;']);
    end
    New_module.opt.Module_settings.folder_out = [MP3_path, 'Tmp'];
	New_module.opt.Module_settings.Table_in = unique(table_in);
    NameMod = [New_module.module_name, num2str(1)];
	pipeline = psom_add_job(pipeline, NameMod, New_module.module_name, Files_in, '', New_module.opt.Module_settings);
    output_database = pipeline.(NameMod).opt.Table_out;
end


output_database = unique(output_database);

  

function MP3_pipeline_execute_button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_execute_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[hObject, eventdata, handles] = Update_MP3_database_After_pipeline_Crash(hObject, eventdata, handles);
if handles.FlagExecutePipe % If the function above did something to the database, dont execute the rest of the function (pipeline execution).
    return
end

% for now we clean all PSOM files  (old execution) in order to decrease
% potential bugs
if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'PSOM'], 'dir') == 7
    rmdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'PSOM'], 's')
end

if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'PSOM'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'PSOM']);
    if status == false
        error('Cannot create the PSOM folder to save the pipeline logs.')
    end
end
opt_pipe.path_logs = [handles.MP3_data.database.Properties.UserData.MP3_data_path,  'PSOM'];
% 
% % set the number of workers as function of the capacity of the machine
% myCluster = parcluster('local');
% opt_pipe.max_queued = myCluster.NumWorkers
% %opt_pipe.mode = 'session';

if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Derived_data'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Derived_data']);
    if status == false
        error('Cannot create the Derived_Data folder to save the results of the computed maps.')
    end
end

if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'ROI_data'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'ROI_data']);
    if status == false
        error('Cannot create the ROI_Data folder to save the results of the computed maps.')
    end
end

if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp']);
    if status == false
        error('Cannot create the Tmp folder to temporarily save the results of the computed maps.')
    end
else
    %  Here we clean all tmp files
    rmdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp'], 's');
    [status, ~, ~] = mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp']);
    if status == false
        error('Cannot create the Tmp folder to temporarily save the results of the computed maps.')
    end
end


%% Create the pipeline from the modules.
if ~isfield(handles, 'MP3_pipeline_ParamsModules') || isempty(fieldnames(handles.MP3_pipeline_ParamsModules))
    warndlg('No pipeline to execute ...')
    return
end
Table_JobDeleted = table();
Names_Mod = fieldnames(handles.MP3_pipeline_ParamsModules);
Pipeline = struct();
for i=1:length(Names_Mod)
    Mod = handles.MP3_pipeline_ParamsModules.(Names_Mod{i});
    ReWritting = CheckReWriting(Mod, handles.MP3_pipeline_TmpDatabase);
    if any(ReWritting) && ~exist('DeleteRewrite','var')
        text = {'WARNING : Some jobs will overwrite existing files. Do you want to :'...
            '     - Remove these jobs or,'...
            '     - Overwrite the existing files ?'};
        answer = questdlg(text, 'Overwritting', 'Remove jobs', 'Overwrite files', 'Cancel', 'Cancel');
        if isempty(answer) || strcmp(answer, 'Cancel')
            return
        elseif strcmp(answer, 'Overwrite files')
            DeleteRewrite = 0;
        elseif strcmp(answer, 'Remove jobs')
            DeleteRewrite = 1;
        end
    end
%<<<<<<< HEAD
    
    JobNames = fieldnames(Mod.Jobs);
    
    for j=1:length(JobNames) % JobNames has the same size as ReWritting
        JobN = JobNames{j};
        if ReWritting(j) && exist('DeleteRewrite','var') && DeleteRewrite
            % Store the Table_out of all jobs we choosed to delete
            Table_JobDeleted = [Table_JobDeleted; Mod.Jobs.(JobN).opt.Table_out];
            Mod.Jobs = rmfield(Mod.Jobs, JobN);
% =======  % Code de Benjamin pour copier les fichiers du dossier Tmp vers
% % celui des données
%     if exist('DeleteRewrite','var') && DeleteRewrite
%         JobNames = fieldnames(Mod.Jobs);
%         for j=1:length(ReWritting)
%             if ReWritting(j)
%                 JobToDelete = JobNames{j};
%                 % copy files to the tmp folder in case another job needs it
%                 % as input
%                 for x = 1: numel( Mod.Jobs.(JobNames{j}).opt.Table_out.SequenceName)
%                     indexes_of_scan_to_copy = find((handles.MP3_data.database.SequenceName == Mod.Jobs.(JobNames{j}).opt.Table_out.SequenceName(x)) & ...
%                         (handles.MP3_data.database.Patient == Mod.Jobs.(JobNames{j}).opt.Table_out.Patient(x)) & ...
%                         (handles.MP3_data.database.Tp == Mod.Jobs.(JobNames{j}).opt.Table_out.Tp(x)))    ;
%                     for y=1:numel(indexes_of_scan_to_copy)
%                         filename_input = strcat(char(handles.MP3_data.database.Path(indexes_of_scan_to_copy(y))),char(handles.MP3_data.database.Filename(indexes_of_scan_to_copy(y))), '.nii');
%                         filename_ouput = strcat(char(handles.MP3_data.database.Properties.UserData.MP3_data_path), 'Tmp', filesep, char(handles.MP3_data.database.Filename(indexes_of_scan_to_copy(y))), '.nii');
%                         copyfile(filename_input, filename_ouput, 'f');
%                         if exist(strrep(filename_input, '.nii', '.json'), 'file') == 2
%                             copyfile(strrep(filename_input, '.nii', '.json'), strrep(filename_ouput, '.nii', '.json'), 'f');
%                         end
%                     end
%                 end
%                 
%                
% 
%                 % then delete the job
%                 Mod.Jobs = rmfield(Mod.Jobs, JobToDelete);
%             end
% >>>>>>> f597f19ac59ff3855e3365cfad49ea101c992f56
        end
    end
    
    %[new_pipeline, output_database] = MP3_pipeline_generate_psom_modules(Mod.Module, Mod.Filters, handles.MP3_pipeline_TmpDatabase, handles.MP3_data.database.Properties.UserData.MP3_data_path);
    Pipeline = smart_pipeline_merge(Pipeline, Mod.Jobs);
    %handles.MP3_pipeline_TmpDatabase = unique([handles.MP3_pipeline_TmpDatabase; Mod.OutputDatabase]);
end



%The files_out of the deleted jobs can be files_in for other jobs. So, if
%the files_out of the deleted jobs exist in the database, we will change
%the path of the relevant files_in from the Tmp folder to the Derived_data,
%Raw_data or ROI_data.

PipelineJobNames = fieldnames(Pipeline);

if isempty(PipelineJobNames)
    warndlg('No pipeline to execute ...')
    return
end


if ~isempty(Table_JobDeleted) % ie if some jobs are deleted
    for i=1:length(PipelineJobNames)
        Mod = Pipeline.(PipelineJobNames{i});
        [Dep, Ia, ~] = intersect(Mod.opt.Table_in, Table_JobDeleted);
        if ~isempty(Dep)
            for l = 1:size(Dep,1)
                Spl = strsplit(char(Mod.opt.Table_in.Path(Ia(l))), filesep);
                if Mod.opt.Table_in.Type(Ia(l)) == 'ROI' || Mod.opt.Table_in.Type(Ia(l)) == 'Cluster'
                    Spl{end-1} = 'ROI_data';
                elseif Mod.opt.Table_in.Type(Ia(l)) == 'Scan' &&  Mod.opt.Table_in.IsRaw(Ia(l)) == '0'
                    Spl{end-1} = 'Derived_data';
                elseif Mod.opt.Table_in.Type(Ia(l)) == 'Scan' &&  Mod.opt.Table_in.IsRaw(Ia(l)) == '1'
                    Spl{end-1} = 'Raw_data';
                end
                Mod.opt.Table_in.Path(Ia(l)) = categorical(cellstr(strjoin(Spl, filesep)));
                InpN = fieldnames(Mod.files_in);
                for k=1:length(InpN)
                    files = Mod.files_in.(InpN{k});
                    for m = 1:length(files)
                        if strcmp(files{m}, [char(Dep.Path(l)), char(Dep.Filename(l)), '.nii'])
                            files{m} = [char(Mod.opt.Table_in.Path(Ia(l))), char(Mod.opt.Table_in.Filename(Ia(l))), '.nii'];
                        end
                    end
                    Mod.files_in.(InpN{k}) = files;
                end
            end
        end
        Pipeline.(PipelineJobNames{i}) = Mod;
    end
end



handles.psom.pipeline = Pipeline;

% display the pipeline
if exist('biograph') == 2
    
    [graph_deps,list_jobs,files_in,files_out,files_clean] = psom_build_dependencies(handles.psom.pipeline);
    bg = biograph(graph_deps,list_jobs);
    
    
    % dolayout(bg);
    %% add editable functions to interact with the biograph
    set(bg, 'NodeCallbacks', @(hObject,eventdata)MP3_pipeline('node_callbacks',hObject));
    set(bg, 'EdgeCallbacks', @(hObject,eventdata)MP3_pipeline('edge_callbacks',hObject));
    view(bg) %, which will bring up the display in a different window.
    set(0, 'ShowHiddenHandles', 'on')
    
    handles.psom.biograph_fig = gcf;
    %set(handles.psom.biograph_ob, 'Name', 'MP3 pipeline manager');
    guidata(hObject, handles);
    %return
    
end

%% Check if some existing files will be rewrited
% 
% Name_Jobs = fieldnames(Pipeline);
% Files_out = {};
% for i=1:length(Name_Jobs)
%     Job = Pipeline.(Name_Jobs{i});
%     for j=1:size(Job.opt.Table_out,1)
%         switch char(Job.opt.Table_out.Type(j))
%             case 'Scan'
%                 Folder = '/Derived_data/';
%             case {'ROI', 'Cluster'}
%                 Folder = '/ROI_data/';
%             case 'Deleted'
%                 continue
%         end
%         FileName = strrep([char(Job.opt.Table_out.Path(j)), char(Job.opt.Table_out.Filename(j)), '.nii'], ['/Tmp/', char(Job.opt.Table_out.Filename(j)), '.nii'], [Folder, char(Job.opt.Table_out.Filename(j)), '.nii']); 
%         Files_out = [Files_out; {FileName}];
%     end
% end


% Overwrite_Files = {};
% for i=1:length(Files_out)
%     if exist(Files_out{i},'file') || exist(strrep(Files_out{i},'.nii', '.nii.gz'),'file')
%         Overwrite_Files = [Overwrite_Files; Files_out{i}];
%     end
% end
% 
% if ~isempty(Overwrite_Files)
%     text = 'WARNING : The execution of this pipeline will overwrite one or several existing files ! Continue ?';
%     answer = questdlg(text, 'Overwritting', 'Yes', 'See those files','No', 'No');
%     if strcmp(answer,'See those files')
%         f = figure;
%         t = uitable(f, 'Position', [30,100,500,300],'Data',Overwrite_Files);
%         %btnDelete = uicontrol('Parent', f, 'Position', [100,50,100,50], 'String', 'Delete All', 'Callback', 'rep = ''Yes''; delete(gcf)');
%         %btnCancel = uicontrol('Parent', f, 'Position', [300,50,100,50], 'String', 'Cancel', 'Callback', 'rep = ''No''; delete(gcf)');
%         uiwait(f)
%         answer = questdlg('So, execute the pipeline ?', 'It''s time to choose', 'Yes', 'No', 'No');
%     end
%     if strcmp(answer, 'No')
%         return
%     end
% end


% Store a database with all the output entries of the pipeline
Output_Table = setdiff(handles.MP3_pipeline_TmpDatabase,handles.MP3_data.database);
if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp/Table_out'], 'dir') ~= 7
    mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp/Table_out']);
end
filename_table_out = [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp/Table_out/Table_out.mat'];
save(filename_table_out,'Output_Table')




% %% execute the pipeline

if isnan(str2double(handles.Select_Number_Workers.String))
    myCluster = parcluster('local');
    opt_pipe.max_queued = myCluster.NumWorkers;
else
    opt_pipe.max_queued = str2double(handles.Select_Number_Workers.String);
end

if handles.MP3_pipeline_radiobuttonPSOM.Value
    if exist(fullfile(opt_pipe.path_logs, 'PIPE_history.txt'), 'file') == 2
        delete(fullfile(opt_pipe.path_logs, 'PIPE_history.txt'))
    end
    
    
    % Check if some files need to be deleted. if yes, replace the
    % DeleteFile jobs by the output ob psom_add_clean.
    
    %pipeline2 = psom_add_clean(pipeline,name_job,files_clean)
    Name_Jobs = fieldnames(handles.psom.pipeline);
    pipe_tmp = handles.psom.pipeline;
    for i=1:length(Name_Jobs)
        Job = handles.psom.pipeline.(Name_Jobs{i});
        if strcmp(Job.command, 'Module_DeleteFile(files_in,files_out,opt);')
            pipe_tmp = rmfield(pipe_tmp, Name_Jobs{i});
            pipe_tmp = psom_add_clean(pipe_tmp, Name_Jobs{i}, Job.files_in);
        end
    end
    handles.psom.pipeline = orderfields(pipe_tmp, Name_Jobs);
    
    psom_run_pipeline(handles.psom.pipeline, opt_pipe);
    % We have to wait until PSOM finish to update every files before
    % loading PIPE_statu.mat
    PSOM_finished = 0;
    while PSOM_finished == 0
        if exist([opt_pipe.path_logs, '/PIPE.exit'], 'file') == 2
            PSOM_finished = 1;
        end
        pause(0.2);   
    end
    Result = load([opt_pipe.path_logs, '/PIPE_status.mat']);

else
    Modules = fieldnames(handles.psom.pipeline);
    Result = struct();
    disp(['Your pipeline (',num2str(length(Modules)), ' jobs) will be executed in Basic Loop :']);
    for i=1:length(Modules)
        Module = getfield(handles.psom.pipeline, Modules{i});
        files_in = Module.files_in;
        files_out = Module.files_out;
        opt = Module.opt;
        disp(['Execution Job ', num2str(i)])
        eval(Module.command);
        Result = setfield(Result, Modules{i}, 'finished');
    end
end

Jobs = fieldnames(handles.psom.pipeline);
update = false;
Deleted_files = table();
for i=1:length(Jobs)
   switch  Result.(Jobs{i}) %getfield(Result, Jobs{i})
       case 'failed'
           disp('FAILED')
       case 'finished' 
           update = true;
           J = getfield(handles.psom.pipeline, Jobs{i});
           if isfield(J, 'files_out')
               A = getfield(J, 'files_out');
           else
               A = [];
           end
           %C = getfield(J, 'files_in');
           % If there is no file_out from a module
           if isempty(A) && isfield(J, 'files_clean')
               A.In1 = J.files_clean;
           end
           if isempty(A)
               continue
           end
           Outputs = fieldnames(A);
           for j=1:length(Outputs)
               B = getfield(A, Outputs{j});
               for k=1:length(B)
                   [path_out, name_out, ~] = fileparts(B{k});
                   path_out = [path_out, filesep];
                   outdb = handles.MP3_pipeline_TmpDatabase(handles.MP3_pipeline_TmpDatabase.Path == path_out, :);
                   outdb = outdb(outdb.Filename == name_out, :);
                   if height(outdb) > 1 % It means we have a Scan/ROI/Cluster, and the same entry with the type Delete.
                       assert(height(outdb) == 2)
                       outdb = outdb(outdb.Type ~= 'Deleted', :);
                       assert(height(outdb) == 1)
                       [~,ia,~] = intersect(handles.MP3_data.database, outdb);
                       handles.MP3_data.database(ia, :) = [];
                       if ~isempty(Deleted_files) && ~isempty(intersect(Deleted_files, outdb))
                           % Check if the files we want to deleted has
                           % already been deleted. If yes, don't do
                           % anything.
                           continue
                       end
                       Deleted_files = [Deleted_files; outdb];
                       if ~handles.MP3_pipeline_radiobuttonPSOM.Value % If basic loop used (not PSOM) 
                           filename_to_delete = [char(outdb.Path), filesep, char(outdb.Filename), '.nii.gz'];
                           if ~exist(filename_to_delete, 'file')
                               filename_to_delete = strrep(filename_to_delete, '.nii.gz', '.nii');
                               if ~exist(filename_to_delete, 'file')
                                   error(['File not found : ', filename_to_delete]);
                               end
                           end
                           delete(filename_to_delete);
                       end
                       switch char(outdb.Type)
                           case 'Scan'
                               json_to_delete = [char(outdb.Path), filesep, char(outdb.Filename), '.json'];
                               delete(json_to_delete);
                           case 'Cluster'
                               mat_to_delete = [char(outdb.Path), filesep, char(outdb.Filename), '.mat'];
                               delete(mat_to_delete);
                       end
                       
                   else
% <<<<<<< HEAD
%                        statusJson = 1;
%                        statusMat = 1;
%                    end
%                    if statusNii && statusJson && statusMat
%                        outdb.Path = NewPath;
%                        % If intersect du triplet patient/Tp/sequencename
%                        % entre mp3_data.database et outdb, on faut
%                        % remplacer l'entrée de la database par celle de
%                        % outdb. A FAIRE
%                        
%                        handles.MP3_data.database = unique([handles.MP3_data.database ; outdb]);
%                        %eval(['delete ' B{k}]);
%                    elseif ~statusNii
%                        warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data''/''ROI_data'' folder.', [name_out, '.nii'])
%                    elseif ~statusJson
%                        warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data''/''ROI_data'' folder.', [name_out, '.json'])
%                    elseif ~statusMat
%                        warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data''/''ROI_data'' folder.', [name_out, '.mat'])
% =======
                       assert(height(outdb) == 1)
                       Folders = strsplit(path_out, filesep);
                       if strcmp(char(outdb.Type), 'Deleted')
                           handles.MP3_data.database = unique([handles.MP3_data.database ; outdb]);
                       end
                       switch char(outdb.Type)
                           case 'Scan'
                               Folders{end-1} = 'Derived_data';
                           case 'ROI'
                               Folders{end-1} = 'ROI_data';
                           case 'Cluster'
                               Folders{end-1} = 'ROI_data';
                       end
                       NewPath = strjoin(Folders, filesep);
                       [statusNii,~] = movefile(B{k}, [NewPath, name_out, '.nii']);

                       if strcmp(char(outdb.Type), 'Scan')
                           [statusJson,~] = movefile(strrep(B{k},'.nii','.json'), [NewPath, name_out, '.json']);
                           statusMat = 1;
                       elseif strcmp(char(outdb.Type), 'Cluster')
                           if exist(strrep(B{k},'.nii','.mat'), 'file')
                               [statusMat,~] = movefile(strrep(B{k},'.nii','.mat'), [NewPath, name_out, '.mat']);
                           end
                           statusJson = 1;
                       else
                           statusJson = 1;
                           statusMat = 1;
                       end
                       if statusNii && statusJson && statusMat
                           outdb.Path = NewPath;
                           handles.MP3_data.database = unique([handles.MP3_data.database ; outdb]);
                           %eval(['delete ' B{k}]);
                       elseif ~statusNii
                           warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data''/''ROI_data'' folder.', [name_out, '.nii'])
                       elseif ~statusJson
                           warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data''/''ROI_data'' folder.', [name_out, '.json'])
                       elseif ~statusMat
                           warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data''/''ROI_data'' folder.', [name_out, '.mat'])
                       end
% >>>>>>> ModuleDeleteFiles
                   end
                   
               end
           end
   end
end


% update MP3 database if needed
if update
    handles2 = guidata(handles.MP3_data.MP3_GUI);
    handles2.database = handles.MP3_data.database;
    guidata(handles.MP3_data.MP3_GUI, handles2);

    MP3('MP3_update_database_display', hObject, eventdata,handles.MP3_data)
    MP3('MP3_menu_save_database_Callback', hObject, eventdata,handles.MP3_data)
    %msgbox('Done', 'Information') ;

    handles.database = handles.MP3_data.database;
    handles = rmfield(handles, 'MP3_pipeline_ParamsModules');
    MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles);
    
end

rmdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep, 'Table_out'], 's');






function [hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles)


Datab = handles.MP3_data.database;
if isfield(handles, 'MP3_pipeline_ParamsModules')
    Names_Mod = fieldnames(handles.MP3_pipeline_ParamsModules);
    for i=1:length(Names_Mod)
        Mod = handles.MP3_pipeline_ParamsModules.(Names_Mod{i});
        Datab = [Datab ; Mod.OutputDatabase];
    end
end
%handles.MP3_pipeline_TmpDatabase = unique(Datab);
handles.MP3_pipeline_TmpDatabase = Datab;



% --- Executes on button press in MP3_pipeline_close_modules_button.
function MP3_pipeline_close_modules_button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_close_modules_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.MP3_pipeline_manager_GUI,'Name') '?'],...
                     ['Close ' get(handles.MP3_pipeline_manager_GUI,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.MP3_pipeline_manager_GUI)

function [Add_list, Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, InitialListing)

if isempty(NewTagListing)
    Add_list = {'NoMoreTags'};
    Remove_list = InitialListing;
elseif sum(contains(InitialListing, NewTagListing)) == length(InitialListing)
    Add_list = NewTagListing;
    Remove_list = {'NoMoreTags'};
else
    Add_list = NewTagListing;
    Remove_list = InitialListing(~contains(InitialListing, NewTagListing));
end




% --- Executes on button press in MP3_pipeline_Add_Tag_Button.
function MP3_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Add_Tag_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if strcmp(handles.Source_selected, 'NoMoreTags')
%    return 
% end
%% Update of the Filtering Table
%if isempty(handles.MP3_pipeline_Filtering_Table.Data)
%    handles.MP3_pipeline_Filtering_Table.Data = table2cell(handles.MP3_data.database);
%end
if length(handles.MP3_pipeline_Unique_Values_Selection) <= 1
    return
end

handles.FilterParameters = [handles.FilterParameters, {handles.MP3_pipeline_Unique_Values_Selection}];

[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module') && ~isempty(fieldnames(handles.new_module))
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI', '1Cluster'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end


%% Update of the add tag popupmenu after an add
Tag_To_Add = handles.Source_selected;
Index = find(contains(handles.Add_list, Tag_To_Add));
NewTagListing = {handles.Add_list{1:Index-1},handles.Add_list{Index+1:end}};
[handles.Add_list, handles.Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);

set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);

guidata(hObject, handles);

% --- Executes on button press in MP3_pipeline_Remove_Tag_Button.
function MP3_pipeline_Remove_Tag_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Remove_Tag_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.Remove_selected, 'NoMoreTags')
    return
end
%% Update of the Filtering table
for i=1:length(handles.FilterParameters)
   if  strcmp(handles.FilterParameters{1,i}{1},handles.Remove_selected)
       index = i;
   end
end
handles.FilterParameters = {handles.FilterParameters{1:index-1}, handles.FilterParameters{index+1:end}};
% 

[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end

if isfield(handles, 'new_module') && ~isempty(fieldnames(handles.new_module))
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI', '1Cluster'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end


%% Update of the remove tag popup menu after a remove
Tag_To_Remove = handles.Remove_selected;
Index = find(contains(handles.Remove_list, Tag_To_Remove));
NewTagListing = {handles.Remove_list{1:Index-1},handles.Remove_list{Index+1:end}};

%% Interchange Remove and Add list to apply the process.
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);

set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
%handles.Source_selected = Tag_To_Remove;
handles.Remove_selected = handles.Remove_list{1};
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
set(handles.MP3_pipeline_remove_tag_popupmenu, 'Value', 1);

guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in MP3_pipeline_Filtering_Table.
function MP3_pipeline_Filtering_Table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Filtering_Table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% NbSelected = length(eventdata.Indices);
% for i = 1:NbSelected
%     NameSelected = eventdata.Source.Data{eventdata.Indices(i,1), eventdata.Indices(i,2)};
%     Tag = handles.MP3_pipeline_Filtering_Table.ColumnName{enventdata.Indices(i,2)};
% end

guidata(hObject, handles);


% --- Executes on button press in MP3_pipeline_pushMP3Selection.
function MP3_pipeline_pushMP3Selection_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_pushMP3Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% define which patient is selected
%data_selected =  MP3('get_data_selected',handles.MP3_data);
data_selected =  finddata_selected(handles.MP3_data);
% add the patient filter
handles.FilterParameters = {};
handles.FilterParameters{1} = {'Patient', char(handles.MP3_pipeline_TmpDatabase.Patient(data_selected(1)))};
handles.FilterParameters{2} = {'Tp', char(handles.MP3_pipeline_TmpDatabase.Tp(data_selected(1)))};

[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module') && isfield(handles.new_module, 'opt')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOr1ROI', 'XScanOrROI', '1ROI', 'XROI', '1Cluster', '1ROIOr1Cluster', 'XROIOrXCluster'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end

NewTagListing = {'Patient', 'Tp'};
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
guidata(hObject, handles);









% --- Executes on selection change in MP3_pipeline_remove_tag_popupmenu.
function MP3_pipeline_remove_tag_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_remove_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Remove_selected = handles.MP3_pipeline_remove_tag_popupmenu.String{handles.MP3_pipeline_remove_tag_popupmenu.Value};
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_remove_tag_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_remove_tag_popupmenu



% --- Executes on button press in MP3_pipeline_pushMP3TPSelection.
function MP3_pipeline_pushMP3TPSelection_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_pushMP3TPSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% define which patient is selected
%data_selected =  MP3('get_data_selected',handles.MP3_data);
data_selected =  finddata_selected(handles.MP3_data);

% add the patient filter
handles.FilterParameters = {};
handles.FilterParameters{1} = {'Patient', char(handles.MP3_pipeline_TmpDatabase.Patient(data_selected(1)))};

[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end

%% Delete a potential scan selected when building a module.
if isfield(handles, 'new_module')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOr1ROI', 'XScanOrROI', '1ROI', 'XROI', '1Cluster', '1ROIOr1Cluster', 'XROIOrXCluster'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end

NewTagListing = {'Patient'};
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
guidata(hObject, handles);


% --- Executes on button press in MP3_pipeline_push_Database.
function MP3_pipeline_push_Database_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_push_Database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.FilterParameters = {};

[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module') && isfield(handles.new_module, 'opt')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOr1ROI', 'XScanOrROI', '1ROI', 'XROI', '1Cluster', '1ROIOr1Cluster', 'XROIOrXCluster'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end

NewTagListing = handles.Add_Tags_listing;
[handles.Add_list, handles.Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);

% 
% %% Update of the add tag popupmenu after an add
% Tag_To_Add = handles.Source_selected;
% Index = find(contains(handles.Add_list, Tag_To_Add));
% NewTagListing = {handles.Add_list{1:Index-1},handles.Add_list{Index+1:end}};
% [handles.Add_list, handles.Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
% 
% set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
% set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
% handles.Source_selected = handles.Add_list{1};
% handles.Remove_selected = handles.Remove_list{1};
% %handles.Remove_selected = Tag_To_Add;
% MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
% set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);


% handles.Tags_listing = handles.MP3_pipeline_TmpDatabase.Properties.VariableNames;
% set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Tags_listing);
% set(handles.MP3_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})
% 
% handles.Source_selected = handles.Tags_listing{1};
% 
% 
% handles.Remove_Tags_listing = {'NoMoreTags'};
% handles.Remove_selected = handles.Remove_Tags_listing{1};
% handles.MP3_pipeline_TmpDatabase.IsRaw = categorical(handles.MP3_pipeline_TmpDatabase.IsRaw);
% 
% handles.MP3_pipeline_Filtered_Table = handles.MP3_pipeline_TmpDatabase;
% handles.MP3_pipeline_Filtering_Table.Data = cellstr(handles.MP3_pipeline_Filtered_Table{:,handles.MP3_pipeline_TagsToPrint});
% handles.MP3_pipeline_Filtering_Table.ColumnName = handles.MP3_pipeline_TagsToPrint;
% 
% 
% 
% handles.FilterParameters = {};
% 
% handles.MP3_pipeline_Unique_Values_Selection = {};
% MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function MP3_pipeline_remove_tag_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_remove_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in MP3_pipeline_Unique_Values_Tag.
function MP3_pipeline_Unique_Values_Tag_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Unique_Values_Tag (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
NbSelect = size(eventdata.Indices,1);
Names = cell(size(eventdata.Indices,1)+1,1);
Names{1,1}=handles.Source_selected;
for i=1:NbSelect
    Names{i+1,1} = eventdata.Source.Data{eventdata.Indices(i,1), eventdata.Indices(i,2)};
end
handles.MP3_pipeline_Unique_Values_Selection = Names;
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);


% --- Executes on button press in MP3_pipeline_Tags_To_Display_Button.
function MP3_pipeline_Tags_To_Display_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Tags_To_Display_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PreviousTagsIndex = find(contains(handles.MP3_pipeline_Filtered_Table.Properties.VariableNames, handles.MP3_pipeline_TagsToPrint));
[selection, ok] = listdlg('PromptString', 'Select the Tags to print among : ', 'InitialValue', PreviousTagsIndex(:), 'ListString', handles.MP3_pipeline_Filtered_Table.Properties.VariableNames);
if ok == 0;
   return 
end
handles.MP3_pipeline_TagsToPrint = handles.MP3_pipeline_Filtered_Table.Properties.VariableNames(selection);


% handles.MP3_pipeline_Filtering_Table.Data = cellstr(handles.MP3_pipeline_Filtered_Table{:,handles.MP3_pipeline_TagsToPrint});
% handles.MP3_pipeline_Filtering_Table.ColumnName = handles.MP3_pipeline_TagsToPrint;

% handles.MP3_pipeline_Filtered_Table = NewTable;
% CellsColoured = DisplayColoredTable(NewTable, handles.MP3_pipeline_TagsToPrint);
% handles.MP3_pipeline_Filtering_Table.Data = CellsColoured;

[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);

guidata(hObject, handles);


% --- Executes on selection change in MP3_pipeline_module_listbox.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_module_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_module_listbox


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MP3_pipeline_pipeline_listbox.
function [hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_pipeline_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelectedIndex = handles.MP3_pipeline_pipeline_listbox.Value;
% if isequal(SelectedIndex, 0) || isempty(SelectedIndex)
%     JobNames = {''};
    
if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    JobNames = {''};
    ColoredJobNames = {''};
else
    %SelectedModule = handles.MP3_pipeline_pipeline_listbox_Raw{SelectedIndex};
    SelectedModule = handles.MP3_pipeline_pipeline_listbox.String{SelectedIndex}; 
    revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
    SelectedModule = revertcolor2(SelectedModule);
    Module = handles.MP3_pipeline_ParamsModules.(SelectedModule);
    if isfield(Module, 'Jobs')
        JobNames = fieldnames(Module.Jobs);
        ColoredJobNames = ColorJobs(Module, handles.MP3_pipeline_TmpDatabase);
    else
        ColoredJobNames = {''};
        JobNames = {''};
    end
end
handles.MP3_pipeline_JobsNames = JobNames;
old_job_seleced = get(handles.MP3_pipeline_JobsList, 'Value');
set(handles.MP3_pipeline_JobsList, 'String',ColoredJobNames);
if old_job_seleced ~=1
    set(handles.MP3_pipeline_JobsList, 'Value', old_job_seleced-1 );
end
MP3_pipeline_JobsList_Callback(hObject, eventdata, handles)
guidata(hObject, handles)


function ColoredJobs = ColorJobs(Module, Tmpdatab)
JobsNames = fieldnames(Module.Jobs);
ColoredJobs = cell(size(JobsNames));

ReWritting = CheckReWriting(Module, Tmpdatab);

color2 = @(color,text) ['<HTML><FONT color="',color,'">',text,'</Font></html>'];
for i=1:length(JobsNames)
    if ReWritting(i)
        ColoredJobs{i} = color2('orange', JobsNames{i});
    else
        ColoredJobs{i} = color2('green', JobsNames{i});
    end
end


% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_pipeline_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_pipeline_listbox


% --- Executes during object creation, after setting all properties.
function MP3_pipeline_pipeline_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_pipeline_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in MP3_pipeline_exectute_pipeline_button.
function MP3_pipeline_exectute_pipeline_button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_exectute_pipeline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in MP3_pipeline_DeleteModule.
function MP3_pipeline_DeleteModule_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_DeleteModule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    return
end


SelectedIndex = handles.MP3_pipeline_pipeline_listbox.Value;
%SelectedModule = handles.MP3_pipeline_pipeline_listbox_Raw{SelectedIndex};
revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
SelectedModule = handles.MP3_pipeline_pipeline_listbox.String{SelectedIndex};
SelectedModule = revertcolor2(SelectedModule);
handles.MP3_pipeline_ParamsModules = rmfield(handles.MP3_pipeline_ParamsModules, SelectedModule);
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);

%set(handles.MP3_pipeline_pipeline_listbox,'String', fieldnames(handles.MP3_pipeline_ParamsModules));
if ~isempty(fieldnames(handles.MP3_pipeline_ParamsModules))
    set(handles.MP3_pipeline_pipeline_listbox, 'Value', 1);
    set(handles.MP3_pipeline_JobsList, 'Value', 1);
    set(handles.MP3_pipeline_JobsParametersFieldsList, 'Value', 1);
    set(handles.MP3_pipeline_JobsParametersValues, 'Value', 1);
    %handles.MP3_pipeline_pipeline_listbox_Raw = fieldnames(handles.MP3_pipeline_ParamsModules);
else 
    handles = rmfield(handles, 'MP3_pipeline_ParamsModules');
    set(handles.MP3_pipeline_pipeline_listbox, 'String', {''});
    set(handles.MP3_pipeline_pipeline_listbox, 'Value', 1);
    set(handles.MP3_pipeline_JobsList, 'String', {''});
    set(handles.MP3_pipeline_JobsList, 'Value', 1);
    set(handles.MP3_pipeline_JobsParametersFieldsList, 'String', {''});
    set(handles.MP3_pipeline_JobsParametersFieldsList, 'Value', 1);
    set(handles.MP3_pipeline_JobsParametersValues, 'String', {''});
    set(handles.MP3_pipeline_JobsParametersValues, 'Value', 1);
    %handles.MP3_pipeline_pipeline_listbox_Raw = {''};
end


%handles.MP3_pipeline_pipeline_listbox_Raw = fieldnames(handles.MP3_pipeline_ParamsModules);
%[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
if isfield(handles, 'MP3_pipeline_ParamsModules')
    [hObject, eventdata, handles] = MP3_pipeline_UpdatePipelineJobs(hObject, eventdata, handles);
end
[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);
[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
%Coloredlistbox = DisplayColoredListbox(revertcolor2(handles.MP3_pipeline_pipeline_listbox.String), handles);
%set(handles.MP3_pipeline_pipeline_listbox,'String', Coloredlistbox);
%[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
MP3_pipeline_JobsList_Callback(hObject, eventdata, handles)
%[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
%Coloredlistbox = DisplayColoredListbox(revertcolor2(handles.MP3_pipeline_pipeline_listbox.String), handles);
%set(handles.MP3_pipeline_pipeline_listbox,'String', Coloredlistbox);

guidata(hObject, handles);


% --- Executes on button press in MP3_pipeline_Edit_Module.
function MP3_pipeline_Edit_Module_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Edit_Module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    return
end



SelectedIndex = handles.MP3_pipeline_pipeline_listbox.Value;
%SelectedModule = handles.MP3_pipeline_pipeline_listbox_Raw{SelectedIndex};
revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
SelectedModule = handles.MP3_pipeline_pipeline_listbox.String{SelectedIndex};
SelectedModule = revertcolor2(SelectedModule);
Module = handles.MP3_pipeline_ParamsModules.(SelectedModule);
handles.MP3_pipeline.EditedModuleName = SelectedModule;
handles.BeforeEditedModuleFilters = handles.FilterParameters;

%% Remove the module we want to edit, and adapt the database
handles.MP3_pipeline_ParamsModules = rmfield(handles.MP3_pipeline_ParamsModules, SelectedModule);
%Stored = handles.MP3_pipeline_pipeline_listbox_Raw;
%handles.MP3_pipeline_pipeline_listbox_Raw = fieldnames(handles.MP3_pipeline_ParamsModules);
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
%[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);

if isfield(handles, 'new_module')
    handles.BeforeEditedModule = handles.new_module;
else
    handles.BeforeEditedModule = struct();
end

handles.new_module = Module.ModuleParams;

module_parameters_string = handles.new_module.opt.table.Names_Display;
module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;

%MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)

handles.module_parameters_string = module_parameters_string;
handles.module_parameters_fields = module_parameters_fields;
    
handles.FilterParameters = [handles.FilterParameters, Module.Filters];


%handles.MP3_pipeline_pipeline_listbox_Raw = Stored;
[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);


%% Add/Remove Popmenu update
TagsUsed = {};
for i=1:length(Module.Filters)
    TagsUsed = [TagsUsed, Module.Filters{i}(1)];
end
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(TagsUsed, handles.Add_Tags_listing);

set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);

handles.MP3_pipeline_module_parameters.Value = 1;
MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)




%% Update figure
set(handles.MP3_pipeline_module_parameters, 'BackgroundColor', [0.5 0.5 0.5]);
set(handles.MP3_pipeline_parameter_setup_table, 'BackgroundColor', [0.6 0.6 0.6;0.4 0.4 0.4]);
set(handles.MP3_pipeline_Unique_Values_Tag, 'BackgroundColor', [0.6 0.6 0.6;0.4 0.4 0.4]);
set(handles.MP3_pipeline_Filtering_Table, 'BackgroundColor', [0.6 0.6 0.6;0.4 0.4 0.4]);

set(handles.MP3_pipeline_clear_pipeline_button, 'Enable', 'off');
set(handles.MP3_pipeline_DeleteModule, 'Enable', 'off');
set(handles.MP3_pipeline_add_module_button, 'Enable', 'off');
set(handles.MP3_pipeline_execute_button, 'Enable', 'off');
set(handles.MP3_pipeline_module_listbox, 'Enable', 'off');
set(handles.MP3_pipeline_Edit_Module, 'Enable', 'off');
set(handles.MP3_pipeline_pipeline_listbox, 'Enable', 'off');
set(handles.MP3_pipeline_JobsList, 'Enable', 'off');
set(handles.MP3_pipeline_JobsParametersFieldsList, 'Enable', 'off');
set(handles.MP3_pipeline_JobsParametersValues, 'Enable', 'off');
set(handles.MP3_pipeline_Delete_Job, 'Enable', 'off');
set(handles.MP3_pipeline_save_pipeline, 'Enable', 'off');
set(handles.MP3_pipeline_load_pipeline, 'Enable', 'off');

set(handles.MP3_pipeline_Save_Module, 'Enable', 'on');

guidata(hObject, handles);
%guidata(hObject, handles);


% --- Executes on button press in MP3_pipeline_Save_Module.
function MP3_pipeline_Save_Module_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Save_Module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%[new_pipeline, output_database] = MP3_pipeline_generate_psom_modules(handles.new_module, handles.FilterParameters, handles.MP3_pipeline_TmpDatabase, handles.MP3_data.database.Properties.UserData.MP3_data_path);

%if isempty(fieldnames(new_pipeline)) && isempty(output_database)
%    return
%end



SaveModule = struct();
SaveModule.Filters = handles.FilterParameters;
SaveModule.ModuleParams = handles.new_module;
%SaveModule.OutputDatabase = output_database;
%SaveModule.Jobs = new_pipeline;
handles.MP3_pipeline_ParamsModules.(handles.MP3_pipeline.EditedModuleName) = SaveModule;
%handles.MP3_pipeline_ParamsModules = orderfields(handles.MP3_pipeline_ParamsModules, handles.MP3_pipeline_pipeline_listbox.String);
revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
fields = revertcolor2(handles.MP3_pipeline_pipeline_listbox.String);
handles.MP3_pipeline_ParamsModules = orderfields(handles.MP3_pipeline_ParamsModules, fields);

[hObject, eventdata, handles] = MP3_pipeline_UpdatePipelineJobs(hObject, eventdata, handles);


if isfield(handles, 'new_module')
    handles = rmfield(handles, 'new_module');
end


[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);

handles.FilterParameters = {};
[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);
%[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);

handles.new_module = handles.BeforeEditedModule;
handles.FilterParameters = handles.BeforeEditedModuleFilters;




%% Display the module we were creating before the Edit (Beginning)
if ~isfield(handles.new_module, 'opt')
    module_parameters_string = {};
    module_parameters_fields = {};
else
    module_parameters_string = handles.new_module.opt.table.Names_Display;
    module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
end
%MP3_pipeline_module_parameters_Callback(hObject, eventdata, handles)

handles.module_parameters_string = module_parameters_string;
handles.module_parameters_fields = module_parameters_fields;
    
[hObject, eventdata, handles]=MP3_pipeline_UpdateTables(hObject, eventdata, handles);

[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);


%% Add/Remove Popmenu update
TagsUsed = {};
for i=1:length(handles.FilterParameters)
    TagsUsed = [TagsUsed, handles.FilterParameters{i}(1)];
end
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(TagsUsed, handles.Add_Tags_listing);

set(handles.MP3_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MP3_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MP3_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MP3_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
%% Display the module we were creating before the Edit (End)




%% Update figure
set(handles.MP3_pipeline_module_parameters, 'BackgroundColor', [0.94 0.94 0.94])
set(handles.MP3_pipeline_parameter_setup_table, 'BackgroundColor', [1 1 1;0.9412 0.9412 0.9412]);
set(handles.MP3_pipeline_Unique_Values_Tag, 'BackgroundColor', [1 1 1;0.9412 0.9412 0.9412]);
set(handles.MP3_pipeline_Filtering_Table, 'BackgroundColor', [1 1 1;0.9412 0.9412 0.9412]);
set(handles.MP3_pipeline_clear_pipeline_button, 'Enable', 'on');
set(handles.MP3_pipeline_DeleteModule, 'Enable', 'on');
set(handles.MP3_pipeline_add_module_button, 'Enable', 'on');
set(handles.MP3_pipeline_execute_button, 'Enable', 'on');
set(handles.MP3_pipeline_module_listbox, 'Enable', 'on');
set(handles.MP3_pipeline_Edit_Module, 'Enable', 'on');
set(handles.MP3_pipeline_pipeline_listbox, 'Enable', 'on');
set(handles.MP3_pipeline_JobsList, 'Enable', 'on');
set(handles.MP3_pipeline_JobsParametersFieldsList, 'Enable', 'on');
set(handles.MP3_pipeline_JobsParametersValues, 'Enable', 'on');
set(handles.MP3_pipeline_Delete_Job, 'Enable', 'on');
set(handles.MP3_pipeline_save_pipeline, 'Enable', 'on');
set(handles.MP3_pipeline_load_pipeline, 'Enable', 'on');

set(handles.MP3_pipeline_Save_Module, 'Enable', 'off');
set(handles.MP3_pipeline_pipeline_listbox, 'Value', 1);
[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
guidata(hObject, handles)

% MP3_pipeline_module_listbox_Callback(hObject, eventdata, handles);


% --- Executes on selection change in MP3_pipeline_JobsList.
function MP3_pipeline_JobsList_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_JobsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelectedModuleIndex = handles.MP3_pipeline_pipeline_listbox.Value;
% if isequal(SelectedModuleIndex, 0) || isempty(SelectedModuleIndex)
%     String = {};
if ~isfield(handles, 'MP3_pipeline_ParamsModules')% || ~isfield(handles.MP3_pipeline_ParamsModules.(handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex}), 'Jobs') || isempty(fieldnames(handles.MP3_pipeline_ParamsModules.(handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex}).Jobs))
    String = {''};
else
    %SelectedModule = handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex};
    revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
    SelectedModule = handles.MP3_pipeline_pipeline_listbox.String{SelectedModuleIndex};
    SelectedModule = revertcolor2(SelectedModule);
    Module = handles.MP3_pipeline_ParamsModules.(SelectedModule);
    if ~isfield(Module, 'Jobs') || isempty(fieldnames(Module.Jobs))
        String = {''};
    else
        SelectedJobIndex = handles.MP3_pipeline_JobsList.Value;
        %SelectedJob = handles.MP3_pipeline_JobsList.String{SelectedJobIndex};
        SelectedJob = handles.MP3_pipeline_JobsNames{SelectedJobIndex};
        Job = Module.Jobs.(SelectedJob);
        %Names = fieldnames(Job);
        String = {};
        
        Inputs = fieldnames(Job.files_in);
        for i=1:length(Inputs)
            String = [String; {['files_in', ' ', Inputs{i}]}];
        end
 
        if ~isempty(Job.files_out)
            Outputs = fieldnames(Job.files_out);
            for i=1:length(Outputs)
                String = [String; {['files_out', ' ', Outputs{i}]}];
            end
        end
        UserFields = Module.ModuleParams.opt.table.PSOM_Fields;
        UserFields = UserFields(~cellfun(@isempty,UserFields));
        for i=1:length(UserFields)
            String = [String; {['opt', ' ', UserFields{i}]}];
        end
    end
end
% for i=2:length(Names)
%     Param = Job.(Names{i});
%     NamesEntries = fieldnames(Param);
%     for j=1:length(NamesEntries)
%         String = [String; {[Names{i}, ' ', NamesEntries{j}]}];
%     end
% end
set(handles.MP3_pipeline_JobsParametersFieldsList, 'String', String)
set(handles.MP3_pipeline_JobsParametersFieldsList, 'Value', 1)
MP3_pipeline_JobsParametersFieldsList_Callback(hObject, eventdata, handles)



% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_JobsList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_JobsList


% --- Executes during object creation, after setting all properties.
function MP3_pipeline_JobsList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_JobsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MP3_pipeline_JobsParametersFieldsList.
function MP3_pipeline_JobsParametersFieldsList_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_JobsParametersFieldsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelectedModuleIndex = handles.MP3_pipeline_pipeline_listbox.Value;
% if isequal(SelectedModuleIndex, 0) || isempty(SelectedModuleIndex)
%     Entrie = {};
if ~isfield(handles, 'MP3_pipeline_ParamsModules')% || ~isfield(handles.MP3_pipeline_ParamsModules.(handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex}), 'Jobs') || isempty(fieldnames(handles.MP3_pipeline_ParamsModules.(handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex}).Jobs))
    Entrie = {''};
else
    %SelectedModule = handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex};
    SelectedModule = handles.MP3_pipeline_pipeline_listbox.String{SelectedModuleIndex};
    revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
    SelectedModule = revertcolor2(SelectedModule);
    Module = handles.MP3_pipeline_ParamsModules.(SelectedModule);
    if ~isfield(Module, 'Jobs') || isempty(fieldnames(Module.Jobs))
        Entrie = {''};
    else
        SelectedJobIndex = handles.MP3_pipeline_JobsList.Value;
        %SelectedJob = handles.MP3_pipeline_JobsList.String{SelectedJobIndex};
        SelectedJob = handles.MP3_pipeline_JobsNames{SelectedJobIndex};
        Job = Module.Jobs.(SelectedJob);
        SelectedParameterFieldIndex = handles.MP3_pipeline_JobsParametersFieldsList.Value;
        SelectedParameterField = handles.MP3_pipeline_JobsParametersFieldsList.String{SelectedParameterFieldIndex};
        Fields = strsplit(SelectedParameterField);
        Param = Job.(Fields{1});
        Entrie = Param.(Fields{2});
        if strcmp(Fields{1}, 'files_in') || strcmp(Fields{1}, 'files_out')
            NewEntrie = [];
            for i=1:length(Entrie)
                [~,name,~] = fileparts(Entrie{i});
                NewEntrie = [NewEntrie; {name}];
            end
            Entrie = NewEntrie;
        end
    end
end
if ~islogical(Entrie) && ~istable(Entrie)
    set(handles.MP3_pipeline_JobsParametersValues, 'Value', 1)
    set(handles.MP3_pipeline_JobsParametersValues, 'String', Entrie)
else
    set(handles.MP3_pipeline_JobsParametersValues, 'String', {'Cannot display this type.'})
end

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_JobsParametersFieldsList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_JobsParametersFieldsList


% --- Executes during object creation, after setting all properties.
function MP3_pipeline_JobsParametersFieldsList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_JobsParametersFieldsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MP3_pipeline_JobsParametersValues.
function MP3_pipeline_JobsParametersValues_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_JobsParametersValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MP3_pipeline_JobsParametersValues contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MP3_pipeline_JobsParametersValues


% --- Executes during object creation, after setting all properties.
function MP3_pipeline_JobsParametersValues_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_JobsParametersValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MP3_pipeline_save_pipeline.
function MP3_pipeline_save_pipeline_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_save_pipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'MP3_pipeline_ParamsModules') || isempty(fieldnames(handles.MP3_pipeline_ParamsModules))
    msgbox('There is no pipeline to save ...', 'No Pipeline');
    return
end

if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Saved_Pipelines'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Saved_Pipelines']);
    if status == false
        error('Cannot create the Saved_Pipelines folder to save the pipelines.')
    end
end


if ~isempty(handles.MP3_pipeline_ParamsModules)
    listing = what([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Saved_Pipelines']);
    pipeline_listing = [listing.mat', 'Other']';
    % Enter the ROI's name
    [pipeline_number,ok] = listdlg('Name', 'Name of the pipeline ?', 'SelectionMode', 'single', 'ListString',  cellstr(pipeline_listing),'ListSize', [200 150],...
        'PromptString', 'Select the name of the pipeline');
    if ok == 0
        return
    end
    
    if strcmp('Other', char(pipeline_listing(pipeline_number))) == 1
        FinalAnswer = inputdlg('Name of the new pipeline ', 'Question?', 1, {''});
        if isempty(FinalAnswer)
            return
        end
        FinalAnswer = {[FinalAnswer{:}, '.mat']};
    else
        FinalAnswer = cellstr(pipeline_listing(pipeline_number));
    end
    FinalAnswer = FinalAnswer{1};
    
    % evaluate if the pipeline name exist already
    % if yes there are 2 options, overwrite or cancel
    Flag = true;
    while any(contains(FinalAnswer, listing.mat)) && Flag
        quest = [FinalAnswer, ' is already used, do you want to overwrite it ?'];
        answer2 = questdlg(quest);
        switch answer2
            case 'Yes'
                Flag = false;
            case {'Cancel', 'No'}
                return
        end
    end
  
    %[file, path] = uiputfile('MyPipeline.mat');
    Pipeline = handles.MP3_pipeline_ParamsModules;
    %selpath = uigetdir(handles.MP3_data.database.Properties.UserData.MP3_data_path,'Please select a file to save your pipeline in.');
    Modules = fieldnames(Pipeline);
    for i=1:length(Modules)
        Module = Pipeline.(Modules{i});
        Module = rmfield(Module, 'Jobs');
        Module = rmfield(Module, 'OutputDatabase');
        Pipeline.(Modules{i}) = Module;
    end
    save([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Saved_Pipelines', filesep, FinalAnswer],'-struct', 'Pipeline');
    msgbox('Pipeline saved!','Done!');
else
    msgbox('There is no pipeline to save ...', 'No Pipeline');
end



% --- Executes on button press in MP3_pipeline_load_pipeline.
function MP3_pipeline_load_pipeline_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_load_pipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'MP3_pipeline_ParamsModules')
    if ~isempty(handles.MP3_pipeline_ParamsModules)
        quest = 'Would you like to replace your current pipeline by the one to load or to merge them ?';
        answer = questdlg(quest, 'What to do with the current pipeline ?' , 'Replace', 'Merge', 'Cancel', 'Cancel');
        if strcmp(answer, 'Cancel') || isempty(answer)
            return
        elseif strcmp(answer, 'Replace')
            [hObject, eventdata, handles] = MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles);
        elseif strcmp(answer, 'Merge')
            %msgbox('Currently not implemented... Sorry !')
            %return
%             while ~isempty(intersect(fieldnames(old_modules), Name_New_Mod))
%                 Name_New_Mod = [handles.new_module.module_name, '_', num2str(j)];
%                 j=j+1;
%             end
        end
    end
end

list = what([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Saved_Pipelines']);
if isempty(list)
    msgbox('There is no saved pipelines');
    return
end
[indx,tf] = listdlg('ListString',list.mat,'PromptString','Select the pipeline to load.', 'SelectionMode','single', 'ListSize',[300,300]);
if tf == 0
    return
end
PipelineName = list.mat{indx};
pipeline_loaded = load([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Saved_Pipelines', filesep, PipelineName]);

if exist('answer', 'var') && strcmp(answer, 'Merge')
    Names_old = fieldnames(handles.MP3_pipeline_ParamsModules); % Existing Pipeline
    Names_loaded = fieldnames(pipeline_loaded); % Loaded Pipeline
    New_pipeline = handles.MP3_pipeline_ParamsModules;
    j=1;
    for i=1:length(Names_loaded)
        Names_loaded_new = Names_loaded{i};
        while ~isempty(intersect(Names_old, Names_loaded_new))
            Names_loaded_new = [Names_loaded{i}, '_', num2str(j)];
            j=j+1;
        end
        New_pipeline.(Names_loaded_new) = pipeline_loaded.(Names_loaded{i});
    end
    pipeline = New_pipeline;
else
    pipeline = pipeline_loaded;
end
    


Modules = fieldnames(pipeline);
for i=1:length(Modules)
    %% Delete all modules filters for now, as we don't really use them efficiently.
    % This allow us to apply on a filtered database an already filtered designed pipeline. 
    pipeline.(Modules{i}).Filters = {};
end

handles.MP3_pipeline_ParamsModules = pipeline;
if ~isequal(handles.MP3_pipeline_TmpDatabase, handles.MP3_pipeline_Filtered_Table)
    quest = 'Would you like to apply the loaded pipeline on the whole database or on the filtered one you defined ?';
    answer = questdlg(quest, 'On which data apply the pipeline ?', 'Whole database', 'Filtered database', 'Whole database');
    if isempty(answer)
        return
    end
    switch answer
        case 'Whole database'
            Tmpdatab = handles.MP3_pipeline_TmpDatabase;
        case 'Filtered database'
            Tmpdatab = handles.MP3_pipeline_Filtered_Table;
    end
else
    Tmpdatab = handles.MP3_pipeline_TmpDatabase;
end

StoredDatab = handles.MP3_pipeline_TmpDatabase;
handles.MP3_pipeline_TmpDatabase = Tmpdatab;
[hObject, eventdata, handles] = MP3_pipeline_UpdatePipelineJobs(hObject, eventdata, handles);
handles.MP3_pipeline_TmpDatabase = StoredDatab;

[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);

% %Tmpdatab = handles.MP3_pipeline_TmpDatabase;
% for i=1:length(Modules)
%     Module = pipeline.(Modules{i});
%     %% Delete all modules filters for now, as we don't really use them efficiently.
%     % This allow us to apply on a filtered database an already filtered designed pipeline. 
%     Module.Filters = {};
%     %%
%     [pipeline_module, output_database_module] = MP3_pipeline_generate_psom_modules(Module.ModuleParams, Module.Filters, Tmpdatab, handles.MP3_data.database.Properties.UserData.MP3_data_path);
% %     if isempty(fieldnames(pipeline_module)) && isempty(output_database_module)
% %         continue
% %     end
%     pipeline.(Modules{i}).Filters = Module.Filters;
%     pipeline.(Modules{i}).Jobs = pipeline_module;
%     pipeline.(Modules{i}).OutputDatabase = output_database_module;
%     Tmpdatab = [Tmpdatab; output_database_module];
% end
% handles.MP3_pipeline_ParamsModules = pipeline;
% 
% 
% set(handles.MP3_pipeline_pipeline_listbox,'String', fieldnames(handles.MP3_pipeline_ParamsModules));
% set(handles.MP3_pipeline_pipeline_listbox,'Value', 1);
% MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
% [hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
% [hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);

guidata(hObject, handles);

function [hObject, eventdata, handles] = MP3_pipeline_UpdatePipelineJobs(hObject, eventdata, handles)

set(handles.MP3_pipeline_manager_GUI, 'pointer', 'watch');
drawnow;

pipeline = handles.MP3_pipeline_ParamsModules;
Modules = fieldnames(pipeline);
%% Attention Code etrange qui distingue le cas d'un load du cas d'un save module.
% * TmpDatabase contient toute la database de MP3 plus les fichiers
%temporaires des précédents modules.
% * Filtered data contient la Tmp database filtrée par les filtres adéquats.
% * La fonction clear_pipeline met à jour la TmpDatabase ET la Filtered_Table
% * Pour le load, je dois avoir une database vierge de tous fichiers
%temporaires.
%  * En réalité, on ne doit pas spécialement distinguer le load du save mais
% surtout le cas "Whole database" et "Filtered database" du load. Et pour le
% gérer, on doit récuperer la database en question (whole ou filtered) dans
% la variable handles.MP3_pipeline_TmpDatabase. Puisque Clear modifie cette
% valeur, il faut la stocker avant. Pour le save, il faut au contraire
% d'abord clear le pipeline précedent, afin de se débarasser des fichiers
% temporaires et ensuite appliquer le pipeline modifié à la table filtrée.
%  * Remarque : la valeur stockée dans handles.MP3_pipeline_TmpDatabase est
% déjà débarassée de tout fichier temporaire dans la fonction load.
if strcmp(eventdata.Source.Tag, 'MP3_pipeline_Save_Module') || strcmp(eventdata.Source.Tag, 'MP3_pipeline_DeleteModule')
    [hObject, eventdata, handles] = MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles);
    Tmpdatab = handles.MP3_pipeline_Filtered_Table;
else
    Tmpdatab = handles.MP3_pipeline_TmpDatabase;
    [hObject, eventdata, handles] = MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles);
end
%% 
% [hObject, eventdata, handles] = MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles);
% [hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
% [hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);
% Tmpdatab = handles.MP3_pipeline_Filtered_Table;


%Tmpdatab = handles.MP3_pipeline_Filtered_Table;
%Tmpdatab = handles.MP3_pipeline_TmpDatabase;

for i=1:length(Modules)
    Module = pipeline.(Modules{i});
%     %% Delete all modules filters for now, as we don't really use them efficiently.
%     % This allow us to apply on a filtered database an already filtered designed pipeline. 
%     Module.Filters = {};
%     %%
    [pipeline_module, output_database_module] = MP3_pipeline_generate_psom_modules(Module.ModuleParams, Module.Filters, Tmpdatab, handles.MP3_data.database.Properties.UserData.MP3_data_path, 0);
    pipeline.(Modules{i}).Filters = Module.Filters;
    pipeline.(Modules{i}).Jobs = pipeline_module;
    %pipeline.(Modules{i}).JobsReWrite = CheckReWriting(pipeline_module, Tmpdatab);
    pipeline.(Modules{i}).OutputDatabase = output_database_module;
    Tmpdatab = [Tmpdatab; output_database_module];
end

handles.MP3_pipeline_ParamsModules = pipeline;
%handles.MP3_pipeline_pipeline_listbox_Raw = fieldnames(handles.MP3_pipeline_ParamsModules);

%[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);

Coloredlistbox = DisplayColoredListbox(handles.MP3_pipeline_pipeline_listbox.String, handles);
set(handles.MP3_pipeline_pipeline_listbox,'String', Coloredlistbox);
%set(handles.MP3_pipeline_pipeline_listbox,'String', fieldnames(handles.MP3_pipeline_ParamsModules));
set(handles.MP3_pipeline_pipeline_listbox,'Value', 1);
[hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);

set(handles.MP3_pipeline_manager_GUI, 'pointer', 'arrow');
guidata(hObject, handles);


function Status = CheckReWriting(Module, Tmpdatab)

JobNames = fieldnames(Module.Jobs);
Status = false(size(JobNames));

for i=1:length(JobNames)
    files_out = Module.Jobs.(JobNames{i}).files_out;
    if isempty(files_out)
        Status(i) = 0;
        continue
    end
    inputNames = fieldnames(files_out);
    for j=1:length(inputNames)
        file = files_out.(inputNames{j});
        [path, filename, ~] = fileparts(file{1});
        Datab = Module.Jobs.(JobNames{i}).opt.Table_out(Module.Jobs.(JobNames{i}).opt.Table_out.Filename == categorical(cellstr(filename)),:);
        assert(height(Datab) == 1);
        
        % Dont check on the path cause the Datab one is still set to Tmp.
        Lia = ismember(Tmpdatab(:,[1:3 5:end]), Datab(:,[1:3 5:end]));
        %UPDATE : on the previous line, we could have the case where 2
        %files have the same triplet Patient/Tp/SequenceName but a
        %different Filename. Then, as its 2 differents entries, the new one
        %will be concatenate to the database and leads to an error in the
        %viewer.
        % So now we are going to check the unicity of this triplet.
        Lib = ismember(Tmpdatab(:,[2 3 8]), Datab(:,[2 3 8]));
        
        % UPDATE : If sum(Lib) is > 1, then an entry with the same triplet will
        % return a warning to the user (orange job), and if sum(Lia) is >
        % 1, then an entry with a different filename (and other tags) will 
        % return a warning to the user (orange job)
        
        
        
        
        if sum(Lia)==1 && sum(Lib)==1
        % Le job ecrit un fichier qui n'existe pas, la seule entree
        % retournee par intersect est l'entree temporaire du job en
        % question
            Status(i) = 0;
        else
        % Le job ecrit un fichier qui existe deja soit reellement soit
        % virtuellement.intersect retourne 2 entree : l'entree temporaire
        % du job en question et une autre, reelle ou virtuelle. Dans les 2
        % cas cela pose un soucis.
            %Folders = strsplit(char(Tmpdatab.Path(ib)),filesep);
            %Fold = Folsers(end-1);
            %if ~strcmp( Fold == 'Tmp')
            
            
            %UPDATE : Depuis la prise en charge du module Delete_file, on peut
            %se retrouver avec 2 entrées qui partagent tout sauf le Type. L'une
            %aura le type Deleted, et l'autre l'un des autres types
            %disponibles.
            %Dans ce cas, on prefere que le job ne s'affiche pas en orange mais
            %en vert.
            ListEntries = Tmpdatab.Type(Lib,:) == categorical({'Deleted'});
            
            if length(ListEntries) == 2 && sum(ListEntries) == 1
                Status(i) = 0;
            else
                Status(i) = 1;
            end
        end        
    end
end

    
    





% --- Executes on button press in MP3_pipeline_Delete_Job.
function MP3_pipeline_Delete_Job_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Delete_Job (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    return
end



SelectedModuleIndex = handles.MP3_pipeline_pipeline_listbox.Value;
%SelectedModule = handles.MP3_pipeline_pipeline_listbox_Raw{SelectedModuleIndex};
SelectedModule = handles.MP3_pipeline_pipeline_listbox.String{SelectedModuleIndex};
revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
SelectedModule = revertcolor2(SelectedModule);
Module = handles.MP3_pipeline_ParamsModules.(SelectedModule);
SelectedJobIndex = handles.MP3_pipeline_JobsList.Value;
%SelectedJob = handles.MP3_pipeline_JobsList.String{SelectedJobIndex};
if isempty(handles.MP3_pipeline_JobsNames)
    return
end
SelectedJob = handles.MP3_pipeline_JobsNames{SelectedJobIndex};
Job = Module.Jobs.(SelectedJob);
Module.Jobs = rmfield(Module.Jobs, SelectedJob);
if isempty(fieldnames(Module.Jobs))
    MP3_pipeline_DeleteModule_Callback(hObject, eventdata, handles)
    handles.MP3_pipeline_ParamsModules = rmfield(handles.MP3_pipeline_ParamsModules, SelectedModule);
    if isempty(fieldnames(handles.MP3_pipeline_ParamsModules))
        handles = rmfield(handles, 'MP3_pipeline_ParamsModules');
    end
else
    
    %% update Output databse of the  selected module if needed.
    if ~isempty(Job.files_out)
        Outputs = fieldnames(Job.files_out);
        for i=1:length(Outputs)
            Files = Job.files_out.(Outputs{i});
            for j=1:length(Files)
                [path, name, ~] = fileparts(Files{j});
                line = Module.OutputDatabase(Module.OutputDatabase.Filename == categorical(cellstr(name)),:);
                assert(line.Path == categorical(cellstr([path, filesep])));
                Module.OutputDatabase(Module.OutputDatabase.Filename == categorical(cellstr(name)),:) = [];
            end
        end
    end
    handles.MP3_pipeline_ParamsModules.(SelectedModule) =  Module;
    [hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
    [hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);
    [hObject, eventdata, handles] = MP3_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles);
    
end
guidata(hObject, handles);

function Coloredlistbox = DisplayColoredListbox(Names, handles)
if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    Coloredlistbox = {''};
    return
end
Names = fieldnames(handles.MP3_pipeline_ParamsModules);
%colergenlistbox = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
color2 = @(color,text) ['<HTML><FONT color="',color,'">',text,'</Font></html>'];
%revertcolor2 = @(string) extractAfter(extractBefore(string,'</Font></html>'), '">');
Coloredlistbox = cell(size(Names));
for i=1:length(Names)
    ReWritting = CheckReWriting(handles.MP3_pipeline_ParamsModules.(Names{i}), handles.MP3_pipeline_TmpDatabase);
%     if isempty(fieldnames(handles.MP3_pipeline_ParamsModules.(Names{i}).Jobs))
%         Coloredlistbox{i} = color2('red', Names{i});
%     else
%         Coloredlistbox{i} = color2('green', Names{i});
%     end
    if isempty(fieldnames(handles.MP3_pipeline_ParamsModules.(Names{i}).Jobs))
        Coloredlistbox{i} = color2('red', Names{i});
    elseif any(ReWritting)
        Coloredlistbox{i} = color2('orange', Names{i});
    else
        Coloredlistbox{i} = color2('green', Names{i});
    end
end



% handles.MP3_pipeline_pipeline_listbox_Raw = fieldnames(handles.MP3_pipeline_ParamsModules);
% Coloredlistbox = DisplayColoredListbox(handles.MP3_pipeline_pipeline_listbox_Raw, handles);
% set(handles.MP3_pipeline_pipeline_listbox,'String', Coloredlistbox);


% --- Executes on button press in MP3_pipeline_Clear_PSOM_history_button.
function MP3_pipeline_Clear_PSOM_history_button_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Clear_PSOM_history_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'PSOM'], 'dir') == 7
    rmdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'PSOM'], 's')
    msgbox('Done', 'Information') ;
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Vectorial_Screenshot_Callback(hObject, eventdata, handles)
% hObject    handle to Vectorial_Screenshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = handles.MP3_pipeline_manager_GUI;

[file,path] = uiputfile('.eps', 'Save the MP3 pipeline manager as eps:');

fig.InvertHardcopy = 'off';  %background
%fig.PaperOrientation = 'landscape';
set(fig,'PaperPositionMode','auto'); % size position

%print(fig, '-bestfit', [path, file],'-dpdf')
%print(fig, '-fillpage', [path, file],'-dpdf')

print(fig, [path, file],'-depsc2')
msgbox('Done', 'Information') ;



function Select_Number_Workers_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Number_Workers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Select_Number_Workers as text
%        str2double(get(hObject,'String')) returns contents of Select_Number_Workers as a double


% --- Executes during object creation, after setting all properties.
function Select_Number_Workers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Select_Number_Workers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [hObject, eventdata, handles] = Update_MP3_database_After_pipeline_Crash(hObject, eventdata, handles)
fname = [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep, 'Table_out', filesep, 'Table_out.mat'];
handles.FlagExecutePipe = 1;
if ~exist(fname)
    handles.FlagExecutePipe = 0;
    return
else
    answer = questdlg('It seems that some files from a previous pipeline execution have not been added to your project database. Want to do it now?', ...
	'Did your previous pipeline crashed?', ...
	'Yes!','No','Cancel','Cancel');
    if strcmp(answer, 'Yes!')
        if isfield(handles, 'MP3_pipeline_pipeline_listbox') && ~isempty(handles.MP3_pipeline_pipeline_listbox.String)
            answer2 = questdlg('This manipulation will erase your current pipeline. Do you want to save it before?', ...
            'It''s now or never', ...
            'Yes!','No','Cancel','Cancel');
            if strcmp(answer2, 'Yes!')
                MP3_pipeline_save_pipeline_Callback(hObject, eventdata, handles);
            end
        end
    else
        handles.FlagExecutePipe = 0;
        return
    end
end
load(fname)
PreviousDB = handles.MP3_data.database;
for i=1:size(Output_Table,1)
    Entry = Output_Table(i,:);
    filename = [char(Entry.Path), char(Entry.Filename), '.nii'];
    if ~exist(filename)
        continue
    end
    switch char(Entry.Type)
        case 'Scan'
            movefile(filename, strrep(filename,[handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep], [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Derived_data', filesep]));
            jsonfilename = strrep(filename, '.nii', '.json');
            if exist(jsonfilename, 'file')
            movefile(jsonfilename, strrep(jsonfilename,[handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep], [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Derived_data', filesep]));
            end
            Entry.Path = categorical(cellstr(handles.MP3_data.database.Properties.UserData.MP3_Derived_data_path));
        case 'ROI'
            movefile(filename, strrep(filename,[handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep], [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'ROI_data', filesep]));
            Entry.Path = categorical(cellstr(handles.MP3_data.database.Properties.UserData.MP3_ROI_path));
        case 'Cluster'
            movefile(filename, strrep(filename,[handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep], [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'ROI_data', filesep]));
            matfilename = strrep(filename, '.nii', '.mat');
            if exist(matfilename, 'file')
             movefile(matfilename, strrep(jsonfilename,[handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep], [handles.MP3_data.database.Properties.UserData.MP3_data_path, 'ROI_data', filesep]));
            Entry.Path = categorical(cellstr(handles.MP3_data.database.Properties.UserData.MP3_ROI_path));
            end
    end
    handles.MP3_data.database = [handles.MP3_data.database; Entry];
    
end

if ~isempty(setdiff(sortrows(handles.MP3_data.database), sortrows(PreviousDB)))
    handles2 = guidata(handles.MP3_data.MP3_GUI);
    handles2.database = handles.MP3_data.database;
    guidata(handles.MP3_data.MP3_GUI, handles2);

    MP3('MP3_update_database_display', hObject, eventdata,handles.MP3_data)
    MP3('MP3_menu_save_database_Callback', hObject, eventdata,handles.MP3_data)
    rmdir([handles.MP3_data.database.Properties.UserData.MP3_data_path, 'Tmp', filesep, 'Table_out'], 's');
    [hObject, eventdata, handles] = MP3_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles);
else
    handles.FlagExecutePipe = 0;
end


% --- Executes on button press in MP3_pipeline_Module_UP.
function MP3_pipeline_Module_UP_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Module_UP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    return
end

module_selected = get(handles.MP3_pipeline_pipeline_listbox, 'Value');
module_names = fieldnames(handles.MP3_pipeline_ParamsModules);
% the first module cannot be moving up
if module_selected == 1
    return
end

module_list = 1:numel(module_names);
module_list(module_selected-1) =module_selected;
module_list(module_selected) =module_selected-1;

handles.MP3_pipeline_ParamsModules = orderfields(handles.MP3_pipeline_ParamsModules,module_list);
guidata(hObject, handles);
MP3_pipeline_update_pipeline_Callback(hObject, eventdata, handles)

% --- Executes on button press in MP3_pipeline_Module_DOWN.
function MP3_pipeline_Module_DOWN_Callback(hObject, eventdata, handles)
% hObject    handle to MP3_pipeline_Module_DOWN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    return
end

module_selected = get(handles.MP3_pipeline_pipeline_listbox, 'Value');
module_names = fieldnames(handles.MP3_pipeline_ParamsModules);
% the last module cannot be moving down
if module_selected == numel(module_names)
    return
end
module_list = 1:numel(module_names);
module_list(module_selected+1) =module_selected;
module_list(module_selected) =module_selected+1;

handles.MP3_pipeline_ParamsModules = orderfields(handles.MP3_pipeline_ParamsModules,module_list);
%save the new order
guidata(hObject, handles);
% update the pipeline
MP3_pipeline_update_pipeline_Callback(hObject, eventdata, handles)

function MP3_pipeline_update_pipeline_Callback(hObject, eventdata, handles)

if ~isfield(handles, 'MP3_pipeline_ParamsModules')
    return
end

if ~isequal(handles.MP3_pipeline_TmpDatabase, handles.MP3_pipeline_Filtered_Table)
    Tmpdatab = handles.MP3_pipeline_Filtered_Table;
else
    Tmpdatab = handles.MP3_pipeline_TmpDatabase;
end

StoredDatab = handles.MP3_pipeline_TmpDatabase;
handles.MP3_pipeline_TmpDatabase = Tmpdatab;
[hObject, eventdata, handles] = MP3_pipeline_UpdatePipelineJobs(hObject, eventdata, handles);
handles.MP3_pipeline_TmpDatabase = StoredDatab;

[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
[hObject, eventdata, handles] = MP3_pipeline_UpdateTables(hObject, eventdata, handles);

guidata(hObject, handles);



% --- Executes on selection change in listbox9.
function listbox9_Callback(hObject, eventdata, handles)
% hObject    handle to listbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox9


% --- Executes during object creation, after setting all properties.
function listbox9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
