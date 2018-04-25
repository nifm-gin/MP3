function varargout = MIA_pipeline(varargin)
% MIA_PIPELINE MATLAB code for MIA_pipeline.fig
%      MIA_PIPELINE, by itself, creates a new MIA_PIPELINE or raises the existing
%      singleton*.
%
%      H = MIA_PIPELINE returns the handle to a new MIA_PIPELINE or the handle to
%      the existing singleton*.
%
%      MIA_PIPELINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIA_PIPELINE.M with the given input arguments.
%
%      MIA_PIPELINE('Property','Value',...) creates a new MIA_PIPELINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MIA_pipeline_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MIA_pipeline_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MIA_pipeline

% Last Modified by GUIDE v2.5 23-Apr-2018 11:35:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MIA_pipeline_OpeningFcn, ...
                   'gui_OutputFcn',  @MIA_pipeline_OutputFcn, ...
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





% --- Executes just before MIA_pipeline is made visible.
function MIA_pipeline_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MIA_pipeline (see VARARGIN)

% Choose default command line output for MIA_pipeline
handles.output = hObject;

handles.MIA_data = varargin{3};
handles.Modules_listing = {'Relaxometry', '   .T1map (Multi Inversion Time)', '   .T1map (Multi Angles)', '   .T2map', '   .Fit_T2_T2star',...
                '   .deltaR2', '   .deltaR2*', '   .MGE2Dfrom3D',...
    'Perfusion', '   .Blood volume fraction (steady-state)', '   .Dynamic Susceptibility Contrast', '   .Vessel Size Imaging (steady-state)', ...
                '   .Vessel Density (steady-state)', '   .Cerebral blood flow (ASL)',  '   .Cerebral blood flow (ASL-Dynamic)',...
                '   .Inversion Efficiency (ASL_InvEff)',...
     'Diffusion', '   .ADCmap',...
     'Permeability', '   .Dynamic Contrast Enhancement (Phenomenology)', '   .Dynamic Contrast Enhancement (Quantitative)',...          
     'Oxygenation', '   .R2prim', '   .SO2map', '   .CMRO2',...
     'MRFingerprint', '   .Vascular MRFingerprint'...
     'SPM', '   .SPM: Coreg (Est)', '   .SPM: Coreg (Est & Res)', '   .SPM: Reslice','   .SPM: Realign', ...
     'Spatial', '   .Smoothing', '   .Arithmetic'...
     };
handles.Module_groups = {'Relaxometry','Perfusion', 'Diffusion', 'Permeability', 'Oxygenation', 'MRFingerprint', 'SPM', 'Spatial' };
 
set(handles.MIA_pipeline_module_listbox, 'String', handles.Modules_listing);
handles.Add_Tags_listing = handles.MIA_data.database.Properties.VariableNames;
handles.Add_list = handles.Add_Tags_listing;
set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_Tags_listing);
set(handles.MIA_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})
handles.Source_selected = handles.Add_Tags_listing{1};
%handles.MIA_pipeline_Unique_Values_Selection = 
handles.MIA_pipeline_TagsToPrint = {'Patient', 'Tp', 'SequenceName'};
handles.Remove_Tags_listing = {'NoMoreTags'};
handles.Remove_selected = handles.Remove_Tags_listing{1};
%handles.MIA_data.database.IsRaw = categorical(handles.MIA_data.database.IsRaw);

handles.MIA_pipeline_TmpDatabase = handles.MIA_data.database;
handles.MIA_pipeline_Filtered_Table = handles.MIA_pipeline_TmpDatabase;
CellsColoured = DisplayColoredTable(handles.MIA_pipeline_TmpDatabase, handles.MIA_pipeline_TagsToPrint);
handles.MIA_pipeline_Filtering_Table.Data = CellsColoured;
% handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_pipeline_Filtered_Table{:,handles.MIA_pipeline_TagsToPrint});
% colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
% for i=1:size(handles.MIA_pipeline_Filtering_Table.Data, 1)
%     handles.MIA_pipeline_Filtering_Table.Data{i,1} = colergen(rgb2hex([1 0 0]), handles.MIA_pipeline_Filtering_Table.Data{i,1});
% end



handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;

% handles.module_parameters_string = {};
% handles.module_parameters_fields = {};

MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
handles.FilterParameters = {};
%MIA_pipeline_Unique_Values_Tag_CellSelectionCallback(hObject, eventdata, handles)
% Update handles structure
handles.MIA_pipeline_Unique_Values_Selection = {};

% update the 'String' of MIA_pipeline_pushMIASelection and MIA_pipeline_pushMIATPSelection push button
data_selected =  MIA2('finddata_selected',handles.MIA_data);
set(handles.MIA_pipeline_pushMIASelection, 'String', [char(handles.MIA_data.database.Patient(data_selected(1))) '-' char(handles.MIA_data.database.Tp(data_selected(1))) ' only'])
set(handles.MIA_pipeline_pushMIATPSelection, 'String', ['All time point of :' char(handles.MIA_data.database.Patient(data_selected(1)))])

%set(findall(gcf,'-property','FontName'), 'FontName', 'Courier')
set(findall(handles.MIA_pipeline_manager_GUI.Children, 'Tag', 'MIA_pipeline_module_parameters'), 'FontName', 'Courier')
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
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex([0 1 1]), CellsColoured{i,j}); %Cyan
            end
        case 'Derived_data'
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex([0 1 0]), CellsColoured{i,j});%Green
            end
        case 'Tmp'
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex([1 1 0]), CellsColoured{i,j});%Yellow
            end
        case 'ROI_data'
            for j=1:size(CellsColoured,2)
                CellsColoured{i,j} = colergen(rgb2hex([1 0 1]), CellsColoured{i,j});%Magenta
            end
            
    end
    
end


% UIWAIT makes MIA_pipeline wait for user response (see UIRESUME)
% uiwait(handles.MIA_pipeline_manager_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = MIA_pipeline_OutputFcn(hObject, eventdata, handles)
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

function [hObject, eventdata, handles] = MIA_pipeline_UpdateTables(hObject, eventdata, handles)
%% Update of Filtering/Filtered Tables
NewTable = handles.MIA_pipeline_TmpDatabase;
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

handles.MIA_pipeline_Filtered_Table = NewTable;
CellsColoured = DisplayColoredTable(NewTable, handles.MIA_pipeline_TagsToPrint);
handles.MIA_pipeline_Filtering_Table.Data = CellsColoured;
handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;
%% Update of Unique Values
if strcmp(handles.Source_selected , 'NoMoreTags')
    TagValues = {};
else
    TagValues = unique(handles.MIA_pipeline_Filtered_Table{:,handles.Source_selected});
end
handles.MIA_pipeline_Unique_Values_Tag.Data = cellstr(TagValues);

%% Update of modules parameters
MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles)
guidata(hObject, handles);



% --- Executes on button press in MIA_pipeline_add_module_button.
function MIA_pipeline_add_module_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_add_module_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[new_pipeline, output_database] = MIA_pipeline_generate_psom_modules(handles.new_module, handles.FilterParameters, handles.MIA_pipeline_TmpDatabase, handles.MIA_data.database.Properties.UserData.MIA_data_path);

if isempty(fieldnames(new_pipeline)) && isempty(output_database)
    return
end


% if isfield(handles, 'tmp_database')
%     handles.tmp_database = [handles.tmp_database ; output_database];
% else
%     handles.tmp_database = output_database;
% end
    
if isfield(handles, 'MIA_pipeline_ParamsModules')
    old_modules = handles.MIA_pipeline_ParamsModules;
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
handles.MIA_pipeline_ParamsModules.(Name_New_Mod) = SaveModule;

%handles.psom.Output_databases = setfield(old_databases, Name_New_Mod, output_database);
%handles.psom.Modules = setfield(old_modules, Name_New_Mod, new_pipeline);


%handles.MIA_pipeline_Filtered_Table = [handles.MIA_pipeline_Filtered_Table ; output_database];
%handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_pipeline_Filtered_Table{:,handles.MIA_pipeline_TagsToPrint});
%handles.psom.pipeline = merged_pipe;
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
%handles.MIA_pipeline_TmpDatabase = unique([handles.MIA_pipeline_TmpDatabase ; output_database]);



[hObject, eventdata, handles] = MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%guidata(hObject, handles);
%MIA_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)
%MIA_pipeline_Remove_Tag_Button_Callback(hObject, eventdata, handles)


%module_listing = get(handles.MIA_pipeline_pipeline_listbox,'String');
%set(handles.MIA_pipeline_pipeline_listbox,'String', [module_listing' {handles.new_module.module_name}]');
set(handles.MIA_pipeline_pipeline_listbox,'String', fieldnames(handles.MIA_pipeline_ParamsModules));

guidata(hObject, handles);

%MIA_pipeline_exectute_module_button_Callback(hObject, eventdata, handles)






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
printdlg(handles.MIA_pipeline_manager_GUI)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.MIA_pipeline_manager_GUI,'Name') '?'],...
                     ['Close ' get(handles.MIA_pipeline_manager_GUI,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.MIA_pipeline_manager_GUI)



% --- Executes on selection change in MIA_pipeline_module_parameters.
function MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_module_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_module_parameters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_module_parameters


parameter_selected = get(handles.MIA_pipeline_module_parameters,'Value');
% display the help associated to the parameter selected
set(handles.MIA_pipeline_parameter_setup_text, 'String', handles.new_module.opt.table.Help{parameter_selected});

switch handles.new_module.opt.table.Type{parameter_selected}
    case 'Text'
        %Text = handles.new_module.opt.table.Default{1}';
        %disp('l MIA_pipeline ligne 225')

        table.ColumnFormat = {'char'};
        table.data = '';
        table.columnName = '';
        table.editable = false;
    case '1Scan1TPXP'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
           
            TP_listing = unique(handles.MIA_pipeline_Filtered_Table.Tp(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(TP_listing),3) = cellstr(TP_listing);
            table.data(1:numel(TP_listing),4) = {false};
            %table.data(1,2) = {'Choose'};
            
            Patients_listing = unique(handles.MIA_pipeline_Filtered_Table.Patient(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(Patients_listing),5) = cellstr(Patients_listing);
            table.data(1:numel(Patients_listing),6) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName', 'Tp', 'Patient'};
        table.columnName = {'SequenceName', 'Select ONE Input','Tp', 'Select ONE Input', 'Patient', 'Select Input'};
        handles.new_module.opt.ColumnNamesInput1Scan1TPXP = table.columnName;
        table.editable = [false true false true false true];
        table.ColumnFormat = {'char'}; 
        
    case '1Scan'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        table.columnName = {'SequenceName', 'Select ONE Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
    case 'XScan'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
        
    case '1ScanOr1ROI'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName);
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            for i=1:numel(SequenceType_listing)
                table.data(i,3) = cellstr(unique(handles.MIA_pipeline_Filtered_Table.Type(handles.MIA_pipeline_Filtered_Table.SequenceName == SequenceType_listing(i))));
            end
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        table.columnName = {'SequenceName', 'Select ONE Input', 'Type'};
        table.editable = [false true false];
        table.ColumnFormat = {'char'};  
    case '1ROI'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'ROI'));
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        table.columnName = {'SequenceName', 'Select ONE Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
    case 'XROI'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'ROI'));
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
     case 'XScanOrXROI'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName);
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
        end
        %handles.new_module.opt.DOF{parameter_selected} = {'SequenceName'};
        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};
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
    otherwise
        table.ColumnFormat = {'char'};
        table.data = '';
        table.columnName = '';
        table.editable = false;
end



%% update the setup table
set(handles.MIA_pipeline_parameter_setup_table, 'ColumnFormat', table.ColumnFormat);
% set names of the columns
set(handles.MIA_pipeline_parameter_setup_table, 'ColumnName', table.columnName);
% set data (default's parameters)
set(handles.MIA_pipeline_parameter_setup_table, 'Data', table.data);
% set each colomn editable
set(handles.MIA_pipeline_parameter_setup_table, 'columnEditable',  table.editable );
if ~isempty(table.data)
    % set ColumnWidth to auto

    merge_Data = [table.columnName; table.data];
    dataSize = size(merge_Data);
    % Create an array to store the max length of data for each column
    maxLen = zeros(1,dataSize(2));
    % Find out the max length of data for each column
    % Iterate over each column
    for i=1:dataSize(2)
        % Iterate over each row
        for j=1:dataSize(1)
            len = length(merge_Data{j,i});
            % Store in maxLen only if its the data is of max length
            if(len > maxLen(1,i))
                maxLen(1,i) = len;
            end
        end
    end
    % Some calibration needed as ColumnWidth is in pixels
    cellMaxLen = num2cell(maxLen*6.5);
    % Set ColumnWidth of UITABLE
    set(handles.MIA_pipeline_parameter_setup_table, 'ColumnWidth', cellMaxLen);
    
end

guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function MIA_pipeline_module_parameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_module_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_setting_windows(hObject, eventdata, handles)




% --- Executes on selection change in MIA_pipeline_parameter_setup.
function MIA_pipeline_parameter_setup_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_parameter_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_parameter_setup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_parameter_setup



% --- Executes during object creation, after setting all properties.
function MIA_pipeline_parameter_setup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_parameter_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MIA_pipeline_clear_pipeline_button.
function MIA_pipeline_clear_pipeline_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_clear_pipeline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'psom')
    set(handles.MIA_pipeline_pipeline_listbox, 'String', '');
    if ~isempty(findobj('Tag', 'BioGraphTool'))
        close(findobj('Tag', 'BioGraphTool'));
    end
        handles = rmfield(handles, 'psom');

end

handles.tmp_database = table();
handles.MIA_pipeline_TmpDatabase = handles.MIA_data.database;
MIA_pipeline_UpdateTables(hObject, eventdata, handles);
guidata(hObject, handles);



function edge_callbacks(hObject, eventdata, handles)
eventdata = [];
handles = guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'));
module_list = get(handles.MIA_pipeline_module_popupmenu, 'String');
get(hObject, 'ID')
% sub_module = strfind(hObject.ID, '_');
% if ~isempty(sub_module)
%     module_name = hObject.ID(1:sub_module-1);
% else
%     module_name =hObject.ID;
% end
% idx = find(ismember(module_list, module_name));
% set(handles.MIA_pipeline_module_popupmenu, 'Value', idx)


function node_callbacks(hObject, ~, handles)
eventdata = [];
handles = guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'));

pipeline_module_names = fieldnames(handles.psom.pipeline);
idx = strcmp(pipeline_module_names,hObject.ID);
module_selected = handles.psom.pipeline.(pipeline_module_names{idx});
module_selected.files_in
%% update handles.MIA_pipeline_module_parameters using the information the module_selected




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
% set(handles.biograph_fig, 'Name', 'MIA pipeline manager');
% guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'), handles);

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


% --- Executes when entered data in editable cell(s) in MIA_pipeline_parameter_setup_table.
function MIA_pipeline_parameter_setup_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_parameter_setup_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
parameter_selected = get(handles.MIA_pipeline_module_parameters,'Value');

%table_data = get(handles.MIA_pipeline_parameter_setup_table, 'Data');
if strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XScan') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XScanOrXROI') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XROI')
    %handles.new_module.opt.parameter_default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data{cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,2)),1};
    handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1Scan') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ScanOr1ROI') || strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ROI')
    if sum(cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,2))) == 1 || sum(cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,2))) == 0
        handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
        
    end
% elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1ScanOr1ROI')
%     if sum(cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,3))) == 1 || sum(cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,3))) == 0
%         handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
%         
%     end
    
elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1Scan1TPXP')
    if isempty(handles.new_module.opt.table.Default{parameter_selected})
        handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
    else
        A = handles.MIA_pipeline_parameter_setup_table.Data(:,2);
        A = A(~cellfun('isempty',A));
        %C = handles.new_module.opt.table.Default{parameter_selected};
        %D = C(:,2);
        %D = D(~cellfun('isempty',D));
        B = handles.MIA_pipeline_parameter_setup_table.Data(:,4);
        B = B(~cellfun('isempty',B));
        %E = C(:,4);
        %E = E(~cellfun('isempty',E));
        %if sum(cell2mat(A)) + sum(cell2mat(D)) == 1 && sum(cell2mat(B)) + sum(cell2mat(E)) == 0
        %    handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
        %elseif sum(cell2mat(A)) + sum(cell2mat(D)) == 0 && sum(cell2mat(B)) + sum(cell2mat(E)) == 1
        %    handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
        %end
        if (sum(cell2mat(A)) == 1 || sum(cell2mat(A)) == 0) && (sum(cell2mat(B)) == 1 || sum(cell2mat(B)) == 0)
            handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
        end
    end
    %handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
else
    handles.new_module.opt.Module_settings = setfield(handles.new_module.opt.Module_settings, handles.new_module.opt.table.PSOM_Fields{parameter_selected},handles.MIA_pipeline_parameter_setup_table.Data{1,1}); 
    %handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data{1,1};
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end
    %patient_listing = table_data(:,1);
% patient_selected = patient_listing(find([table_data{:,2}]' == true));

% % case go back to all patient
% if table_data{1,2} == 1 && handles.new_module.files_in_filter_data{1,2} == 0
%     idex_patient = true(numel(handles.MIA_data.database.Patient),1);
%     table_data(1,2) = {true};
%     table_data(2:numel(unique(handles.MIA_data.database.Patient))+1,2) = {false};
%     % from all patient to a specific patient
% elseif sum(find([table_data{:,2}]' == true)) >1 && table_data{1,2} == 1
%     table_data{1,2} = false;
%     patient_selected = patient_listing(find([table_data{:,2}]' == true));
%     idex_patient = handles.MIA_data.database.Patient == patient_selected;
%     
%     % if all patient is selected
% elseif sum(find([table_data{:,2}]' == true))  == 1 && table_data{1,2} == 1
%     idex_patient = true(numel(handles.MIA_data.database.Patient),1);
%     % 1 or more patient (but not all)
% else
%     idex_patient = false(numel(handles.MIA_data.database.Patient),1);
%     for i=1:numel(patient_selected)
%         idex_patient = idex_patient | handles.MIA_data.database.Patient == patient_selected(i);
%     end
% end
% 
% tp_listing = table_data(:,3);
% tp_selected = tp_listing(find([table_data{:,4}]' == true));
% 
% % case go back to all time point
% if table_data{1,4} == 1 && handles.new_module.files_in_filter_data{1,4} == 0
%     idex_tp = true(numel(handles.MIA_data.database.Tp),1);
%     table_data(1,4) = {true};
%     table_data(2:numel(unique(handles.MIA_data.database.Tp))+1,4) = {false};
%     % from all time point to a specific time point
% elseif sum(find([table_data{:,4}]' == true)) >1 && table_data{1,4} == 1
%     table_data{1,4} = false;
%     tp_selected = tp_listing(find([table_data{:,4}]' == true));
%     idex_tp = handles.MIA_data.database.Tp == tp_selected;
%     
%     % if all time point is selected
% elseif sum(find([table_data{:,4}]' == true))  == 1 && table_data{1,4} == 1
%     idex_tp = true(numel(handles.MIA_data.database.Tp),1);
%     % if 1 or more time point (but not all)
% else
%     idex_tp = false(numel(handles.MIA_data.database.Tp),1);
%     for i=1:numel(tp_selected)
%         idex_tp = idex_tp | handles.MIA_data.database.Tp == tp_selected(i);
%     end
% end
% 
% SequenceName_listing =table_data(:,5);
% SequenceName_selected = SequenceName_listing(find([table_data{:,6}]' == true));
% if isempty(SequenceName_selected)
%     index_SequenceName =true(numel(handles.MIA_data.database.Tp,1));
% else
%      index_SequenceName = false(numel(handles.MIA_data.database.Tp),1);
%     for i=1:numel(SequenceName_selected)
%         index_SequenceName = index_SequenceName | handles.MIA_data.database.SequenceName == SequenceName_selected(i);
%     end
% end
% handles.new_module.files_in = char(handles.MIA_data.database.Filename(idex_patient & idex_tp & index_SequenceName));
% handles.new_module.files_in_index = find(idex_patient & idex_tp & index_SequenceName);
% handles.new_module.files_in_filter_data = table_data;

MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)
%MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
  %             handles.new_module.files_in_filter_name = {'Patient Name', '', 'Time Point','', 'Sequence Name',''};
        %             Patient_listing = unique(handles.MIA_data.database.Patient);
        %             Tp_listing = unique(handles.MIA_data.database.Tp);
        %             SequenceName_listing = unique(handles.MIA_data.database.SequenceName);
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
        %             set(handles.MIA_pipeline_parameter_setup,  'String', cellstr(handles.MIA_data.database.nii));

function output_file_names = MIA_pipeline_generate_file_name(handles, database_indexes, output_extention)

output_file_names = [...
    char(handles.MIA_data.database.Patient(database_indexes)) ...
    repmat('-', [numel(database_indexes),1])...
    char(handles.MIA_data.database.Tp(database_indexes))...
    repmat('-', [numel(database_indexes),1])...
    repmat(output_extention, [numel(database_indexes),1])...
    repmat('_', [numel(database_indexes),1])...
    repmat(datestr(now,'yyyymmdd-HHMMSSFFF'), [numel(database_indexes),1])...
    repmat('.nii', [numel(database_indexes),1])
    ] ;
output_file_names  = strrep(cellstr(output_file_names), ' ', '');


% --- Executes on selection change in MIA_pipeline_add_tag_popupmenu.
function MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)

% hObject    handle to MIA_pipeline_add_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Source_selected = handles.MIA_pipeline_add_tag_popupmenu.String{handles.MIA_pipeline_add_tag_popupmenu.Value};
%[val, ind] = max(contains(handles.MIA_pipeline_Filtered_Table.Properties.VariableNames, handles.Source_selected));
if strcmp(handles.Source_selected, 'NoMoreTags')
    handles.MIA_pipeline_Unique_Values_Tag.Data = {};
    handles.MIA_pipeline_Unique_Values_Tag.ColumnName = {};
else

    TagValues = unique(handles.MIA_pipeline_Filtered_Table{:,handles.Source_selected});
%TagValues = [{'all'} ; unique(handles.MIA_pipeline_Filtering_Table.Data(:,ind))];
%TagValues = [{'all'} ; cellstr(char(unique(handles.MIA_data.database{:,handles.Source_selected})))];

    handles.MIA_pipeline_Unique_Values_Tag.Data = cellstr(TagValues);
    handles.MIA_pipeline_Unique_Values_Tag.ColumnName = {handles.Source_selected};
end


guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_add_tag_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_add_tag_popupmenu


% --- Executes during object creation, after setting all properties.
function MIA_pipeline_add_tag_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_add_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');

end

function nii_json_fullfilename = fullfilename(handles, nii_index, ext)

nii_json_fullfilename = [char(handles.MIA_data.database.Path(nii_index)) char(handles.MIA_data.database.Filename(nii_index)) ext];


% --- Executes on selection change in MIA_pipeline_module_listbox.
function MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_module_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_module_listbox


% set(handles.MIA_pipeline_parameter_setup, 'Value', 1);
set(handles.MIA_pipeline_module_parameters, 'Value', 1);
module_selected = get(handles.MIA_pipeline_module_listbox, 'Value');
if isfield(handles, 'new_module')
    handles = rmfield(handles, 'new_module');
end

ismodule = 0;
switch char(handles.Modules_listing(module_selected))
    case handles.Module_groups	
        module_parameters_string = [char(handles.Modules_listing(module_selected)) ' modules'];
    case '   .SPM: Coreg (Est & Res)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Coreg_Est_Res('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Coreg_Est_Res(files_in,files_out,opt)';
        handles.new_module.module_name = 'Module_Coreg_Est_Res';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
     case '   .SPM: Coreg (Est)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Coreg_Est('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Coreg_Est(files_in,files_out,opt)';
        handles.new_module.module_name = 'Module_Coreg_Est';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;   
    case '   .T2map'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T2map('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_T2map(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_T2map';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;   
    case '   .ADCmap'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_ADCmap('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_ADCmap(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_ADCmap';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Fit_T2_T2star'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Fit_T2_T2star('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Fit_T2_T2star(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Fit_T2_T2star';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Smoothing'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Smoothing('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Smoothing(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Smoothing';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Dynamic Susceptibility Contrast'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Susceptibility('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Susceptibility(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Susceptibility';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .T1map (Multi Angles)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T1map_MultiAngles('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_T1map_MultiAngles(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_T1map_MultiAngles';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Arithmetic'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Arithmetic('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Arithmetic(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Arithmetic';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Inversion Efficiency (ASL_InvEff)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_ASL_InvEff('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_ASL_InvEff(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_ASL_InvEff';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .MGE2Dfrom3D'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_MGE2Dfrom3D('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_MGE2Dfrom3D(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_MGE2Dfrom3D';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .deltaR2'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_DeltaR2('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_DeltaR2(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_DeltaR2';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .deltaR2*'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_DeltaR2Star('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_DeltaR2Star(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_DeltaR2Star';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Dynamic Contrast Enhancement (Phenomenology)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_DCE_phenomeno('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_DCE_phenomeno(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_DCE_phenomeno';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .R2prim'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_R2Prim('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_R2Prim(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_R2Prim';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Blood volume fraction (steady-state)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_BVf('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_BVf(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_BVf';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Vessel Size Imaging (steady-state)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_VSI('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_VSI(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_VSI';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .Vessel Density (steady-state)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Density('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Density(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Density';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .SO2map'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_SO2('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_SO2(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_SO2';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    case '   .CMRO2'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_CMRO2('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_CMRO2(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_CMRO2';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
%     case '   .T1map (Multi Inversion Time)'
%         [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T1map_MIT('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = Module_T1map_MIT(char(files_in),files_out,opt)';
%         handles.new_module.module_name = 'Module_T1map_MIT';
%         module_parameters_string = handles.new_module.opt.table.Names_Display;
%         module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
%         ismodule = 1;
    case '   .Cerebral blood flow (ASL)'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_CBF('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_CBF(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_CBF';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;
        ismodule = 1;
    otherwise
        module_parameters_string = 'Not Implemented yet!!';    
        set(handles.MIA_pipeline_parameter_setup_text, 'String', '');

end

set(handles.MIA_pipeline_module_parameters, 'String', char(module_parameters_string));
if ismodule
    MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)
    %TableSize = handles.MIA_pipeline_module_parameters.Position(3);
    %FigureSizeInPix = handles.MIA_pipeline_manager_GUI.Position(3);
    %TableSizeInPix = TableSize*FigureSizeInPix;
    %DotsPerInch = get(0,'ScreenPixelsPerInch');
    %DotsPerChar = handles.MIA_pipeline_module_parameters.FontSize;
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
    set(handles.MIA_pipeline_parameter_setup_table, 'ColumnName', table.columnName);
    % set data (default's parameters)
    set(handles.MIA_pipeline_parameter_setup_table, 'Data', table.data);
    % set each colomn editable
    set(handles.MIA_pipeline_parameter_setup_table, 'columnEditable',  table.editable );
end
   

%% save the data
guidata(hObject, handles);
guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'), handles);


function [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles)
    ActualValues = cell(size(handles.module_parameters_string));
    StrToDisplay = cell(size(handles.module_parameters_string));
    LString = zeros(size(handles.module_parameters_string));
    for i=1:length(handles.module_parameters_fields)
        if strcmp(handles.module_parameters_fields{i}, '')
            ActualValues{i} = ' ';
        else
            if isnumeric(handles.new_module.opt.Module_settings.(handles.module_parameters_fields{i}))
                ActualValues{i} = num2str(handles.new_module.opt.Module_settings.(handles.module_parameters_fields{i}));
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
    
    
    %set(handles.MIA_pipeline_module_parameters, 'String', char(module_parameters_string));
    set(handles.MIA_pipeline_module_parameters, 'String', char(StrToDisplay));



% --- Executes during object creation, after setting all properties.
function MIA_pipeline_module_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function [pipeline, output_database] = MIA_pipeline_generate_psom_modules(New_module, FilterParameters, TmpDatabase, MIA_path)
% hObject    handle to MIA_pipeline_execute_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pipeline = struct();
Types = New_module.opt.table.Type;
ScanInputs = find(contains(Types, {'Scan', 'ROI'}));
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
            EmptyParams{i} = EmptyParams{i} +1;
        end
    end
    if ~exist('Databtmp', 'var')
        output_database = table();
        pipeline = struct();
        text = ['Please select at least one parameter (Input ',num2str(i),').'];
        warndlg(text)
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
        warndlg('You forgot to select a mandatory scan.', 'Missing Scan');
        pipeline = struct();
        output_database = table();
        return
    end
end



RefInput = New_module.opt.Module_settings.RefInput;
RefDatab = DatabaseInput{RefInput};
RefMat = MatricesInputs{RefInput};

UTag1 = unique(getfield(RefDatab, Tag1));
UTag2 = unique(getfield(RefDatab, Tag2));


FinalMat = cell(NbScanInput,1);
FinalMat{RefInput} = RefMat;

for i=1:NbScanInput
    if i~=RefInput
        Mattmp = cell(length(UTag1), length(UTag2));
        for j=1:length(UTag1)
            for k=1:length(UTag2)
                Databtmp = DatabaseInput{i};
                if ~isempty(Databtmp)
                    Databtmp2 = Databtmp(getfield(Databtmp, Tag1) == UTag1(j), :);
                    Databtmp3 = Databtmp2(getfield(Databtmp2, Tag2) == UTag2(k), :);
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
    InToReshape(find(all(cellfun(@isempty,InToReshape),2)),:) = [];
    
    %if size(InToReshape,1) == 1 && size(InToReshape,2) == 1 && EmptyParams{InputToReshape} == 0
    if EmptyParams{InputToReshape} == 0
        NewIn = repmat(InToReshape, size(RefMat));
        assert(size(NewIn,1) == size(RefMat,1));
        assert(size(NewIn,2) == size(RefMat,2));
        assert(size(NewIn,3) == size(RefMat,3));
        FinalMat{InputToReshape} = NewIn;
    elseif size(InToReshape,1) == 1 && size(InToReshape,2) == 1 && EmptyParams{InputToReshape} == 1
        string = ['Your ' Tag1 ' selection or your ' Tag2 ' selection is empty, but this multiple Tag is useless because it leads to a unique file. So, all the files from the input 2 will be coregistered on this unique file.'];
        choice = questdlg(string, 'Unique file detected','Continue', 'Return', 'Return');
        switch choice
            case 'Continue'
            NewIn = repmat(InToReshape, size(RefMat));
            assert(size(NewIn,1) == size(RefMat,1));
            assert(size(NewIn,2) == size(RefMat,2));
            assert(size(NewIn,3) == size(RefMat,3));
            FinalMat{InputToReshape} = NewIn;
            case 'Return'
                output_database = table();
                pipeline = struct();
                return
        end
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
                    assert(size(databtmp, 1) == 1);
                    InTags = databtmp(1,:);
                    table_in = [table_in; InTags];
                end
            end
        end
        New_module.opt.Module_settings.folder_out = [MIA_path, 'Tmp'];
        New_module.opt.Module_settings.Table_in = unique(table_in);
        pipeline = psom_add_job(pipeline, [New_module.module_name, num2str(i)], New_module.module_name, Files_in, '', New_module.opt.Module_settings);
        Mod_Struct = getfield(pipeline, [New_module.module_name, num2str(i)]);
        output_database = [output_database; Mod_Struct.opt.Table_out];
    end
end
output_database = unique(output_database);

    
    
    
    
    
    
    


function MIA_pipeline_execute_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_execute_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist([handles.MIA_data.database.Properties.UserData.MIA_data_path, 'PSOM'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MIA_data.database.Properties.UserData.MIA_data_path, 'PSOM']);
    if status == false
        error('Cannot create the PSOM folder to save the pipeline logs.')
    end
end
opt_pipe.path_logs = [handles.MIA_data.database.Properties.UserData.MIA_data_path,  'PSOM'];

if exist([handles.MIA_data.database.Properties.UserData.MIA_data_path,  'Derived_data'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MIA_data.database.Properties.UserData.MIA_data_path,  'Derived_data']);
    if status == false
        error('Cannot create the Derived_Data folder to save the results of the computed maps.')
    end
end


if exist([handles.MIA_data.database.Properties.UserData.MIA_data_path,  'Tmp'],'dir') ~= 7
    [status, ~, ~] = mkdir([handles.MIA_data.database.Properties.UserData.MIA_data_path,  'Tmp']);
    if status == false
        error('Cannot create the Tmp folder to temporarily save the results of the computed maps.')
    end
end


%% Create the pipeline from the modules.
Names_Mod = fieldnames(handles.MIA_pipeline_ParamsModules);
%Pipeline = handles.psom.Modules.(Names_Mod{1});
Pipeline = struct();
for i=1:length(Names_Mod)
    Mod = handles.MIA_pipeline_ParamsModules.(Names_Mod{i});
    %[new_pipeline, output_database] = MIA_pipeline_generate_psom_modules(Mod.Module, Mod.Filters, handles.MIA_pipeline_TmpDatabase, handles.MIA_data.database.Properties.UserData.MIA_data_path);
    Pipeline = smart_pipeline_merge(Pipeline, Mod.Jobs);
    %handles.MIA_pipeline_TmpDatabase = unique([handles.MIA_pipeline_TmpDatabase; Mod.OutputDatabase]);
end

handles.psom.pipeline = Pipeline;

% display the pipeline
if exist('biograph') == 2
    
    [graph_deps,list_jobs,files_in,files_out,files_clean] = psom_build_dependencies(handles.psom.pipeline);
    bg = biograph(graph_deps,list_jobs);
    
    
    % dolayout(bg);
    %% add editable functions to interact with the biograph
    set(bg, 'NodeCallbacks', @(hObject,eventdata)MIA_pipeline('node_callbacks',hObject));
    set(bg, 'EdgeCallbacks', @(hObject,eventdata)MIA_pipeline('edge_callbacks',hObject));
    view(bg) %, which will bring up the display in a different window.
    set(0, 'ShowHiddenHandles', 'on')
    
    handles.psom.biograph_fig = gcf;
    %set(handles.psom.biograph_ob, 'Name', 'MIA pipeline manager');
    guidata(hObject, handles);
    %return
    
end

%% exectute the pipeline
if handles.MIA_pipeline_radiobuttonPSOM.Value
    psom_run_pipeline(handles.psom.pipeline, opt_pipe)
    Result = load([opt_pipe.path_logs, '/PIPE_status_backup.mat']);
else
    Modules = fieldnames(handles.psom.pipeline);
    Result = struct();
    for i=1:length(Modules)
        Module = getfield(handles.psom.pipeline, Modules{i});
        files_in = Module.files_in;
        files_out = Module.files_out;
        opt = Module.opt;
        eval(Module.command)
        Result = setfield(Result, Modules{i}, 'finished');
        
    end
end

Jobs = fieldnames(handles.psom.pipeline);
update = false;
for i=1:length(Jobs)
   switch getfield(Result, Jobs{i})
       case 'failed'
           disp('FAILED')
       case 'finished'
           update = true;
           J = getfield(handles.psom.pipeline, Jobs{i});
           A = getfield(J, 'files_out');
           C = getfield(J, 'files_in');
           Outputs = fieldnames(A);
           for j=1:length(Outputs)
               B = getfield(A, Outputs{j});
               for k=1:length(B)
                   [path_out, name_out, ~] = fileparts(B{k});
                   path_out = [path_out, filesep];
                   outdb = handles.MIA_pipeline_TmpDatabase(handles.MIA_pipeline_TmpDatabase.Path == path_out, :);
                   outdb = outdb(outdb.Filename == name_out, :);
                   assert(height(outdb) == 1)
                   Folders = strsplit(path_out, filesep);
                   Folders{end-1} = 'Derived_data';
                   NewPath = strjoin(Folders, filesep);
                   [statusNii,~] = movefile(B{k}, [NewPath, name_out, '.nii']);
                   
                   [statusJson,~] = movefile(strrep(B{k},'.nii','.json'), [NewPath, name_out, '.json']);
                   if statusNii && statusJson
                       outdb.Path = NewPath;
                       %eval(['delete ' B{k}]);
                   elseif ~statusNii
                       warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data'' folder.', [name_out, '.nii'])
                   elseif ~statusJson
                       warning('Cannot move the file %s from the ''Tmp'' to the ''Derived_data'' folder.', [name_out, '.json'])
                   end
                   handles.MIA_data.database = unique([handles.MIA_data.database ; outdb]);
               end
           end
   end
end


% update MIA database if needed
if update
    handles2 = guidata(handles.MIA_data.MIA_GUI);
    handles2.database = handles.MIA_data.database;
    guidata(handles.MIA_data.MIA_GUI, handles2);

%handles2.MIA_GUI, handles
%guidata(handles.MIA_data, handles2);
%handles2 = guidata(handles.MIA_data.MIA_GUI);
%MIA_update_database_display(hObject, eventdata, handles)
%MIA('MIA_update_database_display'
%handles2.database = handles.MIA_data.database;

    MIA2('MIA_update_database_display', hObject, eventdata,handles.MIA_data)
    msgbox('Done', 'Information') ;

   close('MIA pipeline Manager')
end
%handles.MIA_pipeline_Filtered_Table = handles.MIA_data.database;
%MIA_pipeline_OpeningFcn(hObject, eventdata, handles)
%set(handles.MIA_pipeline.Filtering_Table, 'Data', 
%a=0;
function [hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles)

Datab = handles.MIA_data.database;
Names_Mod = fieldnames(handles.MIA_pipeline_ParamsModules);
for i=1:length(Names_Mod)
    Mod = handles.MIA_pipeline_ParamsModules.(Names_Mod{i});
    Datab = [Datab ; Mod.OutputDatabase];
end
handles.MIA_pipeline_TmpDatabase = unique(Datab);



% --- Executes on button press in MIA_pipeline_close_modules_button.
function MIA_pipeline_close_modules_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_close_modules_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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




% --- Executes on button press in MIA_pipeline_Add_Tag_Button.
function MIA_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Add_Tag_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if strcmp(handles.Source_selected, 'NoMoreTags')
%    return 
% end
%% Update of the Filtering Table
%if isempty(handles.MIA_pipeline_Filtering_Table.Data)
%    handles.MIA_pipeline_Filtering_Table.Data = table2cell(handles.MIA_data.database);
%end
if length(handles.MIA_pipeline_Unique_Values_Selection) <= 1
    return
end

handles.FilterParameters = [handles.FilterParameters, {handles.MIA_pipeline_Unique_Values_Selection}];

[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end


%% Update of the add tag popupmenu after an add
Tag_To_Add = handles.Source_selected;
Index = find(contains(handles.Add_list, Tag_To_Add));
NewTagListing = {handles.Add_list{1:Index-1},handles.Add_list{Index+1:end}};
[handles.Add_list, handles.Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);

set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);

guidata(hObject, handles);

% --- Executes on button press in MIA_pipeline_Remove_Tag_Button.
function MIA_pipeline_Remove_Tag_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Remove_Tag_Button (see GCBO)
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

[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end

if isfield(handles, 'new_module')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI'}))
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

set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
%handles.Source_selected = Tag_To_Remove;
handles.Remove_selected = handles.Remove_list{1};
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);


guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in MIA_pipeline_Filtering_Table.
function MIA_pipeline_Filtering_Table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Filtering_Table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% NbSelected = length(eventdata.Indices);
% for i = 1:NbSelected
%     NameSelected = eventdata.Source.Data{eventdata.Indices(i,1), eventdata.Indices(i,2)};
%     Tag = handles.MIA_pipeline_Filtering_Table.ColumnName{enventdata.Indices(i,2)};
% end

guidata(hObject, handles);


% --- Executes on button press in MIA_pipeline_pushMIASelection.
function MIA_pipeline_pushMIASelection_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_pushMIASelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% define which patient is selected
data_selected =  MIA2('get_data_selected',handles.MIA_data);
% add the patient filter
handles.FilterParameters{1} = {'Patient', char(handles.MIA_pipeline_TmpDatabase.Patient(data_selected(1)))};
handles.FilterParameters{2} = {'Tp', char(handles.MIA_pipeline_TmpDatabase.Tp(data_selected(1)))};

[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end

NewTagListing = {'Patient', 'Tp'};
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
guidata(hObject, handles);









% --- Executes on selection change in MIA_pipeline_remove_tag_popupmenu.
function MIA_pipeline_remove_tag_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_remove_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Remove_selected = handles.MIA_pipeline_remove_tag_popupmenu.String{handles.MIA_pipeline_remove_tag_popupmenu.Value};
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_remove_tag_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_remove_tag_popupmenu



% --- Executes on button press in MIA_pipeline_pushMIATPSelection.
function MIA_pipeline_pushMIATPSelection_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_pushMIATPSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% define which patient is selected
data_selected =  MIA2('get_data_selected',handles.MIA_data);
% add the patient filter
handles.FilterParameters{1} = {'Patient', char(handles.MIA_pipeline_TmpDatabase.Patient(data_selected(1)))};

[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end

NewTagListing = {'Patient'};
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
guidata(hObject, handles);


% --- Executes on button press in MIA_pipeline_push_Database.
function MIA_pipeline_push_Database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_push_Database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.FilterParameters = {};

[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles);
if isfield(handles, 'module_parameters_fields') && isfield(handles, 'module_parameters_string')
    [hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);
end


if isfield(handles, 'new_module')
    for i=1:length(handles.new_module.opt.table.Default)
        if any(strcmp(handles.new_module.opt.table.Type{i}, {'1Scan', 'XScan', '1ScanOrROI', 'XScanOrROI', '1ROI', 'XROI'}))
            handles.new_module.opt.table.Default{i} = [];
        end
    end
end

NewTagListing = handles.Add_Tags_listing;
[handles.Add_list, handles.Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);

% 
% %% Update of the add tag popupmenu after an add
% Tag_To_Add = handles.Source_selected;
% Index = find(contains(handles.Add_list, Tag_To_Add));
% NewTagListing = {handles.Add_list{1:Index-1},handles.Add_list{Index+1:end}};
% [handles.Add_list, handles.Remove_list]=UpdateAdd_Remove_Popup(NewTagListing, handles.Add_Tags_listing);
% 
% set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
% set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
% handles.Source_selected = handles.Add_list{1};
% handles.Remove_selected = handles.Remove_list{1};
% %handles.Remove_selected = Tag_To_Add;
% MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
% set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);


% handles.Tags_listing = handles.MIA_pipeline_TmpDatabase.Properties.VariableNames;
% set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Tags_listing);
% set(handles.MIA_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})
% 
% handles.Source_selected = handles.Tags_listing{1};
% 
% 
% handles.Remove_Tags_listing = {'NoMoreTags'};
% handles.Remove_selected = handles.Remove_Tags_listing{1};
% handles.MIA_pipeline_TmpDatabase.IsRaw = categorical(handles.MIA_pipeline_TmpDatabase.IsRaw);
% 
% handles.MIA_pipeline_Filtered_Table = handles.MIA_pipeline_TmpDatabase;
% handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_pipeline_Filtered_Table{:,handles.MIA_pipeline_TagsToPrint});
% handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;
% 
% 
% 
% handles.FilterParameters = {};
% 
% handles.MIA_pipeline_Unique_Values_Selection = {};
% MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function MIA_pipeline_remove_tag_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_remove_tag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in MIA_pipeline_Unique_Values_Tag.
function MIA_pipeline_Unique_Values_Tag_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Unique_Values_Tag (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
NbSelect = size(eventdata.Indices,1);
Names = cell(size(eventdata.Indices,1)+1,1);
Names{1,1}=handles.Source_selected;
for i=1:NbSelect
    Names{i+1,1} = eventdata.Source.Data{eventdata.Indices(i,1), eventdata.Indices(i,2)};
end
handles.MIA_pipeline_Unique_Values_Selection = Names;
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);


% --- Executes on button press in MIA_pipeline_Tags_To_Display_Button.
function MIA_pipeline_Tags_To_Display_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Tags_To_Display_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PreviousTagsIndex = find(contains(handles.MIA_pipeline_Filtered_Table.Properties.VariableNames, handles.MIA_pipeline_TagsToPrint));
[selection, ok] = listdlg('PromptString', 'Select the Tags to print among : ', 'InitialValue', PreviousTagsIndex(:), 'ListString', handles.MIA_pipeline_Filtered_Table.Properties.VariableNames);
if ok == 0;
   return 
end
handles.MIA_pipeline_TagsToPrint = handles.MIA_pipeline_Filtered_Table.Properties.VariableNames(selection);


% handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_pipeline_Filtered_Table{:,handles.MIA_pipeline_TagsToPrint});
% handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;

% handles.MIA_pipeline_Filtered_Table = NewTable;
% CellsColoured = DisplayColoredTable(NewTable, handles.MIA_pipeline_TagsToPrint);
% handles.MIA_pipeline_Filtering_Table.Data = CellsColoured;

[hObject, eventdata, handles] = MIA_pipeline_UpdateTables(hObject, eventdata, handles);

guidata(hObject, handles);


% --- Executes on selection change in MIA_pipeline_module_listbox.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_module_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_module_listbox


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_module_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MIA_pipeline_pipeline_listbox.
function MIA_pipeline_pipeline_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_pipeline_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_pipeline_pipeline_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_pipeline_pipeline_listbox


% --- Executes during object creation, after setting all properties.
function MIA_pipeline_pipeline_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_pipeline_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in MIA_pipeline_exectute_pipeline_button.
function MIA_pipeline_exectute_pipeline_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_exectute_pipeline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in MIA_pipeline_DeleteModule.
function MIA_pipeline_DeleteModule_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_DeleteModule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelectedIndex = handles.MIA_pipeline_pipeline_listbox.Value;
SelectedModule = handles.MIA_pipeline_pipeline_listbox.String{SelectedIndex};
handles.MIA_pipeline_ParamsModules = rmfield(handles.MIA_pipeline_ParamsModules, SelectedModule);
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
[hObject, eventdata, handles] = MIA_pipeline_UpdateTables(hObject, eventdata, handles);
set(handles.MIA_pipeline_pipeline_listbox,'String', fieldnames(handles.MIA_pipeline_ParamsModules));
guidata(hObject, handles);


% --- Executes on button press in MIA_pipeline_Edit_Module.
function MIA_pipeline_Edit_Module_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Edit_Module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelectedIndex = handles.MIA_pipeline_pipeline_listbox.Value;
SelectedModule = handles.MIA_pipeline_pipeline_listbox.String{SelectedIndex};
Module = handles.MIA_pipeline_ParamsModules.(SelectedModule);
handles.MIA_pipeline.EditedModuleName = SelectedModule;
handles.BeforeEditedModuleFilters = handles.FilterParameters;

%% Remove the module we want to edit, and adapt the database
handles.MIA_pipeline_ParamsModules = rmfield(handles.MIA_pipeline_ParamsModules, SelectedModule);
[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);
[hObject, eventdata, handles] = MIA_pipeline_UpdateTables(hObject, eventdata, handles);

handles.BeforeEditedModule = handles.new_module;
handles.new_module = Module.ModuleParams;

module_parameters_string = handles.new_module.opt.table.Names_Display;
module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;

%MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)

handles.module_parameters_string = module_parameters_string;
handles.module_parameters_fields = module_parameters_fields;
    
handles.FilterParameters = Module.Filters;
[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);

[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);


%% Add/Remove Popmenu update
TagsUsed = {};
for i=1:length(Module.Filters)
    TagsUsed = [TagsUsed, Module.Filters{i}(1)];
end
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(TagsUsed, handles.Add_Tags_listing);

set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);

handles.MIA_pipeline_module_parameters.Value = 1;
MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)




%% Update figure
set(handles.MIA_pipeline_module_parameters, 'BackgroundColor', [0.5 0.5 0.5]);
set(handles.MIA_pipeline_parameter_setup_table, 'BackgroundColor', [0.6 0.6 0.6;0.4 0.4 0.4]);
set(handles.MIA_pipeline_Unique_Values_Tag, 'BackgroundColor', [0.6 0.6 0.6;0.4 0.4 0.4]);
set(handles.MIA_pipeline_Filtering_Table, 'BackgroundColor', [0.6 0.6 0.6;0.4 0.4 0.4]);

set(handles.MIA_pipeline_clear_pipeline_button, 'Enable', 'off');
set(handles.MIA_pipeline_DeleteModule, 'Enable', 'off');
set(handles.MIA_pipeline_add_module_button, 'Enable', 'off');
set(handles.MIA_pipeline_execute_button, 'Enable', 'off');
set(handles.MIA_pipeline_module_listbox, 'Enable', 'off');
set(handles.MIA_pipeline_Edit_Module, 'Enable', 'off');
set(handles.MIA_pipeline_pipeline_listbox, 'Enable', 'off');

set(handles.MIA_pipeline_Save_Module, 'Enable', 'on');

guidata(hObject, handles);
%guidata(hObject, handles);


% --- Executes on button press in MIA_pipeline_Save_Module.
function MIA_pipeline_Save_Module_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Save_Module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[new_pipeline, output_database] = MIA_pipeline_generate_psom_modules(handles.new_module, handles.FilterParameters, handles.MIA_pipeline_TmpDatabase, handles.MIA_data.database.Properties.UserData.MIA_data_path);

if isempty(fieldnames(new_pipeline)) && isempty(output_database)
    return
end


SaveModule = struct();
SaveModule.Filters = handles.FilterParameters;
SaveModule.ModuleParams = handles.new_module;
SaveModule.OutputDatabase = output_database;
SaveModule.Jobs = new_pipeline;
handles.MIA_pipeline_ParamsModules.(handles.MIA_pipeline.EditedModuleName) = SaveModule;

if isfield(handles, 'new_module')
    handles = rmfield(handles, 'new_module');
end


[hObject, eventdata, handles] = UpdateTmpDatabase(hObject, eventdata, handles);

handles.FilterParameters = {};
[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);
%[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);

handles.new_module = handles.BeforeEditedModule;
handles.FilterParameters = handles.BeforeEditedModuleFilters;




%% Display the module we were creating before the Edit (Beggining)
module_parameters_string = handles.new_module.opt.table.Names_Display;
module_parameters_fields = handles.new_module.opt.table.PSOM_Fields;

%MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)

handles.module_parameters_string = module_parameters_string;
handles.module_parameters_fields = module_parameters_fields;
    
[hObject, eventdata, handles]=MIA_pipeline_UpdateTables(hObject, eventdata, handles);

[hObject, eventdata, handles] = UpdateParameters_listbox(hObject, eventdata, handles);


%% Add/Remove Popmenu update
TagsUsed = {};
for i=1:length(handles.FilterParameters)
    TagsUsed = [TagsUsed, handles.FilterParameters{i}(1)];
end
[handles.Remove_list, handles.Add_list]=UpdateAdd_Remove_Popup(TagsUsed, handles.Add_Tags_listing);

set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Add_list);
set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
handles.Source_selected = handles.Add_list{1};
handles.Remove_selected = handles.Remove_list{1};
%handles.Remove_selected = Tag_To_Add;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
set(handles.MIA_pipeline_remove_tag_popupmenu,'String',handles.Remove_list);
%% Display the module we were creating before the Edit (End)




%% Update figure
set(handles.MIA_pipeline_module_parameters, 'BackgroundColor', [0.94 0.94 0.94])
set(handles.MIA_pipeline_parameter_setup_table, 'BackgroundColor', [1 1 1;0.9412 0.9412 0.9412]);
set(handles.MIA_pipeline_Unique_Values_Tag, 'BackgroundColor', [1 1 1;0.9412 0.9412 0.9412]);
set(handles.MIA_pipeline_Filtering_Table, 'BackgroundColor', [1 1 1;0.9412 0.9412 0.9412]);
set(handles.MIA_pipeline_clear_pipeline_button, 'Enable', 'on');
set(handles.MIA_pipeline_DeleteModule, 'Enable', 'on');
set(handles.MIA_pipeline_add_module_button, 'Enable', 'on');
set(handles.MIA_pipeline_execute_button, 'Enable', 'on');
set(handles.MIA_pipeline_module_listbox, 'Enable', 'on');
set(handles.MIA_pipeline_Edit_Module, 'Enable', 'on');
set(handles.MIA_pipeline_pipeline_listbox, 'Enable', 'on');

set(handles.MIA_pipeline_Save_Module, 'Enable', 'off');

% MIA_pipeline_module_listbox_Callback(hObject, eventdata, handles);
