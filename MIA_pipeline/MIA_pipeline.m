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

% Last Modified by GUIDE v2.5 24-Jan-2018 11:20:20

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
handles.Modules_listing = {'Relaxometry', '   .T1map (Multi Inversion Time)', '   .T1map (Multi Angles)', '   .T2map', '   .T2*map',...
                '   .deltaR2', '   .deltaR2*',...
    'Perfusion', '   .Blood volume fraction (steady-state)', '   .Dynamic Susceptibility Contrast', '   .Vessel Size Imaging (steady-state)', ...
                '   .Vessel Densisty (steady-state)', '   .Cerebral blood flow (ASL)',  '   .Cerebral blood flow (ASL-Dynamic)',...
     'Permeability', '   .Dynamic Contrast Enhancement (Phenomenology)', '   .Dynamic Contrast Enhancement (Quantitative)',...          
     'Oxygenation', '   .R2prim', '   .SO2map', '   .CMRO2',...
     'MRFingerprint', '   .Vascular MRFingerprint'...
     'SPM', '   .SPM: Coreg', '   .SPM: Reslice','   .SPM: Realign', ...
     'Spatial', '   .Smoothing'...
     };
handles.Module_groups = {'Relaxometry','Perfusion', 'Permeability', 'Oxygenation', 'MRFingerprint', 'SPM', 'Spatial' };
 
 
 
%  
%  {'Arithmetic', 'Mean slices', 'Smooth', 'Add slices', ...
%    'SPM: Realign (Over time)', 'Same registration as', 'Normalization', 'Repair outlier',...
%     'Remove images', 'Shift images', 'Import Atlas (and ROI)','Export to Nifti'};
%     'T2starcorr3D', 'ASL_InvEff',

set(handles.MIA_pipeline_module_listbox, 'String', handles.Modules_listing);
handles.Tags_listing = handles.MIA_data.database.Properties.VariableNames;
set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Tags_listing);
set(handles.MIA_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})
handles.Source_selected = handles.Tags_listing{1};
%handles.MIA_pipeline_Unique_Values_Selection = 
handles.MIA_pipeline_TagsToPrint = {'Patient', 'Tp', 'SequenceName'};
handles.Remove_Tags_listing = {'NoMoreTags'};
handles.Remove_selected = handles.Remove_Tags_listing{1};
%handles.MIA_data.database.IsRaw = categorical(handles.MIA_data.database.IsRaw);

handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_data.database{:,handles.MIA_pipeline_TagsToPrint});
handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;

handles.MIA_pipeline_Filtered_Table = handles.MIA_data.database;
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
handles.FilterParameters = {};
%MIA_pipeline_Unique_Values_Tag_CellSelectionCallback(hObject, eventdata, handles)
% Update handles structure
handles.MIA_pipeline_Unique_Values_Selection = {};

% update the 'String' of MIA_pipeline_pushMIASelection and MIA_pipeline_pushMIATPSelection push button
data_selected =  MIA2('get_data_selected',handles.MIA_data);
set(handles.MIA_pipeline_pushMIASelection, 'String', [char(handles.MIA_data.database.Patient(data_selected)) '-' char(handles.MIA_data.database.Tp(data_selected)) ' only'])
set(handles.MIA_pipeline_pushMIATPSelection, 'String', ['All time point of :' char(handles.MIA_data.database.Patient(data_selected))])
guidata(hObject, handles);




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

% --- Executes on button press in MIA_pipeline_add_module_button.
function MIA_pipeline_add_module_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_add_module_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
module_selected = get(handles.MIA_pipeline_module_popupmenu,'Value');
module_list = get(handles.MIA_pipeline_module_popupmenu,'String');
if isfield(handles, 'biograph_fig')
    if ~isempty(findobj('Tag', 'BioGraphTool'))
        close(handles.biograph_fig);
    end
    handles = rmfield(handles, 'biograph_fig');
    handles = rmfield(handles, 'biograph_obj');
end
if isfield(handles, 'pipeline')
    if sum(strcmp(fieldnames(handles.pipeline), char(module_list(module_selected))))
        module_nbr = strfind(fieldnames(handles.pipeline), char(module_list(module_selected)));
        module_nbr = sum([module_nbr{:}]);
        module_name = [char(module_list(module_selected)), '_' num2str(module_nbr+1)];
        eval(['handles.pipeline.' module_name '= handles.new_module;']);      
    else
        eval(['handles.pipeline.' char(module_list(module_selected)) '= handles.new_module;']);
    end
else
    eval(['handles.pipeline.' char(module_list(module_selected)) '= handles.new_module;']);
    
end

% for i = 1:size(handles.new_module.files_in,1)
%     files_in.(['subject' num2str(i)]).file1 = handles.new_module.files_in(i,:);
% end
% pipeline = struct();
% % Get the list of subjects from files_in
% list_subject = fieldnames(files_in);
% 
% opt.folder_out = '/test/';
% 
% % [files_in,files_out,opt] = brick_name(files_in,opt.folder_out,opt);
% 
% % Loop over subjects
% for num_s = 1:length(list_subject)
%     % Plug the ?fmri? input files of the subjects in the job
%     job_in = files_in.(list_subject{num_s,:}).file1;
%     % Use the default output name
%     job_out = '';
%     % Force a specific folder organization for outputs
%     opt.test_brick.folder_out = [opt.folder_out list_subject{num_s} filesep];
%     % Give a name to the jobs
%     job_name = ['test_brick_' list_subject{num_s}];
%     % The name of the employed brick
%     brick = 'test_brick';
%     % Add the job to the pipeline
%     pipeline = psom_add_job(pipeline, job_name,brick,job_in,job_out,opt.test_brick);
%     % The outputs of this brick are just
%     % intermediate outputs :biograph_obj
%     % clean these up as soon as possible
% %     pipeline = psom_add_clean(pipeline, [job_name ...
% %         '_clean'],pipeline.(job_name).files_out);
% end

% display the pipeline
handles.biograph_obj = psom_visu_dependencies(handles.pipeline);
set(0, 'ShowHiddenHandles', 'on')
handles.biograph_fig = gcf;
set(handles.biograph_fig, 'Name', 'MIA pipeline manager');
guidata(hObject, handles);





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

switch handles.new_module.opt.table.Type{parameter_selected}
    case '1Scan1TPXP'
        if isempty(handles.new_module.opt.table.Default{parameter_selected})
            SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
            table.data(1:numel(SequenceType_listing),2) = {false};
            %SequenceType_listing = cellstr(SequenceType_listing)';
            %table.ColumnFormat = SequenceType_listing;
            %table.data(1,1) = {'Choose'};
            
            TP_listing = unique(handles.MIA_pipeline_Filtered_Table.Tp(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            %TP_listing = cellstr(TP_listing)';
            %table.ColumnFormat = {SequenceType_listing, TP_listing};
            table.data(1:numel(TP_listing),3) = cellstr(TP_listing);
            table.data(1:numel(TP_listing),4) = {false};
            %table.data(1,2) = {'Choose'};
            
            Patients_listing = unique(handles.MIA_pipeline_Filtered_Table.Patient(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
            table.data(1:numel(Patients_listing),5) = cellstr(Patients_listing);
            table.data(1:numel(Patients_listing),6) = {false};
        else
            table.data = handles.new_module.opt.table.Default{parameter_selected};
            %table.data = handles.new_module.opt.parameter_default{parameter_selected};
        end
        
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
            %table.data = handles.new_module.opt.parameter_default{parameter_selected};
        end

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
            %table.data = handles.new_module.opt.parameter_default{parameter_selected};
        end

        table.columnName = {'SequenceName', 'Select Input'};
        table.editable = [false true];
        table.ColumnFormat = {'char'};

    case 'cell'
        %Def_values = handles.new_module.opt.table.Default(parameter_selected);
        table.ColumnFormat = handles.new_module.opt.table.Default(parameter_selected);
        %set(handles.MIA_pipeline_parameter_setup_table, 'ColumnFormat', table.ColumnFormat);
        %Propositions = handles.new_module.opt.table.Default{parameter_selected};
        %table.data = Propositions(1);
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


% if iscell(handles.new_module.opt.parameter_type{parameter_selected})
%    OptionFields = fieldnames(handles.new_module.opt.Module_settings);
%    table.data(1:numel(handles.new_module.opt.parameter_type{parameter_selected}),1) = handles.new_module.opt.parameter_type{parameter_selected};
%    table.columnName = {OptionFields{parameter_selected}};
%    table.editable = false;
% else
% 
% 
%     switch handles.new_module.opt.parameter_type{parameter_selected}
%         case 'Scan'
%             if isempty(handles.new_module.opt.parameter_default{parameter_selected})
%                 SequenceType_listing = unique(handles.MIA_pipeline_Filtered_Table.SequenceName(handles.MIA_pipeline_Filtered_Table.Type == 'Scan'));
%                 table.data(1:numel(SequenceType_listing),1) = cellstr(SequenceType_listing);
%                 table.data(1:numel(SequenceType_listing),2) = {false};
%             else
%                 table.data = handles.new_module.opt.parameter_default{parameter_selected};
%             end
% 
%             table.columnName = {'Scan Type', 'Select Input'};
%             table.editable = [false true];
%             %table.Database = handles.MIA_data.database;
%             %else % case of selecting input data from an output module
%             %end
%         case 'char'
%             OptionFields = fieldnames(handles.new_module.opt.Module_settings);
%             %[val, ind] = max(endsWith(OptionFields, handles.MIA_pipeline_module_parameters.String{parameter_selected}, 'IgnoreCase', true));
%             %table.data(1,1) = cellstr(handles.MIA_pipeline_module_parameters.String{parameter_selected});
%             %table.data(1,2) = cellstr(getfield(handles.new_module.opt.Module_settings, OptionFields{parameter_selected}));
%             table.data(1,1) = cellstr(handles.new_module.opt.parameter_default{parameter_selected});
%             table.columnName = cellstr(handles.MIA_pipeline_module_parameters.String{parameter_selected});
% 
% 
%             table.editable = [true];
%             
%         case 'numeric'
%             OptionFields = fieldnames(handles.new_module.opt.Module_settings);
%             %[val, ind] = max(endsWith(OptionFields, handles.MIA_pipeline_module_parameters.String{parameter_selected}, 'IgnoreCase', true));
%             %table.data(1,1) = {getfield(handles.new_module.opt.Module_settings, OptionFields{parameter_selected})};
%             table.data(1,1) = {handles.new_module.opt.parameter_default{parameter_selected}};
%             %table.columnName = {OptionFields{parameter_selected}};
%             table.columnName = {handles.new_module.opt.parameter_list{parameter_selected}};
%             table.editable = true;
%             
%         otherwise
%             table.data = '';
%             table.columnName = '';
%             table.editable = false;
%     end
% end

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
    cellMaxLen = num2cell(maxLen*6);
    % Set ColumnWidth of UITABLE
    set(handles.MIA_pipeline_parameter_setup_table, 'ColumnWidth', cellMaxLen);
    
end










% parameter_selected = get(handles.MIA_pipeline_module_parameters,'Value');
% parameter_list = get(handles.MIA_pipeline_module_parameters,'String');
% set(handles.MIA_pipeline_parameter_setup, 'Value', 1);
% table_is_present = 0;
% % sub parameter --> find the ascending name
% if strfind(char(parameter_list(parameter_selected)), '  .')
%     for i=parameter_selected-1:-1:1
%         if isempty(strfind(char(parameter_list(i)), '  .'))
%             set(handles.MIA_pipeline_parameter_setup, 'String', eval(['handles.new_module.', char(parameter_list(i)), strrep(char(parameter_list(parameter_selected)), '  .', '.')]));
%             return
%         end
%     end
% elseif  ~isstruct(eval(['handles.new_module.', char(parameter_list(parameter_selected))]))
%     if strcmp(char(parameter_list(parameter_selected)), 'files_in')
%         table_is_present = 1;     
%         set(handles.MIA_pipeline_parameter_setup, 'String', handles.new_module.files_in)
%     elseif strcmp(char(parameter_list(parameter_selected)), 'files_out')
%         output_extention = handles.new_module.SequenceName;
%         handles.new_module.files_out = generate_file_name(handles, handles.new_module.files_in_index,  output_extention);
%         set(handles.MIA_pipeline_parameter_setup, 'String', handles.new_module.files_out)
%     else
%         set(handles.MIA_pipeline_parameter_setup, 'String', eval(['handles.new_module.',char(parameter_list(parameter_selected))]));
%     end
% else
%     set(handles.MIA_pipeline_parameter_setup, 'String', '');
% end
% if isfield(handles.new_module, 'files_in_filter_data') && table_is_present
%     % set names of the columns
%     set(handles.MIA_pipeline_parameter_setup_table, 'ColumnName', handles.new_module.files_in_filter_name);
%     % set data (default's parameters)
%     set(handles.MIA_pipeline_parameter_setup_table, 'Data', handles.new_module.files_in_filter_data);
%     % set columnFormat (option for each parameters)
%     set(handles.MIA_pipeline_parameter_setup_table, 'columnFormat', handles.new_module.files_in_filter_format);
%     % set ColumnWidth to auto
%     width= num2cell(max(cellfun(@length,[handles.new_module.files_in_filter_data(1,:)' handles.new_module.files_in_filter_name']')).*7);
%     set(handles.MIA_pipeline_parameter_setup_table,'ColumnWidth',width)
%     % set each colomn editable
%     set(handles.MIA_pipeline_parameter_setup_table, 'columnEditable', logical(handles.new_module.files_in_filter_editable));
%     
% else
%     set(handles.MIA_pipeline_parameter_setup_table, 'Data', '');
%     set(handles.MIA_pipeline_parameter_setup_table, 'ColumnName', '');
% 
% end


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


% --- Executes on button press in MIA_pipeline_clear_modules_button.
function MIA_pipeline_clear_modules_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_clear_modules_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'pipeline')
    handles = rmfield(handles, 'pipeline');
    if ~isempty(findobj('Tag', 'BioGraphTool'))
        close(handles.biograph_fig);
    end
    handles = rmfield(handles, 'biograph_fig');
    handles = rmfield(handles, 'biograph_obj');
end
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


function node_callbacks(hObject, eventdata, handles)
eventdata = [];
handles = guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'));
module_list = get(handles.MIA_pipeline_module_popupmenu, 'String');

sub_module = strfind(hObject.ID, '_');
if ~isempty(sub_module)
    module_name = hObject.ID(1:sub_module-1);
else
    module_name =hObject.ID;
end
idx = find(ismember(module_list, module_name));
set(handles.MIA_pipeline_module_popupmenu, 'Value', idx)

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
handles.new_module = handles.pipeline.(char(hObject.ID));
update_setting_windows(hObject, eventdata, handles)


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
if strcmp(handles.new_module.opt.table.Type{parameter_selected}, 'XScan')
    %handles.new_module.opt.parameter_default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data{cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,2)),1};
    handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
elseif strcmp(handles.new_module.opt.table.Type{parameter_selected}, '1Scan')
    if sum(cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,2))) == 1 || sum(cell2mat(handles.MIA_pipeline_parameter_setup_table.Data(:,2))) == 0
        handles.new_module.opt.table.Default{parameter_selected} = handles.MIA_pipeline_parameter_setup_table.Data;
    end
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

function output_file_names = generate_file_name(handles, database_indexes, output_extention)

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
%     {'Relaxometry', '   .T1map (Multi Inversion Time)', '   .T1map (Multi Angles)', '   .T2map', '   .T2*map',...
%                 '   .deltaR2', '   .deltaR2*',...
%     'Perfusion', '   .Blood volume fraction (steady-state)', '   .Vessel Size Imaging (steady-state)', ...
%                 '   .Vessel Densisty (steady-state)', '   .Cerebral blood flow (ASL)',  '   .Cerebral blood flow (ASL-Dynamic)',...
%      'Permeability', '   .Dynamic Contrast Enhancement (Phenomenology)', '   .Dynamic Contrast Enhancement (Quantitative)',...          
%      'Oxygenation', '   .R2prim', '   .SO2map', '   .CMRO2',...
%      'MRFingerprint', '   .Vascular MRFingerprint'...
%      'SPM', '   .SPM: Coreg', '   .SPM: Realign&Coreg', '   .SPM: Reslice','   .SPM: Realign', ...
%      };
    case handles.Module_groups	
        module_parameters_string = [char(handles.Modules_listing(module_selected)) ' modules'];
    case '   .SPM: Coreg'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Coreg('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Coreg(files_in,files_out,opt)';
        handles.new_module.module_name = 'Module_Coreg';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        ismodule = 1;
    case '   .T2map'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_T2map('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_T2map(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_T2map';
        %module_parameters_string = handles.new_module.opt.parameter_list;
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        ismodule = 1;
    case '   .Smoothing'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Smoothing('',  '', '');
%         handles.new_module.command = '[files_in,files_out,opt] = module_Smooting(char(pipeline.new_module.files_in),pipeline.new_module.files_out,pipeline.new_module.opt)';
        handles.new_module.command = '[files_in,files_out,opt] = Module_Smoothing(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Smoothing';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        %module_parameters_string = handles.new_module.opt.parameter_list;
        ismodule = 1;
    case '   .Dynamic Susceptibility Contrast'
        [handles.new_module.files_in ,handles.new_module.files_out ,handles.new_module.opt] = Module_Susceptibility('',  '', '');
        handles.new_module.command = '[files_in,files_out,opt] = Module_Susceptibility(char(files_in),files_out,opt)';
        handles.new_module.module_name = 'Module_Susceptibility';
        module_parameters_string = handles.new_module.opt.table.Names_Display;
        ismodule = 1;
        
        %         parameter_list = {'files_in', 'files_out', 'parameters', '  .DSC_parameter1', '  .DSC_parameter2'};
        
%         
%         handles.new_module.command = 'T2_module(files_in,files_out,parameter)';
%         handles.new_module.files_in_index = 1:size(handles.MIA_data.database.Patient);
%         handles.new_module.files_in = cellstr(handles.MIA_data.database.Filename);
%         handles.new_module.files_in_filter_name = {'Patient Name', '', 'Time Point','', 'Sequence Name',''};
%         Patient_listing = unique(handles.MIA_data.database.Patient);
%         Tp_listing = unique(handles.MIA_data.database.Tp);
%         SequenceName_listing = unique(handles.MIA_data.database.SequenceName);
%         handles.new_module.files_in_filter_data = {'all', false, 'all', false, char(SequenceName_listing(1)),false};
%         handles.new_module.files_in_filter_data(1:numel(Patient_listing)+1,1) = ['all' cellstr(Patient_listing)'];
%         handles.new_module.files_in_filter_data(1,2) = {true};
%         handles.new_module.files_in_filter_data(2:numel(Patient_listing)+1,2) = {false};
%         
%         handles.new_module.files_in_filter_data(1:numel(Tp_listing)+1,3) = ['all' cellstr(Tp_listing)'];
%         handles.new_module.files_in_filter_data(1,4) = {true};
%         handles.new_module.files_in_filter_data(2:numel(Patient_listing)+1,4) = {false};
%         
%         handles.new_module.files_in_filter_data(1:numel(SequenceName_listing),5) = cellstr(SequenceName_listing)';
%         handles.new_module.files_in_filter_data(1:numel(SequenceName_listing),6) = {false};
%         
%         handles.new_module.files_in_filter_format = {'char', 'logical','char', 'logical','char', 'logical' };
%         handles.new_module.files_in_filter_editable = [0 1 0 1 0 1];
%         
%         handles.new_module.SequenceName = 'T2map';
%         handles.new_module.files_out = generate_file_name(handles, handles.new_module.files_in_index,  handles.new_module.SequenceName);
%         handles.new_module.opt.threshold = 5;
%         handles.new_module.opt.flag_test =1;
%         parameter_list ={'files_in', 'files_out', 'parameters', '  .threshold'};
        
    otherwise
        module_parameters_string = 'Not Implemented yet!!';     
end
set(handles.MIA_pipeline_module_parameters, 'String', char(module_parameters_string));

if ismodule
    MIA_pipeline_module_parameters_Callback(hObject, eventdata, handles)
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
guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'), handles);






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


% --- Executes on button press in MIA_pipeline_exectute_module_button.
function MIA_pipeline_exectute_module_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_exectute_module_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist([handles.MIA_data.database.Properties.UserData.MIA_data_path, 'PSOM'],'dir') == 7
    opt_pipe.path_logs = [handles.MIA_data.database.Properties.UserData.MIA_data_path, 'PSOM'];
else
    [status, ~, ~] = mkdir([handles.MIA_data.database.Properties.UserData.MIA_data_path, 'PSOM']);
    if status == false
        error('Cannot create the PSOM folder to save the pipeline logs.')
    end
end
Types = handles.new_module.opt.table.Type;
ScanInputs = find(contains(Types, 'Scan'));
NbScanInput = length(ScanInputs);
% opt = struct();
% for i=1:size(handles.new_module.opt.table,1)
%     if isempty(handles.new_module.opt.table.PSOM_Fields{i})
%         
%     else
%         opt = setfield(opt, handles.new_module.opt.table.PSOM_Fields{i}, handles.new_module.opt.table.Default{i});
%     end
% end
%opt.flag_test = 0;
pipeline = struct();

switch NbScanInput
    case 1
        Scans = handles.new_module.opt.table.Default{ScanInputs};
        if isempty(Scans)
            return
        end
        ScansSelected = Scans(cell2mat(Scans(:,2)),1);
        NewTable = table();
        for i = 1:length(ScansSelected)
            NewTable = [NewTable; handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.SequenceName == ScansSelected{i},:)];
        end
        NewTable = unique(NewTable);
        
        NbJobs = size(NewTable,1);
        %pipeline = struct();
        for i=1:NbJobs
            Files_in.In1 = {[char(NewTable.Path(i)), char(NewTable.Filename(i)), '.nii']};
            pipeline = psom_add_job(pipeline, [handles.new_module.module_name, num2str(i)], handles.new_module.module_name, Files_in, '', handles.new_module.opt.Module_settings);
        end
    case 3
        %[Input1 Input2 Input3] = handles.new_module.opt.table.Default{ScanInputs};
        All_selected_Files = cell(3,1);
        for i=1:length(ScanInputs)
            Input = handles.new_module.opt.table.Default{ScanInputs(i)};
            NbParameters = size(Input,2)/2;
            FinalInputTable = {};
            %FinalInputTable = handles.MIA_pipeline_Filtered_Table;
            for j=1:NbParameters
                A = Input(:,2*j);
                A = A(~cellfun('isempty',A));
                ParamsSelected = Input(cell2mat(A),2*j-1);
                if i==1 && j == 3
                    PatientsInput1 = Input(cell2mat(A),5);
                end
                %NewTable = table();
                
                for k = 1:length(ParamsSelected)
                    %TagTableTest = [TagTableTest;NewTable(NewTable{:,Tag}==SelectedValue,:)];
                    %NewTable = [NewTable; FinalInputTable(getfield(FinalInputTable,handles.new_module.opt.ColumnNamesInput1Scan1TPXP{2*j-1}) == categorical(cellstr(ParamsSelected{k})),:)];
                    %NewTable = unique(NewTable);
                    %FinalInputTable{j,k} = ParamsSelected{k};
                    Selection{j,k,i} = ParamsSelected{k};
                end
                %FinalInputTable = NewTable;
                
            end
            %All_selected_Files{i} = FinalInputTable;
        end
        %[Input1 Input2 Input3] = handles.new_module.opt.table.Default{ScanInputs};
        %A = Input1(:,6);
        %A = A(~cellfun('isempty',A));
        %PatientsInput1 = Input1(cell2mat(A),5);
        %A = Input1(:,4);
        %A = A(~cellfun('isempty',A));
        %TPInput1 = Input1(cell2mat(A),3);
        %A = Input1(:,2);
        %A = A(~cellfun('isempty',A));
        %ScanInput1 = Input1(cell2mat(A),1);
        
        Compteur = 0;
        
        if length(PatientsInput1) == 1
            
            Input1Table = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.Patient == Selection(3,1,1),:);
            Input1Table = Input1Table(Input1Table.SequenceName == Selection(1,1,1),:);
            Input1Table = Input1Table(Input1Table.Tp == Selection(2,1,1),:);
            if isempty(Input1Table)
                disp('Impossible to find the file (Input1)');
            end
            File1 = [char(Input1Table.Path), char(Input1Table.Filename), '.nii'];
            %[char(NewTable.Path(i)), char(NewTable.Filename(i)), '.nii']
            
            Input2Table = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.SequenceName == Selection(1,1,2),:);
            PatientsInput2 = Input2Table.Patient;
            for i=1:length(PatientsInput2)
                Input2TableTmp = Input2Table(Input2Table.Patient == PatientsInput2(i), :);
                TpInput2 = Input2TableTmp.Tp;
                for j=1:length(TpInput2)
                    File2Try = [char(Input2TableTmp(Input2TableTmp.Tp == TpInput2(j),:).Path), char(Input2TableTmp(Input2TableTmp.Tp == TpInput2(j),:).Filename), '.nii'];
                    if ~strcmp(File2Try, File1)
                        File2 = File2Try;
                        Input3Table = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.Patient == PatientsInput2(i),:);
                        Input3Table = Input3Table(Input3Table.Tp == TpInput2(j),:);
                        if size(Selection,3) <3
                            File3 = {};
                        else
                            Scans = Selection(1,:,3);
                            Scans = Scans(~cellfun('isempty',Scans));
                            File3 = cell(length(Scans),1);
                            for k=1:length(Scans)
                                File3Test = [char(Input3Table(Input3Table.SequenceName == Scans{k},:).Path), char(Input3Table(Input3Table.SequenceName == Scans{k},:).Filename), '.nii'];
                                if sum(strcmp(File3Test, {File1, File2})) ~= 1
                                    File3{k,1} = File3Test;
                                end
                            end
                        end
                        %%% ADD JOBS HERE
                        Files_in.In1 = {File1};
                        Files_in.In2 = {File2};
                        Files_in.In3 = File3(~cellfun('isempty',File3));
                        pipeline = psom_add_job(pipeline, [handles.new_module.module_name, num2str(Compteur)], handles.new_module.module_name, Files_in, '', handles.new_module.opt.Module_settings);
                        Compteur = Compteur+1;
                    end
                end
            end
            
            
        else
            for i=1:length(PatientsInput1)
                Input1Table = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.Patient == Selection(3,i,1),:);
                Input1Table = Input1Table(Input1Table.SequenceName == Selection(1,1,1),:);
                Input1Table = Input1Table(Input1Table.Tp == Selection(2,1,1),:);
                if isempty(Input1Table)
                    disp('Impossible to find the file (Input1) for');
                    i
                end
                File1 = [char(Input1Table.Path), char(Input1Table.Filename), '.nii'];
                %[char(NewTable.Path(i)), char(NewTable.Filename(i)), '.nii']
                
                Input2Table = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.Patient == Selection(3,i,1),:);
                Input2Table = Input2Table(Input2Table.SequenceName == Selection(1,1,2),:);
                TpInput2 = Input2Table.Tp;
                for j=1:length(TpInput2)
                   File2Try = [char(Input2Table(Input2Table.Tp == TpInput2(j),:).Path), char(Input2Table(Input2Table.Tp == TpInput2(j),:).Filename), '.nii'];
                   if ~strcmp(File2Try, File1)
                       File2 = File2Try;

                       Input3Table = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.Patient == Selection(3,i,1),:);
                       Input3Table = Input3Table(Input3Table.Tp == TpInput2(j),:);
                       if size(Selection,3) <3
                           File3 = {};
                       else
                           Scans = Selection(1,:,3);
                           Scans = Scans(~cellfun('isempty',Scans));
                           File3 = cell(length(Scans),1);
                           for k=1:length(Scans)
                               File3Test = [char(Input3Table(Input3Table.SequenceName == Scans{k},:).Path), char(Input3Table(Input3Table.SequenceName == Scans{k},:).Filename), '.nii'];
                               if sum(strcmp(File3Test, {File1, File2})) ~= 1
                                   File3{k,1} = File3Test;
                               end
                           end
                       end                   
                       %%% ADD JOBS HERE
                       Files_in.In1 = {File1};
                       Files_in.In2 = {File2};
                       Files_in.In3 = File3(~cellfun('isempty',File3));
                       pipeline = psom_add_job(pipeline, [handles.new_module.module_name, num2str(Compteur)], handles.new_module.module_name, Files_in, '', handles.new_module.opt.Module_settings);
                       Compteur = Compteur+1;
                   end
                end
                
                
            end
        end
            
            
%             switch Types{ScanInputs(i)}
%                 case '1Scan'
%                     ScanSelected = Input(cell2mat(Input(:,2)),1);
%                     NewTable = handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.SequenceName == ScanSelected{i},:);
%                 case 'XScan'
%                     ScansSelected = Input(cell2mat(Input(:,2)),1);
%                     NewTable = table();
%                     for i = 1:length(ScansSelected)
%                         NewTable = [NewTable; handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.SequenceName == ScansSelected{i},:)];
%                     end
%                     NewTable = unique(NewTable);
%                 case '1Scan1TPXP'
%                     ScansSelected = Input(cell2mat(Input(:,2)),1);
%                     TpSelected = Input(cell2mat(Input(:,2)),1)
%                 
%             end
end

    
% %answer = inputdlg({'Where do you want to save the reporting files of the pipeline ?'}, 'Save reporting files', 1, {[pwd, '/LogsPSOM']});
% 
% 
% 
% % ScansSelected = Scans(cell2mat(Scans(:,2)),1);
% % NewTable = table();
% % for i = 1:length(ScansSelected)
% %     NewTable = [NewTable; handles.MIA_pipeline_Filtered_Table(handles.MIA_pipeline_Filtered_Table.SequenceName == ScansSelected{i},:)];
% % end
% % NewTable = unique(NewTable);
% 
% %extnii = repmat('.nii', size(NewTable, 1), 1);
% handles.new_module.files_in = struct();
% % for i=1:size(NewTable, 1)
% %     %handles.new_module.files_in{i,1} = [char(NewTable.Path(i)), char(NewTable.Filename(i)), '.nii'];
% %     handles.new_module.files_in = setfield(handles.new_module.files_in, ['Scan', num2str(i)], [char(NewTable.Path(i)), char(NewTable.Filename(i)), '.nii']);
% % end
% handles.new_module.files_out = '';
% 
% handles.new_module.opt2 = handles.new_module.opt;
% handles.new_module.opt = handles.new_module.opt2.Module_settings;
% handles.new_module.opt.flag_test = 0;
% %handles.new_module.opt.restart = {'new_module'};
% 
% 
% 
% for i=1:length(handles.new_module.opt2.parameter_link_psom)
%     NameGUI = handles.new_module.opt2.parameter_link_psom{i, 2};
%     NameField = handles.new_module.opt2.parameter_link_psom{i, 1};
%     Vec = contains(handles.new_module.opt2.parameter_list, NameGUI);
%     DefaultValue = handles.new_module.opt2.parameter_default{Vec};
%     handles.new_module.opt = setfield(handles.new_module.opt, NameField, DefaultValue);
% end
% %pipeline.new_module = handles.new_module;
% %list_fields = fieldnames(handles.new_module.files_in);
% pipeline = struct();
% for i=1:size(NewTable, 1)
%     pipeline = psom_add_job(pipeline, ['job_Scan', num2str(i)], handles.new_module.module_name, [char(NewTable.Path(i)), char(NewTable.Filename(i)), '.nii'], handles.new_module.files_out, handles.new_module.opt);
% end



psom_run_pipeline(pipeline, opt_pipe)

Result = load([opt_pipe.path_logs, '/PIPE_status_backup.mat']);
Jobs = fieldnames(pipeline);
update = false;
for i=1:length(Jobs)
   switch getfield(Result, Jobs{i})
       case 'failed'
           disp('FAILED')
       case 'finished'
           update = true;
           J = getfield(pipeline, Jobs{i});
           A = getfield(J, 'files_out');
           C = getfield(J, 'files_in');
           Outputs = fieldnames(A);
           for j=1:length(Outputs)
               B = getfield(A, Outputs{j});
               D = getfield(C, Outputs{j});
               if length(D) == length(B)
                   for k=1:length(B)
                       [~, name_in, ~] = fileparts(D{k});
                       Vec  = handles.MIA_data.database.Filename == name_in;
                       Tags_in = handles.MIA_data.database(Vec, :);
                       Tags_out = Tags_in;
                       [~, name_out, ~] = fileparts(B{k});
                       Tags_out.Filename = categorical(cellstr(name_out));
                       Tags_out.IsRaw = categorical(0);
                       if isfield(J.opt, 'output_filename_ext')
                           Tags_out.SequenceName = categorical(cellstr([char(Tags_in.SequenceName), '_', J.opt.output_filename_ext]));
                       elseif isfield(J.opt, 'output_filename_prefix')
                           Tags_out.SequenceName = categorical(cellstr([J.opt.output_filename_prefix, '_',char(Tags_in.SequenceName)]));
                       else
                           error('No output_filename_ext or output_filename_prefix')
                       end
                       Tags_out.IsRaw = double(Tags_out.IsRaw);
                       handles.MIA_data.database = unique([handles.MIA_data.database ; Tags_out]);

                   end
               else
                   for k=1:length(B)
                       [~, name_in, ~] = fileparts(D{1});
                       Vec  = handles.MIA_data.database.Filename == name_in;
                       Tags_in = handles.MIA_data.database(Vec, :);
                       Tags_out = Tags_in;
                       [~, name_out, ~] = fileparts(B{k});
                       Tags_out.Filename = categorical(cellstr(name_out));
                       Tags_out.IsRaw = categorical(0);
                       Spl = split(name_out, '_');
                       Type = Spl{end};
                       Tags_out.SequenceName = Type;
%                        if isfield(J.opt, 'output_filename_ext')
%                            Tags_out.SequenceName = categorical(cellstr([char(Tags_in.SequenceName), J.opt.output_filename_ext]));
%                        elseif isfield(J.opt, 'output_filename_prefix')
%                            Tags_out.SequenceName = categorical(cellstr([J.opt.output_filename_prefix, char(Tags_in.SequenceName)]));
%                        else
%                            error('No output_filename_ext or output_filename_prefix')
%                        end
                       Tags_out.IsRaw = double(Tags_out.IsRaw);
                       handles.MIA_data.database = unique([handles.MIA_data.database ; Tags_out]);

                   end
               end
           end
%            [path_in, name_in, ext_in] = fileparts(J.files_in);
%            Vec  = handles.MIA_data.database.Filename == name_in;
%            Tags_in = handles.MIA_data.database(Vec, :);
%            Tags_out = Tags_in;
%            [path_out, name_out, ext_out] = fileparts(J.files_out);
%            Tags_out.Filename = categorical(cellstr(name_out));
%            Tags_out.IsRaw = categorical(0);
%            Tags_out.SequenceName = categorical(cellstr([char(Tags_in.SequenceName), J.opt.output_filename_ext]));
%            handles.MIA_data.database = unique([handles.MIA_data.database ; Tags_out]);
   end
end
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
    close('MIA pipeline Manager')
end
%handles.MIA_pipeline_Filtered_Table = handles.MIA_data.database;
%MIA_pipeline_OpeningFcn(hObject, eventdata, handles)
%set(handles.MIA_pipeline.Filtering_Table, 'Data', 
%a=0;



% --- Executes on button press in MIA_pipeline_close_modules_button.
function MIA_pipeline_close_modules_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_close_modules_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in MIA_pipeline_Add_Tag_Button.
function MIA_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Add_Tag_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.Source_selected, 'NoMoreTags')
   return 
end
%% Update of the Filtering Table
%if isempty(handles.MIA_pipeline_Filtering_Table.Data)
%    handles.MIA_pipeline_Filtering_Table.Data = table2cell(handles.MIA_data.database);
%end
if length(handles.MIA_pipeline_Unique_Values_Selection) <= 1
    return
end

handles.FilterParameters = [handles.FilterParameters, {handles.MIA_pipeline_Unique_Values_Selection}];
NewTable = handles.MIA_data.database;
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
handles.MIA_pipeline_Filtering_Table.Data = cellstr(NewTable{:,handles.MIA_pipeline_TagsToPrint});
handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;
% NewTagValues = [{'all'} ; cellstr(char(unique(handles.MIA_data.database{:,handles.Source_selected})))];
% SizeVert = max(size(NewTagValues,1), size(handles.MIA_pipeline_Filtering_Table.Data,1));
% SizeHor = size(handles.MIA_pipeline_Filtering_Table.Data,2)+1;
% FullTable = cell(SizeVert,SizeHor);
% FullTable(1:size(handles.MIA_pipeline_Filtering_Table.Data,1),1:size(handles.MIA_pipeline_Filtering_Table.Data,2)) = handles.MIA_pipeline_Filtering_Table.Data;
% FullTable(1:size(NewTagValues,1),size(handles.MIA_pipeline_Filtering_Table.Data,2)+1) = NewTagValues;
% handles.MIA_pipeline_Filtering_Table.Data = FullTable;
% handles.MIA_pipeline_Filtering_Table.ColumnName = [handles.MIA_pipeline_Filtering_Table.ColumnName; {handles.Source_selected}];
%% Update of the add tag popupmenu after an add
Tag_To_Add = handles.Source_selected;
Index = find(contains(handles.Tags_listing, Tag_To_Add));
NewTagListing = {handles.Tags_listing{1:Index-1},handles.Tags_listing{Index+1:end}};
if isempty(NewTagListing)
    set(handles.MIA_pipeline_add_tag_popupmenu, 'String', {'NoMoreTags'})
    set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
    handles.Source_selected = {'NoMoreTags'};
    handles.Tags_listing = {'NoMoreTags'};

else
    set(handles.MIA_pipeline_add_tag_popupmenu, 'String', NewTagListing);
    set(handles.MIA_pipeline_add_tag_popupmenu, 'Value', 1);
    handles.Source_selected = NewTagListing{1};
    handles.Tags_listing = NewTagListing;
    MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
end
%% Update of the remove tag popupmenu after an add
if contains(handles.MIA_pipeline_remove_tag_popupmenu.String, 'NoMoreTags')
    set(handles.MIA_pipeline_remove_tag_popupmenu,'String', {Tag_To_Add});
    handles.Remove_selected = Tag_To_Add;
else
%handles.Remove_Tag_Listing = {handles.Remove_Tag_Listing, Tag_To_Add};
    if size(handles.MIA_pipeline_remove_tag_popupmenu.String,1)==1
        set(handles.MIA_pipeline_remove_tag_popupmenu,'String',[handles.MIA_pipeline_remove_tag_popupmenu.String; {Tag_To_Add}]);
        handles.Remove_Tags_listing = handles.MIA_pipeline_remove_tag_popupmenu.String;
    else
        set(handles.MIA_pipeline_remove_tag_popupmenu,'String',[handles.MIA_pipeline_remove_tag_popupmenu.String; Tag_To_Add]);
        handles.Remove_Tags_listing = handles.MIA_pipeline_remove_tag_popupmenu.String;
    end
end
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
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

NewTable = handles.MIA_data.database;
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
handles.MIA_pipeline_Filtering_Table.Data = cellstr(NewTable{:,handles.MIA_pipeline_TagsToPrint});
handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;

% Index = find(contains(handles.MIA_pipeline_Filtering_Table.ColumnName, handles.Remove_selected));
% handles.MIA_pipeline_Filtering_Table.Data(:,Index) = [];
% [val, index]=min(sum(cellfun(@isempty,handles.MIA_pipeline_Filtering_Table.Data)));
% handles.MIA_pipeline_Filtering_Table.Data = handles.MIA_pipeline_Filtering_Table.Data(~cellfun(@isempty, handles.MIA_pipeline_Filtering_Table.Data(:,index)), :);
% %handles.MIA_pipeline_Filtering_Table.Data = handles.MIA_pipeline_Filtering_Table.Data(~cellfun('isempty',handles.MIA_pipeline_Filtering_Table.Data));
% handles.MIA_pipeline_Filtering_Table.ColumnName(Index) = [];
%% Update of the remove tag popup menu after a remove
Tag_To_Add = handles.Remove_selected;
Index = find(contains(handles.Remove_Tags_listing, Tag_To_Add));
NewTagListing = {handles.Remove_Tags_listing{1:Index-1},handles.Remove_Tags_listing{Index+1:end}};
if isempty(NewTagListing)
    set(handles.MIA_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})
    set(handles.MIA_pipeline_remove_tag_popupmenu, 'Value', 1);
    handles.Remove_selected = {'NoMoreTags'};
    handles.Remove_Tags_listing = {'NoMoreTags'};
else
    set(handles.MIA_pipeline_remove_tag_popupmenu, 'String', NewTagListing);
    set(handles.MIA_pipeline_remove_tag_popupmenu, 'Value', 1);
    handles.Remove_selected = NewTagListing{1};
    handles.Remove_Tags_listing = NewTagListing;
    MIA_pipeline_remove_tag_popupmenu_Callback(hObject, eventdata, handles)
end
%% Update of the add tag popup menu after a remove
if contains(handles.MIA_pipeline_add_tag_popupmenu.String, 'NoMoreTags')
    set(handles.MIA_pipeline_add_tag_popupmenu,'String', {Tag_To_Add});
else
%handles.Remove_Tag_Listing = {handles.Remove_Tag_Listing, Tag_To_Add};
    if size(handles.MIA_pipeline_add_tag_popupmenu.String,1)==1
        set(handles.MIA_pipeline_add_tag_popupmenu,'String',[handles.MIA_pipeline_add_tag_popupmenu.String; {Tag_To_Add}]);
        handles.Tags_listing = handles.MIA_pipeline_add_tag_popupmenu.String;
    else
        set(handles.MIA_pipeline_add_tag_popupmenu,'String',[handles.MIA_pipeline_add_tag_popupmenu.String; Tag_To_Add]);
        handles.Tags_listing = handles.MIA_pipeline_add_tag_popupmenu.String;
    end
end
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
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

% If there is an filter loaded remove them all  
if ~isempty(handles.Remove_Tags_listing)
    MIA_pipeline_push_Database_Callback(hObject, eventdata, handles)
    handles = guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'));
end
% define which patient is selected
data_selected =  MIA2('get_data_selected',handles.MIA_data);
% add the patient filter
handles.Source_selected = 'Patient';
handles.MIA_pipeline_Unique_Values_Selection= {'Patient', char(handles.MIA_data.database.Patient(data_selected))}; 
MIA_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)
% retrieve UI data
handles = guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'));

% Add the time point filter
handles.MIA_pipeline_Unique_Values_Selection= {'Tp', char(handles.MIA_data.database.Tp(data_selected))}; 
handles.Source_selected = 'Tp';
MIA_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)

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
if ~isempty(handles.Remove_Tags_listing)
    MIA_pipeline_push_Database_Callback(hObject, eventdata, handles)
    handles = guidata(findobj('Tag', 'MIA_pipeline_manager_GUI'));
end
% define which patient is selected
data_selected =  MIA2('get_data_selected',handles.MIA_data);
% add the patient filter
handles.Source_selected = 'Patient';
handles.MIA_pipeline_Unique_Values_Selection= {'Patient', char(handles.MIA_data.database.Patient(data_selected))}; 
MIA_pipeline_Add_Tag_Button_Callback(hObject, eventdata, handles)


% --- Executes on button press in MIA_pipeline_push_Database.
function MIA_pipeline_push_Database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_push_Database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.Tags_listing = handles.MIA_data.database.Properties.VariableNames;
set(handles.MIA_pipeline_add_tag_popupmenu, 'String', handles.Tags_listing);
set(handles.MIA_pipeline_remove_tag_popupmenu, 'String', {'NoMoreTags'})

handles.Source_selected = handles.Tags_listing{1};


handles.Remove_Tags_listing = {'NoMoreTags'};
handles.Remove_selected = handles.Remove_Tags_listing{1};
handles.MIA_data.database.IsRaw = categorical(handles.MIA_data.database.IsRaw);

handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_data.database{:,handles.MIA_pipeline_TagsToPrint});
handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;

handles.MIA_pipeline_Filtered_Table = handles.MIA_data.database;

handles.FilterParameters = {};

handles.MIA_pipeline_Unique_Values_Selection = {};
MIA_pipeline_add_tag_popupmenu_Callback(hObject, eventdata, handles)
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

handles.MIA_pipeline_Filtering_Table.Data = cellstr(handles.MIA_pipeline_Filtered_Table{:,handles.MIA_pipeline_TagsToPrint});
handles.MIA_pipeline_Filtering_Table.ColumnName = handles.MIA_pipeline_TagsToPrint;
guidata(hObject, handles);
