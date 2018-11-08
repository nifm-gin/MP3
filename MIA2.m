function varargout = MIA2(varargin)
% MIA2 MATLAB code for MIA2.fig
%      MIA2, by itself, creates a new MIA2 or raises the existing
%      singleton*.
%
%      H = MIA2 returns the handle to a new MIA2 or the handle to
%      the existing singleton*.
%
%      MIA2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIA2.M with the given input arguments.
%
%      MIA2('Property','Value',...) creates a new MIA2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MIA2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MIA2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MIA2


% Last Modified by GUIDE v2.5 18-Oct-2018 15:52:27


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MIA2_OpeningFcn, ...
    'gui_OutputFcn',  @MIA2_OutputFcn, ...
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


% --- Executes just before MIA2 is made visible.
function MIA2_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MIA2 (see VARARGIN)

% Choose default command line output for MIA2
handles.output = hObject;

% init stuct/variables
handles.VOIs = {'Other'};
handles.histo = {'Other'}; %{'Pimo', 'ColIV', 'Tc+I_cerveau', 'Tc+I_ref'};
handles.resolution = [1 64 112 128 192 256 384 512 3000];
handles.colors ={'b', 'g', 'm', 'c', 'r', 'k', 'y', 'navy',...
    'u1','turquoise','slateblue',	'springgreen',	'maroon',...
    'purple',	'u2',	'olive',	'u3','chartreuse',	'u4',	'sky',...
    'u5',	'orange',	'u6',	'u7',	'u8',	'gray'};
handles.colors_rgb = [0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 0 0; 1 1 1];
%load('rgb_color_table.mat', 'num');
%handles.colors_rgb = num;
handles.colormap = get(handles.MIA_colormap_popupmenu,'String');
handles.markers ={'o','s', 'd', 'p', 'h', '+', '*', 'x'};
table_data(1,1) = {'Voxel values'};
set(handles.MIA_table_pixel_values, 'Data', table_data);
set(handles.MIA_table1, 'Data', {'', '', '', '', ''});
handles.table1.cluster = [];
handles.table1.cluster_row = [];
handles.mode = 1;
handles.view_mode = 'Axial';

handles.display_option.view_pixel_on_map = 0;
handles.display_option.view_pixel_on_plot = 0;
handles.display_option.view_plot = 1;
handles.display_option.manual_contrast = 0;
set(handles.MIA_menu_view_plot, 'Check', 'on');

for i=1:4
    stri = num2str(i);
    set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
    set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
end


% add MRIManager.jar to the classpath (dynamic classpath)
[filepath,name,ext] = fileparts(which('MRIManager.jar'));
javaclasspath(fullfile(filepath,[name,ext]));
% save the java skin used

handles.original_Java_LookAndFeel = javax.swing.UIManager.getLookAndFeel;
handles.original_Java_LookAndFeel = sprintf('%s',handles.original_Java_LookAndFeel);
handles.original_Java_LookAndFeel = split(handles.original_Java_LookAndFeel,' - ');
handles.original_Java_LookAndFeel = extractBefore(handles.original_Java_LookAndFeel{end},']');
handles.original_Java_LookAndFeel = ['[LookAndFeel] ',handles.original_Java_LookAndFeel];

% A = javax.swing.UIManager.getLookAndFeel;
% char(A.getClass) % 'class com.jgoodies.looks.plastic.Plastic3DLookAndFeel' Replace class by [LookAndFeel]


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MIA2 wait for user response (see UIRESUME)
% uiwait(handles.MIA_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = MIA2_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes when user attempts to close MIA_GUI.
function MIA_GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MIA_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
selection = questdlg('Do you want to close (did you save you database)?',...
    'Warning',...
    'Yes','No','Yes');
if isempty(selection)
    return
end
switch selection
    case 'No'
        return
end

if ~isempty(findobj('type', 'figure', 'name', 'MIA pipeline Manager'))
    close((findobj('type', 'figure', 'name', 'MIA pipeline Manager')));
end



delete(hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MIA_slider_slice.
function MIA_slider_slice_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MIA_slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if ~isfield(handles, 'data_loaded') && ~isfield(handles, 'data_selected_for_PRM')
    return
end
Slice_min = get(hObject,'Min');
Slice_max = get(hObject,'Max');
Position = get(hObject,'Position');

SliderBarWidth = Position(3)/Slice_max;
set(hObject,'UserData',[Slice_min Slice_max SliderBarWidth]);

cp = get(handles.MIA_GUI,'CurrentPoint');
newValue = round((cp(1,1))/SliderBarWidth);
if  newValue ==  get(handles.MIA_slider_slice,'Value')
    return
elseif newValue > Slice_max
    newValue = Slice_max;
elseif newValue < Slice_min
    newValue = Slice_min;
end
set(handles.MIA_slider_slice,'Value',newValue);

MIA_update_axes(hObject, eventdata, handles)
%
% set(handles.MIA_GUI,'WindowButtonMotionFcn',{@MIA_slider_on_move,handles})
% set(handles.MIA_GUI,'WindowButtonUpFcn',{@MIA_slider_release_click,handles})


% --- Executes during object creation, after setting all properties.
function MIA_slider_slice_CreateFcn(hObject, ~, ~)
% hObject    handle to MIA_slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in MIA_name_list.
function MIA_name_list_Callback(hObject,  eventdata, ~)
% hObject    handle to MIA_name_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_name_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_name_list
handles = guidata(hObject);
if ~isfield(handles, 'database')
    return
end
set(handles.MIA_time_points_list, 'Value', 1);
set(handles.MIA_scans_list, 'Value', 1);
set(handles.MIA_file_list, 'Value', 1);
MIA_update_database_display(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function MIA_name_list_CreateFcn(hObject, ~, ~)
% hObject    handle to MIA_name_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in MIA_scans_list.
function MIA_scans_list_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_scans_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_scans_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_scans_list

if ~isfield(handles, 'database')
    return
end
if numel(get(handles.MIA_name_list, 'Value')) >1 || numel(get(handles.MIA_time_points_list, 'Value')) > 1
    return
end
set(handles.MIA_file_list, 'Value', 1);
guidata(hObject, handles);

MIA_update_database_display(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function MIA_scans_list_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
% hObject    handle to MIA_scans_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MIA_update_database_display(hObject, eventdata, handles)
% handles = guidata(gcf, findobj('Tag', 'MIA_GUI'));
handles = guidata(handles.MIA_GUI);
if ~isfield(handles, 'database')
    return
end
if isempty(handles.database)
    % if the table is empty, clear all lists
    set(handles.MIA_name_list, 'String', '', 'Value', 1);
    set(handles.MIA_time_points_list, 'String', '', 'Value', 1);
    set(handles.MIA_scans_list, 'String', '', 'Value', 1);
    set(handles.MIA_file_list, 'String', '', 'Value', 1);
    return
end

patient_id = get(handles.MIA_name_list, 'Value');

id_listing = unique(handles.database.Patient,'stable');
set(handles.MIA_name_list, 'String', char(id_listing));
if numel(patient_id)~= 1
    return
end

Patient_filter = handles.database.Patient== id_listing(patient_id);
tp_listing = unique(handles.database.Tp(Patient_filter),'stable');
% check if the new time point listing is not shorter than the old one. If
% so update MIA_time_points_list 'Value'
if numel(tp_listing) < get(handles.MIA_time_points_list, 'Value')
    set(handles.MIA_time_points_list, 'String', char(tp_listing), 'Value', numel(tp_listing));
else
    set(handles.MIA_time_points_list, 'String', char(tp_listing));
end
time_point = get(handles.MIA_time_points_list, 'Value');

if get(handles.MIA_scan_VOIs_button, 'Value') == 0 %display parameters list
    is_scan =  handles.database.Type == 'Scan';
    tp_filter = handles.database.Tp== tp_listing(time_point);
    sequence_listing = handles.database.SequenceName(Patient_filter & tp_filter & is_scan);
    
    % check if the sequence listing is not shorter than the old one. If
    % so update MIA_scans_list 'Value'
    if numel(sequence_listing) < get(handles.MIA_scans_list, 'Value')
        set(handles.MIA_scans_list, 'String', char(sequence_listing), 'Value', numel(sequence_listing));
    else
        set(handles.MIA_scans_list, 'String', char(sequence_listing));
    end
    scan = get(handles.MIA_scans_list, 'Value');
    
    file_text= cell(1, numel(sequence_listing(scan)));
    for i=1:numel(sequence_listing(scan))
        sequence_filter =  handles.database.SequenceName== sequence_listing(scan(i));
        file_text(i) = cellstr(handles.database.Filename(Patient_filter & tp_filter & sequence_filter & is_scan));
        
    end
    set(handles.MIA_file_list, 'String', file_text);
    
else %display VOIs list
    %is_ROI=  handles.database.Type == 'ROI';
    is_ROI = handles.database.Type == 'ROI' | handles.database.Type == 'Cluster';
    tp_filter = handles.database.Tp== tp_listing(time_point);
    sequence_listing = handles.database.SequenceName(Patient_filter & tp_filter & is_ROI);
    if isempty(sequence_listing)
        set(handles.MIA_scans_list, 'String', '');
        return
    end
    scan = get(handles.MIA_scans_list, 'Value');
    set(handles.MIA_scans_list, 'String', char(sequence_listing));
    file_text= cell(1, numel(sequence_listing(scan)));
    for i=1:numel(sequence_listing(scan))
        sequence_filter =  handles.database.SequenceName== sequence_listing(scan(i));
        file_text(i) = cellstr(handles.database.Filename(Patient_filter & tp_filter & sequence_filter & is_ROI));
    end
    set(handles.MIA_file_list, 'String', file_text);
    
    
end
% if the pipeline Manager is open, update the information : patient selected
% update the 'String' of MIA_pipeline_pushMIASelection and MIA_pipeline_pushMIATPSelection push button




if ~isempty(findobj('type', 'figure', 'name', 'MIA pipeline Manager'))
    % Get the hObject of MIA_pipeline
    h = findobj('Tag', 'MIA_pipeline_manager_GUI');
    % Get the handles of MIA_pipeline
    data = guidata(h);
    % Update the handles of MIA_pipeline by stocking the latest version of
    % MIA handles.
    data.MIA_data = handles;

    % Don't touch the original eventdata, just in case.
    eventdata2 = eventdata;
    %Update the MIA_pipeline tmp_database
    [h, ~, data] = MIA_pipeline('UpdateTmpDatabase', h, eventdata2, data);
    [~, ~, data] = MIA_pipeline('MIA_pipeline_UpdateTables', h, eventdata2, data);
    clear('eventdata2')
    guidata(h, data)
end



if ~isempty(findobj('Tag', 'MIA_pipeline_pushMIASelection'))
    data_selected = finddata_selected(handles);
    if size(char(handles.database.Patient(data_selected)),1) > 1
        return
    else
        set(findobj('Tag', 'MIA_pipeline_pushMIASelection'), 'String', [char(handles.database.Patient(data_selected(1))) '-' char(handles.database.Tp(data_selected(1))) ' only'])
        set(findobj('Tag', 'MIA_pipeline_pushMIATPSelection'), 'String', ['All time point of :' char(handles.database.Patient(data_selected(1)))])
    end
end


% --- Executes on selection change in MIA_time_points_list.
function MIA_time_points_list_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_time_points_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_time_points_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_time_points_list
if ~isfield(handles, 'database')
    return
end
set(handles.MIA_scans_list, 'Value', 1);
set(handles.MIA_file_list, 'Value', 1);
guidata(hObject, handles);
if numel(get(handles.MIA_time_points_list, 'Value')) >1 ||...
        numel(get(handles.MIA_name_list, 'Value')) >1
else
    MIA_update_database_display(hObject, eventdata, handles)
end



% --- Executes during object creation, after setting all properties.
function MIA_time_points_list_CreateFcn(hObject, ~, ~)
% hObject    handle to MIA_time_points_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in MIA_scan_VOIs_button.
function MIA_scan_VOIs_button_Callback(hObject, eventdata, ~)
% hObject    handle to MIA_scan_VOIs_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MIA_scan_VOIs_button

handles = guidata(hObject);
if ~isfield(handles, 'database')
    return
end
patient = get(handles.MIA_name_list, 'Value');
timepoint = get(handles.MIA_time_points_list, 'Value');
if numel(patient) > 1 || numel(timepoint) >1
    set(handles.MIA_scan_VOIs_button, 'Value', 0)
    return
end

value = get(handles.MIA_scan_VOIs_button, 'Value');
set(handles.MIA_scans_list, 'Value', 1);
if value == 0
    set(handles.MIA_scan_VOIs_button, 'String', 'Scans');
else
    set(handles.MIA_scan_VOIs_button, 'String', 'VOIs');
end
MIA_update_database_display(hObject, eventdata, handles);




% --- Executes on selection change in MIA_file_list.
function MIA_file_list_Callback(~, ~, ~)
% hObject    handle to MIA_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_file_list


% --- Executes during object creation, after setting all properties.
function MIA_file_list_CreateFcn(hObject, ~, ~)
% hObject    handle to MIA_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function MIA_rename_name_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_rename_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
data_selected = finddata_selected(handles);
if numel(data_selected) >1
    warndlg('Please select only one Patient', 'Warning');
    return
end

name_option = [cellstr(unique(handles.database.Patient(handles.database.Type == 'Scan')))' 'Other']';

[new_Patient_name, ok1] = listdlg('PromptString','Select the new scan name:',...
    'Name', 'Select a Name',...
    'SelectionMode','single',...
    'ListSize', [400 300],...
    'ListString',name_option);

if ok1 == 0
    return
end
if strcmp('Other',name_option(new_Patient_name)) == 1
    NewPatient = inputdlg('Name of the new Scan ', 'Question?', 1, {''});
else
    NewPatient =name_option(new_Patient_name);
end

%% update the database with the new name
% but first check if the new scan name does not exist for this patient and
% time point

if find(handles.database.Patient == handles.database.Patient(data_selected) &...
        handles.database.Tp == handles.database.Tp(data_selected) & ...
        handles.database.SequenceName == NewPatient) > 0
    msgbox('A Scan with the same name already exist for this patient at this time point') ;
    return
end

idx_scan_to_rename = find(handles.database.Patient == handles.database.Patient(data_selected));
for i=1:numel(idx_scan_to_rename)
    new_nii_filename = strrep(cellstr(handles.database.Filename(idx_scan_to_rename(i))), cellstr(handles.database.Patient(idx_scan_to_rename(i))), NewPatient);
    
    % rename the scan file
    if  exist(fullfilename(handles, idx_scan_to_rename(i), '.nii'), 'file') == 0
        warning_text = sprintf('##$ This file no not exist\n##$ %s',...
            fullfilename(handles, idx_scan_to_rename(i), '.nii'));
        msgbox(warning_text, 'rename file warning') ;
    elseif exist(string(strcat(cellstr(handles.database.Path(idx_scan_to_rename(i))),new_nii_filename{:},'.nii')), 'file') == 2
        msgbox('The new .nii file exist already!!') ;
        
    else
        movefile(fullfilename(handles, idx_scan_to_rename(i), '.nii'), strcat(char(handles.database.Path(idx_scan_to_rename(i))),new_nii_filename{:},'.nii'), 'f')
        if exist(fullfilename(handles, idx_scan_to_rename(i), '.json'), 'file') == 2
            movefile(fullfilename(handles, idx_scan_to_rename(i), '.json'), strcat(char(handles.database.Path(idx_scan_to_rename(i))),new_nii_filename{:},'.json'), 'f');
        end
    end
    
    % update the Filename field in the table
    handles.database.Patient(idx_scan_to_rename(i)) = NewPatient;
    handles.database.Filename(idx_scan_to_rename(i)) = new_nii_filename;
end


% save the structure
guidata(hObject, handles);

set(handles.MIA_name_list, 'Value', 1);

% update graph and display
MIA_update_database_display(hObject, eventdata, handles);


function MIA_remove_name_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_remove_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end
data_selected = finddata_selected(handles);
patient_name = unique(handles.database.Patient(data_selected));
user_response = questdlg(['Do you want to delete every data of ' char(patient_name) '??'], 'Warning', 'Yes', 'No', 'Cancel', 'Cancel');
if strcmp(user_response, 'Cancel') || strcmp(user_response, 'No')
    return
end

nii_index = find(handles.database.Patient == patient_name);

MIA_remove_scan(hObject, eventdata, handles, nii_index)




% --------------------------------------------------------------------
function MIA_name_right_click_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_name_right_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end


% patient_num = get(handles.MIA_name_list, 'Value');
% omit_obj = findobj(handles.MIA_name_right_click, 'Label', 'Omit');
% set(omit_obj, 'Checked', 'off');
% if numel(patient_num) ==1 && ~isempty(handles.database(patient_num).group)
%     group_name = handles.database(patient_num).group;
%     if handles.database(patient_num).omit == 1
%         set(omit_obj, 'Checked', 'on');
%     end
% end
%
% show_menu_obj = findobj(handles.MIA_name_right_click, 'Label', 'Show');
% delete(get(show_menu_obj, 'Children'))  %remove the 'old' show menu
%
% group = reshape({handles.database.group},1,[]);
% group_list = unique(group);
% %group_list = unique(group(1:2:end));
%
% if ~isfield(handles, 'database_all') || numel(handles.database_all) == numel(handles.database)
%     for i = 1:numel(group_list)
%         uimenu(show_menu_obj, 'Label', group_list{i},...
%             'Callback', @(hObject,eventdata)MIA2('MIA_show_group_submenu',hObject,eventdata,guidata(hObject)));
%     end
% else
%     uimenu(show_menu_obj, 'Label', 'all',...
%         'Callback', @(hObject,eventdata)MIA2('MIA_show_group_submenu',hObject,eventdata,guidata(hObject)));
% end
% guidata(hObject, handles);
%
% function MIA_show_group_submenu(hObject, eventdata, handles)
%
% show = get(hObject, 'Label');
% if strcmp('all', show)
%     database_tmp = handles.database_all;
%     for i = 1:numel(handles.database)
%         name_size = numel(handles.database(i).name);
%         match_name = find(strncmp(handles.database(i).name, {database_tmp.name}', name_size) ==1);
%         if numel(match_name) >1
%             for j =1:numel(match_name)
%                 if numel(handles.database(match_name(j)).name) == name_size;
%                     tmp = match_name(j);
%                 end
%             end
%             match_name = tmp;
%         end
%         database_tmp(match_name) = handles.database(i);
%     end
%     handles.database = database_tmp;
%     handles= rmfield(handles, 'database_all');
% else
%     n = 1;
%     for i=1:numel(handles.database)
%         if strncmp(show, handles.database(i).group, numel(show)) == 1 &&...
%                 numel(show) ==  numel(handles.database(i).group)
%             database_tmp(n) = handles.database(i);
%             n=n+1;
%         end
%     end
%     handles.database_all = handles.database;
%     handles.database = database_tmp;
% end
% set(handles.MIA_name_list, 'Value', 1);
% set(handles.MIA_time_points_list, 'Value', 1);
% guidata(handles.MIA_GUI, handles);
% MIA_update_database_display(hObject, eventdata, handles);

% --------------------------------------------------------------------
function MIA_open_database_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to MIA_open_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'database')
    selection = questdlg('Have you saved the present database?',...
        'Warning',...
        'Yes','No','Yes');
    switch selection
        case 'No'
            return
    end
end
path_root=pwd;
%[filename, pathname]=uigetfile('*.mat','Open Mat File','MultiSelect','off');
selpath = uigetdir(path_root,'Select the project''s folder you want to open');
pathname = selpath;
%listfiles = what(selpath)
filename = 'MIA_database.mat';
if pathname == 0
    return
else
    if exist(fullfile(pathname, filename)) ~= 2
        errordlg('The folder you selected might be corrupt. Please select a folder containing a MIA_database.mat file.', 'Cannot open project');
        return
    end
    cd(pathname);
    if ~strcmp(class(filename),'double') %#ok<STISA>
        %reset everything
        handles = MIA_clear_data(hObject, eventdata, handles);
        set(handles.MIA_name_list, 'Value', 1);
        set(handles.MIA_time_points_list, 'Value', 1);
        set(handles.MIA_scan_VOIs_button, 'Value', 0);
        set(handles.MIA_scans_list, 'Value', 1);
        
        database = load(filename);
        handles.database = database.database;

        set(handles.MIA_name_list, 'String', handles.database.Properties.UserData.db_filename)
        
        % update database path (in case the directory has moved) 
        handles.database.Properties.UserData.db_filename = filename;
        new_patient_directory = strcat(pathname, filesep);
        handles.database.Properties.UserData.MIA_data_path  = new_patient_directory;
        handles.database.Properties.UserData.MIA_Raw_data_path = [new_patient_directory, 'Raw_data', filesep];
        handles.database.Properties.UserData.MIA_ROI_path = [new_patient_directory, 'ROI_data', filesep];
        handles.database.Properties.UserData.MIA_Derived_data_path = [new_patient_directory, 'Derived_data', filesep];
        handles.database.Properties.UserData.PSOM_path = [new_patient_directory, 'PSOM', filesep];
        % update the path in the table
        %handles.database.Path(handles.database.Type == 'Scan') = handles.database.Properties.UserData.MIA_Raw_data_path;
        handles.database.Path(handles.database.IsRaw == '0' & handles.database.Type == 'Scan',:) = handles.database.Properties.UserData.MIA_Derived_data_path ;
        handles.database.Path(handles.database.IsRaw == '1' & handles.database.Type == 'Scan',:) = handles.database.Properties.UserData.MIA_Raw_data_path;
        handles.database.Path(handles.database.Type == 'ROI') = handles.database.Properties.UserData.MIA_ROI_path;
        handles.database.Path(handles.database.Type == 'Cluster') = handles.database.Properties.UserData.MIA_ROI_path;
        
        
        guidata(hObject, handles);
    end
    cd(path_root);
end
MIA_update_figureName(hObject, eventdata, handles)

guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);



% --------------------------------------------------------------------
function MIA_time_points_right_click_Callback(~, ~, ~)
% hObject    handle to MIA_time_points_right_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function MIA_rename_time_point_Callback(hObject, eventdata,handles)
% hObject    handle to MIA_rename_time_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
data_selected = finddata_selected(handles);
if numel(data_selected) >1
    warndlg('Please select only one Time point', 'Warning');
    return
end

name_option = [cellstr(unique(handles.database.Tp(handles.database.Type == 'Scan')))' 'Other']';

[new_TP_name, ok1] = listdlg('PromptString','Select the new scan name:',...
    'Name', 'Select a Name',...
    'SelectionMode','single',...
    'ListSize', [400 300],...
    'ListString',name_option);

if ok1 == 0
    return
end
if strcmp('Other',name_option(new_TP_name)) == 1
    NewTp = inputdlg('Name of the new Scan ', 'Question?', 1, {''});
else
    NewTp =name_option(new_TP_name);
end

%% update the database with the new name
% but first check if the new scan name does not exist for this patient and
% time point

% faire le ROI vs SCAN
if find(handles.database.Patient == handles.database.Patient(data_selected) &...
        handles.database.Tp == handles.database.Tp(data_selected) & ...
        handles.database.SequenceName == NewTp) > 0
    msgbox('A Scan with the same name already exist for this patient at this time point') ;
    return
end

idx_scan_to_rename = find(handles.database.Patient == handles.database.Patient(data_selected) & handles.database.Tp == handles.database.Tp(data_selected));
for i=1:numel(idx_scan_to_rename)
    new_nii_filename = strrep(cellstr(handles.database.Filename(idx_scan_to_rename(i))), cellstr(handles.database.Tp(idx_scan_to_rename(i))), NewTp);
    
    % rename the scan file
    if  exist(fullfilename(handles, idx_scan_to_rename(i), '.nii'), 'file') == 0
        warning_text = sprintf('##$ This file no not exist\n##$ %s',...
            fullfilename(handles, idx_scan_to_rename(i), '.nii'));
        msgbox(warning_text, 'rename file warning') ;
    elseif exist(string(strcat(cellstr(handles.database.Path(idx_scan_to_rename(i))),new_nii_filename{:},'.nii')), 'file') == 2
        msgbox('The new .nii file exist already!!') ;
        
    else
        movefile(fullfilename(handles, idx_scan_to_rename(i), '.nii'), strcat(char(handles.database.Path(idx_scan_to_rename(i))),new_nii_filename{:},'.nii'), 'f')
        if exist(fullfilename(handles, idx_scan_to_rename(i), '.json'), 'file') == 2
            movefile(fullfilename(handles, idx_scan_to_rename(i), '.json'), strcat(char(handles.database.Path(idx_scan_to_rename(i))),new_nii_filename{:},'.json'), 'f');
        end
    end
    
    % update the Filename field in the table
    handles.database.Tp(idx_scan_to_rename(i)) = NewTp;
    handles.database.Filename(idx_scan_to_rename(i)) = new_nii_filename;
end


% save the structure
guidata(hObject, handles);

%% update graph and display
MIA_update_database_display(hObject, eventdata, handles);



% --------------------------------------------------------------------
function MIA_rename_scan_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_rename_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
data_selected = finddata_selected(handles);
if numel(data_selected) >1
    warndlg('Please select only one scan', 'Warning');
    return
end

if get(handles.MIA_scan_VOIs_button, 'Value')
    name_option = [cellstr(unique(handles.database.SequenceName(handles.database.Type == 'ROI')))' 'Other']';
else
    name_option = [cellstr(unique(handles.database.SequenceName(handles.database.Type == 'Scan')))' 'Other']';
end

[new_scan_name, ok1] = listdlg('PromptString','Select the new scan name:',...
    'Name', 'Select a Name',...
    'SelectionMode','single',...
    'ListSize', [400 300],...
    'ListString',name_option);

if ok1 == 0
    return
end
if strcmp('Other',name_option(new_scan_name)) == 1
    newparameter = inputdlg('Name of the new Scan ', 'Question?', 1, {''});
else
    newparameter =name_option(new_scan_name);
end

%% update the database with the new name
% but first check if the new scan name does not exist for this patient and
% time point

% faire le ROI vs SCAN
if find(handles.database.Patient == handles.database.Patient(data_selected) &...
        handles.database.Tp == handles.database.Tp(data_selected) & ...
        handles.database.SequenceName == newparameter) > 0
    msgbox('A Scan with the same name already exist for this patient at this time point') ;
    return
end


new_nii_filename = strrep(cellstr(handles.database.Filename(data_selected)), cellstr(handles.database.SequenceName(data_selected)), newparameter);

% rename the scan file
if  exist(fullfilename(handles, data_selected, '.nii'), 'file') == 0
    warning_text = sprintf('##$ This file no not exist\n##$ %s',...
        fullfilename(handles, data_selected, '.nii'));
    msgbox(warning_text, 'rename file warning') ;
elseif exist(string(strcat(cellstr(handles.database.Path(data_selected)),new_nii_filename{:},'.nii')), 'file') == 2
    msgbox('The new .nii file exist already!!') ;
    
else
    movefile(fullfilename(handles, data_selected, '.nii'), strcat(char(handles.database.Path(data_selected)),new_nii_filename{:},'.nii'), 'f')
    if exist(fullfilename(handles, data_selected, '.json'), 'file') == 2
        movefile(fullfilename(handles, data_selected, '.json'), strcat(char(handles.database.Path(data_selected)),new_nii_filename{:},'.json'), 'f');
    end
end

% update the Filename field in the table
handles.database.SequenceName(data_selected) = newparameter;
handles.database.Filename(data_selected) = new_nii_filename;

% save the structure
guidata(hObject, handles);

%% update graph and display
MIA_update_database_display(hObject, eventdata, handles);




% --------------------------------------------------------------------
function MIA_remove_scan_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_remove_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end

nii_index = finddata_selected(handles);
user_response = questdlg('Do you want to delete these data ??', 'Warning', 'Yes', 'No', 'Cancel', 'Cancel');
if strcmp(user_response, 'Cancel') || strcmp(user_response, 'No')
    return
end
MIA_remove_scan(hObject, eventdata, handles, nii_index)

% --------------------------------------------------------------------
function MIA_remove_time_point_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_remove_time_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end
data_selected = finddata_selected(handles);
Time_point_selected = unique(handles.database.Tp(data_selected));
patient_selected = unique(handles.database.Patient(data_selected));
user_response = questdlg(['Do you want to delete every data of ' char(Time_point_selected) '??'], 'Warning', 'Yes', 'No', 'Cancel', 'Cancel');
if strcmp(user_response, 'Cancel') || strcmp(user_response, 'No')
    return
end
nii_index = find(handles.database.Patient == patient_selected & handles.database.Tp == Time_point_selected);

MIA_remove_scan(hObject, eventdata, handles, nii_index)


% --------------------------------------------------------------------
function MIA_ScanVoi_right_click_Callback(~, ~, ~)
% hObject    handle to MIA_ScanVoi_right_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function MIA_update_figureName(~, ~, handles)

num_timepoint = 0;
patient_listing = unique(handles.database.Patient);
for i=1:numel(patient_listing)
    num_timepoint = num_timepoint+ numel(unique(handles.database.Tp(handles.database.Patient == patient_listing(i))));
end
old_title = get(handles.MIA_GUI, 'Name');
title = [old_title(1:strfind(old_title, ';')+1), num2str(numel(patient_listing)),...
    ' patients and ',  num2str(num_timepoint), ' ','time points'];
set(handles.MIA_GUI, 'Name', title);

function handles = MIA_update_handles_parameters_VOIs(handles)


% parameters_list = [];
% Voi_list = [];
% for patient=1:numel(handles.database)
%     for time_point = 1:numel(handles.database(patient).day)
%         parameters_list = [parameters_list handles.database(patient).day(time_point).parameters];
%         Voi_list = [Voi_list  handles.database(patient).day(time_point).VOIs];
%     end
% end
% handles.VOIs = [unique(Voi_list) 'Other'];
VOIs_list = {};
for i=1:numel(handles.database)
    for j = 1:numel(handles.database(i).day)
        VOIs_list = [VOIs_list handles.database(i).day(j).VOIs];
    end
end

handles.VOIs = [unique(VOIs_list), 'Other'];
guidata(handles.MIA_GUI, handles);


function MIA_load_axes_Callback(hObject, eventdata, handles)

if ~isfield(handles, 'database')
    return
end
%tstart = tic;
scan = get(handles.MIA_scans_list, 'Value');
% Load VOIs
if handles.mode == 2 && numel(scan) > 1
    if get(handles.MIA_scan_VOIs_button, 'Value') == 1
        warndlg('Please select only 1 VOI when using the PRM mode','Warning');
    else
        warndlg('Please select only 1 Scan when using the PRM mode','Warning');
    end
    return
end
if get(handles.MIA_scan_VOIs_button, 'Value') && isfield(handles, 'data_loaded') || ...
        get(handles.MIA_scan_VOIs_button, 'Value') && isfield(handles, 'data_selected_for_PRM')
    
    handles = MIA_load_VOIs(hObject, eventdata, handles);
    MIA_update_axes(hObject, eventdata, handles)
    %toc(tstart)
    return
elseif get(handles.MIA_scan_VOIs_button, 'Value') && ~isfield(handles, 'data_loaded')
    warndlg('Please load a scan first','Warning');
    return
end
handles = MIA_clear_data(hObject, eventdata, handles);

% Load Scans
if handles.mode == 1
    handles = MIA_load_axes_single(hObject, eventdata, handles);
else
    handles = MIA_load_axes_PRM(hObject, eventdata, handles);
    % if all conditions are not present --> return
    if ~isfield(handles, 'data_loaded')
        return
    end
    list_day = ['-1', string(handles.data_loaded.info_data_loaded.Tp)'];
    set(handles.MIA_PRM_ref_popupmenu, 'String', list_day', 'Value', 2);
    %set MIA_PRM_slider
    set(handles.MIA_PRM_slider_tp, 'Max', handles.data_loaded.number_of_scan);
    set(handles.MIA_PRM_slider_tp,'Value',1);
    set(handles.MIA_PRM_slider_tp,'Min',1);
    
    set(handles.MIA_PRM_slider_tp,'Visible', 'on');
    set(handles.MIA_PRM_slider_tp,'SliderStep',[1/(handles.data_loaded.number_of_scan-1) min(5/(handles.data_loaded.number_of_scan-1),1)]);
end

% display a waiting symbol
set(handles.MIA_GUI, 'pointer', 'watch');
drawnow;

MIA_update_axes(hObject, eventdata, handles)

set(handles.MIA_GUI, 'pointer', 'arrow');
%toc(tstart)


function handles = MIA_load_VOIs(hObject, ~, handles)

data_selected = finddata_selected(handles);
handles.data_loaded.info_data_loaded(handles.data_loaded.info_data_loaded.Type == 'ROI',:) =[];
if isfield(handles.data_loaded, 'ROI')
    handles.data_loaded= rmfield(handles.data_loaded, 'ROI');
    handles.data_loaded= rmfield(handles.data_loaded, 'number_of_ROI');
end
if isfield(handles.data_displayed, 'ROI')
    handles.data_displayed= rmfield(handles.data_displayed, 'ROI');
end

handles.data_loaded.info_data_loaded(handles.data_loaded.info_data_loaded.Type == 'Cluster',:) =[];
if isfield(handles.data_loaded, 'Cluster')
    handles.data_loaded= rmfield(handles.data_loaded, 'Cluster');
    handles.data_loaded= rmfield(handles.data_loaded, 'number_of_Cluster');
end
if isfield(handles.data_displayed, 'Cluster')
    handles.data_displayed= rmfield(handles.data_displayed, 'Cluster');
end


handles.data_loaded.number_of_ROI = 0;
handles.data_loaded.number_of_Cluster = 0;
for i = 1:numel(data_selected)
    if ~exist(fullfilename(handles, data_selected(i), '.nii'), 'file') && exist(fullfilename(handles, data_selected(i), '.nii.gz'), 'file')
        gunzip(fullfilename(handles, data_selected(i), '.nii.gz'));
        assert(exist(fullfilename(handles, data_selected(i), '.nii'), 'file')==2)
        delete(fullfilename(handles, data_selected(i), '.nii.gz'))
        
    end
    fid_nii=fopen(fullfilename(handles, data_selected(i), '.nii'),'r');
    if fid_nii>0
        fclose(fid_nii);
        % scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');
        if strcmp(char(handles.database(data_selected(i),:).Type), 'ROI')
            %% read and load the nii file
            handles.data_loaded.number_of_ROI = handles.data_loaded.number_of_ROI+1;
            handles.data_loaded.ROI(handles.data_loaded.number_of_ROI).V = spm_vol(char(fullfilename(handles, data_selected(i), '.nii')));
            %           handles.data_loaded.ROI(i).nii = read_volume(handles.data_loaded.ROI(i).V, handles.data_loaded.Scan(scan_of_reference).V,3, handles.view_mode);
            %           handles.data_loaded.ROI(i).nii(handles.data_loaded.ROI(i).nii>0) = 1;
            
            handles.data_loaded.info_data_loaded = [handles.data_loaded.info_data_loaded; handles.database(data_selected(i),:)];
        elseif strcmp(char(handles.database(data_selected(i),:).Type), 'Cluster')
            %% read and load the nii file
            handles.data_loaded.number_of_Cluster = handles.data_loaded.number_of_Cluster+1;
            handles.data_loaded.Cluster(handles.data_loaded.number_of_Cluster).V = spm_vol(char(fullfilename(handles, data_selected(i), '.nii')));
            %           handles.data_loaded.ROI(i).nii = read_volume(handles.data_loaded.ROI(i).V, handles.data_loaded.Scan(scan_of_reference).V,3, handles.view_mode);
            %           handles.data_loaded.ROI(i).nii(handles.data_loaded.ROI(i).nii>0) = 1;
            
            handles.data_loaded.info_data_loaded = [handles.data_loaded.info_data_loaded; handles.database(data_selected(i),:)];
        end
    else
        warndlg('something wrong with the data. Nii or json file is missing','Warning');
        return
    end
    guidata(hObject, handles);
end







function handles = MIA_load_axes_single(hObject, ~, handles)

data_selected = finddata_selected(handles);

if numel(data_selected) > 4  % select only the 4 first scan
    data_selected = data_selected(1:4);
end
for i = 1:numel(data_selected)
    if ~exist(fullfilename(handles, data_selected(i), '.nii'), 'file') && exist(fullfilename(handles, data_selected(i), '.nii.gz'), 'file')
        gunzip(fullfilename(handles, data_selected(i), '.nii.gz'));
        assert(exist(fullfilename(handles, data_selected(i), '.nii'), 'file')==2)
        delete(fullfilename(handles, data_selected(i), '.nii.gz'))
    end
    fid_nii=fopen(fullfilename(handles, data_selected(i), '.nii'),'r');
    fid_json=fopen(fullfilename(handles, data_selected(i), '.json'),'r');
    if fid_nii>0 && fid_json>0
        fclose(fid_nii);
        fclose(fid_json);
        
        %% read and load the json file
        handles.data_loaded.Scan(i).json = spm_jsonread(fullfilename(handles, data_selected(i), '.json'));
        
        %% read and load the nii file
        handles.data_loaded.Scan(i).V =spm_vol(fullfilename(handles, data_selected(i), '.nii'));
    else
        warndlg('something wrong with the data. Nii of json file is missing','Warning');
        return
    end
    clear new
end

set(handles.MIA_patient_information_title, 'String', [char(unique(handles.database.Patient(data_selected))) '_' char(unique(handles.database.Tp(data_selected)))]);
set(handles.MIA_orientation_space_popupmenu, 'String',  char(unique(handles.database.SequenceName(data_selected),'stable')), 'Value', 1);
if numel(data_selected) > 1
    set(handles.MIA_orientation_space_popupmenu, 'Visible', 'on');
    set(handles.MIA_orientation_space_text, 'Visible', 'on');
else
    set(handles.MIA_orientation_space_popupmenu, 'Visible', 'off');
    set(handles.MIA_orientation_space_text, 'Visible', 'off');
end
handles.data_loaded.number_of_scan = numel(data_selected);
handles.data_loaded.info_data_loaded = handles.database(data_selected,:);

guidata(hObject, handles);





if ~isempty(findobj('type', 'figure', 'name', 'FileHistory'))
    % Get the hObject of MIA_pipeline
    h = findobj('type', 'figure', 'name', 'FileHistory');
    % Get the handles of MIA_pipeline
    data = guidata(h);
    % Update the handles of MIA_pipeline by stocking the latest version of
    % MIA handles.
    data.MIA_data = handles;

    % Don't touch the original eventdata, just in case.
    %Update the MIA_pipeline tmp_database
    data.FileHistory_JobsListbox.Value = 1;
    [h,data] = FileHistory('UpdateJobsList', h, data);
    %[~, ~, data] = MIA_pipeline('MIA_pipeline_UpdateTables', h, eventdata2, data);
    guidata(h, data)
end











function handles = MIA_load_axes_PRM(hObject, ~, handles)
% PRM mode i.e. need to open the one parameter (diffusion
% or perfusion or...) for every time point

data_selected = finddata_selected(handles);
if numel(data_selected) ~= 1
    warndlg('In PRM mode you can open only on scan!!', 'Warning');
    return
end


% find indice of the same scan name across each time point to the selected
% patient
data_to_load = find(handles.database.Patient == handles.database.Patient(data_selected) &...
    handles.database.SequenceName == handles.database.SequenceName(data_selected));
[~, idx] =sort(handles.database.Tp(data_to_load));
data_to_load = data_to_load(idx);
if numel(data_to_load) <2
    warndlg(strcat({'Need more than one '},  char(handles.database.SequenceName(data_selected)), ' scan to run the PRM mode') ,'Warning');
    return
end

for i = 1:numel(data_to_load)
    if ~exist(fullfilename(handles, data_to_load(i), '.nii'), 'file') && exist(fullfilename(handles, data_to_load(i), '.nii.gz'), 'file')
        gunzip(fullfilename(handles, data_to_load(i), '.nii.gz'));
        assert(exist(fullfilename(handles, data_to_load(i), '.nii'), 'file')==2)
        delete(fullfilename(handles, data_to_load(i), '.nii.gz'))
    end
    fid_nii=fopen(fullfilename(handles, data_to_load(i), '.nii'),'r');
    fid_json=fopen(fullfilename(handles, data_to_load(i), '.json'),'r');
    if fid_nii>0 && fid_json>0
        fclose(fid_nii);
        fclose(fid_json);
        
        %% read and load the json file
        handles.data_loaded.Scan(i).json = spm_jsonread(fullfilename(handles, data_to_load(i), '.json'));
        
        %% read and load the nii file
        handles.data_loaded.Scan(i).V =spm_vol(fullfilename(handles, data_to_load(i), '.nii'));
    else
        warndlg('something wrong with the data. Nii of json file is missing','Warning');
        return
    end
    clear new
end

set(handles.MIA_patient_information_title, 'String', [char(unique(handles.database.Patient(data_to_load))) '_' char(unique(handles.database.SequenceName(data_to_load)))]);
set(handles.MIA_orientation_space_popupmenu, 'Visible', 'off', 'Value', 1);
handles.data_loaded.number_of_scan = numel(data_to_load);
handles.data_loaded.info_data_loaded = handles.database(data_to_load,:);

%hide used windows
set(handles.MIA_data3, 'Visible', 'off');
set(handles.MIA_data3_title, 'Visible', 'off');
set(handles.MIA_data3, 'HandleVisibility', 'off');

set(handles.MIA_data4, 'Visible', 'off');
set(handles.MIA_data4_title, 'Visible', 'off');
set(handles.MIA_data4, 'HandleVisibility', 'off');

for i=1:2
    stri = num2str(i);
    %  and clear Axes unused
    eval(['cla(handles.MIA_data' stri ');']);
    if i > handles.data_loaded.number_of_scan
        set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off', 'Value', 1);
        set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off', 'Value', 1);
        
    else
        mat_size = handles.data_loaded.Scan(i).V(1).private.dat.dim;
        eval(['set(handles.MIA_data' stri '_title, ''String'', [char(handles.database.SequenceName(data_to_load(i))) ''_'' char(handles.database.Tp(data_to_load(i)))]);']);
        
        if handles.data_loaded.number_of_scan > i-1 && length(mat_size) < 3 % 2D data
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Value', 1);
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
            
        elseif handles.data_loaded.number_of_scan > i-1 && length(mat_size) == 3  % 3D data
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max',  mat_size(3),...
                'SliderStep',[1/(mat_size(3)-1) min(5/(mat_size(3)-1),1)]);
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Value', 1);
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
            
        elseif handles.data_loaded.number_of_scan > i-1 && length(mat_size) == 4 % 4D data
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max',  mat_size(4),...
                'SliderStep',[1/(mat_size(4)-1) min(5/(mat_size(4)-1),1)]);
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
            
        elseif handles.data_loaded.number_of_scan > i-1 && length(handles.data_loaded.Scan(1).V(1).private.dat.dim) == 5 % 5D data
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max',  mat_size(4),...
                'SliderStep',[1/(mat_size(4)-1) min(5/(mat_size(4)-1),1)]);
            
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max', mat_size(5),...
                'SliderStep',[1/(mat_size(5)-1) min(5/(mat_size(5)-1),1)]);
        end
    end
end
% update MIA_table_pixel_values header
% guidata(hObject, handles);
col_header{:,1} = '';
col_header{:,2} = [char(handles.database.SequenceName(data_to_load(1))) '_' char(handles.database.Tp(data_to_load(1)))];
col_header{:,3} = [char(handles.database.SequenceName(data_to_load(2))) '_' char(handles.database.Tp(data_to_load(2)))];

set(handles.MIA_table_pixel_values, 'ColumnName', col_header);
set(handles.MIA_table1, 'ColumnName', {'',char(handles.database.SequenceName(data_to_load(1)))});



% --- Executes during object creation, after setting all properties.
function MIA_resolution_popupmenu_CreateFcn(hObject, ~, ~)
% hObject    handle to MIA_resolution_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = MIA_update_sliders(hObject, eventdata, handles)

%set popupmenu(s) (echo and expt for each axes) and clear Axes if needed
if handles.mode == 1
    for i=1:4 %handles.data_loaded.number_of_scan
        stri = num2str(i);
        %  and clear Axes unused
        eval(['cla(handles.MIA_data' stri ');']);
        if i > handles.data_loaded.number_of_scan
            set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off', 'Value', 1);
            set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off', 'Value', 1);
            
        else
            mat_size = handles.data_loaded.Scan(i).V(1).private.dat.dim;
            eval(['set(handles.MIA_data' stri '_title, ''String'', char(handles.data_loaded.info_data_loaded.SequenceName(i)));']);
            
            if handles.data_loaded.number_of_scan > i-1 && length(mat_size) < 3 % 2D data
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Value', 1);
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
                
            elseif handles.data_loaded.number_of_scan > i-1 && length(mat_size) == 3  % 3D data
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
                %             set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max',  mat_size(3),...
                %                 'SliderStep',[1/(mat_size(3)-1) min(5/(mat_size(3)-1),1)]);
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Value', 1);
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
                
            elseif handles.data_loaded.number_of_scan > i-1 && length(mat_size) == 4 % 4D data
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max',  mat_size(4),...
                    'SliderStep',[1/(mat_size(4)-1) min(5/(mat_size(4)-1),1)]);
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
                
            elseif handles.data_loaded.number_of_scan > i-1 && length(handles.data_loaded.Scan(i).V(1).private.dat.dim) == 5 % 5D data
                set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max',  mat_size(4),...
                    'SliderStep',[1/(mat_size(4)-1) min(5/(mat_size(4)-1),1)]);
                
                set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max', mat_size(5),...
                    'SliderStep',[1/(mat_size(5)-1) min(5/(mat_size(5)-1),1)]);
            end
        end
    end
    
    % resize windows handles.MIA_data1
    switch handles.data_loaded.number_of_scan
        case 1
            set(handles.MIA_data1, 'Position', [0.0188 0.0800 0.5117 0.68]);
            set(handles.MIA_data1_title, 'Visible', 'on');
            set(handles.MIA_data1_echo_slider, 'Position', [0.0334 0.0448 0.3827 0.0168]);
            set(handles.MIA_data1_expt_slider , 'Position', [0.0334 0.0280 0.3827 0.0168]);
            
            set(handles.MIA_data2, 'Visible', 'off');
            set(handles.MIA_data2_title, 'Visible', 'off');
            set(handles.MIA_data2, 'HandleVisibility', 'off');
            
            set(handles.MIA_data3, 'Visible', 'off');
            set(handles.MIA_data3_title, 'Visible', 'off');
            set(handles.MIA_data3, 'HandleVisibility', 'off');
            
            set(handles.MIA_data4, 'Visible', 'off');
            set(handles.MIA_data4_title, 'Visible', 'off');
            set(handles.MIA_data4, 'HandleVisibility', 'off');
            
            set(handles.MIA_slider_slice, 'Position', [0.0334 0.0112 0.3827 0.0168]);
            
        case 2
            set(handles.MIA_data1, 'Position', [0.0188 0.4529 0.2523 0.3140]);
            set(handles.MIA_data1_echo_slider, 'Position', [0.0188 0.4294 0.2518 0.0168]);
            set(handles.MIA_data1_expt_slider , 'Position', [0.0188 0.4126 0.2518 0.0168]);
            
            set(handles.MIA_data2, 'Visible', 'on');
            set(handles.MIA_data2_title, 'Visible', 'on');
            set(handles.MIA_data2, 'HandleVisibility', 'on');
            set(handles.MIA_data3, 'Visible', 'off');
            set(handles.MIA_data3_title, 'Visible', 'off');
            set(handles.MIA_data3, 'HandleVisibility', 'off');
            set(handles.MIA_data4, 'Visible', 'off');
            set(handles.MIA_data4_title, 'Visible', 'off');
            set(handles.MIA_data4, 'HandleVisibility', 'off');
            set(handles.MIA_slider_slice, 'Position', [0.0334 0.395 0.3827 0.0168]);
        case 3
            set(handles.MIA_data1, 'Position', [0.0188 0.4529 0.2523 0.3140]);
            set(handles.MIA_data1_echo_slider, 'Position', [0.0188 0.4294 0.2518 0.0168]);
            set(handles.MIA_data1_expt_slider , 'Position', [0.0188 0.4126 0.2518 0.0168]);
            
            set(handles.MIA_data2, 'Visible', 'on');
            set(handles.MIA_data2_title, 'Visible', 'on');
            set(handles.MIA_data2, 'HandleVisibility', 'on');
            set(handles.MIA_data3, 'Visible', 'on');
            set(handles.MIA_data3_title, 'Visible', 'on');
            set(handles.MIA_data3, 'HandleVisibility', 'on');
            set(handles.MIA_data4, 'Visible', 'off');
            set(handles.MIA_data4_title, 'Visible', 'off');
            set(handles.MIA_data4, 'HandleVisibility', 'off');
            set(handles.MIA_slider_slice, 'Position', [0.0334 0.0112 0.3827 0.0168]);
            
        case 4
            set(handles.MIA_data1, 'Position', [0.0188 0.4529 0.2523 0.3140]);
            set(handles.MIA_data1_echo_slider, 'Position', [0.0188 0.4294 0.2518 0.0168]);
            set(handles.MIA_data1_expt_slider , 'Position', [0.0188 0.4126 0.2518 0.0168]);
            
            set(handles.MIA_data2, 'Visible', 'on');
            set(handles.MIA_data2_title, 'Visible', 'on');
            set(handles.MIA_data2, 'HandleVisibility', 'on');
            set(handles.MIA_data3, 'Visible', 'on');
            set(handles.MIA_data3_title, 'Visible', 'on');
            set(handles.MIA_data3, 'HandleVisibility', 'on');
            set(handles.MIA_data4, 'Visible', 'off');
            set(handles.MIA_data4_title, 'Visible', 'on');
            set(handles.MIA_data4, 'HandleVisibility', 'on');
            set(handles.MIA_slider_slice, 'Position', [0.0334 0.0112 0.3827 0.0168]);
    end
else
    
    %hide used windows
    set(handles.MIA_data3, 'Visible', 'off');
    set(handles.MIA_data3_title, 'Visible', 'off');
    set(handles.MIA_data3, 'HandleVisibility', 'off');
    
    set(handles.MIA_data4, 'Visible', 'off');
    set(handles.MIA_data4_title, 'Visible', 'off');
    set(handles.MIA_data4, 'HandleVisibility', 'off');
end

% set MIA_slider_slice
if  length(handles.data_loaded.Scan(1).V(1).private.dat.dim) == 2  %handles.data_loaded.Scan(1).nii
    set(handles.MIA_slider_slice,'Visible', 'off', 'Value', 1);
    set(handles.MIA_slider_slice,'Max', 1);
    set(handles.MIA_slider_slice,'Value',1);
else
    set(handles.MIA_slider_slice,'Visible', 'on');
    set(handles.MIA_slider_slice,'Min',1);
    set(handles.MIA_slider_slice, 'Max', size(handles.data_displayed.image,3) );
    set(handles.MIA_slider_slice,'Value',1);
    set(handles.MIA_slider_slice,'SliderStep',[1/(size(handles.data_displayed.image,3) -1) min(5/(size(handles.data_displayed.image,3) -1),1)]);
    %set(handles.MIA_slider_slice,'SliderStep',[1/(handles.data_loaded.Scan(1).V(1).private.dat.dim(3) -1) min(5/(handles.data_loaded.Scan(1).V(1).private.dat.dim(3) -1),1)]);
    
end

% update MIA_table_pixel_values header
col_header(:,1) = {'','','','',''};
col_header(2:numel(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'Scan'))+1,1)=cellstr(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'Scan')');
set(handles.MIA_table_pixel_values, 'ColumnName', col_header);
set(handles.MIA_table1, 'ColumnName', col_header);

% reset MIA_plot1
set(get(handles.MIA_plot1, 'XLabel'), 'String', '');
set(get(handles.MIA_plot1, 'YLabel'), 'String', '');
set(get(handles.MIA_plot1, 'ZLabel'), 'String', '');
if ~isempty(findobj('Tag', 'Colorbar'))
    cbfreeze('del');
end




function MIA_update_axes(hObject, eventdata, handles)
%handles = guidata(hObject);

if ~isfield(handles, 'data_loaded')
    return
end


%save data displayed

if ~strcmp(get(hObject, 'Tag'), 'MIA_slider_slice')
    % Update image_displayed matrix
    
    if (isfield(handles, 'data_loaded') && ~(strcmp(get(hObject, 'Tag'), 'MIA_load_axes') && get(handles.MIA_scan_VOIs_button, 'Value'))) && ...
            ~strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
        
        handles = MIA_update_image_displayed(hObject, eventdata, handles);
        
        % Setup every siders, popupmenu when new dataset are loaded
        if strcmp(get(hObject, 'Tag'), 'MIA_load_axes') || strcmp(get(hObject, 'Tag'), 'MIA_Saggital_view_button')...
                || strcmp(get(hObject, 'Tag'), 'MIA_Axial_view_button')  || ...
                strcmp(get(hObject, 'Tag'), 'MIA_Coronal_view_button')|| ...
                strcmp(get(hObject, 'Tag'), 'MIA_orientation_space_popupmenu')
            handles =MIA_update_sliders(hObject, eventdata, handles);
        end
        
        guidata(hObject, handles);
    end
    
    % update the ROI matrix (new ROI, resized...)
    if isfield(handles.data_loaded, 'ROI')
        handles = MIA_update_VOI_displayed(hObject, eventdata, handles);
        
        %     %update MIA_plot1 if needed
        if (strcmp(get(hObject, 'Tag'), 'MIA_load_axes') || strcmp(get(hObject, 'Tag'), 'MIA_PRM_slider_tp') || ...
                (strcmp(get(hObject, 'Tag'), 'MIA_PRM_CI')) || strcmp(get(hObject, 'Tag'), 'MIA_PRM_ref_popupmenu') ...
               || strcmp(get(hObject, 'Tag'), 'MIA_orientation_space_popupmenu'))  && ...
                handles.display_option.view_plot == 1
            if handles.mode == 1
                handles =MIA_update_plot1_single(hObject,handles);
            else
                handles = MIA_update_plot1_PRM(hObject, handles);
            end
            
        end
    end
    
    
    if isfield(handles.data_loaded, 'Cluster')
        handles = MIA_update_VOI_displayed(hObject, eventdata, handles);
        %% In the future, we will have the possibility to update MIA_plot1 here.
    end
    % Update the VOI_cluster matrix (new cluster, resized...)
    if isfield(handles, 'ROI_cluster_resized')
        handles = MIA_update_VOI_cluster_displayed(hObject,handles);
    end
end
slice_nbr = get(handles.MIA_slider_slice, 'Value');
% is zommed?
if numel(get(handles.MIA_data1, 'Children')) ~=0 && ~strcmp(get(hObject, 'Tag'), 'MIA_orientation_space_popupmenu')
    origInfo = getappdata(handles.MIA_data1, 'matlab_graphics_resetplotview');
    if isempty(origInfo)
        isZoomed = false;
    elseif isequal(get(handles.MIA_data1,'XLim'), origInfo.XLim) && ...
            isequal(get(handles.MIA_data1,'YLim'), origInfo.YLim)% && ...
        isZoomed = false;
    else
        isZoomed = true;
        XLim_zoomed = get(handles.MIA_data1,'XLim');
        YLim_zoomed = get(handles.MIA_data1,'YLim');
    end
else
    isZoomed = false;
end

if isfield(handles.data_loaded, 'Cluster')
    if length(handles.data_loaded.Cluster) > 1
        msgbox('Please select only one cluster');
        return
    end
end


% display every data available (image, ROI, cluster...)
if isfield(handles, 'data_displayed')
    number_of_data_to_displayed = numel(fieldnames(handles.data_displayed));
    axe = zeros(1,size( handles.data_displayed.image,4));
    colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
    for i=1:size(handles.data_displayed.image,4)
        stri = num2str(i);
        current_data = ['MIA_data' stri];
        % store current contrast
        current_contrast = get(handles.(current_data), 'Clim');
        if     ~isempty(get(handles.(current_data), 'Children'))
            delete(get(handles.(current_data), 'Children'));
        end
        %         switch number_of_data_to_displayed
        %             case 1 % image only
        if isfield(handles.data_displayed, 'image')
            image_to_display =squeeze(handles.data_displayed.image(:,:,slice_nbr,i));
            
            
            % update image and apply contrast using the values stored (if it is not a new scan
            % loaded)
            if handles.display_option.manual_contrast == 1 && (strcmp(get(hObject, 'Tag'), 'MIA_slider_slice') || ...
                    strcmp(get(hObject, 'Tag'), 'MIA_new_roi') || strcmp(get(hObject, 'Tag'), 'MIA_PRM_slider_tp') || ...
                    strcmp(get(hObject, 'Tag'), 'MIA_load_axes') || strcmp(get(hObject, 'Tag'), 'MIA_PRM_ref_popupmenu') || ...
                     strcmp(get(hObject, 'Tag'), 'MIA_PRM_slider_trans'))
                
                image(image_to_display,'CDataMapping','Scaled','Parent', handles.(current_data),'Tag',current_data);
                set(handles.(current_data), 'Clim', current_contrast );
            else
                % Exlude the 2 extremum (min, max) of the Clim
                tmp = image_to_display;
                tmp(tmp ==0 ) = NaN;

                min = prctile_copy(tmp(:),1);
                max = prctile_copy(tmp(:),99);
                if  ~isnan(min*max) && sum([min max] ~= [0 0]) ~= 0 &&  min ~= max
                    image(image_to_display,'CDataMapping','Scaled','Parent', handles.(current_data),'Tag',current_data);
                    set(handles.(current_data), 'Clim', [min max]);
                end
                
            end
            
            
            
            
            % apply the colormap selected
            colormap(handles.(current_data),  colormap_selected{:})
            if isZoomed == true
                set(handles.(current_data),'XLim', XLim_zoomed);
                set(handles.(current_data),'YLim', YLim_zoomed);
                setappdata(handles.(current_data),'matlab_graphics_resetplotview', origInfo);
            else
                set(handles.(current_data), 'XLim', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 2)]);
                set(handles.(current_data), 'YLim', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 1)]);
            end
            set(handles.(current_data), 'Visible', 'on', 'XTick' , [], 'YTick', []);
            
        end
        
        if isfield(handles.data_displayed, 'ROI')
            hold(handles.(current_data), 'on');
            if strcmp(get(handles.MIA_menu_roi_fill, 'Checked'), 'on')
                fillroi = true;
                trans = get(handles.MIA_PRM_slider_trans, 'Value')/100;
            else
                fillroi = false;
            end
            % if ROI on the slice
            if handles.mode == 1
                %  ROI_indices = find(handles.data_loaded.info_data_loaded.Type == 'ROI');
                for x = 1:numel(handles.data_loaded.ROI)
                    if handles.data_displayed.ROI.on_slice(x,slice_nbr) == 1
                        
                        roi_a_appliquer=handles.data_loaded.ROI(x).nii(:,:,slice_nbr);
                        if fillroi
                            roiRGB = repmat(roi_a_appliquer,[1 1 3]) .* permute(repmat(rgb(handles.colors{x}),[size(roi_a_appliquer,1) 1 size(roi_a_appliquer,1)]),[1 3 2]);
                            image(roiRGB,'AlphaData',roi_a_appliquer*trans,'CDataMapping','Scaled','Parent', handles.(current_data),'Tag',current_data);
                        else
                            contour(handles.(current_data), roi_a_appliquer, 1, 'Color',rgb(handles.colors{x}),...
                                'Visible', 'on',...
                                'tag','ROI_contour');
                        end
                    end
                end
            else
                if handles.data_displayed.ROI.on_slice(1,slice_nbr) == 1
                    roi_a_appliquer=handles.data_loaded.ROI(1).nii(:,:,slice_nbr);
                    if fillroi
                        roiRGB = repmat(roi_a_appliquer,[1 1 3]) .* permute(repmat(rgb(handles.colors{1}),[size(roi_a_appliquer,1) 1 size(roi_a_appliquer,1)]),[1 3 2]);
                        image(roiRGB,'AlphaData',roi_a_appliquer*trans,'CDataMapping','Scaled','Parent', handles.(sprintf('MIA_data%d', i)),'Tag',sprintf('data%d', i));
                    else
                        contour(handles.(sprintf('MIA_data%d', i)), roi_a_appliquer, 1, 'Color',rgb(handles.colors{1}),...
                            'Visible', 'on',...
                            'tag','ROI_contour');
                    end
                    eval(['hold(handles.MIA_data' stri ', ''on'');']);
                    eval(['image(squeeze(handles.data_loaded.PRM.map(:,:,slice_nbr,:)), ''CDataMapping'',''Scaled'', ''parent'', handles.MIA_data' stri ', ''AlphaData'',handles.data_loaded.PRM.trans(:,:,slice_nbr), ''Tag'', ''data' stri '_ROI_cluster'');']);
                    eval(['hold(handles.MIA_data' stri ', ''off'');']);
                end
            end
            set(handles.(current_data), 'Visible', 'on', 'XTick' , [], 'YTick', []);
            hold(handles.(current_data), 'off');
            
        end
        
        if isfield(handles.data_displayed, 'Cluster')
            if handles.data_displayed.Cluster.on_slice(slice_nbr)
                Im_binary = handles.data_displayed.Cluster.data{slice_nbr} >0;
                Im_binary = Im_binary * handles.data_displayed.Cluster.trans;
                %Im_To_Display = zeros(size(handles.data_displayed.Cluster.data{slice_nbr},1), size(handles.data_displayed.Cluster.data{slice_nbr},2),3);
                Im_To_Display = label2rgb(handles.data_displayed.Cluster.data{slice_nbr}, handles.colors_rgb);
                hold(handles.(current_data), 'on');
                image(squeeze(Im_To_Display), 'CDataMapping','Scaled', 'parent',  handles.(current_data), 'AlphaData',Im_binary)%, 'Tag'', ''data' stri '_ROI_cluster')
                hold(handles.(current_data), 'off');
                %                 eval(['hold(handles.MIA_data' stri ', ''on'');']);
                %                 eval(['image(squeeze(handles.data_displayed.VOI_cluster.data(:,:,slice_nbr,:)), ''CDataMapping'',''Scaled'', ''parent'', handles.MIA_data' stri ', ''AlphaData'',handles.data_displayed.VOI_cluster.trans(:,:,slice_nbr), ''Tag'', ''data' stri '_ROI_cluster'');']);    %Affichage du cluster selectionn??
                %                 eval(['hold(handles.MIA_data' stri ', ''off'');']);
            end
        end
        
             
        axe(i) = handles.(current_data);
        %%%%%%%% activate clic on graph MIA_dataXX_ButtonDownFcn
        set(get(handles.(current_data), 'Children'), 'HitTest', 'off');
        set(handles.(current_data),'ButtonDownFcn', @MIA_clic_on_image);
        set(get(handles.(current_data), 'Children'), 'ButtonDownFcn', @MIA_clic_on_image);
        
        
    end
    
end

%Si on affiche un cluster, on affiche les statistiques calcul??es
%pr??cedemment
if exist('Informations', 'var')
    if number_of_data_to_displayed == 3
        set(handles.MIA_table1, 'Data', [StatVol' Espace' StatTranche']');
        set(handles.MIA_table1, 'ColumnName',NomsCartes);
    end
end
linkaxes(axe, 'xy');

guidata(hObject, handles);





function handles = MIA_update_image_displayed(hObject, eventdata, handles)

scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');

switch get(hObject, 'Tag')
    case {'MIA_new_roi', 'MIA_PRM_slider_trans', 'MIA_PRM_CI'}
        return
    case {'MIA_data1_echo_slider', 'MIA_data1_expt_slider'}
        
        data1_echo_nbr = round(get(handles.MIA_data1_echo_slider, 'Value'));
        data1_expt_nbr = round(get(handles.MIA_data1_expt_slider, 'Value'));
        if handles.mode == 1
            handles.data_displayed.image(:,:,:,1) = read_slice(handles.data_loaded.Scan(1).V, handles.data_loaded.Scan(scan_of_reference).V, data1_echo_nbr, data1_expt_nbr, handles.view_mode);
        else
            scan_number = get(handles.MIA_PRM_ref_popupmenu, 'Value');
            if scan_number > 1
                scan_number = scan_number-1;
            else %scan number dynamic 1 prior post-time point
                scan_number = get(handles.MIA_PRM_slider_tp, 'Value') -1;
                if scan_number == 0  % case PRM_ref = -1 and slider = 1
                    scan_number = 1;
                end
            end
            handles.data_displayed.image(:,:,:,1) = read_slice(handles.data_loaded.Scan(scan_number).V, handles.data_loaded.Scan(scan_of_reference).V, data1_echo_nbr, data1_expt_nbr, handles.view_mode);
        end
    case {'MIA_data2_echo_slider', 'MIA_data2_expt_slider'}
        
        data2_echo_nbr = round(get(handles.MIA_data2_echo_slider, 'Value'));
        data2_expt_nbr = round(get(handles.MIA_data2_expt_slider, 'Value'));
        if handles.mode == 1
            handles.data_displayed.image(:,:,:,2) = read_slice(handles.data_loaded.Scan(2).V, handles.data_loaded.Scan(scan_of_reference).V, data2_echo_nbr, data2_expt_nbr, handles.view_mode);
        else
            scan_number = get(handles.MIA_PRM_slider_tp, 'Value');
            handles.data_displayed.image(:,:,:,2) = read_slice(handles.data_loaded.Scan(scan_number).V, handles.data_loaded.Scan(scan_of_reference).V, data2_echo_nbr, data2_expt_nbr, handles.view_mode);
        end
    case {'MIA_data3_echo_slider', 'MIA_data3_expt_slider'}
        data3_echo_nbr = round(get(handles.MIA_data3_echo_slider, 'Value'));
        data3_expt_nbr = round(get(handles.MIA_data3_expt_slider, 'Value'));
        handles.data_displayed.image(:,:,:,3) = read_slice(handles.data_loaded.Scan(3).V, handles.data_loaded.Scan(scan_of_reference).V, data3_echo_nbr, data3_expt_nbr, handles.view_mode);
    case {'MIA_data4_echo_slider', 'MIA_data4_expt_slider'}
        data4_echo_nbr = round(get(handles.MIA_data4_echo_slider, 'Value'));
        data4_expt_nbr = round(get(handles.MIA_data4_expt_slider, 'Value'));
        handles.data_displayed.image(:,:,:,4) = read_slice(handles.data_loaded.Scan(4).V, handles.data_loaded.Scan(scan_of_reference).V, data4_echo_nbr, data4_expt_nbr, handles.view_mode);
    otherwise
        if isfield(handles, 'data_displayed')
            handles = rmfield(handles, 'data_displayed');
        end
        if handles.mode == 1
            
            for i=1:handles.data_loaded.number_of_scan
                stri = num2str(i);
                eval(['data' stri '_echo_nbr = round(get(handles.MIA_data' stri '_echo_slider, ''Value''));']);
                eval(['data' stri '_expt_nbr = round(get(handles.MIA_data' stri '_expt_slider, ''Value''));']);
                
                eval(['ima' stri '= read_slice(handles.data_loaded.Scan(i).V, handles.data_loaded.Scan(scan_of_reference).V, data' stri '_echo_nbr, data' stri '_expt_nbr, handles.view_mode);']);
                
                handles.data_displayed.image(:,:,:,i) = eval(['ima' num2str(i)]);
            end
            
        else
            for i=1:2
                if i == 1 % select pre-scan
                    scan_number = get(handles.MIA_PRM_ref_popupmenu, 'Value');
                    if scan_number > 1
                        scan_number = scan_number-1;
                    else %scan number dynamic 1 prior post-time point
                        scan_number = get(handles.MIA_PRM_slider_tp, 'Value') -1;
                        if scan_number == 0  % case PRM_ref = -1 and slider = 1
                            scan_number = 1;
                        end
                    end
                else % select post-scan
                    scan_number = get(handles.MIA_PRM_slider_tp, 'Value');
                    
                end
                stri = num2str(i);
                eval(['data' stri '_echo_nbr = round(get(handles.MIA_data' stri '_echo_slider, ''Value''));']);
                eval(['data' stri '_expt_nbr = round(get(handles.MIA_data' stri '_expt_slider, ''Value''));']);
                eval(['ima' stri '= read_slice(handles.data_loaded.Scan(scan_number).V, handles.data_loaded.Scan(scan_of_reference).V, data' stri '_echo_nbr, data' stri '_expt_nbr, handles.view_mode);']);
                
                handles.data_displayed.image(:,:,:,i) = eval(['ima' num2str(i)]);
                
                % update title
                eval(['set(handles.MIA_data' stri '_title, ''String'', [char(handles.data_loaded.info_data_loaded.SequenceName(scan_number)) ''_'' char(handles.data_loaded.info_data_loaded.Tp(scan_number))], ''Visible'', ''on'');']);
                
            end
            
        end
        %         if ~strcmp(get(hObject, 'Tag'), 'MIA_PRM_slider_tp')
        %             if ~size(handles.data_displayed.image, 3) == 1
        %                 set(handles.MIA_slider_slice,'Visible', 'off');
        %             else
        %                 set(handles.MIA_slider_slice,'Visible', 'on');
        %                 set(handles.MIA_slider_slice,'Value', 1);
        %                 set(handles.MIA_slider_slice, 'max', size(handles.data_displayed.image, 3), 'Value', 1, 'SliderStep',[1/(size(handles.data_displayed.image, 3) -1) min(5/(size(handles.data_displayed.image, 3) -1),1)]);
        %             end
        %         end
end


guidata(hObject, handles);



function handles = MIA_update_VOI_displayed(hObject, eventdata, handles)

scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');

if isfield(handles.data_loaded, 'ROI')
    for i = 1:numel(handles.data_loaded.ROI)
        switch get(hObject, 'Tag')
            case {'MIA_load_axes', 'MIA_Axial_view_button', 'MIA_Saggital_view_button', 'MIA_Coronal_view_button'}
                handles.data_loaded.ROI(i).nii = read_volume(handles.data_loaded.ROI(i).V, handles.data_loaded.Scan(scan_of_reference).V,'auto', handles.view_mode);
                handles.data_loaded.ROI(i).nii(handles.data_loaded.ROI(i).nii>0) = 1;
        end
        for slice_nbr=1:get(handles.MIA_slider_slice, 'Max')
            roi_a_appliquer=handles.data_loaded.ROI(i).nii(:,:,slice_nbr);
            roi_contour=contourc(double(roi_a_appliquer),1);
            handles.data_displayed.ROI.data(i,slice_nbr) = {roi_contour};
            if sum(sum(handles.data_loaded.ROI(i).nii(:,:,slice_nbr))) > 1
                handles.data_displayed.ROI.on_slice(i,slice_nbr) = 1;
            else
                handles.data_displayed.ROI.on_slice(i,slice_nbr) = 0;
            end
        end
    end
end

if isfield(handles.data_loaded, 'Cluster')
    for i = 1:numel(handles.data_loaded.Cluster)
        switch get(hObject, 'Tag')
            case {'MIA_load_axes', 'MIA_Axial_view_button', 'MIA_Saggital_view_button', 'MIA_Coronal_view_button'}
                handles.data_loaded.Cluster(i).nii = read_volume(handles.data_loaded.Cluster(i).V, handles.data_loaded.Scan(scan_of_reference).V,'auto', handles.view_mode);
                
                handles.data_loaded.Cluster(i).nii(0<handles.data_loaded.Cluster(i).nii & handles.data_loaded.Cluster(i).nii<1) = 1;
                handles.data_loaded.Cluster(i).nii = round(handles.data_loaded.Cluster(i).nii);
        end
        for slice_nbr=1:get(handles.MIA_slider_slice, 'Max')
            
            roi_a_appliquer=handles.data_loaded.Cluster(i).nii(:,:,slice_nbr);
            roi_contour = roi_a_appliquer;
            handles.data_displayed.Cluster.data(i,slice_nbr) = {roi_contour};
            trans = round(handles.MIA_PRM_slider_trans.Value)/100;
            if trans == 0
                trans = 0.01;
            end
            handles.data_displayed.Cluster.trans = trans;
            if sum(sum(handles.data_loaded.Cluster(i).nii(:,:,slice_nbr))) > 1
                handles.data_displayed.Cluster.on_slice(i,slice_nbr) = 1;
            else
                handles.data_displayed.Cluster.on_slice(i,slice_nbr) = 0;
            end
        end
    end
end



function handles = MIA_find_VOI_coordonates(~,handles)

handles.data_ploted.coordonates = [];
if handles.mode == 1
    scan_number = handles.data_loaded(1).number_of_scan;
else
    scan_number = 2;
end
for ii = 1:numel(handles.ROI_selected_resized)
    strii=num2str(ii); %#ok<NASGU>
    dim = size(handles.data_displayed.image);
    if numel(dim) <4
        tmp =  zeros(dim(1:numel(dim)));
    else
        tmp = zeros(dim(1:4));
    end
    for i = 1:scan_number
        stri = num2str(i);
        %check if the data belong to the histo data
        if handles.mode == 1
            if sum(strcmp(handles.data_selected.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected.image(i).reco.color_map, 'Color')
                if size(handles.data_displayed.image,3)==1
                    test = squeeze(handles.data_selected_resized.image(i).reco.data(:,:,1,:,:)).*repmat(handles.ROI_selected_resized(ii).data.value, [1 1 3]);
                    test(test == 0) = NaN;
                    tmp(:,:,:,i) = nanmean(test,3);
                else
                    test = squeeze(handles.data_selected_resized.image(i).reco.data(:,:,1,:,:)).*repmat(handles.ROI_selected_resized(ii).data.value, [1 1 1 3]);
                    test(test == 0) = NaN;
                    tmp(:,:,:,i) = nanmean(test,4);
                end
            else
                if size(handles.data_displayed.image,3)==1 %only 1 slice
                    tmp(:,:,:,i) = squeeze(handles.data_displayed.image(:,:,1,i,1)).*handles.ROI_selected_resized(ii).data.value;
                else
                    tmp(:,:,:,i) = squeeze(handles.data_displayed.image(:,:,:,i,1)).*handles.ROI_selected_resized(ii).data.value;
                end
            end
        else
            if sum(strcmp(handles.data_selected_for_PRM.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected_for_PRM.image(i).reco.color_map, 'Color')
                if size(handles.data_displayed.image,3)==1
                    test = squeeze(handles.data_selected_for_PRM_resized.image(i).reco.data(:,:,1,:,:)).*repmat(handles.ROI_selected_resized(ii).data.value, [1 1 3]);
                    test(test == 0) = NaN;
                    tmp(:,:,:,i) = nanmean(test,3);
                else
                    test = squeeze(handles.data_selected_resized.image(i).reco.data(:,:,1,:,:)).*repmat(handles.ROI_selected_resized(ii).data.value, [1 1 1 3]);
                    test(test == 0) = NaN;
                    tmp(:,:,:,i) = nanmean(test,4);
                end
            else
                if size(handles.data_displayed.image,3)==1 %only 1 slice
                    tmp(:,:,:,i) = squeeze(handles.data_displayed.image(:,:,1,i,1)).*handles.ROI_selected_resized(ii).data.value;
                else
                    tmp(:,:,:,i) = squeeze(handles.data_displayed.image(:,:,:,i,1)).*handles.ROI_selected_resized(ii).data.value;
                end
            end
        end
    end
    % assuming no value = 0 and no NaN in any parameters
    switch scan_number
        case 1
            %             test_result = tmp(:,:,:,1);
            test_result = tmp(:,:,:,1) ~=0 & ~isnan(tmp(:,:,:,1));
        case 2
            test_result = tmp(:,:,:,1) ~=0 & ~isnan(tmp(:,:,:,1)) &...
                tmp(:,:,:,2) ~=0 & ~isnan(tmp(:,:,:,2));
        case 3
            test_result = tmp(:,:,:,1) ~=0 & ~isnan(tmp(:,:,:,1)) &...
                tmp(:,:,:,2) ~=0 & ~isnan(tmp(:,:,:,2))&...
                tmp(:,:,:,3) ~=0 & ~isnan(tmp(:,:,:,3));
        case 4
            test_result = tmp(:,:,:,1) ~=0 & ~isnan(tmp(:,:,:,1)) &...
                tmp(:,:,:,2) ~=0 & ~isnan(tmp(:,:,:,2)) &...
                tmp(:,:,:,3) ~=0 & ~isnan(tmp(:,:,:,3)) &...
                tmp(:,:,:,4) ~=0 & ~isnan(tmp(:,:,:,4));
    end
    tmp2 = findn(test_result(:,:,:,1) ~=0);
    if isempty(tmp2)
        return
    end
    coordonates = zeros(size(tmp2,1), 8);
    if size(tmp2,2) == 2 % if one slice only
        coordonates(:,1:2) = tmp2;
        coordonates(:,3) = ones(size(tmp2,1),1);
    else
        coordonates(:,1:3) = tmp2;
    end
    coordonates(:,4) = ii;
    for j = 1:scan_number
        for i = 1:size(coordonates,1)
            coordonates(i,4+j)=tmp(coordonates(i,1),coordonates(i,2),coordonates(i,3),j);
        end
    end
    if ~isempty(coordonates)
        % remove the NaN pixels
        indices = findn(isnan(coordonates));
        if ~isempty(indices)
            indices = unique(indices(:,1));
            for j = 1:numel(indices)
                indice = indices(numel(indices)+1-j);
                for i = indice:size(coordonates(:,1),1)-1
                    coordonates(i,:)=coordonates(i+1,:);
                end
            end
            tmp = coordonates(1:end-numel(indices),:);
            coordonates = tmp;
        end
    end
    handles.data_ploted.coordonates = [handles.data_ploted.coordonates' coordonates']';
    clear coordonates indices tmp tmp2
end




function handles = MIA_update_plot1_single(hObject, handles)
if ~isempty(get(handles.MIA_plot1, 'Children'))
    delete(get(handles.MIA_plot1, 'Children'));
    legend(handles.MIA_plot1,'off');
    hold(handles.MIA_plot1, 'off');
end
%coordonates = handles.data_ploted.coordonates;
ROI_names = char(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'ROI'));
%Cluster_names = char(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'Cluster'));
Scan_names = char(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'Scan'));

table_data = cell(2*handles.data_loaded.number_of_ROI,handles.data_loaded(1).number_of_scan+1);
ROI_indices = find(handles.data_loaded.info_data_loaded.Type == 'ROI');
%Cluster_indices = find(handles.data_loaded.info_data_loaded.Type == 'Cluster');
for ii = 1:handles.data_loaded.number_of_ROI
    voi_empty = 0;
    strii=num2str(ii);
    if ~isempty(get(handles.MIA_plot1, 'Children'))
        hold(handles.MIA_plot1, 'on');
    end
    if voi_empty == 0
        ROI_binary = handles.data_loaded.ROI(ii).nii;
        ROI_binary(abs(ROI_binary)>0) =1;
        
        switch handles.data_loaded(1).number_of_scan
            case 1
                VOI_data  = handles.data_displayed.image .* ROI_binary;
                VOI_data = reshape(VOI_data, [size(VOI_data,1)*size(VOI_data,2)*size(VOI_data,3),1]);
                VOI_data(VOI_data == 0) =[];
                %nbin = numel(voi_data(:,5))/(numel(voi_data(:,5))/15);
                nbin = 150;
                if ii > 1
                    hold(handles.MIA_plot1, 'on');
                end
                [f, xi] = histnorm(VOI_data,nbin);
                plot(handles.MIA_plot1,xi,f,...
                    'Color',rgb(handles.colors{ii}),...
                    'Tag', strcat('MIA_plot1_1d', strii));
                legende_txt(ii,:) = [Scan_names, '-', ROI_names(ii,:)];
                clear f xi
            case 2
                % find x data
                VOI_data_x  = squeeze(handles.data_displayed.image(:,:,:,1)) .* ROI_binary;
                VOI_data(:,1) = reshape(VOI_data_x, [size(VOI_data_x,1)*size(VOI_data_x,2)*size(VOI_data_x,3),1]);
                %VOI_data_x(VOI_data_x == 0) =[];
                % find y data
                VOI_data_y  = squeeze(handles.data_displayed.image(:,:,:,2)) .* ROI_binary;
                VOI_data(:,2) = reshape(VOI_data_y, [size(VOI_data_y,1)*size(VOI_data_y,2)*size(VOI_data_y,3),1]);
                
                % keep only voxel which has x and y values
                VOI_data((VOI_data(:,1).*VOI_data(:,2)) == 0,:) =[];
                
                scatter(handles.MIA_plot1, VOI_data(:,1), VOI_data(:,2), 'filled',...
                    'SizeData', 20,...
                    'MarkerFaceColor',rgb(handles.colors{ii}),...
                    'MarkerEdgeColor',rgb(handles.colors{ii}),...
                    'Visible', 'on',...
                    'Tag', strcat('MIA_plot1_2d', strii));
                %uistack(findobj('Tag', strcat('MIA_plot1_2d', strii)), 'bottom');
                legende_txt(ii,:) = ROI_names(ii,:);
                
            case 3
                % find x data
                VOI_data_x  = squeeze(handles.data_displayed.image(:,:,:,1)) .* ROI_binary;
                VOI_data(:,1) = reshape(VOI_data_x, [size(VOI_data_x,1)*size(VOI_data_x,2)*size(VOI_data_x,3),1]);
                % find y data
                VOI_data_y  = squeeze(handles.data_displayed.image(:,:,:,2)) .* ROI_binary;
                VOI_data(:,2) = reshape(VOI_data_y, [size(VOI_data_y,1)*size(VOI_data_y,2)*size(VOI_data_y,3),1]);
                % find z data
                VOI_data_z  = squeeze(handles.data_displayed.image(:,:,:,3)) .* ROI_binary;
                VOI_data(:,3) = reshape(VOI_data_z, [size(VOI_data_z,1)*size(VOI_data_z,2)*size(VOI_data_z,3),1]);
                
                
                % keep only voxel which has x and y values
                VOI_data((VOI_data(:,1).*VOI_data(:,2).*VOI_data(:,3)) == 0,:) =[];
                
                scatter3(handles.MIA_plot1, VOI_data(:,1), VOI_data(:,2), VOI_data(:,3),...
                    'filled',...
                    'SizeData', 20,...
                    'MarkerFaceColor',rgb(handles.colors{ii}),...
                    'MarkerEdgeColor',rgb(handles.colors{ii}),...
                    'Tag', strcat('MIA_plot1_3d', strii));
                legende_txt(ii,:) = ROI_names(ii,:);
            case 4
                % find x data
                VOI_data_x  = squeeze(handles.data_displayed.image(:,:,:,1)) .* ROI_binary;
                VOI_data(:,1) = reshape(VOI_data_x, [size(VOI_data_x,1)*size(VOI_data_x,2)*size(VOI_data_x,3),1]);
                % find y data
                VOI_data_y  = squeeze(handles.data_displayed.image(:,:,:,2)) .* ROI_binary;
                VOI_data(:,2) = reshape(VOI_data_y, [size(VOI_data_y,1)*size(VOI_data_y,2)*size(VOI_data_y,3),1]);
                % find z data
                VOI_data_z  = squeeze(handles.data_displayed.image(:,:,:,3)) .* ROI_binary;
                VOI_data(:,3) = reshape(VOI_data_z, [size(VOI_data_z,1)*size(VOI_data_z,2)*size(VOI_data_z,3),1]);
                % find t data
                VOI_data_t  = squeeze(handles.data_displayed.image(:,:,:,4)) .* ROI_binary;
                VOI_data(:,4) = reshape(VOI_data_t, [size(VOI_data_t,1)*size(VOI_data_t,2)*size(VOI_data_t,3),1]);
                
                % keep only voxel which has x and y values
                VOI_data((VOI_data(:,1).*VOI_data(:,2).*VOI_data(:,3).*VOI_data(:,4)) == 0,:) =[];
                % remove nan
                %                 VOI_data(isnan(VOI_data(:,1)),:) = [];
                VOI_data(isnan(mean(VOI_data,2)),:) = [];
                
                color4d = zeros(size(VOI_data(:,4),1),3);
                tmp = jet(256);
                mini=min(VOI_data(:,4));
                maxi=max(VOI_data(:,4));
                for i = 1:size(VOI_data(:,4),1)
                    color4d(i,:) = tmp(round((VOI_data(i,4)-mini)/(maxi-mini)*255)+1,:);
                end
                cbfreeze('del');
                scatter3(handles.MIA_plot1, VOI_data(:,1), VOI_data(:,2), VOI_data(:,3), 20, color4d, 'filled',...
                    'Marker', handles.markers{ii},...
                    'Tag', strcat('MIA_plot1_4d', strii));
                if ~strcmp(get(gco, 'Tag'), 'speedy_run_button')
                    colormap(jet);
                    colorbar('peer', handles.MIA_plot1);
                end
                legende_txt(ii,:) = ROI_names(ii,:);
        end
        %update table
        % first compute the volume of one voxel
        scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');
        nifti_info = niftiinfo(handles.data_loaded.Scan(scan_of_reference).V(1).fname);
        voxel_volume = prod(nifti_info.raw.pixdim(2:4));
        
        table_data(ii*3-2,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(ii))) '-Volume (mm3)']};
        table_data(ii*3-1,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(ii))) '-mean']};
        table_data(ii*3,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(ii))) '-SD']};
        VOI_data(VOI_data==0)=nan;
        
        table_data(ii*3-2,2:2+size(VOI_data,2)-1) = num2cell(numel(VOI_data)*voxel_volume);
        table_data(ii*3-1,2:2+size(VOI_data,2)-1) = num2cell(nanmean(VOI_data));
        table_data(ii*3,2:2+size(VOI_data,2)-1) = num2cell(nanstd(VOI_data));
    end
    clear VOI_data
end

clear index voi_data
if ii == 1
    hold(handles.MIA_plot1, 'on');
end

if handles.data_loaded(1).number_of_scan == 1
    set(get(handles.MIA_plot1, 'YLabel'), 'String', 'Frequency Normalized');
else
    set(get(handles.MIA_plot1, 'YLabel'), 'String', get(handles.MIA_data2_title, 'String'));
end
legend(handles.MIA_plot1,legende_txt, 'Location','NorthEast');
set(get(handles.MIA_plot1, 'XLabel'), 'String', get(handles.MIA_data1_title, 'String'));

set(get(handles.MIA_plot1, 'ZLabel'), 'String', get(handles.MIA_data3_title, 'String'));
hold(handles.MIA_plot1, 'off');



%update table1

set(handles.MIA_table1, 'Data', table_data);

if ~isempty(table_data)
    % set ColumnWidth to auto

    merge_Data = [get(handles.MIA_table1, 'ColumnName')'; table_data];
    %merge_Data = table_data;
    dataSize = size(merge_Data);
    % Create an array to store the max length of data for each column
    maxLen = zeros(1,dataSize(2));
    % Find out the max length of data for each column
    % Iterate over each column
    for i=1:dataSize(2)
        % Iterate over each row
        for j=1:dataSize(1)
            if strcmp(class(merge_Data{j,i}), 'char') %#ok<ISCHR>
                len = length(merge_Data{j,i});
            elseif strcmp(class(merge_Data{j,i}), 'double') %#ok<STISA>
                len = length(num2str(merge_Data{j,i}));
            end
            
            % Store in maxLen only if its the data is of max length
            if(len > maxLen(1,i))
                maxLen(1,i) = len;
            end
        end
    end
    % Some calibration needed as ColumnWidth is in pixels
    cellMaxLen = num2cell(maxLen*7.5);
    % Set ColumnWidth of UITABLE
    set(handles.MIA_table1, 'ColumnWidth', cellMaxLen);
    
end




% --------------------------------------------------------------------
function MIA_plot1_right_click_Callback(hObject, eventdata, ~)
% hObject    handle to MIA_plot1_right_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% patient('figure_patient_ButtonDownFcn',hObject,eventdata,guidata(hObject))

% --------------------------------------------------------------------
function MIA_plot1_3d_view_Callback(~, ~, handles)
% hObject    handle to MIA_plot1_3d_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.MIA_plot1_3d_view, 'Checked'), 'off')
    set(handles.MIA_plot1_3d_view, 'Checked', 'on')
    set(handles.MIA_plot1_3d_view, 'Label', 'Desactivate 3d view');
    rotate3d(handles.MIA_plot1, 'on');
    
else
    set(handles.MIA_plot1_3d_view, 'Checked', 'off')
    set(handles.MIA_plot1_3d_view, 'Label', 'Activate 3d view');
    rotate3d(handles.MIA_plot1, 'off');
end


% --- Executes on mouse press over axes background.
function MIA_plot1_ButtonDownFcn(~, ~, ~)
% hObject    handle to MIA_plot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on mouse motion over figure - except title and menu.
function MIA_GUI_WindowButtonMotionFcn(hObject, ~, handles)
% hObject    handle to MIA_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'data_loaded') && ~isfield(handles, 'data_selected_for_PRM')
    return
end
slice_nbre = get(handles.MIA_slider_slice, 'Value');
% Current position of each axes in percentage of the MIA_GUI size
position_plot1 = get(handles.MIA_plot1, 'Position');
position_data1 = get(handles.MIA_data1, 'Position');
position_data2 = get(handles.MIA_data2, 'Position');
position_data3 = get(handles.MIA_data3, 'Position');
position_data4 = get(handles.MIA_data4, 'Position');

% currPt in percentage of the MIA_GUI size
currPt=get(handles.MIA_GUI,'CurrentPoint');

if ~isempty(findobj('Tag', 'Pixel_contour'))
    delete(findobj('Tag', 'Pixel_contour'))
end
if ~isempty(findobj('Tag', 'CurrentDot'))
    delete(findobj('Tag', 'CurrentDot'))
end
if currPt(1) > position_data1(1) && currPt(1) < position_data1(1)+position_data1(3) && ...
        currPt(2) > position_data1(2) && currPt(2) < position_data1(2)+position_data1(4)
    
    currPt_on_axe=get(handles.MIA_data1,'CurrentPoint');
    currPt_on_axe(:,3)=slice_nbre;
    if isfield(handles, 'data_ploted') && ~isempty(handles.data_ploted.coordonates)
        MIA_draw_pixel(hObject,handles,currPt_on_axe);
    end
    MIA_table1_add_pixel_value(hObject,handles,currPt_on_axe);
    
elseif currPt(1) > position_data2(1) && currPt(1) < position_data2(1)+position_data2(3) && ...
        currPt(2) > position_data2(2) && currPt(2) < position_data2(2)+position_data2(4)
    currPt_on_axe=get(handles.MIA_data2,'CurrentPoint');
    currPt_on_axe(:,3)=slice_nbre;
    if handles.mode == 2 || handles.data_loaded(1).number_of_scan >=2
        MIA_table1_add_pixel_value(hObject,handles,currPt_on_axe);
        if isfield(handles, 'data_ploted') && ~isempty(handles.data_ploted.coordonates)
            MIA_draw_pixel(hObject,handles,currPt_on_axe);
        end
    end
elseif currPt(1) > position_data3(1) && currPt(1) < position_data3(1)+position_data3(3) && ...
        currPt(2) > position_data3(2) && currPt(2) < position_data3(2)+position_data3(4)
    if handles.mode == 1 && handles.data_loaded(1).number_of_scan >=3
        currPt_on_axe=get(handles.MIA_data3,'CurrentPoint');
        currPt_on_axe(:,3)=slice_nbre;
        if isfield(handles, 'data_ploted') && ~isempty(handles.data_ploted.coordonates)
            MIA_draw_pixel(hObject,handles,currPt_on_axe);
        end
        MIA_table1_add_pixel_value(hObject,handles,currPt_on_axe);
    end
    
elseif currPt(1) > position_data4(1) && currPt(1) < position_data4(1)+position_data4(3) && ...
        currPt(2) > position_data4(2) && currPt(2) < position_data4(2)+position_data4(4)
    if handles.mode == 1 && handles.data_loaded(1).number_of_scan >=4
        currPt_on_axe=get(handles.MIA_data4,'CurrentPoint');
        currPt_on_axe(:,3)=slice_nbre;
        if isfield(handles, 'data_ploted') && ~isempty(handles.data_ploted.coordonates)
            MIA_draw_pixel(hObject,handles,currPt_on_axe);
        end
        MIA_table1_add_pixel_value(hObject,handles,currPt_on_axe);
    end
    
elseif currPt(1) > position_plot1(1) && currPt(1) < position_plot1(1)+position_plot1(3) && ...
        currPt(2) > position_plot1(2) && currPt(2) < position_plot1(2)+position_plot1(4)
    currPt_on_axe=get(handles.MIA_plot1,'CurrentPoint');
    MIA_plot1_curr_postion(hObject, handles, currPt_on_axe);  %%%  cursor on plot 1
else
    % Clear  values in MIA_table1
    
    set(handles.MIA_table_pixel_values, 'Data', {'Voxel values', '', '', '',''});
end


% --- Executes on mouse press over axes background.
function MIA_clic_on_image(hObject, eventdata, handles)
% hObject    handle to patient_graph1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = guidata(hObject);

if ~strcmp(get(handles.MIA_GUI,'SelectionType'),'normal')
    G.initpnt=get(gca,'currentpoint');
    G.initClim = get(gca,'Clim');
    set(handles.MIA_GUI,'userdata',G);
    set(handles.MIA_GUI, 'WindowButtonMotionFcn',@MIA_AdjWL);
    
    return
end
% cannot plot anything in PRM mode (yet)
if handles.mode == 2
    return
end
% check if a 4d data is loaded
dimension_size = zeros([size(handles.data_displayed.image,4),7]);
for i=1:size(handles.data_displayed.image,4)
    dimension_size(i,1:numel(handles.data_loaded.Scan(i).V(1).private.dat.dim)) = handles.data_loaded.Scan(i).V(1).private.dat.dim;
end
%remove x, y and Z dimensions
dimension_size(:,1:3) =[];
% find only the 4D data (exclude < and > of 4d data)
fourD_data = find(sum(dimension_size,2)  > 0 & sum(dimension_size > 0,2) == 1);
if isempty(fourD_data)
    return
end

% get pixel coordonates
slice_nbre = get(handles.MIA_slider_slice, 'Value');
tag = get(get(hObject, 'Children'), 'Tag');
if size(tag,1)>1
    currPt_on_axe=eval(['get(handles.' tag{end} ',''CurrentPoint'');']);
else
    currPt_on_axe=eval(['get(handles.' tag ',''CurrentPoint'');']);
end
[pixel_coordinates_2d] = [round(currPt_on_axe(1,1)) round(currPt_on_axe(1,2)) round(currPt_on_axe(1,3))];
voxel = pixel_coordinates_2d(1:2);
if voxel(1) == 0 || voxel(2)==0
    return  %bug somewhere
end
% clean old plot (if needed)
if ~isempty(get(handles.MIA_plot1, 'Children'))
    delete(get(handles.MIA_plot1, 'Children'));
    legend(handles.MIA_plot1,'off');
    hold(handles.MIA_plot1, 'off');
end
% display a waiting symbol
set(handles.MIA_GUI, 'pointer', 'watch');
drawnow;

scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');
legende_txt = cell(numel(fourD_data),1);
for i = 1:numel(fourD_data)
    strii = num2str(i);
    tmp  = read_volume(handles.data_loaded.Scan(fourD_data(i)).V, handles.data_loaded.Scan(scan_of_reference).V,'auto', handles.view_mode);
    y_data = squeeze(tmp(voxel(2), voxel(1), slice_nbre,:));
    x_data = 1:size(tmp,4);
    plot(handles.MIA_plot1,x_data,y_data,...
        'Color',rgb(handles.colors{i}),...
        'Tag', strcat('MIA_plot1_1d', strii));
    hold(handles.MIA_plot1, 'on');
    legende_txt{i,1} = char(handles.data_loaded.info_data_loaded.SequenceName(fourD_data(i)));
    
end
% add the legend
if ~isempty(legende_txt)
    legend(handles.MIA_plot1, legende_txt, 'Location','NorthEast');
end

%% Code pour extraire les courbes de bolus d'un pixel de perf.
% CBV_G_Norm  = read_volume(handles.data_loaded.Scan(3).V, handles.data_loaded.Scan(scan_of_reference).V,3, handles.view_mode);
% CBV_Norm  = read_volume(handles.data_loaded.Scan(4).V, handles.data_loaded.Scan(scan_of_reference).V,3, handles.view_mode);
% 
% Voxel3D = [voxel, slice_nbre];
% 
% A = struct('File', handles.data_loaded.Scan(fourD_data(i)).V(1).fname, 'Curve', y_data, 'Voxel', Voxel3D, 'CBV_G_Norm', CBV_G_Norm(Voxel3D(2), Voxel3D(1), Voxel3D(3)), 'CBV_Norm', CBV_Norm(Voxel3D(2), Voxel3D(1), Voxel3D(3)))
% 

%save('/home/cbrossard/Bureau/Comparaison_Courbes_Perf/Type1_2.mat','A')

set(handles.MIA_GUI, 'pointer', 'arrow');





function MIA_table1_add_pixel_value(~,handles,pixel_coordinates)

if ~isfield(handles, 'data_displayed')
    return
end
[pixel_coordinates_2d] = [round(pixel_coordinates(1,1)) round(pixel_coordinates(1,2)) round(pixel_coordinates(1,3))];
if pixel_coordinates_2d(1)<=0 || pixel_coordinates_2d(1) > size(handles.data_displayed.image,2) || ...
        pixel_coordinates_2d(2)<=0 || pixel_coordinates_2d(2) > size(handles.data_displayed.image,1) || ...
        pixel_coordinates_2d(3)<=0 || pixel_coordinates_2d(3) > size(handles.data_displayed.image,3)
    return
end
% if pixel_coordinates_2d(1) > handles.resolution_selected || pixel_coordinates_2d(2) > handles.resolution_selected ||...
%         pixel_coordinates_2d(1) == 0 || pixel_coordinates_2d(2) == 0
%     return % bug somewhere !!!
% end
if ~isfield(handles, 'data_displayed')
    return
end
values = zeros(1,size(handles.data_displayed.image, 4));
Sizes = zeros(1,size(handles.data_displayed.image, 4));
for i = 1:size(handles.data_displayed.image, 4)
    if numel(size(handles.data_displayed.image))>4 && max(max(handles.data_displayed.image(:, :,pixel_coordinates_2d(3),i, 2))) > 0
        histo_data_in_gray = rgb2gray(uint8(squeeze(handles.data_displayed.image(:, :,pixel_coordinates_2d(3),i,:))));
        pixel_value = 255-histo_data_in_gray(pixel_coordinates_2d(2), pixel_coordinates_2d(1));
        if pixel_value == 255
            pixel_value = 0;
        end
        values(i) = pixel_value;
        Sizes(i) = length(num2str(values(i)))+100;
        
    else
        % rgb2ind(uint8(squeeze(handles.data_displayed.image(pixel_coordinates_2d(2), pixel_coordinates_2d(1),pixel_coordinates_2d(3),2, :))), 32)
        values(i) = handles.data_displayed.image(pixel_coordinates_2d(2), pixel_coordinates_2d(1),pixel_coordinates_2d(3),i,1);
        %[pixel_coordinates_2d(2),pixel_coordinates_2d(1), pixel_coordinates_2d(3)]
        %values(i)
        Sizes(i) = length(num2str(values(i)))+100;
    end
end
table_data = get(handles.MIA_table_pixel_values, 'data');
table_data(1,2:1+size(handles.data_displayed.image, 4)) = num2cell(values);

set(handles.MIA_table_pixel_values, 'ColumnWidth', [{'auto'},num2cell(Sizes)]);
set(handles.MIA_table_pixel_values, 'Data', table_data);


function  MIA_draw_pixel(hObject,handles, pixel_coordinates)
% hObject    handle to MIA_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% currPt=get(handles.MIA_data1,'CurrentPoint');
[pixel_coordinates_2d] = [round(pixel_coordinates(1,1)) round(pixel_coordinates(1,2)) round(pixel_coordinates(1,3))];
if pixel_coordinates_2d(3) ~= get(handles.MIA_slider_slice, 'Value'); % pixel-of_interest not on this slice
    return
end
voxel = pixel_coordinates_2d(1:2);
if voxel(1) == 0 || voxel(2)==0
    return  %bug somewhere
end

% check the display option
if  handles.display_option.view_pixel_on_map == 1
    tmp.map = zeros([size(handles.data_displayed.image, 1),size(handles.data_displayed.image, 2),3]);
    tmp.map(voxel(2),voxel(1),:) = [1 0 0];
    tmp.trans = zeros(size(handles.data_displayed.image, 1),size(handles.data_displayed.image, 2));
    tmp.trans(voxel(2),voxel(1)) = 1;  %#ok<STRNU>
    for i = 1:size(handles.data_displayed.image, 4)
        stri = num2str(i);
        eval(['hold(handles.MIA_data' stri ', ''on'');'])
        eval(['image(tmp.map,''parent'', handles.MIA_data' stri ', ''AlphaData'',tmp.trans,  ''Tag'', ''Pixel_contour'');'])
        eval(['hold(handles.MIA_data' stri ', ''off'');'])
    end
end
if  handles.display_option.view_pixel_on_plot == 1 && size(handles.data_displayed.image, 4) ~= 1
    if isfield(handles, 'ROI_selected_resized')
        if isfield(handles, 'data_ploted')
            x = findn(handles.data_ploted.coordonates(:,1) == pixel_coordinates_2d(2));
            if ~isempty(x)
                y = findn(handles.data_ploted.coordonates(x(:,1),2) == pixel_coordinates_2d(1));
                if ~isempty(y)
                    z = findn(handles.data_ploted.coordonates(x(y(:,1),1),3) == pixel_coordinates_2d(3));
                    if ~isempty(z)
                        closest_Pt= x(y(z(:,1),1),1);
                        hold(handles.MIA_plot1, 'on');
                        scatter3(handles.MIA_plot1, handles.data_ploted.coordonates(closest_Pt,5), handles.data_ploted.coordonates(closest_Pt,6), handles.data_ploted.coordonates(closest_Pt,7),...
                            'Tag', 'CurrentDot', 'MarkerFaceColor', [1 0 0]);
                        uistack(findobj('Tag', 'CurrentDot'), 'top');
                        hold(handles.MIA_plot1, 'off');
                        guidata(hObject, handles);
                    end
                end
            else
                return
            end
        else
            return
        end
    end
end

function  MIA_plot1_curr_postion(hObject, handles, currPt)
% hObject    handle to MIA_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_ploted') || size(handles.data_displayed.image, 4) == 1
    return
end
if handles.display_option.view_pixel_on_plot == 1
    x_scale = get(handles.MIA_plot1, 'XTick');
    y_scale = get(handles.MIA_plot1, 'YTick');
    z_scale = get(handles.MIA_plot1, 'ZTick');
    
    currPt1(1,1) = (currPt(1,1)-min(x_scale))/(max(x_scale)-min(x_scale));
    currPt1(1,2) = (currPt(1,2)-min(y_scale))/(max(y_scale)-min(y_scale));
    currPt1(1,3) = (currPt(1,3)-min(z_scale))/(max(z_scale)-min(z_scale));
    
    currPt2(1,1) = (currPt(2,1)-min(x_scale))/(max(x_scale)-min(x_scale));
    currPt2(1,2) = (currPt(2,2)-min(y_scale))/(max(y_scale)-min(y_scale));
    currPt2(1,3) = (currPt(2,3)-min(z_scale))/(max(z_scale)-min(z_scale));
    
    dist=zeros(numel(handles.data_ploted.coordonates(:,5)),1);
    handles.data_ploted.coordonates_normalized = handles.data_ploted.coordonates;
    for i=1:numel(handles.data_ploted.coordonates(:,5))
        handles.data_ploted.coordonates_normalized(i,5) = (handles.data_ploted.coordonates(i,5)-min(x_scale))/(max(x_scale)-min(x_scale));
        handles.data_ploted.coordonates_normalized(i,6) = (handles.data_ploted.coordonates(i,6)-min(y_scale))/(max(y_scale)-min(y_scale));
        handles.data_ploted.coordonates_normalized(i,7) = (handles.data_ploted.coordonates(i,7)-min(z_scale))/(max(z_scale)-min(z_scale));
        dist(i)=norm(cross(currPt2-currPt1,handles.data_ploted.coordonates_normalized(i,5:7)-currPt1))/norm(currPt2-currPt1);
    end
    
    [~, closest_Pt] =  min(dist);
    
    if isfield(handles, 'data_ploted')
        hold(handles.MIA_plot1, 'on');
        scatter3(handles.MIA_plot1, handles.data_ploted.coordonates(closest_Pt,5), handles.data_ploted.coordonates(closest_Pt,6), handles.data_ploted.coordonates(closest_Pt,7),...
            'Tag', 'CurrentDot', 'MarkerFaceColor', [1 0 0]);
        uistack(findobj('Tag', 'CurrentDot'), 'bottom');
        guidata(hObject, handles);
    end
    MIA_table1_add_pixel_value(hObject,handles,handles.data_ploted.coordonates(closest_Pt,1:3))
end
% draw pixel on the images
if handles.display_option.view_pixel_on_map == 1
    if handles.display_option.view_pixel_on_plot == 0
        x_scale = get(handles.MIA_plot1, 'XTick');
        y_scale = get(handles.MIA_plot1, 'YTick');
        z_scale = get(handles.MIA_plot1, 'ZTick');
        
        currPt1(1,1) = (currPt(1,1)-min(x_scale))/(max(x_scale)-min(x_scale));
        currPt1(1,2) = (currPt(1,2)-min(y_scale))/(max(y_scale)-min(y_scale));
        currPt1(1,3) = (currPt(1,3)-min(z_scale))/(max(z_scale)-min(z_scale));
        
        currPt2(1,1) = (currPt(2,1)-min(x_scale))/(max(x_scale)-min(x_scale));
        currPt2(1,2) = (currPt(2,2)-min(y_scale))/(max(y_scale)-min(y_scale));
        currPt2(1,3) = (currPt(2,3)-min(z_scale))/(max(z_scale)-min(z_scale));
        
        dist=zeros(numel(handles.data_ploted.coordonates(:,5)),1);
        handles.data_ploted.coordonates_normalized = handles.data_ploted.coordonates;
        for i=1:numel(handles.data_ploted.coordonates(:,5))
            handles.data_ploted.coordonates_normalized(i,5) = (handles.data_ploted.coordonates(i,5)-min(x_scale))/(max(x_scale)-min(x_scale));
            handles.data_ploted.coordonates_normalized(i,6) = (handles.data_ploted.coordonates(i,6)-min(y_scale))/(max(y_scale)-min(y_scale));
            handles.data_ploted.coordonates_normalized(i,7) = (handles.data_ploted.coordonates(i,7)-min(z_scale))/(max(z_scale)-min(z_scale));
            dist(i)=norm(cross(currPt2-currPt1,handles.data_ploted.coordonates_normalized(i,5:7)-currPt1))/norm(currPt2-currPt1);
        end
        
        [~, closest_Pt] =  min(dist);
    end
    pixel_coordinates =[handles.data_ploted.coordonates(closest_Pt,2),handles.data_ploted.coordonates(closest_Pt,1),handles.data_ploted.coordonates(closest_Pt,3)];
    [pixel_coordinates_2d] = [round(pixel_coordinates(1,1)) round(pixel_coordinates(1,2)) round(pixel_coordinates(1,3))];
    MIA_table1_add_pixel_value(hObject,handles,pixel_coordinates_2d)
    
    if pixel_coordinates_2d(3) ~= get(handles.MIA_slider_slice, 'Value'); % pixel-of_interest not on this slice
        return
    end
    
    [pixel_coordinates_2d] = [round(pixel_coordinates(1,1)) round(pixel_coordinates(1,2)) round(pixel_coordinates(1,3))];
    if pixel_coordinates_2d(3) ~= get(handles.MIA_slider_slice, 'Value'); % pixel-of_interest not on this slice
        return
    end
    
    voxel = pixel_coordinates_2d(1:2);
    tmp.map = zeros([size(handles.data_displayed.image, 1),size(handles.data_displayed.image, 2),3]);
    tmp.map(voxel(2),voxel(1),:) = [1 0 0];
    tmp.trans = zeros(size(handles.data_displayed.image, 1),size(handles.data_displayed.image, 2));
    tmp.trans(voxel(2),voxel(1)) = 1; %#ok<STRNU>
    
    for i = 1:size(handles.data_displayed.image, 4)
        stri = num2str(i);
        eval(['hold(handles.MIA_data' stri ', ''on'');'])
        eval(['image(tmp.map, ''parent'', handles.MIA_data' stri ', ''AlphaData'',tmp.trans, ''Tag'', ''Pixel_contour'');'])
        eval(['hold(handles.MIA_data' stri ', ''off'');'])
    end
end





% --- Executes on button press in MIA_mode_single_button.
function MIA_mode_single_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_mode_single_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update display mode PRM (multi-time point) --> single time point


set(handles.MIA_mode_single_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_mode_multi_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.mode = 1;
if ~isfield(handles, 'database')
    return
end
handles = MIA_clear_data(hObject, eventdata, handles);


% PRM option off
set(handles.MIA_PRM_CI, 'Visible', 'off');
set(handles.MIA_PRM_CI_title, 'Visible', 'off');
set(handles.MIA_PRM_ref_title, 'Visible', 'off');
set(handles.MIA_PRM_ref_popupmenu, 'Visible', 'off');
set(handles.MIA_PRM_slider_tp, 'Visible', 'off');

% single time point option on
set(handles.MIA_data1, 'Visible', 'on', 'Position', [0.0188 0.4529 0.2523 0.3140], 'XTick', [], 'YTick', []);
set(handles.MIA_data2, 'Visible', 'on', 'XTick', [], 'YTick', []);
set(handles.MIA_data2_title, 'Visible', 'on');
set(handles.MIA_data3, 'Visible', 'on', 'XTick', [], 'YTick', []);
set(handles.MIA_data3_title, 'Visible', 'on');
set(handles.MIA_data4, 'Visible', 'on', 'XTick', [], 'YTick', []);
set(handles.MIA_data4_title, 'Visible', 'on');
set(handles.MIA_slider_slice, 'Position', [0.0334 0.0434 0.3827 0.0168]);

guidata(hObject, handles)

function handles = MIA_clear_data(hObject, eventdata, handles)

if isfield(handles, 'data_loaded')
    handles = rmfield(handles, 'data_loaded');
    delete(get(handles.MIA_data1, 'Children'));
    delete(get(handles.MIA_data2, 'Children'));
    delete(get(handles.MIA_data3, 'Children'));
    delete(get(handles.MIA_data4, 'Children'));
    set(handles.MIA_data1_echo_slider, 'Visible', 'off');
    set(handles.MIA_data1_expt_slider, 'Visible', 'off');
    set(handles.MIA_data2_echo_slider, 'Visible', 'off');
    set(handles.MIA_data2_expt_slider, 'Visible', 'off');
    set(handles.MIA_data1_title, 'String', '');
    set(handles.MIA_data2_title, 'String', '');
    set(handles.MIA_data3_title, 'String', '');
    set(handles.MIA_data4_title, 'String', '');
end
% reset the pointer
set(handles.MIA_GUI, 'pointer', 'arrow');

% reset every cluster for now
handles.table1.cluster = [];
handles.table1.cluster_row = [];

%uncheck the mask option
%set(handles.MIA_menu_define_mask, 'Checked', 'off');

% restore the manual contrast to 0;
handles.display_option.manual_contrast = 0;

if isfield(handles, 'data_selected_resized')
    handles = rmfield(handles, 'data_selected_resized');
end
if isfield(handles, 'data_displayed')
    handles = rmfield(handles, 'data_displayed');
end


% close clip table if open
if ~isempty(findobj('Tag', 'MIA_clip_table'))
    delete(findobj('Tag', 'MIA_clip_table'));
end

if isfield(handles, 'data_selected_for_PRM')
    handles = rmfield(handles, 'data_selected_for_PRM');
    delete(get(handles.MIA_data1, 'Children'));
    delete(get(handles.MIA_data2, 'Children'));
    set(handles.MIA_data1_echo_slider, 'Visible', 'off');
    set(handles.MIA_data1_expt_slider, 'Visible', 'off');
    set(handles.MIA_data1_title, 'String', '');
    set(handles.MIA_data2_title, 'String', '');
end
if isfield(handles, 'data_selected_for_PRM_resized')
    handles = rmfield(handles, 'data_selected_for_PRM_resized');
end

if isfield(handles, 'ROI_selected')
    handles = rmfield(handles, 'ROI_selected');
end
if isfield(handles, 'ROI_selected_resized')
    handles = rmfield(handles, 'ROI_selected_resized');
end
if isfield(handles, 'ROI_PRM_resized')
    handles = rmfield(handles, 'ROI_PRM_resized');
end

if isfield(handles, 'ROI_cluster_resized')
    handles = rmfield(handles, 'ROI_cluster_resized');
end

if ~isempty(get(handles.MIA_plot1, 'Children'))
    delete(get(handles.MIA_plot1, 'Children'));
    legend(handles.MIA_plot1,'off');
    hold(handles.MIA_plot1, 'off');
end

if isfield(handles, 'data_ploted')
    handles = rmfield(handles, 'data_ploted');
end

% clear voxel table ColumName
set(handles.MIA_table_pixel_values, 'ColumnName', {''});

% clear table_1
set(handles.MIA_table1, 'ColumnName', {''});
set(handles.MIA_table1, 'Data', {''});
set(handles.MIA_patient_information_title, 'String', 'No images');

set(handles.MIA_orientation_space_popupmenu, 'String', 'Select orientation', 'Value', 1)

if ~isempty(findobj('Tag', 'legend'))
    delete(findobj('Tag', 'legend'));
end

% set every slider value to 1
set(handles.MIA_data1_echo_slider, 'Value', 1);
set(handles.MIA_data2_echo_slider, 'Value', 1);
set(handles.MIA_data3_echo_slider, 'Value', 1);
set(handles.MIA_data4_echo_slider, 'Value', 1);
set(handles.MIA_data1_expt_slider, 'Value', 1);
set(handles.MIA_data2_expt_slider, 'Value', 1);
set(handles.MIA_data3_expt_slider, 'Value', 1);
set(handles.MIA_data4_expt_slider, 'Value', 1);



% --- Executes on button press in MIA_mode_multi_button.
function MIA_mode_multi_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_mode_multi_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update display mode single time point --> multi-time point


set(handles.MIA_mode_multi_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_mode_single_button,'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.mode = 2;
if ~isfield(handles, 'database')
    return
end
handles = MIA_clear_data(hObject, eventdata, handles);


% single time point option off
set(handles.MIA_data2, 'Visible', 'off');
set(handles.MIA_data2_echo_slider, 'Visible', 'off');
set(handles.MIA_data2_expt_slider, 'Visible', 'off');
set(handles.MIA_data3, 'Visible', 'off');
set(handles.MIA_data3_echo_slider, 'Visible', 'off');
set(handles.MIA_data3_expt_slider, 'Visible', 'off');
set(handles.MIA_data3_title, 'Visible', 'off');
set(handles.MIA_data4, 'Visible', 'off');
set(handles.MIA_data4_echo_slider, 'Visible', 'off');
set(handles.MIA_data4_expt_slider, 'Visible', 'off');
set(handles.MIA_data4_title, 'Visible', 'off');


% PRM option on


set(handles.MIA_slider_slice, 'Position', [0.0334 0.340 0.3827 0.0168]);
set(handles.MIA_data1, 'Visible', 'on', 'Position', [0.0188 0.4529 0.2523 0.3140],  'XTick', [], 'YTick', []);
set(handles.MIA_data1_echo_slider, 'Position', [0.0188 0.4294 0.2518 0.0168]);
set(handles.MIA_data1_expt_slider, 'Position', [0.0188 0.4294 0.2518 0.0168])
set(handles.MIA_data2, 'Visible', 'on', 'XTick', [], 'YTick', []);
set(handles.MIA_PRM_CI, 'Visible', 'on');
set(handles.MIA_PRM_CI_title, 'Visible', 'on');
set(handles.MIA_PRM_slider_trans, 'Value', 50)
set(handles.MIA_PRM_slider_trans, 'min', 0)
set(handles.MIA_PRM_slider_trans, 'max', 100)
set(handles.MIA_PRM_ref_title, 'Visible', 'on', 'Position', [0.0229 0.36 0.0542 0.0266]);
set(handles.MIA_PRM_ref_popupmenu, 'Visible', 'on', 'String', '-1', 'Value', 1,'Position', [0.0792 0.366 0.0803 0.0224]);
set(handles.MIA_PRM_slider_tp, 'Visible', 'on', 'Position', [ 0.2815 0.365 0.2398 0.0168]);


guidata(hObject, handles)




function MIA_PRM_ref_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MIA_PRM_ref as text
%        str2double(get(hObject,'String')) returns contents of MIA_PRM_ref as a double


% --- Executes during object creation, after setting all properties.
function MIA_PRM_ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MIA_PRM_ref_popupmenu.
function MIA_PRM_ref_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_ref_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_PRM_ref_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_PRM_ref_popupmenu
MIA_update_axes(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MIA_PRM_ref_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_ref_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function MIA_PRM_slider_tp_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_slider_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.MIA_PRM_slider_tp, 'Value', round(get(handles.MIA_PRM_slider_tp, 'Value')));

MIA_update_axes(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MIA_PRM_slider_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_slider_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function MIA_PRM_CI_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MIA_PRM_CI as text
%        str2double(get(hObject,'String')) returns contents of MIA_PRM_CI as a double
MIA_update_axes(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MIA_PRM_CI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = MIA_update_plot1_PRM(hObject, handles)

if ~isempty(get(handles.MIA_plot1, 'Children'))
    delete(get(handles.MIA_plot1, 'Children'));
    legend(handles.MIA_plot1,'off');
    hold(handles.MIA_plot1, 'off');
end

%coordonates = handles.data_ploted.coordonates;
ROI_binary = handles.data_loaded.ROI.nii;
ROI_binary(abs(ROI_binary)>0) = true;

% create a matrice which contains the coordonates (x,y,z) of each voxel)
VOI_data = findn(ones(size(ROI_binary)));
% find x data
VOI_data_x  = squeeze(handles.data_displayed.image(:,:,:,1)) .* ROI_binary;
VOI_data(:,4) = reshape(VOI_data_x, [size(VOI_data_x,1)*size(VOI_data_x,2)*size(VOI_data_x,3),1]);

% find y data
VOI_data_y  = squeeze(handles.data_displayed.image(:,:,:,2)) .* ROI_binary;
VOI_data(:,5) = reshape(VOI_data_y, [size(VOI_data_y,1)*size(VOI_data_y,2)*size(VOI_data_y,3),1]);


% keep only voxel which has x and y values
VOI_data((VOI_data(:,4).*VOI_data(:,5)) == 0,:) =[];
VOI_data(isnan(VOI_data(:,4)),:)=[];
VOI_data(isnan(VOI_data(:,5)),:)=[];

% determine PRM+, PRM0 and PRM-
CI = str2double(get(handles.MIA_PRM_CI, 'String'));
PRMr = VOI_data(:,5) > VOI_data(:,4)+CI;
PRMb = VOI_data(:,5) < VOI_data(:,4)-CI;
PRMg = ~or(PRMr, PRMb);

scatter(handles.MIA_plot1,VOI_data(PRMr,4), VOI_data(PRMr,5), 20, [1 0 0]);%,'r')%...
% 'MarkerFaceColor','r');
%     'filled',...
%     'SizeData', 20,...
%     'MarkerEdgeColor','r',...
%     'Visible', 'on',...
%     'Tag', 'MIA_plot1_PRM_r');
hold(handles.MIA_plot1, 'on');
scatter(handles.MIA_plot1,VOI_data(PRMb,4), VOI_data(PRMb,5), 20, [0 0 1]);
%     'filled',...
%     'SizeData', 20,...
%     'MarkerFaceColor','b',...
%     'MarkerEdgeColor','b',...
%     'Visible', 'on',...
%     'Tag', 'MIA_plot1_PRM_b');

scatter(handles.MIA_plot1,VOI_data(PRMg,4), VOI_data(PRMg,5), 20, [0 1 0]);
%     'filled',...
%     'SizeData', 20,...
%     'MarkerFaceColor','g',...
%     'MarkerEdgeColor','g',...
%     'Visible', 'on',...
%     'Tag', 'MIA_plot1_PRM_g');

% plot CI line
x_values = get(handles.MIA_plot1, 'XTick');
y_values = get(handles.MIA_plot1, 'YTick');
min_value = min([x_values y_values]);
max_value = max([x_values y_values]);

plot(handles.MIA_plot1, [min_value, max_value], [min_value, max_value]);
plot(handles.MIA_plot1, [min_value, max_value], [min_value-CI, max_value-CI]);
plot(handles.MIA_plot1, [min_value, max_value], [min_value+CI, max_value+CI]);

set(get(handles.MIA_plot1, 'XLabel'), 'String', get(handles.MIA_data1_title, 'String'));
set(get(handles.MIA_plot1, 'YLabel'), 'String', get(handles.MIA_data2_title, 'String'));

hold(handles.MIA_plot1, 'off');

% save data
if numel(size(handles.data_loaded.ROI.nii)) == 2
    handles.data_loaded.PRM.map=zeros([size(handles.data_loaded.ROI.nii),1,3]);
else
    handles.data_loaded.PRM.map=zeros([size(handles.data_loaded.ROI.nii),3]);
end

handles.data_loaded.PRM.trans=zeros(size(handles.data_loaded.ROI.nii));
trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
for i =1:numel(PRMr)
    if PRMr(i) == 1
        handles.data_loaded.PRM.map(VOI_data(i,1), VOI_data(i,2), VOI_data(i,3),:) = [1 0 0];
    elseif PRMg(i) == 1
        handles.data_loaded.PRM.map(VOI_data(i,1), VOI_data(i,2), VOI_data(i,3),:) = [0 1 0];
    elseif PRMb(i) ==1
        handles.data_loaded.PRM.map(VOI_data(i,1), VOI_data(i,2), VOI_data(i,3),:) = [0 0 1];
    end
    handles.data_loaded.PRM.trans(VOI_data(i,1), VOI_data(i,2), VOI_data(i,3)) = trans;
end

%update table1
table_data = cell(8,2);
table_data(1:8,1) = {'PRMr (%)', 'PRMg (%)', 'PRMb (%)', 'Mean_pre', 'SD_pre', 'Mean_post', 'SD_post', 'Vol (mm3)'};

table_data(1,2) = {sum(PRMr)/numel(VOI_data(:,1))*100};
table_data(2,2) = {sum(PRMg)/numel(VOI_data(:,1))*100};
table_data(3,2) = {sum(PRMb)/numel(VOI_data(:,1))*100};
table_data(4,2) = {mean(VOI_data(:,4))};
table_data(5,2) = {std(VOI_data(:,4))};
table_data(6,2) = {mean(VOI_data(:,5))};
table_data(7,2) = {std(VOI_data(:,5))};

voxel_volume = abs(handles.data_loaded.Scan(1).V(1).mat(1,1)*...
    handles.data_loaded.Scan(1).V(1).mat(2,2)*...
    handles.data_loaded.Scan(1).V(1).mat(3,3));

table_data(8,2) = {numel(VOI_data(:,1))*voxel_volume};
set(handles.MIA_table1, 'Data', table_data);

guidata(hObject, handles);

function handles = MIA_update_PRM_Overlay_map(hObject,handles)

slice_nbr = get(handles.MIA_slider_slice, 'Value');
for slice_nbr=1:get(handles.MIA_slider_slice, 'Max')
    handles.data_displayed.PRM.data(:,:,slice_nbr,:)=squeeze(handles.ROI_PRM_resized.map(:,:,slice_nbr,:));
    handles.data_displayed.PRM.trans(:,:,slice_nbr) = handles.ROI_PRM_resized.trans(:,:,slice_nbr);
end


function handles = MIA_update_VOI_cluster_displayed(hObject,handles)
for slice_nbr=1:get(handles.MIA_slider_slice, 'Max')
    handles.data_displayed.VOI_cluster.data(:,:,slice_nbr,:)=squeeze(handles.ROI_cluster_resized.map(:,:,slice_nbr,:));
    handles.data_displayed.VOI_cluster.trans(:,:,slice_nbr) = handles.ROI_cluster_resized.trans(:,:,slice_nbr);
end


% --- Executes on slider movement.
function MIA_PRM_slider_trans_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_slider_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if strcmp(get(handles.MIA_menu_roi_fill, 'Checked'), 'on')
    MIA_update_axes(hObject, eventdata, handles)
elseif isfield(handles, 'brain_extraction_ROI')
    MIA_Brain_Extraction(hObject, eventdata, handles)
else
    if isfield(handles.data_loaded, 'PRM')
        trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
        index = findn(handles.data_loaded.PRM.map ~=0);
        for i = 1: size(index,1)
            handles.data_loaded.PRM.trans(index(i,1),index(i,2),index(i,3)) = trans;
        end
        guidata(hObject,handles)
        MIA_update_axes(hObject, eventdata, handles)
    elseif isfield(handles, 'data_loaded') && isfield(handles, 'data_displayed')
        if isfield(handles.data_displayed, 'Cluster')
            trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
            if trans == 0
                trans = 0.01;
            end
            handles.data_displayed.Cluster.trans = trans;
            guidata(hObject,handles)
            MIA_update_axes(hObject, eventdata, handles)
        end
    end
end


% --- Executes during object creation, after setting all properties.
function MIA_PRM_slider_trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_PRM_slider_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function MIA_menu_display_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_menu_view_voxel_on_map_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_view_voxel_on_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.MIA_menu_view_voxel_on_map, 'Check'), 'off')
    set(handles.MIA_menu_view_voxel_on_map, 'Check', 'on');
    handles.display_option.view_pixel_on_map = 1;
    
else
    set(handles.MIA_menu_view_voxel_on_map, 'Check', 'off');
    handles.display_option.view_pixel_on_map = 0;
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function MIA_menu_view_voxel_on_plot_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_view_voxel_on_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.MIA_menu_view_voxel_on_plot, 'Check'), 'off')
    set(handles.MIA_menu_view_voxel_on_plot, 'Check', 'on');
    handles.display_option.view_pixel_on_plot = 1;
    
else
    set(handles.MIA_menu_view_voxel_on_plot, 'Check', 'off');
    handles.display_option.view_pixel_on_plot = 0;
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function MIA_menu_view_plot_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_view_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.MIA_menu_view_plot, 'Check'), 'off')
    set(handles.MIA_menu_view_plot, 'Check', 'on');
    handles.display_option.view_plot = 1;
    guidata(hObject, handles)
    if isfield(handles,'ROI_selected_resized')
        
        if handles.mode == 1
            handles =MIA_update_plot1_single(hObject,handles);
        else
            handles =MIA_update_plot1_PRM(hObject,handles);
        end
    end
else
    set(handles.MIA_menu_view_plot, 'Check', 'off');
    handles.display_option.view_plot = 0;
    guidata(hObject, handles)
    if ~isempty(get(handles.MIA_plot1, 'Children'))
        delete(get(handles.MIA_plot1, 'Children'));
        legend(handles.MIA_plot1,'off');
        hold(handles.MIA_plot1, 'off');
        set(handles.MIA_plot1, 'XTick', [], 'YTick', []);
        set(get(handles.MIA_plot1, 'XLabel'), 'String', '');
        set(get(handles.MIA_plot1, 'YLabel'), 'String', '');
        set(get(handles.MIA_plot1, 'ZLabel'), 'String', '');
    end
    if isfield(handles, 'data_ploted')
        handles = rmfield(handles, 'data_ploted');
    end
    
end



% --------------------------------------------------------------------
function MIA_menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function MIA_plot1_plot_all_voxels_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_plot1_plot_all_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'data_loaded')
    return
end
% clear data
if isfield(handles, 'ROI_selected')
    handles = rmfield(handles, 'ROI_selected');
    handles = rmfield(handles, 'ROI_selected_resized');
end
if isfield(handles, 'ROI_PRM_resized')
    handles = rmfield(handles, 'ROI_PRM_resized');
end
if isfield(handles, 'ROI_cluster_resized')
    handles = rmfield(handles, 'ROI_cluster_resized');
end
if ~isempty(get(handles.MIA_plot1, 'Children'))
    delete(get(handles.MIA_plot1, 'Children'));
    legend(handles.MIA_plot1,'off');
    hold(handles.MIA_plot1, 'off');
end
if isfield(handles, 'data_ploted')
    handles = rmfield(handles, 'data_ploted');
end
set(handles.MIA_table1, 'Data', {'', '', '', '', ''});
%%%

for i=1:get(handles.MIA_slider_slice, 'Max')
    handles.ROI_selected_resized.data.value = ones([handles.resolution_selected handles.resolution_selected get(handles.MIA_slider_slice, 'Max')]);
    handles.ROI_selected_resized.name = 'all';
end

guidata(hObject, handles);
MIA_update_axes(hObject, eventdata, handles)


% --------------------------------------------------------------------
function MIA_menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MIA_menu_load_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in MIA_new_ROI_dyn.
function MIA_new_ROI_dyn_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_new_ROI_dyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 MIA_new_roi_Callback(hObject, eventdata, handles)

% --- Executes on button press in MIA_new_ROI_dyn.
function MIA_new_roi_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_new_ROI_dyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% return if no data loaded
if ~isfield(handles, 'data_displayed')
    return
end

slice_nbr = get(handles.MIA_slider_slice, 'Value');


% select on which image will be use to draw the VOI
if handles.data_loaded(1).number_of_scan > 1
    if handles.mode == 1
        [which_image,ok] = listdlg('Name', 'Bip', 'ListString', unique(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'Scan')),'ListSize', [200 150], 'PromptString', 'Which image you want to use?', 'SelectionMode', 'single');
        if ok == 0
            return
        end
        image_number = num2str(which_image);
    else
        image_number = '1';
        which_image = 1;
    end
else
    which_image = 1;
    image_number = '1';
    
end
if handles.mode == 1
    [x, y, z] = size(handles.data_displayed.image(:,:,:,which_image));
else
    warndlg('Please draw your ROI in the single mode (not PRM)', 'Warning');
    return
end


ROI_listing = [unique(handles.database.SequenceName(handles.database.Type == 'ROI'))', 'Other']';
% Enter the ROI's name
[VOI_number,ok] = listdlg('Name', 'Bip', 'SelectionMode', 'single', 'ListString',  ROI_listing,'ListSize', [200 150],...
    'PromptString', 'Select the Roi''s name');
if ok == 0
    return
end

if strcmp('Other', char(ROI_listing(VOI_number))) == 1
    newVOI_name = inputdlg('Name of the new VOI ', 'Question?', 1, {''});
    if isempty(newVOI_name)
        return
    end
    if sum(strcmp(newVOI_name, cellstr(ROI_listing)) > 0)
        warning_text = sprintf('Can not draw this "New" ROI, because %s ROI exist already,\nplease try again and select %s ROI in the ROIs list proposed',...
            newVOI_name{:}, newVOI_name{:});
        msgbox(warning_text, 'ROI warning') ;
        return
    end
    ROI_case = 'New ROI';
    ROI_loaded_idex = 0; % ROI does not exist to it is not loaded
    ROI_idex = [];  %ROI does not exist to in the database
    
else
    %% this ROI name exist already
    % check if this patient at this time point has this ROI
    newVOI_name = cellstr(ROI_listing(VOI_number));
    ROI_idex = find(handles.database.Patient == unique(handles.data_loaded.info_data_loaded.Patient) &...
        handles.database.Tp == unique(handles.data_loaded.info_data_loaded.Tp) &...
        handles.database.Type == 'ROI' & ...
        handles.database.SequenceName == newVOI_name, 1);
    %
    if isempty(ROI_idex) % this ROI does not exist for this patient/time point
        ROI_case = 'New ROI';
        ROI_loaded_idex = 0; % ROI does not exist to it is not loaded
    else % The Voi selected exist for this patient/time_point
        %% check if the ROI is already loaded, if not --> load the ROI
        if isfield(handles.data_loaded, 'ROI') && ...3D_
                sum(strcmp(strcat(cellstr(handles.data_loaded.info_data_loaded.Path(handles.data_loaded.info_data_loaded.Type == 'ROI')), cellstr(handles.data_loaded.info_data_loaded.Filename(handles.data_loaded.info_data_loaded.Type == 'ROI'))),...
                fullfilename(handles, ROI_idex, '')))
            %% the ROI is already loaded
        else
            % load the ROI
            if isfield(handles.data_loaded, 'ROI')
                handles.data_loaded.ROI(numel(handles.data_loaded.ROI)+1).V = spm_vol(fullfilename(handles, ROI_idex, '.nii'));
                handles.data_loaded.ROI(end).nii = read_volume(handles.data_loaded.ROI(end).V , handles.data_loaded.Scan(which_image).V(1), 'auto', handles.view_mode);
            else
                handles.data_loaded.ROI.V = spm_vol(fullfilename(handles, ROI_idex, '.nii'));
                handles.data_loaded.ROI.nii = read_volume(handles.data_loaded.ROI.V , handles.data_loaded.Scan(which_image).V(1),'auto', handles.view_mode);
            end
            % update the data_loaded structure with the new ROI
            handles.data_loaded.number_of_ROI = size(handles.data_loaded.ROI,1);
            handles.data_loaded.info_data_loaded = [handles.data_loaded.info_data_loaded; handles.database(ROI_idex,:)];
            % save the handles
            guidata(hObject, handles)
            % update the fiure
            MIA_update_axes(hObject, eventdata, handles)
        end
        ROI_loaded_listing = handles.data_loaded.info_data_loaded.Filename(handles.data_loaded.info_data_loaded.Type == 'ROI');
        ROI_loaded_idex  = ROI_loaded_listing == handles.database.Filename(ROI_idex);
        
        
        ROI_data_loaded_idex = handles.data_loaded.info_data_loaded.Type == 'ROI' & handles.data_loaded.info_data_loaded.Path == handles.database.Path(ROI_idex) & handles.data_loaded.info_data_loaded.Filename == handles.database.Filename(ROI_idex);
        if sum(sum(handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr))) > 0
            list_of_choices = {'Union', 'Instersection'};
            [user_response,ok] = listdlg('Name', 'Bip', 'ListString', list_of_choices,'ListSize', [200 150], 'PromptString', 'What do you want to do?', 'SelectionMode', 'single');
            if isempty(user_response) || ok == 0
                return
            end
            ROI_case = list_of_choices{user_response};
        else
            ROI_case = 'New slice';
        end
    end
end

switch ROI_case
    case {'New ROI', 'Union', 'Instersection', 'New slice'}
        if handles.display_option.view_pixel_on_map
            % unselect the option 'view_pixel_on_map' 
            set(handles.MIA_menu_view_voxel_on_map, 'Check', 'off');
            handles.display_option.view_pixel_on_map = 0;
            guidata(hObject, handles);
            
            if strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
                eval(['hroi=impoly(handles.MIA_data' image_number ',[]);']);
                position = getPosition(hroi);
            elseif strcmp(get(hObject, 'Tag'), 'MIA_new_ROI_dyn')
                Img=squeeze(handles.data_displayed.image(:,:,get(handles.MIA_slider_slice, 'Value'), str2double(image_number)));
                eval(['[Drawn_ROI_matrice, position] = MIA_new_ROI_dyn(hObject, eventdata, handles, handles.MIA_data' image_number ', Img, str2double(image_number));']);
                position = position'; %#ok<NODEF>
            end
            
            % select the option 'view_pixel_on_map' 
            set(handles.MIA_menu_view_voxel_on_map, 'Check', 'on');
            handles.display_option.view_pixel_on_map = 1;
            guidata(hObject, handles);
        else
            if strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
                eval(['hroi=impoly(handles.MIA_data' image_number ',[]);']);
                position = getPosition(hroi);
            elseif strcmp(get(hObject, 'Tag'), 'MIA_new_ROI_dyn')
                Img=squeeze(handles.data_displayed.image(:,:,get(handles.MIA_slider_slice, 'Value'), str2double(image_number)));
                eval(['[Drawn_ROI_matrice, position] = MIA_new_ROI_dyn(hObject, eventdata, handles, handles.MIA_data' image_number ', Img, str2double(image_number));']);             
                if isempty(Drawn_ROI_matrice)
                    return
                end
                position = position'; %#ok<NODEF>
            end
            
        end
    
end





Scan_of_reference_selected = get(handles.MIA_orientation_space_popupmenu, 'Value');
% Scan_of_reference_listing = get(handles.MIA_orientation_space_popupmenu, 'String');
% Scan_of_reference = Scan_of_reference_listing(Scan_of_reference_selected,:);

if isempty(position) || size(position, 1) ==1
    % case of delete?
    return
end
switch ROI_case    
    case 'New slice'  % new ROI slice added to an existing ROI
        if strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
            handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = createMask(hroi, findobj('Tag', ['MIA_data' image_number]));
        elseif strcmp(get(hObject, 'Tag'), 'MIA_new_ROI_dyn')
            handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = Drawn_ROI_matrice;
        end
        
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(handles.data_loaded.ROI(ROI_loaded_idex).nii, handles.data_loaded.Scan(Scan_of_reference_selected).V, handles.view_mode);
        
    case 'New ROI'

        % tested code
        %create the VOI matrix
        ROI_matrix=zeros([x y z]);
        
         if strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
                ROI_matrix(:,:,slice_nbr)=createMask(hroi, findobj('Tag', ['MIA_data' image_number]));
            elseif strcmp(get(hObject, 'Tag'), 'MIA_new_ROI_dyn')
                ROI_matrix(:,:,slice_nbr)=Drawn_ROI_matrice;
         end
        ROI_matrix = write_volume(ROI_matrix, handles.data_loaded.Scan(Scan_of_reference_selected).V, handles.view_mode );
        
    case 'Union'
          if strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
            handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) + createMask(hroi, findobj('Tag', ['MIA_data' image_number]));
        elseif strcmp(get(hObject, 'Tag'), 'MIA_new_ROI_dyn')
             handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) + Drawn_ROI_matrice;
          end
       
       
        %         handles.data_loaded.ROI(ROI_loaded_idex).nii(handles.data_loaded.ROI(ROI_loaded_idex).nii > 0 ) = 1;
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(handles.data_loaded.ROI(ROI_loaded_idex).nii, handles.data_loaded.Scan(Scan_of_reference_selected).V, handles.view_mode);
    case 'Instersection'
        if strcmp(get(hObject, 'Tag'), 'MIA_new_roi')
            handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = or(handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr), createMask(hroi, findobj('Tag', ['MIA_data' image_number]))) - createMask(hroi, findobj('Tag', ['MIA_data' image_number]));
        elseif strcmp(get(hObject, 'Tag'), 'MIA_new_ROI_dyn')
            handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = or(handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr), Drawn_ROI_matrice) - Drawn_ROI_matrice;        
        end
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(handles.data_loaded.ROI(ROI_loaded_idex).nii, handles.data_loaded.Scan(Scan_of_reference_selected).V, handles.view_mode);
end
ROI_matrix = logical(ROI_matrix);
% create a new V structure based on the image used to drow the VOI
% Choose any of the spm data types here. Get list with spm_type(spm_type)
% To save space with binary images, choose 'int8' or 'uint8'
outputDatatype = 'int8';
V_ROI = handles.data_loaded.Scan(Scan_of_reference_selected).V(1);                     % Copy existing volume description from reference image


file_name = strcat(char(handles.data_loaded.info_data_loaded.Patient(Scan_of_reference_selected)),...
    '-', char(handles.data_loaded.info_data_loaded.Tp(Scan_of_reference_selected)),...
    '-ROI',...
    '-', newVOI_name{:} ,...
    '_',datestr(now,'yyyymmdd-HHMMSSFFF'));          %Specify file name to write to
% Create the ROI folder if needed
if exist(handles.database.Properties.UserData.MIA_ROI_path, 'dir') ~= 7
    status = mkdir([handles.database.Properties.UserData.MIA_data_path, 'ROI_data', filesep]);
    if status == 0
        warndlg('You do not have the right to write in the folder!', 'Warning');
        return
    end
end
V_ROI.fname =  [handles.database.Properties.UserData.MIA_ROI_path, file_name, '.nii'] ;

V_ROI = rmfield(V_ROI,'private'); % Delete old nifti header. Will be recreated to match new image properties
V_ROI.dt(1) = spm_type(outputDatatype); % save images in specified format
if spm_type(outputDatatype,'intt')
    V_ROI = rmfield(V_ROI,'pinfo'); % integer datatype : let spm_write_vol decide on scaling
else
    V_ROI.pinfo(1:2) = [1;0];       % do not apply any scaling when saving as float data
end
% save the ROI in nii file (could be a new ROI or and old but updated)
[ROI_matrix, FinalMat] = CropROI(ROI_matrix, V_ROI.mat);
V_ROI.dim = [size(ROI_matrix,1), size(ROI_matrix,2), size(ROI_matrix,3)];
% V_ROI.mat = adaptedstruct(1).mat;
V_ROI.mat = FinalMat;
% indextomove = 3; %% Not X, not Y, but Z
% V_ROI.mat(indextomove,4) = V_ROI.mat(indextomove,4)+(slice_nbr-1)*V_ROI.mat(indextomove,indextomove);

spm_write_vol(V_ROI,ROI_matrix);
% spm_jsonwrite(fullfilename(handles, roi_index, '.json'), hroi_data);
% load the new ROI
if isfield(handles.data_loaded, 'ROI')
    if sum(ROI_loaded_idex) ~= 0 % update the structre for an updated ROI
        handles.data_loaded.ROI(ROI_loaded_idex).V = spm_vol(V_ROI.fname);
        handles.data_loaded.ROI(ROI_loaded_idex).nii = read_volume(handles.data_loaded.ROI(ROI_loaded_idex).V , handles.data_loaded.Scan(Scan_of_reference_selected).V,'auto', handles.view_mode);
    else %% add new ROI to the data_loaded_ROI structure (another ROI is already loaded)
        handles.data_loaded.ROI(numel(handles.data_loaded.ROI)+1).V = spm_vol(V_ROI.fname);
        handles.data_loaded.ROI(end).nii = read_volume(handles.data_loaded.ROI(end).V , handles.data_loaded.Scan(Scan_of_reference_selected).V,'auto', handles.view_mode);
    end
else %% add the new ROI as load ROI
    handles.data_loaded.ROI.V = spm_vol(V_ROI.fname);
    handles.data_loaded.ROI.nii = read_volume(handles.data_loaded.ROI.V , handles.data_loaded.Scan(Scan_of_reference_selected).V,'auto', handles.view_mode);
end
%% if an ROI has been updated --> delete the old nii file and update the database
if ~isempty(ROI_idex)
    delete(fullfilename(handles, ROI_idex, '.nii'))
    handles.database.Filename(ROI_idex) = file_name;
    handles.data_loaded.info_data_loaded.Filename(ROI_data_loaded_idex) = file_name;
else
    %% add new ROI data to the database
    new_data = table(categorical(cellstr('Undefined')), categorical(handles.data_loaded.info_data_loaded.Patient(which_image)),...
        categorical(handles.data_loaded.info_data_loaded.Tp(which_image)),...
        categorical(cellstr([handles.database.Properties.UserData.MIA_data_path, 'ROI_data', filesep])), categorical(cellstr(file_name)),...
        categorical(cellstr('ROI')), categorical(0), categorical(newVOI_name),...
        'VariableNames', {'Group','Patient', 'Tp', 'Path', 'Filename', 'Type', 'IsRaw', 'SequenceName'});
    
    handles.database = [handles.database;new_data ];
    
    % update the data_loaded structure with the new ROI
    handles.data_loaded.number_of_ROI = size(handles.data_loaded.ROI,1);
    handles.data_loaded.info_data_loaded = [handles.data_loaded.info_data_loaded; new_data];
    
end

%% if an ROI has been updated --> delete the old nii file

guidata(hObject, handles)
MIA_update_axes(hObject, eventdata, handles)



% --- Executes on button press in MIA_test_button.
function [ROI_matrice, position] = MIA_new_ROI_dyn(hObject, eventdata, handles, axe, Img, image_selected)
% hObject    handle to MIA_test_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% this function is adapted from 
% Reference: <Li Wang, Lei He, Arabinda Mishra, Chunming Li. 
% Active Contours Driven by Local Gaussian Distribution Fitting Energy.
% Signal Processing, 89(12), 2009,p. 2435-2447>


NumIter = 5000; %iterations
timestep=0.1; %time step
mu=0.1/timestep;% level set regularization term, please refer to "Chunming Li and et al. Level Set Evolution Without Re-initialization: A New Variational Formulation, CVPR 2005"
sigma = 5;%size of kernel

epsilon = 1;
c0 = 2; % the constant value 
lambda1=1.0;%outer weight, please refer to "Chunming Li and et al,  Minimization of Region-Scalable Fitting Energy for Image Segmentation, IEEE Trans. Image Processing, vol. 17 (10), pp. 1940-1949, 2008"
lambda2=1.0;%inner weight
%if lambda1>lambda2; tend to inflate
%if lambda1<lambda2; tend to deflate
nu = 0.0005*255*255;%length term

% terme important, li� � la d�finition de la limite de la ROI 
% alf = 30 by default
alf = 30;%data term weight 


h = imellipse(axe);
phi =createMask(h);
phi=double(phi);
phi(phi==0) = -1;
phi = sign(phi).*c0.*-1;

first_ellipse = phi;
tmp = h.getPosition;
first_centroid = [tmp(1) + tmp(3)/2, tmp(2) + tmp(4)/2];
delete(h);

Ksigma=fspecial('gaussian',round(2*sigma)*2 + 1,sigma); %  kernel


ONE=ones(size(Img));
KONE = imfilter(ONE,Ksigma,'replicate');  
KI = imfilter(Img,Ksigma,'replicate');  
KI2 = imfilter(Img.^2,Ksigma,'replicate'); 

hold on,[~,~] = contour(axe, phi,[0 0],'r','linewidth',1, 'tag','ROI_active'); hold off
pause(0.5)

edges_found = 0;
slice_selected = get(handles.MIA_slider_slice, 'Value');
nbr_of_slice = get(handles.MIA_slider_slice, 'Max');

indice=0;
ROI_size = nan([NumIter, 1]);
for_loop_stopped = 0;
for iter = 1:NumIter
    
    phi =evolution_LGD(Img,phi,epsilon,Ksigma,KONE,KI,KI2,mu,nu,lambda1,lambda2,timestep,alf);
    
    if(mod(iter,10) == 0)
        indice = indice+1;
        if ~isempty(findobj('Tag', 'ROI_active'))
            delete(findobj('Tag', 'ROI_active'))
        end
        hold on,[~,~] = contour(axe, phi,[0 0],'r','linewidth',1, 'tag','ROI_active'); hold off
        
        BW = ~imbinarize(phi, 0);
        if sum(BW(:)) == 0
            
            ROI_matrice = [];
            position = [];
            return
        end
        test = regionprops(BW);
        test = struct2table(test);
        
        idx = findClosestCentroids(first_centroid, test.Centroid);
        ROI_size(indice)=test.Area(idx);
        
        if indice > 3 && ...
                isequal(ROI_size(indice),ROI_size(indice-1), ROI_size(indice-2)) %, ROI_size(indice-3), ROI_size(indice-4), ROI_size(indice-5))
            break
        end
        pause(0.01);
    end
    
end


ROI_matrice = bwlabel(BW);
blob_to_select = ROI_matrice(round(first_centroid(2)),round(first_centroid(1)));
ROI_matrice(ROI_matrice~=blob_to_select)=0;
ROI_matrice = imbinarize(ROI_matrice);

if ~isempty(findobj('Tag', 'ROI_active'))
    delete(findobj('Tag', 'ROI_active'))
end

if sum(ROI_matrice(:)) == 0
    ROI_matrice = [];
    position = [];
    
    return
end
hold on,[position, ~] = contour(axe, ROI_matrice,1,'r','linewidth',1, 'tag','ROI_active'); hold off

if ~exist('position', 'var')
    position = [];
end




function idx = findClosestCentroids(X, centroids)

% findClosestCentroids computes the closest centroid for each point based
% on the Euclidean distance between the point and the centroid

% Initialize variables
K = size(centroids, 1); 
idx = zeros(size(X,1), 1); % returns index of closest centroid

for i=1:size(X,1)
    temp = X(i,:);
    [~,idx(i,1)] = min(sum(((bsxfun(@minus,temp,centroids)).^2),2));
end



% --------------------------------------------------------------------
function MIA_cloneScanVoi_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_cloneScanVoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
data_selected = finddata_selected(handles);
if numel(data_selected) >1
    warndlg('Please select only one scan', 'Warning');
    return
end
new_parameter_name = strcat(char(handles.database.SequenceName(data_selected)), '-cloned');
if ~exist(fullfilename(handles, data_selected, '.nii'), 'file') && exist(fullfilename(handles, data_selected, '.nii.gz'), 'file')
    gunzip(fullfilename(handles, data_selected, '.nii.gz'));
    assert(exist(fullfilename(handles, data_selected, '.nii'), 'file')==2)
    delete(fullfilename(handles, data_selected, '.nii.gz'))
end
fid_nii=fopen(fullfilename(handles, data_selected, '.nii'),'r');
if fid_nii>0
    fclose(fid_nii);
    % update the database
    new_scan = handles.database(data_selected,:);
    new_scan.SequenceName = new_parameter_name;
    new_scan.Filename =  strrep(char(new_scan.Filename), char(handles.database.SequenceName(data_selected)), new_parameter_name);
    handles.database(size(handles.database,1)+1,:) = new_scan;
    % create the new files
    copyfile(fullfilename(handles, data_selected, '.nii'), fullfilename(handles, size(handles.database,1), '.nii'), 'f');
    if exist(fullfilename(handles, data_selected, '.json'), 'file') == 2
        copyfile(fullfilename(handles, data_selected, '.json'), fullfilename(handles, size(handles.database,1), '.json'), 'f');
    end
else
    warndlg('something is wrong this the data',  'Warning');
end
% save handles and update display
guidata(hObject, handles)
MIA_update_database_display(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MIA_save_database_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_save_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% check if the database has been saved already. If to goto saveas function
if ~isfield(handles, 'database')
    return
end

MIA_menu_save_database_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MIA_menu_save_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_save_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles.database.Properties.UserData , 'db_filename')
    %     [filename,pathname] = uiputfile('.mat','Please enter the name of the database', handles.database.Properties.UserData.MIA_data_path );%pathstr(1:file_sep(end)));
    %     if pathname == 0
    %         return
    %     end
    handles.database.Properties.UserData.db_filename = 'MIA_database.mat';
    %handles.database.Properties.UserData.MIA_data_path = pathname;
end

handles.database.Properties.UserData.VOIs = handles.VOIs;
handles.database.Properties.UserData.histo = handles.histo;
database = handles.database; %#ok<NASGU>
save(fullfile(handles.database.Properties.UserData.MIA_data_path, handles.database.Properties.UserData.db_filename), 'database');

% save handles
guidata(hObject, handles)

msgbox('Done', 'Information') ;


% --------------------------------------------------------------------
function MIA_show_group_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_show_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_name_properties_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_name_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_add_info_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_add_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
patient = get(handles.MIA_name_list, 'Value');
old_group_name = {handles.database(patient).group};
group_name= inputdlg('Enter the group name', 'Add group name',1, old_group_name(1));
if isempty(group_name)
    return
end
for i = 1:numel(patient)
    handles.database(patient(i)).group = group_name{:};
end


if isfield(handles, 'database_all') && numel(handles.database) < numel(handles.database_all)
    for i = 1:numel(patient)
        for j = 1:numel(handles.database_all)
            if strcmp(handles.database(patient(i)).name, handles.database_all(j).name) ==1
                handles.database_all(j).group = group_name;
            end
        end
    end
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function MIA_menu_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_menu_delete_from_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_delete_from_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end

selection = questdlg('What do you want to delete?',...
    'Warning',...
    'Parameters','VOIs', 'Cancel', 'Cancel');

switch selection
    case 'Parameters'
        name_option = cellstr(unique(handles.database.SequenceName(handles.database.Type == 'Scan')));
        if isempty(name_option)
            warndlg('There are no Scan in the database', 'Warning');
            return
        end
    case 'VOIs'
        name_option = cellstr(unique(handles.database.SequenceName(handles.database.Type == 'ROI')));
        if isempty(name_option)
            warndlg('There are no ROI in the database', 'Warning');
            return
        end
    case 'Cancel'
        return
end

[ScansToRemove, ok1] = listdlg('PromptString','Select the new scan name:',...
    'Name', 'Select a Name',...
    'SelectionMode','Multiple',...
    'ListSize', [400 300],...
    'ListString',name_option);

if ok1 == 0
    return
end

user_response = questdlg('Do you really want to delete every data ?', 'Warning', 'Yes', 'No', 'Cancel', 'Cancel');
if strcmp(user_response, 'Cancel') || strcmp(user_response, 'No')
    return
end

ScansToRemove =name_option(ScansToRemove);
nii_index = [];
for i=1:numel(ScansToRemove)
    nii_index = [nii_index' find(handles.database.SequenceName == ScansToRemove(i))']';
end


MIA_remove_scan(hObject, eventdata, handles, nii_index)

msgbox('Done', 'Message') ;


% --------------------------------------------------------------------
function MIA_sort_name_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_sort_name_up_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_name_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
% store the information of the new order
handles.database.Properties.UserData.Order_data_display{1} =  'ascend';
%update the database
handles.database = sortrows(handles.database,{'Patient', 'Tp', 'SequenceName'},handles.database.Properties.UserData.Order_data_display);
% store the newdatabase
guidata(hObject, handles);
%update the dispay
MIA_update_database_display(hObject, eventdata, handles);


% --------------------------------------------------------------------
function MIA_sort_name_down_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_name_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
% store the information of the new order
handles.database.Properties.UserData.Order_data_display{1} =  'descend';
%update the database
handles.database = sortrows(handles.database,{'Patient', 'Tp', 'SequenceName'},handles.database.Properties.UserData.Order_data_display);
% store the newdatabase
guidata(hObject, handles);
%update the dispay
MIA_update_database_display(hObject, eventdata, handles);

% --------------------------------------------------------------------
function MIA_sort_scan_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_sort_time_point_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_time_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_sort_time_point_up_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_time_point_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end

% store the information of the new order
handles.database.Properties.UserData.Order_data_display{2} =  'ascend';
%update the database
handles.database = sortrows(handles.database,{'Patient', 'Tp', 'SequenceName'}, handles.database.Properties.UserData.Order_data_display);
% store the newdatabase
guidata(hObject, handles);
%update the dispay
MIA_update_database_display(hObject, eventdata, handles);


% --------------------------------------------------------------------
function MIA_sort_time_point_down_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_time_point_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end

% store the information of the new order
handles.database.Properties.UserData.Order_data_display{2} =  'descend';
%update the database
handles.database = sortrows(handles.database,{'Patient', 'Tp', 'SequenceName'},handles.database.Properties.UserData.Order_data_display);
% store the newdatabase
guidata(hObject, handles);
%update the dispay
MIA_update_database_display(hObject, eventdata, handles);


% --------------------------------------------------------------------
function MIA_sort_scan_up_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_scan_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
% store the information of the new order
handles.database.Properties.UserData.Order_data_display{3} =  'ascend';
%update the database
%% New Code
fields = fieldnames(handles.database);
% [~, indPat] = max(strcmp(fields, {'Patient'}));
% [~, indTp] = max(strcmp(fields, {'Tp'}));
% [~, indSeq] = max(strcmp(fields, {'SequenceName'}));
% Indexs = [indPat, indTp, indSeq];
[~, indSeq] = max(strcmp(fields, {'SequenceName'}));
Indexs = indSeq;
[~, index] = sortrows(cellstr(handles.database{:,:}), Indexs, handles.database.Properties.UserData.Order_data_display{3});
handles.database = handles.database(index,:);

%% Old Code
%handles.database = sortrows(handles.database,{'Patient', 'Tp', 'SequenceName'},handles.database.Properties.UserData.Order_data_display);
%%
% store the newdatabase
guidata(hObject, handles);
%update the dispay
MIA_update_database_display(hObject, eventdata, handles);



% --------------------------------------------------------------------
function MIA_sort_scan_down_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_scan_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = guidata(hObject);
if ~isfield(handles, 'database')
    return
end

% store the information of the new order
handles.database.Properties.UserData.Order_data_display{3} =  'descend';
%update the database
handles.database = sortrows(handles.database,{'Patient', 'Tp', 'SequenceName'},handles.database.Properties.UserData.Order_data_display);
% store the newdatabase
guidata(hObject, handles);
%update the dispay
MIA_update_database_display(hObject, eventdata, handles);



% --------------------------------------------------------------------
function MIA_menu_roi_fill_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_roi_fill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.MIA_menu_roi_fill, 'Checked'), 'on')
    set(handles.MIA_menu_roi_fill, 'Checked', 'off')
else
    set(handles.MIA_menu_roi_fill, 'Checked', 'on')
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MIA_scans_list.
function MIA_scans_list_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MIA_scans_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~strcmp(get(gcf,'SelectionType'),'open')
    %wenn normal click is happened, nothing to do
    return;
end
fprintf(1,'\nI am doing a single-click.\n\n');


% --------------------------------------------------------------------
function MIA_menu_rename_from_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_rename_from_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if ~isfield(handles, 'database')
    return
end
selection = questdlg('What do you want to rename?',...
    'Warning',...
    'Parameters','VOIs', 'Cancel', 'Cancel');

switch selection
    case 'Parameters'
        old_name_listing = cellstr(unique(handles.database.SequenceName(handles.database.Type == 'Scan')));
        new_name_listing = [cellstr(unique(handles.database.SequenceName(handles.database.Type == 'Scan')))' 'Other']';
    case 'VOIs'
        old_name_listing = cellstr(unique(handles.database.SequenceName(handles.database.Type == 'ROI')));
        new_name_listing = [cellstr(unique(handles.database.SequenceName(handles.database.Type == 'ROI')))' 'Other']';
    case 'Cancel'
        return
end
[old_scan_name,ok] = listdlg('Name', 'Question?', 'ListString', old_name_listing,'SelectionMode','single','ListSize', [500 500],...
    'PromptString', 'Select the parameters/VOIs you want to rename');
if ok == 0 || isempty(old_scan_name)
    return
end
old_scan_name = categorical(old_name_listing(old_scan_name));


[new_scan_name, ok1] = listdlg('PromptString','Select the new scan name:',...
    'Name', 'Select a Name',...
    'SelectionMode','single',...
    'ListSize', [400 300],...
    'ListString',new_name_listing);

if ok1 == 0
    return
end
if strcmp('Other',new_name_listing(new_scan_name)) == 1
    newparameter_name = inputdlg('Name of the new Scan ', 'Question?', 1, {''});
    newparameter_name = clean_variable_name(newparameter_name, '');
    newparameter_name =categorical(newparameter_name);
else
    newparameter_name =categorical(new_name_listing(new_scan_name));
end
idx_to_update = find(handles.database.SequenceName == old_scan_name);
%% update the database with the new name
% but first check if the new scan name does not exist for this patient and
% time point
for i = 1:numel(idx_to_update)
    % faire le ROI vs SCAN
    if find(handles.database.Patient == handles.database.Patient(idx_to_update(i)) &...
            handles.database.Tp == handles.database.Tp(idx_to_update(i)) & ...
            handles.database.SequenceName == newparameter_name) > 0
        msgbox('A Scan with the same name already exist for this patient at this time point') ;
        return
    end
    
    
    new_nii_filename = strrep(cellstr(handles.database.Filename(idx_to_update(i))), cellstr(handles.database.SequenceName(idx_to_update(i))), cellstr(newparameter_name));
    
    % rename the scan file .nii.gz
    if  exist(fullfilename(handles, idx_to_update(i), '.nii.gz'), 'file') == 0
        warning_text = sprintf('##$ This file no not exist\n##$ %s',...
            fullfilename(handles, idx_to_update(i), '.nii.gz'));
        msgbox(warning_text, 'rename file warning') ;
    elseif exist(string(strcat(cellstr(handles.database.Path(idx_to_update(i))),new_nii_filename{:},'.nii.gz')), 'file') == 2
        msgbox('The new .nii.gz file exist already!!') ;
        
    else
        movefile(fullfilename(handles, idx_to_update(i), '.nii.gz'), strcat(char(handles.database.Path(idx_to_update(i))),new_nii_filename{:},'.nii.gz'), 'f')
        if exist(fullfilename(handles, idx_to_update(i), '.json'), 'file') == 2
            movefile(fullfilename(handles, idx_to_update(i), '.json'), strcat(char(handles.database.Path(idx_to_update(i))),new_nii_filename{:},'.json'), 'f');
        end
    end
    
    % rename the scan file .nii
    if  exist(fullfilename(handles, idx_to_update(i), '.nii'), 'file') == 0
        warning_text = sprintf('##$ This file no not exist\n##$ %s',...
            fullfilename(handles, idx_to_update(i), '.nii'));
        msgbox(warning_text, 'rename file warning') ;
    elseif exist(string(strcat(cellstr(handles.database.Path(idx_to_update(i))),new_nii_filename{:},'.nii')), 'file') == 2
        msgbox('The new .nii file exist already!!') ;
        
    else
        movefile(fullfilename(handles, idx_to_update(i), '.nii'), strcat(char(handles.database.Path(idx_to_update(i))),new_nii_filename{:},'.nii'), 'f')
        if exist(fullfilename(handles, idx_to_update(i), '.json'), 'file') == 2
            movefile(fullfilename(handles, idx_to_update(i), '.json'), strcat(char(handles.database.Path(idx_to_update(i))),new_nii_filename{:},'.json'), 'f');
        end
    end
    
    
    % update the Filename field in the table
    handles.database.SequenceName(idx_to_update(i)) = newparameter_name;
    handles.database.Filename(idx_to_update(i)) = new_nii_filename;
    
    % save the structure
    guidata(hObject, handles);
end

%% update graph and display
MIA_update_database_display(hObject, eventdata, handles);

msgbox('Done', 'Message') ;



function MIA_AdjWL(varargin)
fh=varargin{1};
G=get(fh,'userdata');
G.cp=get(gca,'currentpoint');

G.x=G.cp(1,1);
G.y=G.cp(1,2);
G.xinit = G.initpnt(1,1);
G.yinit = G.initpnt(1,2);
G.dx = G.x-G.xinit;
G.dy = G.y-G.yinit;
G.clim = G.initClim+G.initClim(2).*[G.dx G.dy]./128;
try
    switch get(fh,'SelectionType')
        case 'extend' % Mid-button, shft+left button,
            %             'extend'
            set(findobj(fh,'Type','axes'),'Clim',G.clim);
        case 'alt' %right-click,ctrl+left button,
            %             'alt'
            set(gca,'Clim',G.clim);
    end
catch
    
end



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function MIA_GUI_WindowButtonUpFcn(hObject, eventdata, handles)

% hObject    handle to MIA_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if ~strcmp(get(handles.MIA_GUI,'SelectionType'),'normal')
    set(handles.MIA_GUI,'WindowButtonMotionFcn',   @(hObject,eventdata)MIA2('MIA_GUI_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject)));
    % save contrast
    handles.display_option.manual_contrast = 1;
    guidata(hObject, handles)
end



% --------------------------------------------------------------------
function MIA_copy_ScanVoi_to_other_tp_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_copy_ScanVoi_to_other_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_selected = finddata_selected(handles);
Tp_listing = unique(handles.database.Tp(handles.database.Patient == handles.database.Patient(data_selected)));
[time_point,ok] = listdlg('PromptString', 'Select 1 or several time point',...
    'Name', 'Question?',...
    'ListSize', [200 300],...
    'ListString', Tp_listing);
if ok == 0
    return
end
for i=1:numel(time_point)
    new_entry = handles.database(data_selected,:);
    new_entry.Tp = Tp_listing(time_point(i));
    new_entry.Filename  = categorical(cellstr(strcat(char(new_entry.Patient), '-', char(new_entry.Tp),'-', char(new_entry.SequenceName))));
    % check if the scan already exist for this time point
    if sum(ismember(handles.database, new_entry)) == 0
        if ~exist(fullfilename(handles, data_selected, '.nii'), 'file') && exist(fullfilename(handles, data_selected, '.nii.gz'), 'file')
            gunzip(fullfilename(handles, data_selected, '.nii.gz'));
            assert(exist(fullfilename(handles, data_selected, '.nii'), 'file')==2)
            delete(fullfilename(handles, data_selected, '.nii.gz'))
        end
        fid_nii=fopen(fullfilename(handles, data_selected, '.nii'),'r');
        if fid_nii>0
            fclose(fid_nii);
            % update the database
            handles.database(size(handles.database,1)+1,:) = new_entry;
            % create the new files
            copyfile(fullfilename(handles, data_selected, '.nii'), fullfilename(handles, size(handles.database,1), '.nii'), 'f');
            if exist(fullfilename(handles, data_selected, '.json'), 'file') == 2
                copyfile(fullfilename(handles, data_selected, '.json'), fullfilename(handles, size(handles.database,1), '.json'), 'f');
            end
        else
            warndlg(['something is wrong this the data : ' fullfilename(handles, data_selected, '.nii')],  'Warning');
            return
        end
    else
        warndlg(['This scan already exist : ' char(new_entry.Filename)],  'Warning');
        return
    end
    
end
%update handes
guidata(hObject, handles)

msgbox('Done', 'logbook') ;


function MIA_warning_duplicate_scan_fcn(handles, parameters,patient_nbr,time_point)

warndlg(['Several scan have the same name: ' parameters{:} ' For the patient ' handles.database(patient_nbr).name ' at day '...
    handles.database(patient_nbr).day(time_point).date] , 'Warning');


function MIA_warning_duplicate_VOI_fcn(handles, VOI,patient_nbr,time_point)

warndlg(['Several VOI have the same name: ' VOI{:} ' For the patient ' handles.database(patient_nbr).name ' at day '...
    handles.database(patient_nbr).day(time_point).date] , 'Warning');



% --------------------------------------------------------------------
function MIA_menu_load_data_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'database')
    MIA_data_path = handles.database.Properties.UserData.MIA_data_path;
    %     MIA_data_path = handles.database.Properties.UserData.MIA_data_path ;
else
    MIA_root_path = uigetdir(pwd, 'Select the directory to save the your new projet');
    if sum(MIA_root_path) == 0
        return
    end
    
    MIA_root_path = [MIA_root_path filesep];
    %% create the output folder ('MIA_data')
    MIA_data_path = MIA_root_path;%[MIA_root_path, 'MIA_data', filesep];
    % Create the output folder if needed
    if exist(MIA_data_path, 'dir') ~= 7
        status = mkdir(MIA_data_path);
        if status == 0
            warndlg('You do not have the right to write in the folder!', 'Warning');
            return
        end
    end
    % Create the RAW data folder if needed
    if exist([MIA_data_path, 'Raw_data', filesep], 'dir') ~= 7
        status = mkdir([MIA_data_path, 'Raw_data', filesep]);
        if status == 0
            warndlg('You do not have the right to write in the folder!', 'Warning');
            return
        end
    end
    
    
    %% create the database structure
    handles.database =  table;%
    handles.database = cell2table(cell(0,8));
    handles.database.Properties.VariableNames = {'Group','Patient', 'Tp', 'Path', 'Filename', 'Type', 'IsRaw', 'SequenceName'};
    handles.database.Properties.UserData.MIA_data_path = MIA_data_path;
    handles.database.Properties.UserData.MIA_Raw_data_path = [MIA_data_path, 'Raw_data', filesep];
    handles.database.Properties.UserData.MIA_Derived_data_path = [MIA_data_path, 'Derived_data', filesep];
    handles.database.Properties.UserData.MIA_ROI_path = [MIA_data_path, 'ROI_data', filesep];
    handles.database.Properties.UserData.Order_data_display = {'ascend','ascend','ascend'};
    
end
%% create the tmp folder
MIA_tmp_folder = [MIA_data_path, 'tmp'];
if exist(MIA_tmp_folder, 'dir') ~= 7
    status = mkdir(MIA_tmp_folder);
    if status == 0
        warndlg('You do not have the right to write in the folder!', 'Warning');
        return
    end
end

%update handes
guidata(hObject, handles)

% [folder, ~, ~] = fileparts(which('MIA2.m'));
% MRIFileManager_path = [folder,filesep,'MRIFileManager', filesep];
% repProject='[ExportNifti] C:\\Users\\omonti\\Desktop\\tmpOMPLOUI';
% namExport='[ExportToMIA] PatientName-StudyName-CreationDate-SeqNumber-Protocol-SequenceName';
% MRIFileManager.FileManagerFrame.main({MIA_tmp_folder_for_java, namExport}

MIA_tmp_folder_for_java = strrep(MIA_tmp_folder, filesep, [filesep filesep]);
MIA_tmp_folder_for_java = ['[ExportNifti]' MIA_tmp_folder_for_java];
% save hObject, hObject and handles in global variable
global hObjectb;
global eventdatab;
global handlesb;
hObjectb=hObject;
eventdatab=eventdata;
handlesb=handles;

namExport =  '[ExportToMIA]PatientName-StudyName-CreationDate-SeqNumber-Protocol-SequenceName-AcquisitionTime';
MRIFileManager.FileManagerFrame.main({MIA_tmp_folder_for_java, namExport handles.original_Java_LookAndFeel, 'NoExitSystem'}) %add lookAndFeel option
%
% system(char(strcat([MRIFileManager_path 'jre/bin/java -jar '], [' ' MRIFileManager_path],...
%     'MRIManager.jar [ExportNifti] ', MIA_tmp_folder_for_java, {' '}, namExport)));

function back_from_MRIManager(data_loaded)
% retrived hObject, hObject and handles in global variable
global hObjectb;
global eventdatab;
global handlesb;
hObject = hObjectb;
eventdata=eventdatab;
handles=handlesb;
clear hObjectb
clear eventdatab
clear handlesb

data_loaded{1} = strrep(data_loaded{1}, 'null', '"''''"');
data_loaded{1} = strrep(data_loaded{1}, '""', '"''''"');
%javax.swing.UIManager.setLookAndFeel(newLnF);
handles = guidata(handles.MIA_GUI);
MIA_tmp_folder = [handles.database.Properties.UserData.MIA_data_path, 'tmp'];

try
    log_file =struct2table(jsondecode(char(data_loaded)));
    
catch ME
    if strcmp(ME.message,'Fields in a scalar structure must have the same number of rows.')
        warndlg('Something is wrong with MRImanager: the output table does not have the same number of rows');
        return
        
    end
end


% log_file_to_update = log_file;
log_file.StudyName = categorical(cellstr(log_file.StudyName));
log_file.CreationDate = categorical(cellstr(log_file.CreationDate));

StudyName_listing = unique(log_file.StudyName);
for i = 1:numel(unique(log_file.StudyName))
    patient_filter = log_file.StudyName == char(StudyName_listing(i));
    CreationDate_listing = unique(log_file.CreationDate(patient_filter));
    for j = 1:numel(CreationDate_listing)
        CreationDate_filter = log_file.CreationDate == char(CreationDate_listing(j));
        index_data_to_import =findn(patient_filter & CreationDate_filter);
        index_data_to_import = index_data_to_import(:,1);
        filesch =  [{char(unique(log_file.StudyName(index_data_to_import)))}, {char(unique(log_file.CreationDate(index_data_to_import)))},...
            {char(log_file.NameFile(index_data_to_import(1),:))}];
        [name_selected, tp_selected]  = Add_Tag_to_new_patient_GUI(handles,filesch);
        if ~isempty(name_selected)
            for m=1:numel(index_data_to_import)
                NAME = char(log_file.NameFile(index_data_to_import(m),:));
                if exist(fullfile(MIA_tmp_folder, [NAME, '.json']), 'file')
                    json_data = spm_jsonread(fullfile(MIA_tmp_folder, [NAME, '.json']));
                    if ~isfield(json_data, 'ProtocolName')
                        json_data.ProtocolName.value = {'Undefined'};
                    elseif isempty(char(json_data.ProtocolName.value))
                        json_data.ProtocolName.value = {'Undefined'};
                    end
                    if ~isempty(handles.database)
                        %% check if a scan with the same SequenceName exist for this patient at this time point. If so, add suffix to the SequenceName (ie. SequenceName(X)
                        if sum(handles.database.Patient ==  char(name_selected) & handles.database.Tp ==  char(tp_selected) &  handles.database.SequenceName == char(clean_variable_name(char(json_data.ProtocolName.value), ''))) == 1
                            nbr_of_seq = sum(handles.database.Patient ==  char(name_selected) &...
                                handles.database.Tp ==  char(tp_selected) &...
                                strncmp(cellstr(handles.database.SequenceName), char(clean_variable_name(char(json_data.ProtocolName.value), '')), length(char(clean_variable_name(char(json_data.ProtocolName.value), '')))-1));
                            seq_name = [char(clean_variable_name(char(json_data.ProtocolName.value), '')) '(' num2str(nbr_of_seq+1) ')'];
                            file_name = strcat(name_selected, '-', tp_selected,'-',seq_name,'_',datestr(now,'yyyymmdd-HHMMSSFFF'));
                        else
                            
                            seq_name = clean_variable_name(char(json_data.ProtocolName.value), '');
                            file_name = strcat(name_selected, '-', tp_selected,'-',seq_name,'_',datestr(now,'yyyymmdd-HHMMSSFFF'));
                        end
                    else
                        seq_name = clean_variable_name(char(json_data.ProtocolName.value), '');
                        file_name = strcat(name_selected , '-', tp_selected,'-',seq_name,'_',datestr(now,'yyyymmdd-HHMMSSFFF'));
                        
                    end
                    new_data = table(categorical(cellstr('Undefined')), categorical(cellstr(name_selected)), categorical(cellstr(tp_selected)), categorical(cellstr(handles.database.Properties.UserData.MIA_Raw_data_path)), categorical(cellstr(file_name)),...
                        categorical(cellstr('Scan')), categorical(1), categorical(cellstr(seq_name)),...
                        'VariableNames', {'Group','Patient', 'Tp', 'Path', 'Filename', 'Type', 'IsRaw', 'SequenceName'});
                    
                    %% add data to the database
                    handles.database = [handles.database; new_data];
                    
                    %%save the files (nii + json)
                    movefile(fullfile(MIA_tmp_folder, [NAME, '.nii']), [char(handles.database.Path(end)), char(handles.database.Filename(end)), '.nii']);
                    movefile(fullfile(MIA_tmp_folder, [NAME, '.json']), [char(handles.database.Path(end)), char(handles.database.Filename(end)), '.json']);
                    
                    guidata(hObject, handles);
                end
            end
        end
    end
end

% delete temps files and folder
rmdir(MIA_tmp_folder, 's')

guidata(hObject, handles);

MIA_update_database_display(hObject, eventdata, handles);

MIA_menu_save_database_Callback(hObject, eventdata, handles)
msgbox('Import Done and database saved!', 'Message') ;


% --------------------------------------------------------------------
function MIA_menu_Clean_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Clean_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%% create a listing of all file present in the database


Data_path = handles.database.Properties.UserData.MIA_data_path;

%% Check if all entries are correctly linked to valid files.

ValidEntries = [];
InvalidEntries = [];
ValidNiftiFiles = {};
InvalidNiftiFiles = {};
ValidJsonFiles = {};
InvalidJsonFiles = {};
ValidMatFiles = {};
InvalidMatFiles = {};
for i=1:height(handles.database)
    if handles.database(i,:).Type == 'Scan'
        Nifti_file = [char(handles.database(i,:).Path), char(handles.database(i,:).Filename), '.nii'];
        Nifti_file_compressed = strrep(Nifti_file, '.nii', '.nii.gz');
        Json_file = [char(handles.database(i,:).Path), char(handles.database(i,:).Filename), '.json'];
        
        if (exist(Nifti_file, 'file')==2 || exist(Nifti_file_compressed, 'file')==2) && exist(Json_file, 'file') == 2
            ValidEntries = [ValidEntries, i];
            ValidNiftiFiles = [ValidNiftiFiles, {Nifti_file}];
            ValidJsonFiles = [ValidJsonFiles, {Json_file}];
        else
            InvalidEntries = [InvalidEntries, i];
            if exist(Nifti_file, 'file')~=2
                InvalidNiftiFiles = [ValidNiftiFiles, {Nifti_file}];
            elseif exist(Json_file, 'file')~=2
                InvalidJsonFiles = [InvalidJsonFiles, {Json_file}];
            end
        end
    elseif handles.database(i,:).Type == 'ROI'
        Nifti_file = [char(handles.database(i,:).Path), char(handles.database(i,:).Filename), '.nii'];
        Nifti_file_compressed = strrep(Nifti_file, '.nii', '.nii.gz');
        if exist(Nifti_file, 'file')==2 || exist(Nifti_file_compressed, 'file')==2
            ValidEntries = [ValidEntries, i];
            ValidNiftiFiles = [ValidNiftiFiles, {Nifti_file}];
        else
            InvalidEntries = [InvalidEntries, i];
            InvalidNiftiFiles = [ValidNiftiFiles, {Nifti_file}];
        end
    elseif handles.database(i,:).Type == 'Cluster'
        Nifti_file = [char(handles.database(i,:).Path), char(handles.database(i,:).Filename), '.nii'];
        Nifti_file_compressed = strrep(Nifti_file, '.nii', '.nii.gz');
        Mat_file = [char(handles.database(i,:).Path), char(handles.database(i,:).Filename), '.mat'];
        
        if (exist(Nifti_file, 'file')==2 || exist(Nifti_file_compressed, 'file')==2) && exist(Mat_file, 'file') == 2
            ValidEntries = [ValidEntries, i];
            ValidNiftiFiles = [ValidNiftiFiles, {Nifti_file}];
            ValidMatFiles = [ValidMatFiles, {Mat_file}];
        else
            InvalidEntries = [InvalidEntries, i];
            if exist(Nifti_file, 'file')~=2
                InvalidNiftiFiles = [ValidNiftiFiles, {Nifti_file}];
            elseif exist(Json_file, 'file')~=2
                InvalidMatFiles = [InvalidMatFiles, {Mat_file}];
            end
        end
        
    else
        warning('Types allowed : Scan, ROI, Cluster. We found an unknown type : %s', handles.database(i,:).Type)
    end
end

%% If there is some invalid entries, the user can choose to delete all of them.
if ~isempty(InvalidEntries)
    text = ['We found some entries in the database for which file(s) are missing in ', Data_path, '. Do you want to delete them from the table ?'];
    answer = questdlg(text, 'Wrong entries', 'Yes', 'See wrong entries','No', 'No');
    %uiwait(gcf)
    %function answer = DeleteWrongEntries(PushButton, EventData)
    if strcmp(answer,'See wrong entries')
        f = figure;
        t = uitable(f, 'Position', [30,100,500,300],'ColumnName', handles.database.Properties.VariableNames,'Data',cellstr(handles.database{InvalidEntries, :}));
        %btnDelete = uicontrol('Parent', f, 'Position', [100,50,100,50], 'String', 'Delete All', 'Callback', 'rep = ''Yes''; delete(gcf)');
        %btnCancel = uicontrol('Parent', f, 'Position', [300,50,100,50], 'String', 'Cancel', 'Callback', 'rep = ''No''; delete(gcf)');
        uiwait(f)
        answer = questdlg('So, should I delete them ?', 'It''s time to choose', 'Yes', 'No', 'No');
    end
    if strcmp(answer, 'Yes')
        handles.database(InvalidEntries,:) = [];
        ValidFiles = [ValidNiftiFiles, ValidJsonFiles, ValidMatFiles];
        msgbox('Done')
    else
        ValidFiles = [ValidNiftiFiles, ValidJsonFiles, InvalidNiftiFiles, InvalidJsonFiles, ValidMatFiles, InvalidMatFiles];
    end
else
    msgbox('All your database''s entries are correctly linked to the files !')
    ValidFiles = [ValidNiftiFiles, ValidJsonFiles, ValidMatFiles];
end

% ValidFiles is the list of all files linked to the database. Other files
% in DirectoriesToProcess are useless.


DirectoriesToProcess = {'Derived_data', 'Raw_data', 'ROI_data'};


if exist([Data_path, 'To_Trash'],'dir') ~= 7
    [status, ~, ~] = mkdir([Data_path, 'To_Trash']);
    if status == false
        error('Cannot create the To_Thrash folder to move the useless files in.')
    end
end

MovedFiles = {};
for i=1:length(DirectoriesToProcess)
    Files = dir([Data_path, DirectoriesToProcess{i}]);
    Files = struct2cell(Files);
    for j=1:size(Files,2)
        file = Files(:,j);
        if ~file{5}
            if ~contains(strrep([file{2}, filesep, file{1}], '.nii.gz', '.nii'), ValidFiles) 
                [status, msg] = movefile([file{2}, filesep, file{1}], [Data_path,  'To_Trash']);
                if status
                    MovedFiles = [MovedFiles; {[file{2}, filesep, file{1}]}];
                else
                    warning('Error when moving the useless files. %s', msg)
                end
            end
        end
    end
end


if isempty(MovedFiles)
    msgbox('The folders containing your database''s files are already pretty clean !')
    %rmdir([Data_path, 'To_Trash']);
else
    text = [num2str(length(MovedFiles)), ' files were just moved from your data folders to the To_Trash folder.'];
    answer = questdlg(text, 'Moved Files', 'OK !', 'See files moved','OK !');
    if strcmp(answer,'See files moved')
        f = figure;
        t = uitable(f, 'Position', [30,100,500,300],'ColumnName', {'Moved Files'},'Data',cellstr(MovedFiles), 'ColumnWidth',{800});
        btnCancel = uicontrol('Parent', f, 'Position', [225,40,100,50], 'String', 'OK', 'Callback', 'delete(gcf)');
    end
end
guidata(hObject, handles)



% --------------------------------------------------------------------
function MIA_menu_Importing_form_another_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Importing_form_another_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dir = uigetdir(pwd, 'Select the project to import');
if dir == 0
    return
end
DatabPath = [dir, filesep, 'MIA_database.mat'];
other_database = load(DatabPath);
other_database = other_database.database;
UPatients = unique(other_database.Patient);
data_listing = {};
linked_table = {};
for i=1:length(UPatients)
    tmp_database = other_database(other_database.Patient == UPatients(i),:);
    UTp = unique(tmp_database.Tp);
    for j=1:length(UTp)
        tmp2_database = table();
        data_listing = [data_listing; {[char(UPatients(i)), ' - ', char(UTp(j))]}];
        tmp2_database = tmp_database(tmp_database.Tp == UTp(j),:);
        linked_table = [linked_table; {tmp2_database}];
    end
end




[data_to_import, ok] = listdlg('PromptString','Select the data you want to import',...
    'Name', 'Select a data',...
    'ListSize', [400 300],...
    'ListString',data_listing);

if ok == 0
    return
end

database_to_import = table();
for i = 1:length(data_to_import)
    database_to_import = [database_to_import; linked_table{data_to_import(i)}];
end

for i=1:size(database_to_import,1)
    switch char(database_to_import.Type(i))
        case 'Scan'
            if char(database_to_import.IsRaw(i)) == '1'
                outfolder = handles.database.Properties.UserData.MIA_Raw_data_path;
                
            else
                outfolder = handles.database.Properties.UserData.MIA_Derived_data_path;
            end
            % Copy nifti file 
            extension = '.nii.gz';
            filename = [char(database_to_import.Path(i)), char(database_to_import.Filename(i)), extension];
            if ~exist(filename, 'file')
                extension = '.nii';
                filename = [char(database_to_import.Path(i)), char(database_to_import.Filename(i)), extension];
            end
            database_to_import.Filename(i) = categorical(cellstr([char(database_to_import.Patient(i)), '_', char(database_to_import.Tp(i)), '_', char(database_to_import.SequenceName(i))]));
            outfilename = [outfolder, char(database_to_import.Filename(i)), extension];
            original_filename = char(database_to_import.Filename(i));
            idx2=2;
            while exist(outfilename)
                database_to_import.Filename(i) = categorical(cellstr([original_filename, '_', num2str(idx2)]));
                outfilename = [outfolder, char(database_to_import.Filename(i)), extension];
                idx2=idx2+1;
            end
            status1 = copyfile(filename, outfilename);
%             if ~exist([outfolder, char(database_to_import.Filename(i)), '.nii'])
%                 status1 = copyfile(filename, outfolder);
%             else
%                 text = ['The file: ', filename, ' hasn''t been imported because there is another file with the same name in the folder: ', outfolder];
%                 msgbox(text);
%                 status1 = 0;
%             end
            % Copy json file
            filename = strrep(filename, extension, '.json');
            outfilename = strrep(outfilename, extension, '.json');
            while exist(outfilename)
                outfilename = [outfolder, char(database_to_import.Filename(i)), '.json'];
                idx2=idx2+1;
            end
            status2 = copyfile(filename, outfilename);
%             if ~exist([outfolder, char(database_to_import.Filename(i)), '.json'])
%                 status2 = copyfile(filename, outfolder);
%             else
%                 text = ['The file: ', filename, ' hasn''t been imported because there is another file with the same name in the folder: ', outfolder];
%                 msgbox(text);
%                 status2 = 0;
%             end
            

            
        case {'ROI', 'Cluster'}
            extension = '.nii.gz';
            outfolder = handles.database.Properties.UserData.MIA_ROI_path;
            filename = [char(database_to_import.Path(i)), char(database_to_import.Filename(i)), extension];
            if ~exist(filename, 'file')
                extension = '.nii';
                filename = [char(database_to_import.Path(i)), char(database_to_import.Filename(i)), extension];
            end
            database_to_import.Filename(i) = categorical(cellstr([char(database_to_import.Patient(i)), '_', char(database_to_import.Tp(i)), '_', char(database_to_import.SequenceName(i))]));
            original_filename = char(database_to_import.Filename(i));
            outfilename = [outfolder, char(database_to_import.Filename(i)), extension];
            idx2=2;
            while exist(outfilename)
                database_to_import.Filename(i) = categorical(cellstr([original_filename, '_', num2str(idx2)]));
                outfilename = [outfolder, char(database_to_import.Filename(i)), extension];
                idx2=idx2+1;
            end
            status1 = copyfile(filename, outfilename);
           
            
    end
    if status1 && status2
        tags = database_to_import(i,:);
        tags.Path = categorical(cellstr(outfolder));
        flag = 1;
        idx = 2;
        while flag
            tmpdatab = handles.database(handles.database.Patient == tags.Patient,:);
            tmpdatab = tmpdatab(tmpdatab.Tp == tags.Tp, :);
            %tmpdatab = tmpdatab(tmpdatab.Type == tags.Type, :); % ???????
            tmpdatab = tmpdatab(tmpdatab.SequenceName == tags.SequenceName, :);
            if size(tmpdatab,1) == 0
                flag = 0;
            else
                tags.SequenceName = [char(tags.SequenceName), '_', num2str(idx)];
                idx = idx+1;
            end
        end
        handles.database = [handles.database; tags];
    end
    
end


title = ['database name: ',filename,'; ', num2str(numel(handles.database)), ' files.'];
set(handles.MIA_GUI, 'Name', title);
set(handles.MIA_name_list, 'Value', 1);


guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);
msgbox('Done', 'Message') ;





% --- Executes on button press in MIA_test_button.
function MIA_test_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_test_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end


% %% code to copy plotted data to clipboard
% data(:,1) = get(get(handles.MIA_plot1, 'Children'), 'XData');
% data(:,2) = get(get(handles.MIA_plot1, 'Children'), 'YData');
% data(isnan(data(:,1)),:) = [];
% data(isnan(data(:,2)),:) = [];
% mat2clipboard(data)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% handles.database.Type(handles.database.SequenceName == 'Mask') = 'ROI';
% guidata(hObject, handles)
%mat2clipboard(get(handles.MIA_table1, 'Data'));
%msgbox('Data copied!', 'logbook') ;


%
% %% code to extact fingerprint from cluster-ROI
% color_RGB = ...
%     [ 1 0 0; ...                  Red
%     0 1 0; ...                  Lime
%     0 0 1; ...                  Blue
%     1 1 0;  ...                 Yellow
%     0 1 1; ...                  Aqua
%     1 0 1; ...                  Fuchsia
%     0 0.5 0.5; ...              Teal
%     1 0.6471 0; ...             Orange
%     0.5 0 0.5; ...              Purple
%     0 0.5 0; ...                Green
%     0.5 0 0; ...                Maroon
%     0 0 0.5; ...                Navy
%     0.5 0.5 0; ...              Olive
%     0.7529 0.7529 0.7529; ...   Silver
%     0.5 0.5 0.5; ...            Gray
%     0 0 0; ...                 Black
%     1 1 1];                     %White
%
% cluster_name = 'Tumeur-clusterGMM';
% number_of_class = 3;
% fingerprint = dataset;
% row = 1;
% for i = 1:numel(handles.database)
%     for j = 1:numel(handles.database(i).day)
%         ROI_nbr = strcmp(cluster_name, handles.database(i).day(j).VOIs) == 1;
%         ROI_filename = strcat(handles.database(i).path, handles.database(i).day(j).VOIs_file(ROI_nbr));
%         fid=fopen(ROI_filename{:} ,'r');
%         if fid>0
%             fclose(fid);
%             load(ROI_filename{:});
%             fingerprint.patient{row,1} = handles.database(i).name;
%             fingerprint.tp{row,1} = handles.database(i).day(j).date;
%             tmp_voxels = [];
%             for sl = 1:numel(uvascroi)
%                 tmp_voxels = cat(1, tmp_voxels, reshape(uvascroi(sl).value, [size(uvascroi(sl).value, 1) * size(uvascroi(sl).value, 2) size(uvascroi(sl).value, 3)]));
%             end
%             roi_size = sum(~isnan(tmp_voxels(:,1)));
%             for k = 1 : number_of_class
%                  tmp = tmp_voxels == repmat(color_RGB(k,:), [size(tmp_voxels, 1), 1]);
%                  eval(['fingerprint.cluster' num2str(k) '(row,1) = sum(tmp(:,1) == 1 & tmp(:,2) == 1 & tmp(:,3) == 1) / roi_size * 100;']);
%             end
%
%            row = row+1;
%         end
%     end
% end
% fingerprint

%% apply GMM model (clustering)

% %% texture analysis
% if isfield(handles, 'data_selected_resized') && isfield(handles, 'ROI_selected_resized')
% vox_size = handles.data_selected_resized.image.reco.fov(1)/handles.data_selected_resized.image.reco.no_samples;
% %[ROIonly,levels,ROIbox,maskBox] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,textType,quantAlgo,Ng)
% [ROIonly,~,~,~] = prepareVolume(...
%     squeeze(handles.data_selected_resized.image.reco.data),...
%     squeeze(handles.ROI_selected_resized.data.value),...
%     'MRscan',...
%     vox_size,...
%     0.8,...
%     1,...
%     handles.data_selected_resized.image.reco.thickness,...
%     'Global');
% [Metric_global] = getGlobalTextures(ROIonly,100)
%
% clear ROIonly
%
% [ROIonly,levels,~,~] = prepareVolume(...
%     squeeze(handles.data_selected_resized.image.reco.data),...
%     squeeze(handles.ROI_selected_resized.data.value),...
%     'MRscan',...
%     vox_size,...
%     0.8,...
%     1,...
%     handles.data_selected_resized.image.reco.thickness,...
%     'Matrix',...
%     'Equal',...
%     32);
%
% [GLCM] = getGLCM(ROIonly,levels);
% [Metric_GLCM] = getGLCMtextures(GLCM)
% [GLRLM] = getGLRLM(ROIonly,levels);
% [Metric_GLRLM] = getGLRLMtextures(GLRLM)
% [GLSZM] = getGLSZM(ROIonly,levels);
% [Metric_GLSZM] = getGLSZMtextures(GLSZM)
% [NGTDM,countValid] = getNGTDM(ROIonly,levels);
% [Metric_NGTDM] = getNGTDMtextures(NGTDM,countValid)
%
% % disp(struct2dataset(catstruct(Metric_global,Metric_GLCM,Metric_GLRLM,Metric_GLSZM,Metric_NGTDM)))
% %  msgbox(struct2table(Metric_global), 'texture Analysis') ;
%
% end

% for i = 1:numel(handles.database)
%     set(handles.MIA_name_list, 'Value', i);
%     for j=1:numel(handles.database(i).day)
%         set(handles.MIA_time_points_list, 'Value', j)
%
%
%         % set and open scan
%         set(handles.MIA_scan_VOIs_button, 'Value', 0)
%         MIA_scan_VOIs_button_Callback(hObject, eventdata, handles);
%         set(handles.MIA_scans_list, 'Value', find(strcmp(handles.database(i).day(j).parameters', {'09_MGEFIDSEpre-60ms'}) | ...
%             strcmp(handles.database(i).day(j).parameters', {'12_MGEFIDSEpost-60ms'}) > 0))
%         MIA_load_axes_Callback(hObject, eventdata, handles)
%
%
%         % set and open VOI
%         set(handles.MIA_scan_VOIs_button, 'Value', 1)
%         MIA_scan_VOIs_button_Callback(hObject, eventdata, handles);
%
%         set(handles.MIA_scans_list, 'Value', find(strcmp(handles.database(i).day(j).VOIs', {'Tumor'}) | ...
%             strcmp(handles.database(i).day(j).VOIs', {'Striat'}) > 0))
%         MIA_load_axes_Callback(hObject, eventdata, handles)
%
%
%         % Move to slice 3
%         set(handles.MIA_slider_slice, 'Value', 3)
%         MIA_slider_slice_Callback(hObject, eventdata, handles)
%
%         % clic on image
%         MIA_clic_on_image(get(hObject, 'Title'), eventdata, handles)
%
%     end
% end

% % code to erode an ROI
% se = strel('disk',4);
% % if 256??? --> se = strel('disk',5);
% for i=1:numel(handles.ROI_selected.data)
%     handles.ROI_selected.data(i).value = imerode(handles.ROI_selected.data(i).value,se);
%     %          handles.ROI_selected.data(i).value = imdilate(handles.ROI_selected.data(i).value,se);
% end
%
% uvascroi = handles.ROI_selected.data;
% save(char(handles.data_selected.roifile),'uvascroi');
% return
% % code pour ligia
% if numel(handles.ROI_selected) ~= 2
%     warndlg('Please open 2 ROI');
%     return
% end
% ROI1_name = handles.ROI_selected(1).data(1).name;
% ROI2_name = handles.ROI_selected(2).data(1).name;
%
% if ischar(ROI1_name)
%     ROI1_name = {ROI1_name};
% end
% if ischar(ROI2_name)
%     ROI2_name = {ROI2_name};
% end
% [Option,ok] = listdlg('Name', 'What do you want to do?', 'ListSize', [300 200], 'ListString', [[ROI1_name{:}, ' - ', ROI2_name{:}]; [ROI2_name{:}, ' - ', ROI1_name{:}]]);
%
% if ok ~= 1
%     return
% end
%
% switch Option
%     % ROI1 - ROI2
%     case 1
%         new_ROI = handles.ROI_selected_resized(1).data.value;
%         new_ROI(handles.ROI_selected_resized(2).data.value==1) = 0;
%         new_ROI_name = [ROI1_name{:}, '-', ROI2_name{:}];
%         % ROI2 - ROI1
%     case 2
%         new_ROI = handles.ROI_selected_resized(2).data.value;
%         new_ROI(handles.ROI_selected_resized(1).data.value==1) = 0;
%         new_ROI_name = [ROI2_name{:}, '-', ROI1_name{:}];
% end
%
% % save new ROI
%
% for i=1:get(handles.MIA_slider_slice, 'Max')
%     uvascroi(i).value = squeeze(new_ROI(:,:,i));
%     uvascroi(i).area_in_pixel=sum(sum(uvascroi(i).value));
%     uvascroi(i).position=[0,0; 0, handles.resolution_selected; handles.resolution_selected, handles.resolution_selected; handles.resolution_selected, 0];
%     uvascroi(i).name={new_ROI_name};
%     uvascroi(i).colorindex=1;
%     %%% case roi drown on histo data
%     uvascroi(i).displayedecho=handles.data_selected.image(1).reco.displayedecho;
%     uvascroi(i).displayedslice=handles.data_selected.image(1).reco.displayedslice;
%     uvascroi(i).displayedexpt=handles.data_selected.image(1).reco.displayedexpt;
%     uvascroi(i).fov_orientation=handles.data_selected.image(1).reco.fov_orientation(1:9,handles.data_selected.image(1).reco.displayedecho,i,handles.data_selected.image(1).reco.displayedexpt);
%     uvascroi(i).no_samples=handles.data_selected.image(1).reco.no_samples;
%     if isfield(handles.data_selected.image(1).reco, 'no_views')
%         uvascroi(i).no_views=handles.data_selected.image(1).reco.no_views;
%     end
%     uvascroi(i).fov=handles.data_selected.image(1).acq.fov;
%     if isfield(handles.data_selected.image(1).reco, 'no_views')
%         uvascroi(i).area_in_mmxmm=(uvascroi(i).fov(1)/uvascroi(i).no_samples)*(uvascroi(i).fov(2)/uvascroi(i).no_views)*uvascroi(i).area_in_pixel;
%     end
%     if isfield(handles.data_selected.image(1), 'reco_number')
%         uvascroi(i).reco_number=handles.data_selected.image(1).reco_number;
%     end
%     if isfield(handles.data_selected.image(1), 'scan_number')
%         uvascroi(i).scan_number=handles.data_selected.image(1).scan_number;
%     end
%     uvascroi(i).filename=handles.data_selected.image(1).filename;
%     uvascroi(i).filenames=sprintf('%s\n',handles.data_selected.image(1).filename);
%     uvascroi(i).fov_offsets=handles.data_selected.image(1).reco.fov_offsets(1:3,handles.data_selected.image(1).reco.displayedecho,i,handles.data_selected.image(1).reco.displayedexpt);
%
%     %     end
%     uvascroi(i).thickness = handles.data_selected.image(1).reco.thickness;
%     uvascroi(i).date=date;
%     uvascroi(i).affiche=0;
% end
% if ~isempty(uvascroi)
%     pathname = [handles.database(handles.data_selected.patient_info.name_nbr).path handles.database(handles.data_selected.patient_info.name_nbr).name handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).date '-VOI-' new_ROI_name '.mat'];
%     save(pathname,'uvascroi');
%     handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs = [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs {new_ROI_name}];
%     voi_file = [handles.database(handles.data_selected.patient_info.name_nbr).name handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).date '-VOI-' new_ROI_name '.mat'];
%     handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs_file = [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs_file {voi_file}];
% end
% % add cluster to the database
% handles.VOIs(end) = {new_ROI_name};
% handles.VOIs = unique(handles.VOIs);
% handles.VOIs(end+1) = {'Other'};
% %%
%
%
% guidata(hObject, handles);
% MIA_update_database_display(hObject, eventdata, handles)




% --- Executes on selection change in MIA_colormap_popupmenu.
function MIA_colormap_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_colormap_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_colormap_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_colormap_popupmenu
if ~isfield(handles, 'data_loaded') && ~isfield(handles, 'data_selected_for_PRM')
    return
end
% get(handles.MIA_slider_slice)
% set(handles.MIA_slider_slice, 'Value', round(get(handles.MIA_slider_slice, 'Value')));


if isfield(handles, 'data_displayed')
    colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
    for i=1:size(handles.data_displayed.image,4)
        colormap(handles.(sprintf('MIA_data%d', i)), colormap_selected{:});
    end
end

% --- Executes during object creation, after setting all properties.
function MIA_colormap_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_colormap_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MIA_data1_expt_slider.
function MIA_slider_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MIA_data1_expt_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
if ~isfield(handles, 'data_loaded') && ~isfield(handles, 'data_selected_for_PRM')
    return
end

Slider_min = get(hObject,'Min');
Slider_max = get(hObject,'Max');
Position = get(hObject,'Position');

SliderBarWidth = Position(3)/Slider_max;
set(hObject,'UserData',[Slider_min Slider_max SliderBarWidth Position]);

cp = get(handles.MIA_GUI,'CurrentPoint');
newValue = round((cp(1,1)-Position(1))/ SliderBarWidth); %0.01 = width of the arrow
if  newValue ==  get(hObject,'Value')
    return
elseif newValue > Slider_max
    newValue = Slider_max;
elseif newValue < Slider_min
    newValue = Slider_min;
end

set(hObject,'Value',newValue);

MIA_update_axes(hObject, eventdata, handles)
%
% set(handles.MIA_GUI,'WindowButtonMotionFcn',{@MIA_slider_on_move,handles})
% set(handles.MIA_GUI,'WindowButtonUpFcn',{@MIA_slider_release_click,handles})
% set(handles.MIA_GUI,'WindowButtonMotionFcn', @(hObject,eventdata)MIA2('MIA_slider_on_move',hObject,eventdata,guidata(hObject)))
% set(handles.MIA_GUI,'WindowButtonUpFcn',
% @(hObject,eventdata)MIA2('MIA_slider_release_click',hObject,eventdata,guidata(hObject)));b
%
%
% function MIA_slider_on_move(hObject, eventdata, handles)
% UserData = get(gco,'UserData');
% if ~isempty(UserData),
%     cp = get(hObject,'CurrentPoint');
%     newValue = round((cp(1,1)-UserData(4))/UserData(3));
%     if  newValue ==  get(gco,'Value')
%         return
%     elseif newValue > UserData(2)
%         newValue = UserData(2);
%     elseif newValue < UserData(1)
%         newValue = UserData(1);
%     end
%     set(gco,'Value',newValue);
%     MIA_update_axes(hObject, eventdata, handles)
% end

%
% function MIA_slider_release_click(hObject, eventdata, handles)
%
% set(gco,'UserData',[]);
set(handles.MIA_GUI,'WindowButtonMotionFcn',   @(hObject,eventdata)MIA2('MIA_GUI_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject)));
% set(handles.MIA_GUI,'WindowButtonUpFcn',   @(hObject,eventdata)MIA2('MIA_GUI_WindowButtonUpFcn',hObject,eventdata,guidata(hObject)));




% --------------------------------------------------------------------
function MIA_stats_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_stats_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if exist(handles.data_selected.roifile{:},'file') == 0
if ~isfield(handles,'data_displayed') || ~isfield(handles.data_displayed, 'Cluster')
    [filename,pathname] = uigetfile('.mat',...
        'Select a cluster to print its statistics',...
        strcat(pwd,'/ClustersGMM'));
    if filename == 0
        return
    end
    load(strcat(pathname,filename),'Informations','Statistiques');
else
    if exist(strrep(handles.data_loaded.Cluster.V.fname, '.nii', '.mat'), 'file')
        load(strrep(handles.data_loaded.Cluster.V.fname, '.nii', '.mat'),'Informations', 'Statistiques');
    else
        msg = ['Cannot open ', strrep(handles.data_loaded.Cluster.V.fname, '.nii', '.mat'), ' : No such file.'];
        warndlg(msg)
        return
    end
end
Signatures(Informations, Statistiques, handles.colors_rgb)


% --- Executes on slider movement.
function MIA_slider_slice_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if ~isfield(handles, 'data_loaded') && ~isfield(handles, 'data_selected_for_PRM')
    return
end
% get(handles.IA_slider_slice)
set(handles.MIA_slider_slice, 'Value', round(get(handles.MIA_slider_slice, 'Value')));
MIA_update_axes(hObject, eventdata, handles)


% --- Executes on selection change in MIA_orientation_space_popupmenu.
function MIA_orientation_space_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_orientation_space_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIA_orientation_space_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIA_orientation_space_popupmenu
scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');
% for i=1:handles.data_loaded.number_of_scan
%      handles.data_loaded.Scan(i).nii = read_volume(handles.data_loaded.Scan(i).V, handles.data_loaded.Scan(scan_of_reference).V,3, handles.view_mode);
% end
if ~isfield(handles, 'data_loaded')
    return
end
if isfield(handles.data_loaded, 'ROI')
    for i=1:numel(handles.data_loaded.ROI)
        handles.data_loaded.ROI(i).nii = read_volume(handles.data_loaded.ROI(i).V, handles.data_loaded.Scan(scan_of_reference).V,'auto', handles.view_mode);
    end
end
% % update slider_slice
% set(handles.MIA_slider_slice, 'max', size(handles.data_loaded.Scan(1).nii, 3), 'Value', 1, 'SliderStep',[1/(size(handles.data_loaded.Scan(1).nii, 3) -1) min(5/(size(handles.data_loaded.Scan(1).nii, 3) -1),1)]);
%
% %update handes
handles.display_option.manual_contrast = 0;
guidata(hObject, handles)
MIA_update_axes(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MIA_orientation_space_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MIA_orientation_space_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MIA_pipeline_Manager.
function MIA_pipeline_Manager_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_Manager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MIA_pipeline_Manager
if ~isfield(handles, 'database')
    return
end
MIA_pipeline(hObject, eventdata, handles)


function nii_json_fullfilename = fullfilename(handles, nii_index, ext)

nii_json_fullfilename = [char(handles.database.Path(nii_index)) char(handles.database.Filename(nii_index)) ext];



function MIA_remove_scan(hObject, eventdata, handles, nii_index)
% delete the file (nii/json) from the hard drive
for i=1:numel(nii_index)
    % delete the nii file
    if exist(fullfilename(handles, nii_index(i), '.nii'), 'file')
        delete(fullfilename(handles, nii_index(i), '.nii'))
    end
    if exist(fullfilename(handles, nii_index(i), '.nii.gz'), 'file')
        delete(fullfilename(handles, nii_index(i), '.nii.gz'))
    end
    % delete the json file if it exist
    if exist(fullfilename(handles, nii_index(i), '.json'), 'file') == 2
        delete(fullfilename(handles, nii_index(i), '.json'))
    end
end


% remove the scan from the database
handles.database(nii_index,:) = [];
switch get(hObject, 'Tag')
    case {'MIA_remove_scan', 'MIA_menu_delete_from_database'}
        set(handles.MIA_scans_list, 'Value', 1);
    case 'MIA_remove_time_point'
        set(handles.MIA_scans_list, 'Value', 1);
        set(handles.MIA_time_points_list, 'Value', 1);
    case 'MIA_remove_name'
        set(handles.MIA_scans_list, 'Value', 1);
        set(handles.MIA_time_points_list, 'Value', 1);
        set(handles.MIA_name_list, 'Value', 1);
end
guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);


% --- Executes on button press in MIA_Coronal_view_button.
function MIA_Coronal_view_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_Coronal_view_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change the colour of the button
set(handles.MIA_Coronal_view_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_Saggital_view_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
set(handles.MIA_Axial_view_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.view_mode = 'Coronal';
handles.display_option.manual_contrast = 0;

guidata(hObject, handles)
MIA_update_axes(hObject, eventdata, handles)

% --- Executes on button press in MIA_Saggital_view_button.
function MIA_Saggital_view_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_Saggital_view_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change the colour of the button
set(handles.MIA_Saggital_view_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_Axial_view_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
set(handles.MIA_Coronal_view_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.view_mode = 'Saggital';
handles.display_option.manual_contrast = 0;

guidata(hObject, handles)
MIA_update_axes(hObject, eventdata, handles)

% --- Executes on button press in MIA_Axial_view_button.
function MIA_Axial_view_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_Axial_view_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change the colour of the button
set(handles.MIA_Axial_view_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_Saggital_view_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
set(handles.MIA_Coronal_view_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.view_mode = 'Axial';
handles.display_option.manual_contrast = 0;

guidata(hObject, handles)
MIA_update_axes(hObject, eventdata, handles)


% --------------------------------------------------------------------
function MIA_plot1_Texture_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_plot1_Texture_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end
%% code to generate texture values for 1 scan and X ROIs
ROI_names = char(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'ROI'));
Scan_names = char(handles.data_loaded.info_data_loaded.SequenceName(handles.data_loaded.info_data_loaded.Type == 'Scan'));
scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');
texture_values = table;
warning('off')
for i = 1:handles.data_loaded.number_of_scan
    for j = 1:handles.data_loaded.number_of_ROI
        volume = squeeze(handles.data_displayed.image(:,:,:,i));
        mask = handles.data_loaded.ROI(j).nii;
        matrix_tmp = volume .* mask;
        matrix_tmp(matrix_tmp==0)=nan;
        
        texture_values.SequenceName(size(texture_values,1)+1) = nominal(Scan_names(i,:));
        texture_values.ROI(size(texture_values,1)) = nominal(ROI_names(j,:));
        texture_values.ROI_Size_mm(size(texture_values,1)) = prod(handles.data_loaded.Scan(scan_of_reference).json.Grid_spacings__X_Y_Z_T_____.value(1:3)) * sum(mask(:));
        %  texture_values.Entropy(size(texture_values,1)) = entropy(matrix_tmp);
        texture_values.Mean(size(texture_values,1)) = nanmean(matrix_tmp(:));
        texture_values.Median(size(texture_values,1)) = nanmedian(matrix_tmp(:));
        texture_values.Percentile_1(size(texture_values,1)) = prctile_copy(matrix_tmp(:),1);
        texture_values.Percentile_5(size(texture_values,1)) = prctile_copy(matrix_tmp(:),5);
        texture_values.Percentile_95(size(texture_values,1)) = prctile_copy(matrix_tmp(:),95);
        texture_values.Percentile_99(size(texture_values,1)) = prctile_copy(matrix_tmp(:),99);
        
        [ROIonly,~,~,~] = prepareVolume(volume,mask,'MRscan',...
            sum(handles.data_loaded.Scan(scan_of_reference).json.Grid_spacings__X_Y_Z_T_____.value(1:2))/2 , ...
            handles.data_loaded.Scan(scan_of_reference).json.Grid_spacings__X_Y_Z_T_____.value(3), 1,...
            'pixelW','Global','Equal',128);
        [texture] =  getGlobalTextures(ROIonly,100);
        texture_values.Kurtosis(size(texture_values,1)) = texture.Kurtosis;
        texture_values.Skewness(size(texture_values,1)) = texture.Skewness;
        texture_values.Variance(size(texture_values,1)) = texture.Variance;
        clear matrix_tmp;
        
    end
end
disp(texture_values);
warning('on')


% --------------------------------------------------------------------
function MIA_Import_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_Import_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ROI_listing = [unique(handles.database.SequenceName(handles.database.Type == 'ROI'))', 'Other']';
% Enter the ROI's name
[VOI_number,ok] = listdlg('Name', 'Bip', 'SelectionMode', 'single', 'ListString',  ROI_listing,'ListSize', [200 150],...
    'PromptString', 'Select the Roi''s name');
if ok == 0
    return
end

if strcmp('Other', char(ROI_listing(VOI_number))) == 1
    newVOI_name = inputdlg('Name of the new VOI ', 'Question?', 1, {''});
    if isempty(newVOI_name)
        return
    end
    if sum(strcmp(newVOI_name, cellstr(ROI_listing)) > 0)
        warning_text = sprintf('Can not import this "New" ROI, because %s ROI exist already,\nplease try again and select %s ROI in the ROIs list proposed',...
            newVOI_name{:}, newVOI_name{:});
        msgbox(warning_text, 'ROI warning') ;
        return
    end
    
else
    newVOI_name = {char(ROI_listing(VOI_number))};
end


[file,path] = uigetfile('*.nii', 'Select your ROI:');

if path == 0
    return
end

Index = get_data_selected(handles);
Button = get(handles.MIA_scan_VOIs_button, 'Value');
if isempty(Index)
    set(handles.MIA_scan_VOIs_button, 'Value',0);
    Index = get_data_selected(handles);
    set(handles.MIA_scan_VOIs_button, 'Value',Button);
end
%set(handles.MIA_scan_VOIs_button, 'Value',1);
Entry_Selected = handles.database(Index,:);

if size(Entry_Selected,1)>1
    Entry_Selected = Entry_Selected(1,:);
end

file_name = strcat(char(Entry_Selected.Patient),...
    '-', char(Entry_Selected.Tp),...
    '-ROI',...
    '-', newVOI_name{1} ,...
    '_',datestr(now,'yyyymmdd-HHMMSSFFF'));



ROI_dest = handles.database.Properties.UserData.MIA_ROI_path;
if ~exist(ROI_dest)
    mkdir(ROI_dest)
end
%status = copyfile([path, file],ROI_dest);

if endsWith(file, '.nii.gz')
    gunzip([path, file]);
    file = strrep(file, '.nii.gz', '.nii');
    status = copyfile([path, file], [handles.database.Properties.UserData.MIA_ROI_path, file_name, '.nii']);
    delete([path, file]);
else
    status = copyfile([path, file], [handles.database.Properties.UserData.MIA_ROI_path, file_name, '.nii']);
end

if ~status
    msgbox('An error occured when trying to copy your ROI file into the MIA ROI_data folder of your project.', 'Copy error') ;
    return
end

New_Entry = Entry_Selected;
New_Entry.Type = categorical(cellstr('ROI'));
New_Entry.IsRaw = categorical(0);
New_Entry.SequenceName = categorical(cellstr(newVOI_name{1}));
New_Entry.Path = categorical(cellstr(handles.database.Properties.UserData.MIA_ROI_path));
New_Entry.Filename = categorical(cellstr(file_name));
handles.database = [handles.database; New_Entry];
set(handles.MIA_scans_list, 'Value',1);
guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles)



% --------------------------------------------------------------------
function MIA_menu_Display_Transform_Matrix_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Display_Transform_Matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'data_loaded')
    return
end
table = handles.data_loaded.info_data_loaded;
table2 = table(table.Type == categorical(cellstr('Scan')),:);
Mats = cell(size(table2, 1),1);
for i=1:size(table2,1)
    file = [char(table2.Path(i)), char(table2.Filename(i)), '.nii'];
    info = niftiinfo(file);
    Mats{i} = info.Transform.T;
    f = figure('Name', char(table2.Filename(i)));
    %title(uf,[char(table2.Filename(i)), '.nii'])
    t = uitable(f,'Data',Mats{i});%,'Position',[20 20 260 204]);
end


% --------------------------------------------------------------------
function MIA_menu_Import_BIDS_data_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Import_BIDS_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'database')
    answer = questdlg('Please be sure that the merged database will not contain some files with the same Patient - Tp - SequenceName tags. This case is currently not implemented. Continue ?');
    if any(strcmp(answer, {'No', 'Cancel'}))
        return
    end
    MIA_data_path = handles.database.Properties.UserData.MIA_data_path;
    %     MIA_data_path = handles.database.Properties.UserData.MIA_data_path ;
    %handles.database = table();
else
    MIA_root_path = uigetdir(pwd, 'Select the directory to save your new projet');
    if sum(MIA_root_path) == 0
        return
    end
    
    MIA_root_path = [MIA_root_path filesep];
    %% create the output folder ('MIA_data')
    MIA_data_path = MIA_root_path;%[MIA_root_path, 'MIA_data', filesep];
    % Create the output folder if needed
    if exist(MIA_data_path, 'dir') ~= 7
        status = mkdir(MIA_data_path);
        if status == 0
            warndlg('You do not have the right to write in the folder!', 'Warning');
            return
        end
    end
    % Create the RAW data folder if needed
    if exist([MIA_data_path, 'Raw_data', filesep], 'dir') ~= 7
        status = mkdir([MIA_data_path, 'Raw_data', filesep]);
        if status == 0
            warndlg('You do not have the right to write in the folder!', 'Warning');
            return
        end
    end
    
    
    %% create the database structure
    handles.database =  table;%
    handles.database = cell2table(cell(0,8));
    handles.database.Properties.VariableNames = {'Group','Patient', 'Tp', 'Path', 'Filename', 'Type', 'IsRaw', 'SequenceName'};
    handles.database.Properties.UserData.MIA_data_path = MIA_data_path;
    handles.database.Properties.UserData.MIA_Raw_data_path = [MIA_data_path, 'Raw_data', filesep];
    handles.database.Properties.UserData.MIA_Derived_data_path = [MIA_data_path, 'Derived_data', filesep];
    handles.database.Properties.UserData.MIA_ROI_path = [MIA_data_path, 'ROI_data', filesep];
    handles.database.Properties.UserData.Order_data_display = {'ascend','ascend','ascend'};
    
end

BIDS_Path = uigetdir(pwd, 'Select the directory containing the BIDS data you want to import');
if sum(BIDS_Path) == 0
    return
end
if ~strcmp(BIDS_Path(end), filesep)
    BIDS_Path = [BIDS_Path, filesep];
end
if ~strcmp(MIA_data_path(end), filesep)
    MIA_data_path = [MIA_data_path, filesep];
end
database = CreateTableFromBIDS(BIDS_Path, MIA_data_path,0);
if ~isempty(database)
    handles.database = [handles.database; database];
end
guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles)


% --------------------------------------------------------------------
function MIA_menu_Compress_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Compress_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data_path = handles.database.Properties.UserData.MIA_data_path;
DirectoriesToProcess = {'Derived_data', 'Raw_data', 'ROI_data'};


if exist([Data_path, 'Uncompressed_Trash'],'dir') ~= 7
    [status, ~, ~] = mkdir([Data_path, 'Uncompressed_Trash']);
    if status == false
        error('Cannot create the Uncompressed_Thrash folder to move the uncompressed files in after compression.')
    end
end

Files_To_Compress = {};
NbCompressions = 0;
for i=1:length(DirectoriesToProcess)
    Files = dir([Data_path, DirectoriesToProcess{i}, '/*.nii']);
    for j=1:length(Files)
        filename = [Files(j).folder, filesep, Files(j).name];
        Files_To_Compress = [Files_To_Compress, {filename}];
        gzip(filename);
        NbCompressions = NbCompressions+1;
        movefile(filename, strrep(filename, ['/', DirectoriesToProcess{i}, '/'], ['/', 'Uncompressed_Trash', '/']));
    end
end
Mess = ['Done! ', num2str(NbCompressions), ' files were compressed and the useless uncompressed files were moved to the ''Uncompressed_Thrash'' folder.'];
msgbox(Mess)


% --- Executes on button press in MIA2_exit.
function MIA2_exit_Callback(hObject, eventdata, handles)
% hObject    handle to MIA2_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 selection = questdlg('Do you want to close (did you save you database)?',...
    'Warning',...
    'Yes','No','Yes');
if isempty(selection)
    return
end
switch selection
    case 'No'
        return
end

if ~isempty(findobj('type', 'figure', 'name', 'MIA pipeline Manager'))
    close((findobj('type', 'figure', 'name', 'MIA pipeline Manager')));
end

delete(handles.MIA_GUI);


% --------------------------------------------------------------------
function MIA_tools_file_history_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_tools_file_history (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isfield(handles, 'data_loaded')
    warndlg('Please load the scan you want to display the history.');
    return
end

if length(handles.data_loaded.Scan) ~= 1
    warndlg('Please load only one scan to display its history.');
    return
end

% JSON = handles.data_loaded.Scan.json;
% 
% if ~isfield(JSON, 'Bricks')
%     warndlg('There is no history available for this file. It might be raw data of the processing has not be done by the Pipeline Manager.');
%     return
% end
% 
% 
% Jobs = JSON.Bricks;


FileHistory(hObject, eventdata, handles)




