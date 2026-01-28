function varargout = Add_Tag_to_new_patient_GUI(varargin)
% ADD_TAG_TO_NEW_PATIENT_GUI MATLAB code for Add_Tag_to_new_patient_GUI.fig
%      ADD_TAG_TO_NEW_PATIENT_GUI, by itself, creates a new ADD_TAG_TO_NEW_PATIENT_GUI or raises the existing
%      singleton*.
%
%      H = ADD_TAG_TO_NEW_PATIENT_GUI returns the handle to a new ADD_TAG_TO_NEW_PATIENT_GUI or the handle to
%      the existing singleton*.
%
%      ADD_TAG_TO_NEW_PATIENT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADD_TAG_TO_NEW_PATIENT_GUI.M with the given input arguments.
%
%      ADD_TAG_TO_NEW_PATIENT_GUI('Property','Value',...) creates a new ADD_TAG_TO_NEW_PATIENT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Add_Tag_to_new_patient_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Add_Tag_to_new_patient_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Add_Tag_to_new_patient_GUI

% Last Modified by GUIDE v2.5 23-Oct-2017 18:25:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Add_Tag_to_new_patient_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Add_Tag_to_new_patient_GUI_OutputFcn, ...
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


% --- Executes just before Add_Tag_to_new_patient_GUI is made visible.
function Add_Tag_to_new_patient_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Add_Tag_to_new_patient_GUI (see VARARGIN)

% Choose default command line output for Add_Tag_to_new_patient_GUI
handles.output = hObject;
handles.output_filename_data = varargin{1};
info_patient = varargin{2};
set(handles.IA_load_brucker_GUI_patient_name, 'String', char(info_patient{3}));
% set popup menu
id_list = {''}; tp_list = {''};
if ~isempty(handles.output_filename_data.database)
    handles.id_list = [id_list cellstr(sort(unique(handles.output_filename_data.database.Patient)))'];
    handles.tp_list = [tp_list cellstr(sort(unique(handles.output_filename_data.database.Tp)))'];
else
    handles.id_list = {''};
    handles.tp_list = {''};
end
set(handles.IA_load_brucker_GUI_id_popupmenu, 'String', handles.id_list);
set(handles.IA_load_brucker_GUI_time_point_popupmenu, 'String', handles.tp_list);
if sum(strcmp(handles.id_list, char(info_patient{1}))) ==0
    id_name = clean_variable_name(char(info_patient{1}), '');
    set(handles.IA_load_brucker_GUI_id_other, 'String', id_name);
else
    set(handles.IA_load_brucker_GUI_id_popupmenu, 'Value', find(strcmp(handles.id_list, info_patient{1}) == 1))
end
if sum(strcmp(handles.tp_list, char(info_patient{2}))) == 0
    tp_name = clean_variable_name(char(info_patient{2}), '');
    set(handles.IA_load_brucker_GUI_time_point_other, 'String', tp_name);
else
    set(handles.IA_load_brucker_GUI_time_point_popupmenu, 'Value', find(strcmp(handles.tp_list, info_patient{2}) == 1))
end



% defind output
handles.name_selected = '';
handles.tp_selected = '';
handles.output1 = [{handles.name_selected}, {handles.tp_selected }];
handles.output2 = 0; 


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Add_Tag_to_new_patient_GUI wait for user response (see UIRESUME)
% uiwait(handles.IA_load_brucker_info_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = Add_Tag_to_new_patient_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
uiwait;
handles = guidata(hObject);

% varargout{1} = handles.output1;
% varargout{2} = handles.output2;
varargout{1} = handles.name_selected;
varargout{2} = handles.tp_selected;
delete(handles.IA_load_brucker_info_GUI)


% --- Executes on selection change in IA_load_brucker_GUI_id_popupmenu.
function IA_load_brucker_GUI_id_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_id_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IA_load_brucker_GUI_id_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IA_load_brucker_GUI_id_popupmenu


% --- Executes during object creation, after setting all properties.
function IA_load_brucker_GUI_id_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_id_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in IA_load_brucker_GUI_time_point_popupmenu.
function IA_load_brucker_GUI_time_point_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_time_point_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IA_load_brucker_GUI_time_point_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IA_load_brucker_GUI_time_point_popupmenu


% --- Executes during object creation, after setting all properties.
function IA_load_brucker_GUI_time_point_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_time_point_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IA_load_brucker_GUI_id_other_Callback(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_id_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IA_load_brucker_GUI_id_other as text
%        str2double(get(hObject,'String')) returns contents of IA_load_brucker_GUI_id_other as a double


% --- Executes during object creation, after setting all properties.
function IA_load_brucker_GUI_id_other_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_id_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IA_load_brucker_GUI_time_point_other_Callback(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_time_point_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IA_load_brucker_GUI_time_point_other as text
%        str2double(get(hObject,'String')) returns contents of IA_load_brucker_GUI_time_point_other as a double


% --- Executes during object creation, after setting all properties.
function IA_load_brucker_GUI_time_point_other_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_time_point_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IA_load_brucker_GUI_OK_button.
function IA_load_brucker_GUI_OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% id's name
if  get(handles.IA_load_brucker_GUI_id_popupmenu, 'Value')> 1 %strcmp('Other?', get(handles.IA_load_brucker_GUI_id_other, 'String'))
    id_nbr = get(handles.IA_load_brucker_GUI_id_popupmenu, 'Value');
    handles.name_selected = handles.id_list{id_nbr};
else
    handles.name_selected = get(handles.IA_load_brucker_GUI_id_other, 'String');
end

% time point's name
if get(handles.IA_load_brucker_GUI_time_point_popupmenu, 'Value') > 1% strcmp('Other?', get(handles.IA_load_brucker_GUI_time_point_other, 'String'))
    id_nbr = get(handles.IA_load_brucker_GUI_time_point_popupmenu, 'Value');
    handles.tp_selected = handles.tp_list{id_nbr};
else
    handles.tp_selected = get(handles.IA_load_brucker_GUI_time_point_other, 'String');
end

handles.output1 = [handles.name_selected, handles.tp_selected]; 
handles.output2 = 1;

% Update handles structure
guidata(hObject, handles);

uiresume;



% --- Executes on button press in IA_load_brucker_GUI_Cancel_button.
function IA_load_brucker_GUI_Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to IA_load_brucker_GUI_Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output2 = 0;
handles.output1 = 0;
% handles.output.statut = 0;
% handles.output.statut = [{handles.name_selected}, {handles.tp_selected}];

% Update handles structure
guidata(hObject, handles);

uiresume;
