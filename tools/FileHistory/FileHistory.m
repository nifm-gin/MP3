function varargout = FileHistory(varargin)
% FILEHISTORY MATLAB code for FileHistory.fig
%      FILEHISTORY, by itself, creates a new FILEHISTORY or raises the existing
%      singleton*.
%
%      H = FILEHISTORY returns the handle to a new FILEHISTORY or the handle to
%      the existing singleton*.
%
%      FILEHISTORY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILEHISTORY.M with the given input arguments.
%
%      FILEHISTORY('Property','Value',...) creates a new FILEHISTORY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FileHistory_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FileHistory_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FileHistory

% Last Modified by GUIDE v2.5 18-Oct-2018 18:21:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FileHistory_OpeningFcn, ...
                   'gui_OutputFcn',  @FileHistory_OutputFcn, ...
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


% --- Executes just before FileHistory is made visible.
function FileHistory_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FileHistory (see VARARGIN)

% Choose default command line output for FileHistory
handles.output = hObject;
handles.MP3_data = varargin{3};

[hObject, handles] = UpdateJobsList(hObject, handles);




% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FileHistory wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FileHistory_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in FileHistory_JobsListbox.
function FileHistory_JobsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to FileHistory_JobsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[hObject, handles] = UpdateJobsList(hObject, handles);
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns FileHistory_JobsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileHistory_JobsListbox


% --- Executes during object creation, after setting all properties.
function FileHistory_JobsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileHistory_JobsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FileHistory_JobsFieldsListbox.
function FileHistory_JobsFieldsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to FileHistory_JobsFieldsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[hObject, handles] = UpdateJobsFieldsList(hObject, handles);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns FileHistory_JobsFieldsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileHistory_JobsFieldsListbox


% --- Executes during object creation, after setting all properties.
function FileHistory_JobsFieldsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileHistory_JobsFieldsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FileHistory_JobsValuesListbox.
function FileHistory_JobsValuesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to FileHistory_JobsValuesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FileHistory_JobsValuesListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileHistory_JobsValuesListbox


% --- Executes during object creation, after setting all properties.
function FileHistory_JobsValuesListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileHistory_JobsValuesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function [hObject, handles] = UpdateJobsList(hObject, handles)
    
    handles.JSON = handles.MP3_data.data_loaded.Scan.json;
    if ~isfield(handles.JSON, 'Bricks')
        set(handles.Filename, 'String', {'No History available for this file. Might be raw data or the processing has not be made in the Pipeline Manager.'});
        set(handles.FileHistory_JobsListbox, 'String', {''});
        set(handles.FileHistory_JobsFieldsListbox, 'String', {''});
        set(handles.FileHistory_JobsValuesListbox, 'String', {''});
        set(handles.FileHistory_JobsListbox, 'Value', 1);
        set(handles.FileHistory_JobsFieldsListbox, 'Value', 1);
        set(handles.FileHistory_JobsValuesListbox, 'Value', 1);
        return
    end
    Filename = handles.MP3_data.data_loaded.Scan.V.fname;
    Parts = strsplit(Filename, filesep);
    ShortFilename = strjoin(Parts(end-1:end), filesep);
    set(handles.Filename, 'String', {ShortFilename});
    Jobs = fieldnames(handles.JSON.Bricks);
    set(handles.FileHistory_JobsListbox, 'String', Jobs);
    set(handles.FileHistory_JobsFieldsListbox, 'Value', 1);
    [hObject, handles] = UpdateJobsFieldsList(hObject, handles);
    set(handles.FileHistory_JobsValuesListbox, 'Value', 1);
    guidata(hObject, handles);
    
    
   
    

function [hObject, handles] = UpdateJobsFieldsList(hObject, handles)
    SelectedField = handles.FileHistory_JobsListbox.String{handles.FileHistory_JobsListbox.Value};
    if isempty(SelectedField)
        return
    end
    Fields = fieldnames(handles.JSON.Bricks.(SelectedField));
    FieldsToDisplay = {};
    for i=1:length(Fields)
       if strcmp(Fields{i}, 'files_in') || strcmp(Fields{i}, 'files_out') || strcmp(Fields{i}, 'opt')
           Fields2 = fieldnames(handles.JSON.Bricks.(SelectedField).(Fields{i}));
           for j=1:length(Fields2)
              FieldsToDisplay = [FieldsToDisplay, {[Fields{i}, ' ', Fields2{j}]}];
           end
       else
           FieldsToDisplay = [FieldsToDisplay, {Fields{i}}];
       end
    end

    set(handles.FileHistory_JobsFieldsListbox, 'String', FieldsToDisplay);
    set(handles.FileHistory_JobsValuesListbox, 'Value', 1);
    [hObject, handles] = UpdateJobsValuesList(hObject, handles);
    guidata(hObject, handles);
    
    
function [hObject, handles] = UpdateJobsValuesList(hObject, handles)
    SelectedField = handles.FileHistory_JobsListbox.String{handles.FileHistory_JobsListbox.Value};
    SelectedParameterField = handles.FileHistory_JobsFieldsListbox.String{handles.FileHistory_JobsFieldsListbox.Value};
    Fields = strsplit(SelectedParameterField);
    Param = handles.JSON.Bricks.(SelectedField).(Fields{1});
    if length(Fields) == 1
        Entrie = Param;
    else
        Entrie = Param.(Fields{2});
    end
    if strcmp(Fields{1}, 'files_in') || strcmp(Fields{1}, 'files_out')
        NewEntrie = [];
        for i=1:length(Entrie)
            EntrieSpl = split(Entrie{i}, filesep);
            name = strjoin(EntrieSpl(end-1:end), filesep);
            %[~,name,~] = fileparts(Entrie{i});
            NewEntrie = [NewEntrie; {name}];
        end
        Entrie = NewEntrie;
    end
    
    if islogical(Entrie)
        Entrie = num2str(Entrie);
    elseif isstruct(Entrie)
        Entrie = {'This type of variable cannot be displayed here. Please open the corresponding json file to explore it.'};
    end
    
    set(handles.FileHistory_JobsValuesListbox, 'String', Entrie);
    set(handles.FileHistory_JobsValuesListbox, 'Value', 1);
    guidata(hObject, handles);

    
    
    
