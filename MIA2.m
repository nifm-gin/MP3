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


% Last Modified by GUIDE v2.5 06-Feb-2018 09:50:47


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

% Load data/ROI from uvasc
if ~isempty(varargin)
    handles.uvascdata= varargin{1};
    handles.data = varargin{2};
    handles.roi = varargin{3};
end
% init stuct/variables
handles.VOIs = {'Other'};
handles.histo = {'Other'}; %{'Pimo', 'ColIV', 'Tc+I_cerveau', 'Tc+I_ref'};
handles.resolution = [1 64 112 128 192 256 384 512 3000];
handles.colors ={'b', 'g', 'y', 'm', 'c', 'r', 'k', 'w', 'navy',...
    'u1','turquoise','slateblue',	'springgreen',	'maroon',...
    'purple',	'u2',	'olive',	'u3','chartreuse',	'u4',	'sky',...
    'u5',	'orange',	'u6',	'u7',	'u8',	'gray'};
handles.colors_rgb = [0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 0 0; 1 1 1];
handles.colormap = get(handles.MIA_colormap_popupmenu,'String');
handles.markers ={'o','s', 'd', 'p', 'h', '+', '*', 'x'};
table_data(1,1) = {'Voxel values'};
set(handles.MIA_table_pixel_values, 'Data', table_data);
set(handles.MIA_table1, 'Data', {'', '', '', '', ''});
handles.table1.cluster = [];
handles.table1.cluster_row = [];
handles.mode = 1;
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
[dpath, ~, ~] = fileparts(which('sendList.m'));
javaclasspath(strcat(dpath, filesep, 'MRIManager.jar'))
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
switch selection,
    case 'No'
        return
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

patient_id = get(handles.MIA_name_list, 'Value');
time_point = get(handles.MIA_time_points_list, 'Value');
scan = get(handles.MIA_scans_list, 'Value');

id_listing = unique(handles.database.Patient);
set(handles.MIA_name_list, 'String', char(id_listing));
if numel(patient_id)~= 1
    return
end

Patient_filter = handles.database.Patient== id_listing(patient_id);
tp_listing = unique(handles.database.Tp(Patient_filter));
set(handles.MIA_time_points_list, 'String', char(tp_listing));

if get(handles.MIA_scan_VOIs_button, 'Value') == 0 %display parameters list
    is_scan =  handles.database.Type == 'Scan';
    tp_filter = handles.database.Tp== tp_listing(time_point);
    sequence_listing = handles.database.SequenceName(Patient_filter & tp_filter & is_scan);
    set(handles.MIA_scans_list, 'String', char(sequence_listing));
    file_text= cell(1, numel(sequence_listing(scan)));
    for i=1:numel(sequence_listing(scan))
        sequence_filter =  handles.database.SequenceName== sequence_listing(scan(i));
        file_text(i) = cellstr(handles.database.Filename(Patient_filter & tp_filter & sequence_filter & is_scan));
        
    end
    set(handles.MIA_file_list, 'String', file_text);
    
else %display VOIs list
    is_ROI=  handles.database.Type == 'ROI';
    tp_filter = handles.database.Tp== tp_listing(time_point);
    sequence_listing = handles.database.SequenceName(Patient_filter & tp_filter & is_ROI);
    if isempty(sequence_listing)
        set(handles.MIA_scans_list, 'String', '');
        return
    end
    set(handles.MIA_scans_list, 'String', char(sequence_listing));
    file_text= cell(1, numel(sequence_listing(scan)));
    for i=1:numel(sequence_listing(scan))
        sequence_filter =  handles.database.SequenceName== sequence_listing(scan(i));
        file_text(i) = cellstr(handles.database.Filename(Patient_filter & tp_filter & sequence_filter & is_ROI));
    end
    set(handles.MIA_file_list, 'String', file_text);
    
    
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


% --- Executes on button press in MIA_path_file.
function MIA_path_file_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_path_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end
data_selected = get_data_selected(handles);
% contents = cellstr(get(handles.MIA_file_list,'String'));
% FileNameLoad = contents{get(hObject,'Value')};
fprintf('\nload(''%s'')\n',handles.database.nii(data_selected));



% --------------------------------------------------------------------
function MIA_add_name_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_add_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global uvascim uvascroi

if strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker') == 1 ||...
        strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker_OM') == 1 ||...
        strcmp(get(hObject, 'Tag'), 'MIA_menu_load_data') == 1 ||...
        strcmp(get(hObject, 'Tag'), 'MIA_menu_load_DICOM') == 1
    filename = handles.load_bruker_tmp.filename;
    directory = handles.load_bruker_tmp.pathname;
    
else
    [filename, directory]=uigetfile({'*.mat;*.nii;','Image Files (*.mat,*.nii...)'}, 'Selection une structure uvasc','MultiSelect','on');
    if isequal(filename,0)
        return
    end
    if ischar(filename)
        filename = {filename};
    end
end

% Can load several structure
for x = 1:numel(filename)
    handles = guidata(hObject);
    fid=fopen(fullfile(directory,filename{x}),'r');
    if fid>0
        fclose(fid);
        [~, ~, ext] = fileparts(filename{x});
        if strcmp(ext, '.nii')
            data_loaded.uvascim= convert_nii2uvascim(fullfile(directory,filename{x}));
        else
            data_loaded = load(fullfile(directory,filename{x}));
        end
        
        if ~or(isfield(data_loaded, 'uvascim'), isfield(data_loaded, 'so2struct'))
            warndlg('Please select a correct uvasc structure');
            clear data_loaded
            return
        end
    else
        warndlg('Please select a correct uvasc structure');
        return
    end
    %     end
    id_nbr = numel(get(handles.MIA_name_list, 'String'));
    % First patient of the database
    if id_nbr == 0
        id_nbr =1;
    end
    if isfield(data_loaded, 'uvascim') && strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker') ~= 1 &&...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker_OM') ~= 1 &&...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_data') ~= 1 &&...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_DICOM') ~= 1
        
        handles.data = data_loaded.uvascim;
        id_info =filename{x};
        sep_position = strfind(id_info, '-');
        %try to figure out the patient's name and time point
        if ~isempty(sep_position) && numel(sep_position) == 2
            id = id_info(1:sep_position(1)-1);
            tp = id_info(sep_position(1)+1:sep_position(2)-1);
        else
            id = id_info;
            tp = id_info;
        end
    elseif strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker') ~= 1 || ...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker_OM') ~= 1 ||...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_data') ~= 1 ||...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_DICOM') ~= 1
        
        handles.data = data_loaded.uvascim;
    elseif isfield(data_loaded, 'so2struct')
        
        % load computed maps
        maps_names = fieldnames(data_loaded.so2struct.map);
        for i = 1:numel(maps_names)
            handles.data.image(i) = eval(['data_loaded.so2struct.map.' maps_names{i}]);
        end
        % correct the bug on the z_offset in the so2stuct
        % the correct z_offset is located in the T2map
        match = strcmp(maps_names, 'T2map');
        scan_ref = find(match==1);
        for i = 1:numel(maps_names)
            if match(i) == 0
                handles.data.image(i).reco = rmfield(handles.data.image(i).reco,'fov_offsets');
                handles.data.image(i).reco.fov_offsets =  handles.data.image(scan_ref).reco.fov_offsets;
                handles.data.image(i).reco.thickness = handles.data.image(scan_ref).reco.thickness;
            end
        end
        for i = 1:numel(maps_names)
            if match(i) == 0
                handles.data.image(i).reco = rmfield(handles.data.image(i).reco,'fov_offsets');
                handles.data.image(i).reco.fov_offsets =  handles.data.image(scan_ref).reco.fov_offsets;
                handles.data.image(i).reco.thickness = handles.data.image(scan_ref).reco.thickness;
            end
        end
        
        id_info =filename{x};
        id_beg = strfind(handles.data.image(1).texte,'Patient name');
        id_end = strfind(handles.data.image(1).texte,'Examination name');
        id = strtrim(handles.data.image(1).texte(id_beg+39:id_end-6));
        tp_beg = strfind(handles.data.image(1).texte,'Examination date/time              :');
        td_end = strfind(handles.data.image(1).texte,'.    Series Type');
        tp = strtrim(handles.data.image(1).texte(tp_beg+36:td_end-14));
        
    end
    clear data_loaded
    handles.data.roifile = '';
    % The soft already know that!
    if strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker') ~=1 &&...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker_OM') ~=1  &&...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_data') ~=1  &&...
            strcmp(get(hObject, 'Tag'), 'MIA_menu_load_DICOM') ~= 1
        dlgbox_title =  strcat('File name:', {' '}, id_info,{' '}');
        answer = inputdlg({'Is the id''s name correct?', 'Is the time point correct?'},dlgbox_title{:}, [1 size(dlgbox_title{:},2)+25], [{id} {tp}]);
        if isempty(answer)
            return
        end
    else
        answer(1) = handles.load_bruker_tmp.answer(1);
        answer(2) = handles.load_bruker_tmp.answer(2);
    end
    % find if the patient exist already
    if ~isempty(strmatch(answer(1), get(handles.MIA_name_list, 'String'), 'exact')) %#ok<*MATCH3>
        id_nbr = strmatch(answer(1), get(handles.MIA_name_list, 'String'), 'exact');
        tmp = handles.database(1:id_nbr);
        % find if the time point exist already
        if ~isempty(strmatch(answer(2), {handles.database(id_nbr).day.date}', 'exact'))
            tp_nbr = strmatch(answer(2), {handles.database(id_nbr).day.date}', 'exact');
        else
            tp_nbr = numel(handles.database(id_nbr).day)+1;
            tmp(id_nbr).day(tp_nbr).date = '';
            tmp(id_nbr).day(tp_nbr).comment = '';
            tmp(id_nbr).day(tp_nbr).VOIs = '';
            tmp(id_nbr).day(tp_nbr).VOIs_file = '';
            tmp(id_nbr).day(tp_nbr).parameters = '';
            tmp(id_nbr).day(tp_nbr).scans_file = '';
        end
    else
        if id_nbr > 1 || id_nbr == 1 && isfield(handles, 'database')
            tmp = handles.database(1:id_nbr);
            id_nbr = id_nbr+1;
            tmp(id_nbr).day.date = '';
            tmp(id_nbr).day.comment = '';
            tmp(id_nbr).day.VOIs = '';
            tmp(id_nbr).day.VOIs_file = '';
            tmp(id_nbr).day.parameters = '';
            tmp(id_nbr).day.scans_file = '';
            tmp(id_nbr).group = 'group name';
            tmp(id_nbr).omit = 0;
        else % First patient of the database
            tmp(1).name = '';
            tmp(1).path = '';
            tmp(1).group = 'group name';
            tmp(1).omit = 0;
            tmp(1).day.date = '';
            tmp(1).day.comment = '';
            tmp(1).day.VOIs = '';
            tmp(1).day.VOIs_file = '';
            tmp(1).day.parameters = '';
            tmp(1).day.scans_file = '';
        end
        tp_nbr = 1;
    end
    % add patient info to the database
    tmp(id_nbr).name = answer{1};
    if isfield(handles,'database')
        tmp(id_nbr).path = handles.database(1).path;
    else
        if isempty(strfind(directory, 'Image_Analyses_data'))
            if strcmp(directory(end), filesep) ~= 1
                tmp(id_nbr).path = [directory 'Image_Analyses_data' filesep];
            else
                tmp(id_nbr).path = [directory filesep 'Image_Analyses_data' filesep];
            end
            % Create a new folder if needed
            %         if exist(tmp(id_nbr).path, 'dir') ~= 7
            status = mkdir(tmp(id_nbr).path);
            if status == 0
                warndlg('You do not the right to write in the folder!', 'Warning');
                %                 clear uvascim
                return
            end
            %         end
        else
            if strcmp(directory(end), filesep) ~= 1
                tmp(id_nbr).path = [directory filesep];
            else
                tmp(id_nbr).path = directory;
            end
        end
    end
    %     end
    
    tmp(id_nbr).day(tp_nbr).date = answer{2};
    
    % find scan info
    size_data = zeros(numel(handles.data.image),5);
    list_scan = cell(numel(handles.data.image),1);
    for i=1:numel(handles.data.image)
        if strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker') ==1 ||...
                strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker_OM') ==1 ||...
                strcmp(get(hObject, 'Tag'), 'MIA_menu_load_data') ==1 ||...
                strcmp(get(hObject, 'Tag'), 'MIA_menu_load_DICOM') == 1
            list_scan{i} = sprintf(handles.load_bruker_tmp.scanname{x});
        else
            size_data(i,1:size(size(handles.data.image(i).reco.data),2))=size(handles.data.image(i).reco.data);
            if isfield(handles.data.image(i).reco,'texte')
                list_scan{i} = sprintf(handles.data.image(i).reco.texte);
            else
                list_scan{i}= sprintf(handles.data.image(i).acq.ppl_name);
            end
        end
    end
    % Split and save 1 structure per parameter
    for i=1:numel(list_scan)
        % Clean the futur file name by removing weird characters: ( ) [
        % ]....
        if strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker') ==1 ||...
                strcmp(get(hObject, 'Tag'), 'MIA_menu_load_bruker_OM') == 1 ||...
                strcmp(get(hObject, 'Tag'), 'MIA_menu_load_data') == 1 ||...
                strcmp(get(hObject, 'Tag'), 'MIA_menu_load_DICOM') == 1
            file_name =  char(filename(x));
            file_name = cellstr(file_name(1:end-4));
        else
            file_name = strcat(answer{1},'_',answer{2},'_',num2str(i), '-', list_scan(i));
        end
        
        file_name = regexprep(file_name,'[: ]','-');
        file_name = regexprep(file_name,'*','star');
        file_name = regexprep(file_name,'(','');
        file_name = regexprep(file_name,')','');
        file_name = regexprep(file_name,',','');
        file_name = regexprep(file_name,'<','');
        file_name = regexprep(file_name,'>','');
        file_name = regexprep(file_name,'/','-');
        
        
        % Split Echoes if for particular scan such as (DTI-EPI)
        if  ~isempty(regexpi(handles.data.image(i).acq.ppl_name,'\w*Dti\w*'))  && ...
                strcmp(handles.data.image(i).reco.echo_label(1), 'signal intensity')
            %save only the ADC map !!
            for ii=3%   %% 1:size(handles.data.image(i).reco.data,3)
                % If data imported from uvasc --> save uvascim
                %tmp_uvascim = uvascim;
                uvascim = handles.data;
                uvascim.image = handles.data.image(i);
                uvascim.image.reco.data = uvascim.image.reco.data(:,:,ii,:);
                uvascim.image.reco.fov_offsets = uvascim.image.reco.fov_offsets(:,ii,:,:);
                uvascim.image.reco.fov_orientation = uvascim.image.reco.fov_orientation(:,ii,:,:);
                uvascim.image.reco.fov_phase_orientation= uvascim.image.reco.fov_phase_orientation(ii,:,:);
                uvascim.image.reco.label= uvascim.image.reco.label(ii,:,:);
                uvascim.image.reco.no_echoes = 1;
                uvascim.image.reco.phaselabel = uvascim.image.reco.phaselabel(ii,:,:);
                uvascim.image.reco.scaling_factor =uvascim.image.reco.scaling_factor(ii,:,:);
                uvascim.image.reco.scaling_offset = uvascim.image.reco.scaling_offset(ii,:,:);
                
                echo_label = handles.data.image(i).reco.echo_label{ii};
                switch echo_label
                    case 'signal intensity'
                        sub_file_name=strcat(file_name{:},'-Diff-SI.mat');
                        parameter_name= {'Diff-SI'};
                    case 'std dev of signal intensity'
                        sub_file_name=strcat(file_name{:},'-Diff-std_of_SI.mat');
                        parameter_name= {'Diff-std_of_SI'};
                    case 'diffusion constant'
                        sub_file_name=strcat(file_name{:},'-ADC.mat');
                        parameter_name= {'ADC'};
                    case 'std dev of diffusion constant'
                        sub_file_name=strcat(file_name{:},'-Diff-std_of_diff_ct.mat');
                        parameter_name= {'Diff-std_of_diff_ct'};
                    case 'std dev of the fit'
                        sub_file_name=strcat(file_name{:},'-Diff-std_of_fit.mat');
                        parameter_name= {'Diff-std_of_fit'};
                    otherwise
                        sub_file_name=strcat(file_name{:},'-',handles.data.image(i).reco.echo_label{ii}, '.mat');
                        parameter_name=strcat(list_scan(i),'-',handles.data.image(i).reco.echo_label{ii});
                end
                if exist(fullfile(strcat(tmp(id_nbr).path,sub_file_name)) , 'file') ~= 2
                    save(fullfile(strcat(tmp(id_nbr).path,sub_file_name)),'uvascim');
                end
                
                if ischar(sub_file_name)
                    sub_file_name = {sub_file_name};
                end
                tmp(id_nbr).day(tp_nbr).parameters = [tmp(id_nbr).day(tp_nbr).parameters parameter_name];
                tmp(id_nbr).day(tp_nbr).scans_file = [tmp(id_nbr).day(tp_nbr).scans_file sub_file_name];
                
                %updade handles.clips
                if isempty(strmatch(parameter_name, handles.clips(:,1), 'exact'))
                    handles.clips(size(handles.clips,1)+1,1) = parameter_name;
                    
                    if( max(~isreal(uvascim.image.reco.data(:)))==1 )  % complex data
                        handles.clips(size(handles.clips,1),2) = {min(abs(uvascim.image.reco.data(:)))};
                        handles.clips(size(handles.clips,1),3) =  {max(abs(uvascim.image.reco.data(:)))};
                    else    % real data
                        handles.clips(size(handles.clips,1),2) = {min(uvascim.image.reco.data(:))};
                        handles.clips(size(handles.clips,1),3) =  {max(uvascim.image.reco.data(:))};
                    end
                end
            end
            
        else
            
            % If data imported from uvasc --> save uvascim
            %tmp_uvascim = uvascim;
            uvascim = handles.data;
            uvascim.image = handles.data.image(i);
            %             uvascfile=fullfile(strcat(tmp(id_nbr).path,file_name{:},'.mat'));
            %             if exist(uvascfile, 'file') ~= 2
            %                 save(uvascfile,'uvascim');
            %             else
            %                 uvascfile=fullfile(strcat(tmp(id_nbr).path,file_name{:},'_',datestr(now,'yyyymmdd-HHMMSS'),'.mat'));
            %                 save(uvascfile,'uvascim');
            %
            %             end
            
            tmp(id_nbr).day(tp_nbr).parameters = [tmp(id_nbr).day(tp_nbr).parameters list_scan(i)];
            tmp(id_nbr).day(tp_nbr).scans_file = [tmp(id_nbr).day(tp_nbr).scans_file strcat(file_name, '.mat')];
            
            %updadehandles.clips
            if ~isfield(handles, 'clips') || isempty(strmatch(list_scan(i), handles.clips(:,1), 'exact'))
                if ~isfield(handles, 'clips')
                    handles.clips(1,1) = list_scan(i);
                    
                    if( max(~isreal(uvascim.image.reco.data(:)))==1 )    % complex data
                        handles.clips(1,2) =  {min(abs(uvascim.image.reco.data(:)))};
                        handles.clips(1,3) =  {max(abs(uvascim.image.reco.data(:)))};
                    else    % real data
                        handles.clips(1,2) =  {min(uvascim.image.reco.data(:))};
                        handles.clips(1,3) =  {max(uvascim.image.reco.data(:))};
                    end
                end
                handles.clips(size(handles.clips,1)+1,1) = list_scan(i);
                if( max(~isreal(uvascim.image.reco.data(:)))==1 )    % complex data
                    handles.clips(size(handles.clips,1),2) = {min(abs(uvascim.image.reco.data(:)))};
                    handles.clips(size(handles.clips,1),3) =  {max(abs(uvascim.image.reco.data(:)))};
                else    % real data
                    handles.clips(size(handles.clips,1),2) = {min(uvascim.image.reco.data(:))};
                    handles.clips(size(handles.clips,1),3) =  {max(uvascim.image.reco.data(:))};
                end
            end
        end
    end
    
    if isfield(handles, 'database')
        % add this patient to the database
        if id_nbr~=numel(handles.database)
            for i = id_nbr+1:numel(handles.database)
                tmp(i)=handles.database(i);
            end
        end
    end
    handles.database = tmp;
    clear tmp
    guidata(hObject, handles);
    MIA_update_database_display(hObject, eventdata, handles);
end
%update Figure name
MIA_update_figureName(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MIA_add_roi_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_add_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, directory]=uigetfile({'*.mat'}, 'Selection of ROI structure -uvascroi-','MultiSelect','on');
if isequal(filename,0)
    return
end
if ischar(filename)
    filename = {filename};
end

% Can load several structure
for x = 1:numel(filename)
    fid=fopen(fullfile(directory,filename{x}),'r');
    if fid>0
        fclose(fid);
        roi_loaded = load(fullfile(directory,filename{x}));
        if ~isfield(roi_loaded, 'uvascroi')
            warndlg(strcat('something wrong with the ROI file:', {' '},filename{x}),'Warning');
            return
        end
    else
        warndlg(strcat('something wrong with the ROI file:', {' '},filename{x}),'Warning');
        return
    end
    dlgbox_title =  strcat('File name:', {' '}, filename{x});
    id_info = filename{x};
    sep_position = strfind(id_info, '_');
    %try to figure out the patient's name and time point
    if ~isempty(sep_position) && numel(sep_position) > 2
        id = id_info(1:sep_position(1)-1);
        tp = id_info(sep_position(1)+1:sep_position(2)-1);
    else
        id = id_info;
        tp = id_info;
    end
    answer = inputdlg({'Is the id''s name correct?', 'Is the time point correct?'},dlgbox_title{:}, [1 size(dlgbox_title{:},2)+25], [{id} {tp}]);
    if isempty(answer)
        return
    end
    % find if the patient exist already
    if ~isempty(strmatch(answer(1), get(handles.MIA_name_list, 'String'), 'exact')) %#ok<*MATCH3>
        id_nbr = strmatch(answer(1), get(handles.MIA_name_list, 'String'), 'exact');
        if ~isempty(strmatch(answer(2), {handles.database(id_nbr).day.date}', 'exact')) %#ok<*MATCH3>
            tp_nbr = strmatch(answer(2), {handles.database(id_nbr).day.date}', 'exact');
        else
            warndlg('No time point corresponds to this ROI','Warning');
            return
        end
    else
        warndlg('No Patient''s id corresponds to this ROI','Warning');
        return
    end
    roi_names = unique({roi_loaded.uvascroi.name}');
    for i =1:numel(roi_names)
        clear roitmp
        aa = 1;
        roi_num = strcmp(roi_names(i), {roi_loaded.uvascroi.name}');
        for z = 1:numel(roi_num)
            if roi_num(z)
                roitmp(aa) = roi_loaded.uvascroi(z); %#ok<AGROW>
                aa = aa+1;
            end
        end
        uvascroi = roitmp;  %#ok<NASGU>
        %save ROI's file -- 1 file per ROI's name
        roifile=fullfile(strcat(handles.database(id_nbr).path,id_info,'-ROI-', roi_names{i},'.mat'));
        save(roifile,'uvascroi');
        %update database
        handles.database(id_nbr).day(tp_nbr).VOIs = [handles.database(id_nbr).day(tp_nbr).VOIs roi_names(i)];
        handles.database(id_nbr).day(tp_nbr).VOIs_file = [handles.database(id_nbr).day(tp_nbr).VOIs_file strcat(id_info,'-ROI-', roi_names(i),'.mat')];
        %updade handles.VOIs
        if isempty(strmatch(roi_names(i), handles.VOIs, 'exact'))
            handles.VOIs = [handles.VOIs(1:end-1), roi_names(i), handles.VOIs(end-1:end)];
        end
    end
    clear roi_name;
    guidata(hObject, handles);
    MIA_update_database_display(hObject, eventdata, handles)
end



% --------------------------------------------------------------------
function MIA_rename_name_Callback(hObject, ~, handles)
% hObject    handle to MIA_rename_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end

handles = guidata(hObject);

patient = get(handles.MIA_name_list, 'Value');
old_patient_name = handles.database(patient).name;
if numel(patient) >1
    warndlg('Please select only one timepoint', 'Warning');
    return
end

new_patient_name = inputdlg('Enter the new patient name',...
    old_patient_name, 1, {old_patient_name});
if isempty(new_patient_name)
    return
end

if isfield(handles, 'database_all') && numel(handles.database) ~= numel(handles.database_all)
    size = length(handles.database(patient).name);
    %     num_patient = find(strncmp(handles.database(patient).name, {handles.database_all.name}', size) == 1);
    num_patient = strncmp(handles.database(patient).name, {handles.database_all.name}', size) == 1;
    handles.database_all(num_patient).name = new_patient_name{:};
end

handles.database(patient).name = new_patient_name{:};

% rename each parameter files
for timepoint = 1:numel(handles.database(patient).day)
    for i = 1:numel(handles.database(patient).day(timepoint).scans_file)
        [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).scans_file{i}]);
        new_name = [new_patient_name{:}  '-' handles.database(patient).day(timepoint).date '-' handles.database(patient).day(timepoint).parameters{i}];
        if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
            warning_text = sprintf('##$ Can not rename this file because it do not exist\n##$ %s',...
                fullfile(PATHSTR,[NAME,EXT]));
            msgbox(warning_text, 'rename file warning') ;
        elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
            user_response = questdlg('The new file name for this scan exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
            if strcmp(user_response, 'No')
                return
            end
            user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
            % overwrite the file if requested
            if strcmp(user_response, 'Yes')
                if ~strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]))
                    movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                    handles.database(patient).day(timepoint).scans_file{i} = [new_name,EXT];
                end
            end
            
        else
            movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
            handles.database(patient).day(timepoint).scans_file{i} = [new_name,EXT];
        end
    end
end
% rename each ROI files
for timepoint = 1:numel(handles.database(patient).day)
    for i = 1:numel(handles.database(patient).day(timepoint).VOIs_file)
        [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).VOIs_file{i}]);
        new_name = [new_patient_name{:} '-' handles.database(patient).day(timepoint).date '-' handles.database(patient).day(timepoint).VOIs{i}];
        if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
            warning_text = sprintf('##$ Can not rename this file because it do not exist\n##$ %s',...
                fullfile(PATHSTR,[NAME,EXT]));
            msgbox(warning_text, 'rename file warning') ;
        elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
            user_response = questdlg('The new file name for this scan exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
            if strcmp(user_response, 'No')
                return
            end
            user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
            % Overwrite the file if requested
            if strcmp(user_response, 'Yes')
                if ~strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]))
                    movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                    handles.database(patient).day(timepoint).VOIs_file{i} = [new_name,EXT];
                end
            end
            
        else
            movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
            handles.database(patient).day(timepoint).VOIs_file{i} = [new_name,EXT];
        end
    end
end






guidata(hObject, handles);

%%% update graph and display
MIA_update_database_display(hObject, eventdata, handles);


% --------------------------------------------------------------------
function MIA_remove_name_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_remove_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end
data_selected = get_data_selected(handles);
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
    switch selection,
        case 'No'
            return
    end
end
path_root=pwd;
[filename, pathname]=uigetfile('*.mat','Open Mat File','MultiSelect','off');

if pathname == 0
    return
else
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
        
        %         VOIs_list = {};
        %         for i=1:numel(handles.database)
        %             for j = 1:numel(handles.database(i).day)
        %                 VOIs_list = [VOIs_list handles.database(i).day(j).VOIs];
        %             end
        %         end
        %
        %         handles.VOIs = [unique(VOIs_list), 'Automatic', 'Other'];
        %update the database ... because there is a bug somewhere...
        %         database.database(1).databaseinfo.VOIs = handles.VOIs;
        
        %         handles.clips = database.database(1).databaseinfo.clips;
        %         if isfield(database.database(1).databaseinfo, 'histo');
        %             handles.histo = database.database(1).databaseinfo.histo;
        %         end
        %         if ~isfield(handles.database(1).databaseinfo, 'voxels_database_filename')
        %             handles.database(1).databaseinfo.voxels_database_filename = [];
        %             handles.database(1).databaseinfo.voxels_database_need_to_update = 0;
        %             handles.database(1).databaseinfo.voxels_database_voxels_to_update = [];
        %         end
        set(handles.MIA_name_list, 'String', handles.database.Properties.UserData.db_filename)
        guidata(hObject, handles);
    end
    cd(path_root);
end
MIA_update_figureName(hObject, eventdata, handles)
% update database.mat path (required if this file has been moved
% handles.database(1).databaseinfo.pathname = pathname;
% handles.database(1).databaseinfo.filename = filename;

guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);



% --------------------------------------------------------------------
function MIA_time_points_right_click_Callback(~, ~, ~)
% hObject    handle to MIA_time_points_right_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_add_time_points_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_add_time_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(get(handles.MIA_name_list, 'String'))
    return
end
patient = get(handles.MIA_name_list, 'Value');
new_time_point = inputdlg({'New time point'}, 'Add time point');
if  isempty(new_time_point)
    return
end
numero = numel(handles.database(patient).day);
if numero == 1 && isempty(handles.database(patient).day(1).date)
    handles.database(patient).day(numero).date = new_time_point{:};
else
    handles.database(patient).day(numero+1).date = new_time_point{:};
end
%update display
guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);


% --------------------------------------------------------------------
function MIA_rename_time_point_Callback(hObject, eventdata,handles)
% hObject    handle to MIA_rename_time_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
patient = get(handles.MIA_name_list, 'Value');
if ~isfield(handles, 'database') || isempty(handles.database(patient).day(1).date)
    return
end


timepoint = get(handles.MIA_time_points_list, 'Value');
old_timepoint_name = handles.database(patient).day(timepoint).date;

if numel(timepoint) >1
    warndlg('Please select only one timepoint', 'Warning');
    return
end

new_timepoint_name = inputdlg('Enter the new time point name',...
    old_timepoint_name, 1, {old_timepoint_name});
if isempty(new_timepoint_name)
    return
end
handles.database(patient).day(timepoint).date = new_timepoint_name{:};
% rename each parameter files
warning_text = {};

for i = 1:numel(handles.database(patient).day(timepoint).scans_file)
    [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).scans_file{i}]);
    cd(PATHSTR);
    new_name = [handles.database(patient).name '-' new_timepoint_name{:} '-' handles.database(patient).day(timepoint).parameters{i}];
    new_name = strrep(new_name, '*', 'star');
    if strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT])) == 1
        continue
    end
    %     if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
    %         warning_text{numel(warning_text)+1} = sprintf('##$ Can not rename this file because it do not exist\n##$ %s',...
    %             fullfile(PATHSTR,[NAME,EXT]));
    %
    %     elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
    %         user_response = questdlg('The new file name for this scan exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
    %         if strcmp(user_response, 'No')
    %             return
    %         end
    %         user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
    %         % Overwrite the file if requested
    %         if strcmp(user_response, 'Yes')
    % %             movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
    %             dos(['rename "' [NAME EXT] '" "' [new_name EXT] '"']);
    %             handles.database(patient).day(timepoint).scans_file{i} = [new_name,EXT];
    %         end
    %
    %     else
    %         movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
    
    if ispc
        dos(['rename "' [NAME EXT] '" "' [new_name EXT] '"']);
    elseif isunix
        movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f')
    end
    handles.database(patient).day(timepoint).scans_file{i} = [new_name,EXT];
    %     end
    % save database after each file renamed
    path_root = pwd;
    cd(handles.database(1).databaseinfo.pathname);
    database = handles.database; %#ok<NASGU>
    save(handles.database(1).databaseinfo.filename, 'database');
    cd(path_root);
    
    guidata(hObject, handles);
    %%% update graph and display
    MIA_update_database_display(hObject, eventdata, handles);
end
% rename each ROI files
for i = 1:numel(handles.database(patient).day(timepoint).VOIs_file)
    [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).VOIs_file{i}]);
    cd(PATHSTR);
    new_name = [handles.database(patient).name '-' new_timepoint_name{:} '-' handles.database(patient).day(timepoint).VOIs{i}];
    if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
        warning_text{numel(warning_text)+1}  = sprintf('##$ Can not rename this file because it do not exist\n##$ %s',...
            fullfile(PATHSTR,[NAME,EXT]));
        %         msgbox(warning_text, 'rename file warning') ;
    elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
        user_response = questdlg('The new file name for this scan exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
        if strcmp(user_response, 'No')
            return
        end
        user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
        % overwrite the file if requested
        if strcmp(user_response, 'Yes')
            %             movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
            if ispc
                dos(['rename "' [NAME EXT] '" "' [new_name EXT] '"']);
            elseif isunix
                movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f')
            end
            handles.database(patient).day(timepoint).VOIs_file{i} = [new_name,EXT];
        end
        
    else
        
        %         movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
        if ispc
            dos(['rename "' [NAME EXT] '" "' [new_name EXT] '"']);
        elseif isunix
            movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f')
        end
        handles.database(patient).day(timepoint).VOIs_file{i} = [new_name,EXT];
        
    end
    % save database after each file renamed
    
    path_root = pwd;
    cd(handles.database(1).databaseinfo.pathname);
    database = handles.database; %#ok<NASGU>
    save(handles.database(1).databaseinfo.filename, 'database');
    cd(path_root);
    
    guidata(hObject, handles);
    %%% update graph and display
    MIA_update_database_display(hObject, eventdata, handles);
end
msgbox(warning_text, 'rename file warning') ;
% msgbox('Done', 'Message') ;


% --------------------------------------------------------------------
function MIA_add_histo_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_add_histo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
patient = get(handles.MIA_name_list, 'Value');
time_point = get(handles.MIA_time_points_list, 'Value');
% path = handles.database(patient).path;
% cd(path)
% scan_number = 1;

if get(handles.MIA_scan_VOIs_button, 'Value') == 0  % add news Scan
    [scan, ok1] = listdlg('PromptString','Select a staining type:',...
        'Name', 'Selecting...',...
        'SelectionMode','single',...
        'ListSize', [400 300],...
        'ListString',handles.histo');
    if ok1 == 0
        return
    end
    if strcmp('Other', handles.histo(scan)) == 1
        newparameter = inputdlg('Name of the new histo scan ', 'Question?', 1, {''});
        if isempty(newparameter)
            return
        end
        handles.histo = [handles.histo(1:end-1), newparameter, handles.histo(end)];
        scan = numel(handles.histo)-1;
        
    end
    % Grey of color scale
    color_map = questdlg('What is the color map?', 'Bip', 'Gray', 'Color', 'Cancel', 'Cancel');
    if isempty(color_map) || strcmp(color_map, 'Cancel')
        return
    end
    % load Pimo or H&E images
    %     if scan == 1 || scan == 2
    [file, pathname]=uigetfile({'*.jpg;*.tif;*.tiff; *.inf; *.wfml;','Image Files (*.jpg,*.tif...)'}, 'Please select an image image','MultiSelect','off');
    if isa(file, 'double')
        return
    end
    [~, ~, ext] =fileparts([pathname file]);
    if file == 0
        return
    end
end
% Load the histology image
if strcmp(ext, '.inf') % case of binary images
    info = importdata(fullfile(pathname,file));
    if info.data(3) == 16
        G = 65535;
        histo_image = multibandread(fullfile(pathname,strcat(info.textdata{2}, '.img')), [info.data(5) info.data(4) 1], 'ubit16=>uint16', 0, 'bsq', 'b');
    else
        histo_image = multibandread(fullfile(pathname,strcat(info.textdata{2}, '.img')), [info.data(5) info.data(4) 1], 'ubit8=>uint8', 0, 'bsq', 'b');
        G= 255;
    end
elseif strcmp(ext, '.wfml') % need to explore sub directories
    fid=fopen(fullfile(pathname,file) ,'r');
    if fid>0
        fclose(fid);
        info_txt = fileread(fullfile(pathname,file));
        char_beg = strfind(info_txt, '<image-path>');
        char_end = strfind(info_txt, '</image-path>');
        for i=1:numel(char_beg)
            file_name= info_txt(char_beg(i)+12:char_end(i)-1);
            if numel(strfind(file_name, '/')) >0
                file_name = strrep(file_name, '/', '\');
            end
            file_path{i} = [pathname file_name];
        end
        file_path = sort(file_path);
        file_path_svs = file_path(((numel(file_path)/4)*3)+1:end);
        [image_selected, ok] = listdlg('PromptString','Select one image:',...
            'Name', 'Select a file',...
            'SelectionMode','single',...
            'ListSize', [400 300],...
            'ListString',file_path_svs');
        if ok == 0
            warndlg('something wrong with the data','Warning');
            return
        end
        info = imfinfo(file_path_svs{image_selected});
        histo_image = importdata(file_path_svs{image_selected});
    else
        warndlg('something wrong with the data','Warning');
        return
    end
    
    
else
    fid=fopen(fullfile(pathname,file) ,'r');
    if fid>0
        fclose(fid);
        info = imfinfo(fullfile(pathname,file));
        histo_image = imread(fullfile(pathname,file));
    else
        warndlg('something wrong with the data','Warning');
        return
    end
end
factor = questdlg('Do you want to multiply your histo by a factor?', 'Bip', 'Yes', 'No', 'Cancel', 'Cancel');
if strcmp(factor, 'Yes')
    factor = inputdlg('Enter the factor', 'bip', 1);
    if isempty(factor)
        return
    end
    factor = str2double(factor);
else
    factor = 1;
end
histo_image = histo_image*factor;

registration = questdlg('Do you want to register you image to an MRI scan or do you want to appply a transformation?', 'Bip', 'New registration', 'Apply transformation', 'Cancel', 'Cancel');
if isempty(registration) || strcmp(registration, 'Cancel')
    return
end
if strcmp(registration, 'New registration') == 1
    % Load the reference MRI scan
    [new_time_point, ok3] = listdlg('PromptString','Select a the time point for the scan of reference:',...
        'Name', 'Select a time point',...
        'SelectionMode','single',...
        'ListSize', [400 300],...
        'ListString',{handles.database(patient).day.date}');
    if ok3 == 0
        return
    end
    
    [scan_ref_nbr, ok4] = listdlg('PromptString','Select a scan of reference:',...
        'Name', 'Select a reference',...
        'SelectionMode','single',...
        'ListSize', [400 300],...
        'ListString',handles.database(patient).day(new_time_point).parameters');
    if ok4 == 0
        return
    end
    fid=fopen(fullfile(handles.database(patient).path,handles.database(patient).day(new_time_point).scans_file{scan_ref_nbr}) ,'r');
    if fid>0
        fclose(fid);
        scan_ref_loaded = load(fullfile(handles.database(patient).path,handles.database(patient).day(new_time_point).scans_file{scan_ref_nbr}));
        scan_ref = squeeze(scan_ref_loaded.uvascim.image.reco.data(:,:,1,:));
    else
        warndlg('something wrong with the data','Warning');
        return
    end
end


guidata(hObject, handles);

% Run the registration module if needed
switch registration
    case 'New registration'
        histo_registred = IRM_Histo_reg_module(hObject, eventdata, scan_ref, histo_image);%(:,:,1));
        while ~isempty(findobj('Tag', 'IHr_GUI'))
            pause(1)
        end
        global output_data %#ok<TLEV>
        histo_registred = output_data;
        clear global
        clear output_data
    case 'Apply transformation'
        
        % Load the histological image already registered
        [new_time_point_histo, ok3] = listdlg('PromptString','Select a the time point for the histological image already registered :',...
            'Name', 'Select a time point',...
            'SelectionMode','single',...
            'ListSize', [400 300],...
            'ListString',{handles.database(patient).day.date}');
        if ok3 == 0
            return
        end
        
        [scan_ref_nbr_histo, ok4] = listdlg('PromptString','Select the histological image already registered:',...
            'Name', 'Select a reference',...
            'SelectionMode','single',...
            'ListSize', [400 300],...
            'ListString',handles.database(patient).day(new_time_point_histo).parameters');
        if ok4 == 0
            return
        end
        fid=fopen(fullfile(handles.database(patient).path,handles.database(patient).day(new_time_point_histo).scans_file{scan_ref_nbr_histo}) ,'r');
        if fid>0
            fclose(fid);
            histo_registered_loaded = load(fullfile(handles.database(patient).path,handles.database(patient).day(new_time_point_histo).scans_file{scan_ref_nbr_histo}));
        else
            warndlg('something wrong with the data','Warning');
            return
        end
        histo_registred.mytform = histo_registered_loaded.uvascim.image.reco.mytform;
        tmp = imread(fullfile(pathname,file), 1);
        if strcmp(color_map, 'Gray')
            histo_registred.image(:,:,1,1) = imtransform(tmp*factor, histo_registred.mytform, 'XData', [1 size(tmp,2)], 'YData',[1 size(tmp,1)]);
        else
            histo_registred.image(:,:,1,1) = imtransform(tmp(:,:,1)*factor, histo_registred.mytform, 'XData', [1 size(tmp,2)], 'YData',[1 size(tmp,1)]);
            histo_registred.image(:,:,1,2) = imtransform(tmp(:,:,2)*factor, histo_registred.mytform, 'XData', [1 size(tmp,2)], 'YData',[1 size(tmp,1)]);
            histo_registred.image(:,:,1,3) = imtransform(tmp(:,:,3)*factor, histo_registred.mytform, 'XData', [1 size(tmp,2)], 'YData',[1 size(tmp,1)]);
        end
end

if isempty(histo_registred)
    return
else
    % if stack
    %     histo_registred = output_data;
    if numel(info) > 2
        for i = 2:numel(info)
            tmp = imread(fullfile(pathname,file), i);
            if isfield(histo_registred, 'image_crop')
                warndlg('A stack registration of croped images is not coded yet', 'Warning');
                return
            else
                histo_registred.image(:,:,1,i) = imtransform(tmp*factor, histo_registred.mytform, 'XData', [1 size(tmp,2)], 'YData',[1 size(tmp,1)]);
            end
        end
    end
end

% Check if a same kind of histo is already in the database
% if yes add the new histo data in the same .mat
if sum(strcmp(handles.histo(scan), handles.database(patient).day(time_point).parameters')) > 0
    
else
    if isfield(histo_registred, 'image_crop')
        % check if the registered image is in color or B/W
        if strcmp(color_map, 'Gray')
            uvascim.image.reco.data(:,:,1,:)= histo_registred.image_crop(:,:);
        else
            uvascim.image.reco.data(:,:,1,1,1:3)=histo_registred.image_crop(:,:,:);
        end
        
    else
        if strcmp(color_map, 'Gray')
            uvascim.image.reco.data(:,:,1,:)= histo_registred.image;
        else
            %uvascim.image.reco.data(:,:,1,1,1:3)=histo_registred.image(:,:,:);
            uvascim.image.reco.data(:,:,1,1,1:3)=  imtransform(histo_image*factor, histo_registred.mytform, 'XData', [1 size(histo_image,2)], 'YData',[1 size(histo_image,1)]); %
        end
    end
    % save the slice thickness, color map and Row file
    slice_thickness = inputdlg('Enter the new slice thickness in micrometer!!', 'bip', 1);
    uvascim.image.reco.thickness = str2double(slice_thickness)/1000;
    uvascim.image.reco.color_map = color_map;
    uvascim.image.reco.no_slices = 1;
    uvascim.image.reco.raw_filename = file;
    uvascim.image.reco.raw_pathname = pathname;
    uvascim.image.reco.mytform = histo_registred.mytform;
    switch registration
        case 'New registration'
            ParamConfig=sprintf('##$Method=Image registered onto \n##$filename=%s\n##$On slice z=%d\n',...
                fullfile(handles.database(patient).path,handles.database(patient).day(new_time_point).scans_file{scan_ref_nbr}),...
                histo_registred.z_offset);
            uvascim.image.reco.paramQuantif = ParamConfig;
            % save the slice offset package (the offset is located in the middle of the
            % referece slice
            if mod(numel(info),2)
                % even
                for i = 1:numel(info)
                    uvascim.image.reco.fov_offsets(1:3,1,i) = [scan_ref_loaded.uvascim.image.reco.fov_offsets(1:2,1,histo_registred.z_offset)' ...
                        scan_ref_loaded.uvascim.image.reco.fov_offsets(3,1,histo_registred.z_offset) - uvascim.image.reco.thickness/2 - ...
                        (uvascim.image.reco.thickness*numel(info)/2) + ...
                        i * uvascim.image.reco.thickness]';
                    uvascim.image.reco.no_slices = uvascim.image.reco.no_slices+1;
                end
            else
                % odd/impair
                for i = 1:numel(info)
                    uvascim.image.reco.fov_offsets(1:3,1,i) = [scan_ref_loaded.uvascim.image.reco.fov_offsets(1:2,1,histo_registred.z_offset)' ...
                        scan_ref_loaded.uvascim.image.reco.fov_offsets(3,1,histo_registred.z_offset) - ...
                        (uvascim.image.reco.thickness*numel(info)/2) + ...
                        i * uvascim.image.reco.thickness]';
                    uvascim.image.reco.no_slices = uvascim.image.reco.no_slices+1;
                end
                
            end
        case 'Apply transformation'
            uvascim.image.reco.paramQuantif = histo_registered_loaded.uvascim.image.reco.paramQuantif;
            uvascim.image.reco.fov_offsets =  histo_registered_loaded.uvascim.image.reco.fov_offsets;
    end
    
    %set clip
    uvascim.image.clip(1) = min(uvascim.image.reco.data(:));
    uvascim.image.clip(2) = max(uvascim.image.reco.data(:));
end
% Save data
dot = strfind(file,'.');
histo_file_name_img = strcat(file(1:dot-1), '_R', file(dot: end));
histo_file_name_mat = strcat(file(1:dot-1), '_R.mat');

% Convert data to PSL (only for scan 3 and 4 because
% autoradiography data)
% if strcmp(handles.histo(scan), 'Tc+I_cerveau') || strcmp(handles.histo(scan), 'Tc+I_ref')
%     if isfield(histo_registred, 'image_crop')
%
%         % for both the crop data registered and unchanged
%         histo_registred.image_crop_double = double(histo_registred.image_crop);
%         ct1 = double(repmat((info.data(1)/100)^2,size(histo_registred.image_crop_double,1), size(histo_registred.image_crop_double,2)));
%         ct2 = double(repmat(4000/info.data(6),size(histo_registred.image_crop_double,1), size(histo_registred.image_crop_double,2)));
%         cte_exp = double(repmat(0.5,size(histo_registred.image_crop_double,1), size(histo_registred.image_crop_double,2)));
%         L = double(repmat(info.data(7),size(histo_registred.image_crop_double,1),size(histo_registred.image_crop_double,2)));
%         exp = L.*(histo_registred.image_crop_double./G)-cte_exp;
%         histo_registred.image_crop = (ct1.*ct2.*10.^exp)./100;
%         % Convert Original data croped --> PSL
%         histo_registred.image_crop_unchanged_double = double(histo_registred.image_crop_unchanged);
%         ct1 = double(repmat((info.data(1)/100)^2,size(histo_registred.image_crop_unchanged_double,1), size(histo_registred.image_crop_unchanged_double,2)));
%         ct2 = double(repmat(4000/info.data(6),size(histo_registred.image_crop_unchanged_double,1), size(histo_registred.image_crop_unchanged_double,2)));
%         cte_exp = double(repmat(0.5,size(histo_registred.image_crop_unchanged_double,1), size(histo_registred.image_crop_unchanged_double,2)));
%         L = double(repmat(info.data(7),size(histo_registred.image_crop_unchanged_double,1),size(histo_registred.image_crop_unchanged_double,2)));
%         exp = L.*(histo_registred.image_crop_unchanged_double./G)-cte_exp;
%         histo_registred.image_crop_unchanged = (ct1.*ct2.*10.^exp)./100;
%         histo_file_name_tiff
%         % define the new file names
%         size_img = size(histo_registred.image_crop);
%         histo_file_name_img = strcat(file(1:dot-1), '-',num2str(size_img(1)), 'x', num2str(size_img(2)),'-crop_PSL-R.tiff');
%         histo_file_name_tiff_unchanged = strcat(file(1:dot-1),'-', num2str(size_img(1)), 'x', num2str(size_img(2)),'-crop_PSL.tiff');
%         histo_file_name_mat = strcat(histo_file_name_img(1:end-4), 'mat');
%         histo_path_name_tiff_unchanged = fullfile(handles.database(patient).path,histo_file_name_tiff_unchanged);
%         histo_path_name_mat_unchanged = fullfile(handles.database(patient).path ,strcat(histo_file_name_tiff_unchanged(1:end-4), 'mat'));
%     else
%         % Convert data to PSL (only for scan 3 and 4 because
%         % autoradiography data)
%         histo_registred.image_double = double(histo_registred.image);
%         ct1 = double(repmat((info.data(1)/100)^2,size(histo_registred.image_double,1), size(histo_registred.image_double,2)));
%         ct2 = double(repmat(4000/info.data(6),size(histo_registred.image_double,1), size(histo_registred.image_double,2)));
%         cte_exp = double(repmat(0.5,size(histo_registred.image_double,1), size(histo_registred.image_double,2)));
%         L = double(repmat(info.data(7),size(histo_registred.image_double,1),size(histo_registred.image_double,2)));
%         exp = L.*(histo_registred.image_double./G)-cte_exp;
%         histo_registred.image = (ct1.*ct2.*10.^exp)./100;
%
%         histo_file_name_img = strcat(file(1:dot-1),'PSL_R');
%         histo_file_name_mat = strcat(histo_file_name_img(1:end-4), 'mat');
%     end
% end

histo_path_name_img = fullfile(handles.database(patient).path,histo_file_name_img);
histo_path_name_mat = fullfile(handles.database(patient).path, histo_file_name_mat);

% save data .mat and .Tiff (.jpg, ....)
if isfield(histo_registred, 'image_crop')
    % save original file name
    uvascim.image.reco.original_file_name = file;
    % Registered data
    imwrite(histo_registred.image_crop ,histo_path_name_img, 'tiff');
    save(histo_path_name_mat,'uvascim');
    % Original data
    imwrite(histo_registred.image_crop_unchanged ,histo_path_name_tiff_unchanged, 'tiff');
    save(histo_path_name_mat_unchanged,'uvascim');
    
else
    % save original file name
    uvascim.image.reco.original_file_name = file;
    save(histo_path_name_mat,'uvascim');
    for i = 1:numel(info)
        if strcmp(file(dot+1: end), 'jpg')
            imwrite(histo_registred.image(:,:,1,i) ,histo_path_name_img, file(dot+1: end));
        elseif strcmp(file(dot+1: end), 'tif') || strcmp(file(dot+1: end), 'tiff')
            % set the tiff format
            tagstruct.ImageLength = size(histo_registred.image,1);
            tagstruct.ImageWidth = size(histo_registred.image,2);
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample = 32;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.SampleFormat = 3;
            tagstruct.SampleFormat = 1;
            %tagstruct.RowsPerStrip = 16;
            tagstruct.PlanarConfiguration = 1;
            tagstruct.Software = 'MATLAB';
            if i == 1
                t = Tiff(histo_path_name_img,'w');
            else
                t = Tiff(histo_path_name_img,'a');
            end
            t.setTag(tagstruct);
            t.write(uint32(histo_registred.image(:,:,1,i)));
            t.close();
        else
            imwrite(histo_registred.image(:,:,1,i) ,histo_path_name_img, file(dot+1: end),  'WriteMode', 'append');
            
        end
        
    end
end

if ~isfield(handles.database(patient).day(time_point), 'parameters')
    handles.database(patient).day(time_point).parameters = handles.histo(scan);
    handles.database(patient).day(time_point).scans_file = {histo_file_name_mat};
else
    handles.database(patient).day(time_point).parameters = [handles.database(patient).day(time_point).parameters handles.histo(scan)];
    handles.database(patient).day(time_point).scans_file = [handles.database(patient).day(time_point).scans_file {histo_file_name_mat}];
end

if isempty(strmatch(handles.histo(scan), handles.clips(:,1), 'exact'))
    handles.clips(numel(handles.clips(:,1))+1,1) =  handles.histo(scan);
    handles.clips(numel(handles.clips(:,1)),2) =  num2cell(uvascim.image.clip(1));
    handles.clips(numel(handles.clips(:,1)),3) =  num2cell(uvascim.image.clip(2));
end

guidata(hObject, handles);
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

% if get(handles.MIA_scan_VOIs_button, 'Value')
%     msgbox('Not coded yet warning') ;
%     %     name_option = handles.VOIs';
%     %     old_scan_name = handles.database(patient).day(timepoint).VOIs_file(scan_name);
%     return
% else
%     name_option = [cellstr(unique(handles.database.SequenceName))' 'Other']';
%     old_scan_name = handles.database.SequenceName(data_selected);
% end
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

% if get(handles.MIA_scan_VOIs_button, 'Value')
%     if strcmp('Other', handles.VOIs(new_scan_name)) == 1
%         newparameter = inputdlg('Name of the new VOI ', 'Question?', 1, {''});
%         handles.VOIs=[handles.VOIs(1:end-1), newparameter, handles.VOIs(end-1:end)];
%         new_scan_name = numel(handles.VOIs)-2;
%         handles.database(patient).day(timepoint).VOIs(scan_name) = handles.VOIs(new_scan_name);
%     else
%         handles.database(patient).day(timepoint).VOIs(scan_name) = handles.VOIs(new_scan_name);
% %         newparameter = name_option(new_scan_name);
% %     end
% % else
% %     if strcmp('Other',name_option(new_scan_name)) == 1
% %         newparameter = inputdlg('Name of the new Scan ', 'Question?', 1, {''});
% %         handles.database(patient).day(timepoint).parameters(scan_name) = newparameter;
% %         if sum(strcmp(newparameter, handles.clips(:,1))) == 0
% %             handles.clips(size(handles.clips,1)+1,1) =newparameter;
% %             handles.clips(size(handles.clips,1),2) = handles.clips(strcmp(old_scan_name, handles.clips(:,1)),2);
% %             handles.clips(size(handles.clips,1),3) = handles.clips(strcmp(old_scan_name, handles.clips(:,1)),3);
% %         end
% %     else
% %         handles.database(patient).day(timepoint).parameters(scan_name) = name_option(new_scan_name);
% %         newparameter = name_option(new_scan_name);
% %     end
% % end
%
% % rename the scan file
% if get(handles.MIA_scan_VOIs_button, 'Value')
%     [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).VOIs_file{scan_name}]);
%     new_name = [handles.database(patient).name '-' handles.database(patient).day(timepoint).date '-' newparameter{:}];
%     if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
%         warning_text = sprintf('##$ This file no not exist\n##$ %s',...
%             fullfile(PATHSTR,[NAME,EXT]));
%         msgbox(warning_text, 'rename file warning') ;
%     elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
%         user_response = questdlg('The new file name for this ROI exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
%         if strcmp(user_response, 'No')
%             return
%         end
%         user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
%         % overwrite the file if requested
%         if strcmp(user_response, 'Yes')
%             if ~strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]))
%                 movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
%                 handles.database(patient).day(timepoint).VOIs_file{scan_name} = [new_name,EXT];
%             end
%             %% Add for modify name in ROI
%             load([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{scan_name}])
%             for itRoi = 1 : numel(uvascroi)
%                 uvascroi(itRoi).name = newparameter{1};
%             end
%             save([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{scan_name}],'uvascroi')
%         end
%
%     else
%         movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
%         handles.database(patient).day(timepoint).VOIs_file{scan_name} = [new_name,EXT];
%         %% Add for modify name in ROI
%         load([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{scan_name}])
%         for itRoi = 1 : numel(uvascroi)
%             uvascroi(itRoi).name = newparameter{1};
%         end
%         save([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{scan_name}],'uvascroi')
%     end
% else
%     [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).scans_file{scan_name}]);
%     new_name = [handles.database(patient).name '-' handles.database(patient).day(timepoint).date '-' newparameter{:}];
%     new_name = strrep(new_name, '*', 'star');
%     if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
%         warning_text = sprintf('##$ This file no not exist\n##$ %s',...
%             fullfile(PATHSTR,[NAME,EXT]));
%         msgbox(warning_text, 'rename file warning') ;
%     elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
%         user_response = questdlg('The new file name exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
%         if strcmp(user_response, 'No')
%             return
%         end
%         user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
%         % overwrite the file if requested
%         if strcmp(user_response, 'Yes')
%             if ~strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]))
%                 movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
%                 handles.database(patient).day(timepoint).scans_file{scan_name} = [new_name,EXT];
%             end
%         end
%
%     else
%         movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
%         handles.database(patient).day(timepoint).scans_file{scan_name} = [new_name,EXT];
%     end
% end
%
% guidata(hObject, handles);
%
% % check if the old parameter exist somewhere else, if not remove it
% parameters_list = [];
% for i=1:numel(handles.database)
%     for j = 1:numel(handles.database(i).day)
%         parameters_list = [parameters_list handles.database(i).day(j).parameters]; %#ok<AGROW>
%
%     end
% end
% parameters_list= unique(parameters_list);
%
% if sum(strcmp(old_scan_name, parameters_list')) == 0
%     match = find(strcmp(old_scan_name,handles.clips(:,1)), 1);
%     handles.clips(match,:) = [];
% end




% --------------------------------------------------------------------
function MIA_remove_scan_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_remove_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end

nii_index = get_data_selected(handles);
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
data_selected = get_data_selected(handles);
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

% --------------------------------------------------------------------
function MIA_change_directory_patient_Callback(hObject, ~, handles)
% hObject    handle to MIA_change_directory_patient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end

%handles = guidata(hObject);
% patient = get(handles.MIA_name_list, 'Value');
new_patient_directory = uigetdir(handles.database.Properties.UserData.MIA_root_path, 'Select Directory');
new_patient_directory = strcat(new_patient_directory, filesep);
if sum(new_patient_directory) == 0
    return
else
    handles.database.Path =categorical(strrep(cellstr(handles.database.Path ),...
        handles.database.Properties.UserData.MIA_root_path, new_patient_directory));
    handles.database.Properties.UserData.MIA_root_path = new_patient_directory;
    handles.database.Properties.UserData.MIA_data_path = [new_patient_directory 'Raw_data' filesep];
end
guidata(hObject, handles);


function MIA_load_axes_Callback(hObject, eventdata, handles)

if ~isfield(handles, 'database')
    return
end
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
    
    list_day = ['-1', string(unique(handles.data_loaded.info_data_loaded.Tp))'];
    set(handles.MIA_PRM_ref_popupmenu, 'String', list_day', 'Value', 2);
    %set MIA_PRM_slider
    set(handles.MIA_PRM_slider_tp, 'Max', handles.data_loaded.number_of_scan);
    set(handles.MIA_PRM_slider_tp,'Value',1);
    set(handles.MIA_PRM_slider_tp,'Min',1);
    if handles.data_loaded.number_of_scan == 2
        set(handles.MIA_PRM_slider_tp,'Visible', 'off');
    else
        set(handles.MIA_PRM_slider_tp,'Visible', 'on');
        set(handles.MIA_PRM_slider_tp,'SliderStep',[1/(handles.data_loaded.number_of_scan-1) min(5/(handles.data_loaded.number_of_scan-1),1)]);
    end
end

MIA_update_axes(hObject, eventdata, handles)


function handles = MIA_load_VOIs(hObject, ~, handles)

data_selected = get_data_selected(handles);
% get patient's information (already loaded)
% if handles.mode == 1
%     patient_loaded = handles.data_selected.patient_info.name_nbr;
%     timepoint_loaded = handles.data_selected.patient_info.timepoint_nbr;
% else
%     patient_loaded = handles.data_selected_for_PRM.patient_info.name_nbr;
%     timepoint_loaded = handles.data_selected_for_PRM.patient_info.timepoint_nbr;
% end
% if patient_loaded ~= patient
%     user_response = questdlg('Wrong PATIENT selected. Do you still want to open the ROIs?', 'Warning', 'Yes', 'No', 'No');
%     if strcmp(user_response, 'No')
%         return
%     end
% end
% if timepoint_loaded ~= time_point
%     if strcmp(get(hObject, 'Tag'), 'MIA_load_axes')
%         user_response = questdlg('Wrong TIME POINT selected. Do you still want to open the ROIs?', 'Warning', 'Yes', 'No', 'No');
%         if strcmp(user_response, 'No')
%             return
%         end
%     end
% end
handles.data_loaded.info_data_loaded(handles.data_loaded.info_data_loaded.Type == 'ROI',:) =[];
if isfield(handles.data_loaded, 'ROI')
    handles.data_loaded= rmfield(handles.data_loaded, 'ROI');
    handles.data_loaded= rmfield(handles.data_loaded, 'number_of_ROI');
end
for i = 1:numel(data_selected)
    
    fid_nii=fopen(fullfilename(handles, data_selected(i), '.nii'),'r');
    if fid_nii>0
        fclose(fid_nii);
        scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');
        
        %% read and load the nii file
        handles.data_loaded.ROI(i).V = spm_vol(char(fullfilename(handles, data_selected(i), '.nii')));
        handles.data_loaded.ROI(i).nii = read_volume(handles.data_loaded.ROI(i).V, handles.data_loaded.Scan(scan_of_reference).V);
        handles.data_loaded.number_of_ROI = size(handles.data_loaded.ROI,1);
        handles.data_loaded.info_data_loaded = [handles.data_loaded.info_data_loaded; handles.database(data_selected(i),:)];
    else
        warndlg('something wrong with the data. Nii of json file is missing','Warning');
        return
    end
    guidata(hObject, handles);
end




function handles = MIA_load_axes_single(hObject, ~, handles)

data_selected = get_data_selected(handles);


if numel(data_selected) > 4  % select only the 4 first scan
    data_selected = data_selected(1:4);
end
for i = 1:numel(data_selected)
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
handles.data_loaded.number_of_scan = numel(data_selected);
handles.data_loaded.info_data_loaded = handles.database(data_selected,:);

%set popupmenu(s) (echo and expt for each axes) and clear Axes if needed
for i=1:4 %handles.data_loaded.number_of_scan
    stri = num2str(i);
    %  and clear Axes unused
    eval(['cla(handles.MIA_data' stri ');']);
    if i > handles.data_loaded.number_of_scan
        set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off', 'Value', 1);
        set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off', 'Value', 1);
        
    else
        mat_size = handles.data_loaded.Scan(i).V(1).private.dat.dim;
        eval(['set(handles.MIA_data' stri '_title, ''String'', char(handles.database.SequenceName(data_selected(i))));']);
        
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
        set(handles.MIA_set_clip, 'Position', [0.4275 0.0042 0.0719 0.025]);
        
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
        set(handles.MIA_set_clip, 'Position', [0.4275 0.385 0.0719 0.025]);
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
        set(handles.MIA_set_clip, 'Position', [0.4275 0.0042 0.0719 0.025]);
        
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
        set(handles.MIA_set_clip, 'Position', [0.4275 0.0042 0.0719 0.025]);
end

% set MIA_slider_slice
if  length(handles.data_loaded.Scan(1).V(1).private.dat.dim) == 2  %handles.data_loaded.Scan(1).nii
    set(handles.MIA_slider_slice,'Visible', 'off', 'Value', 1);
    set(handles.MIA_slider_slice,'Max', 1);
    set(handles.MIA_slider_slice,'Value',1);
else
    set(handles.MIA_slider_slice,'Visible', 'on');
    set(handles.MIA_slider_slice,'Min',1);
    set(handles.MIA_slider_slice, 'Max', handles.data_loaded.Scan(1).V(1).private.dat.dim(3) );
    set(handles.MIA_slider_slice,'Value',1);
    set(handles.MIA_slider_slice,'SliderStep',[1/(handles.data_loaded.Scan(1).V(1).private.dat.dim(3) -1) min(5/(handles.data_loaded.Scan(1).V(1).private.dat.dim(3) -1),1)]);
end

% update MIA_table_pixel_values header
% guidata(hObject, handles);
%col_header(:,1) = {'','','','',''};
% if size(handles.database(patient).day(time_point).parameters(scan),2)>=4
%     col_header(2:5,1) = handles.database(patient).day(time_point).parameters(scan(1:4))';
% else
%col_header(2:size(handles.database(patient).day(time_point).parameters(scan),2)+1,1) = handles.database(patient).day(time_point).parameters(scan)';
% end
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
guidata(hObject, handles);

% Update table1 with the scan loaded and the clip values
% for i=1:handles.data_loaded(1).number_of_scan
%     scan_name= strcmp(handles.clips(:,1)', handles.database(patient).day(time_point).parameters(scan(i)));
%     scan_name = find(scan_name ==1);
%     if ~isempty(scan_name)
%         if numel(scan_name) == 1
%             handles.table1.clips(i,:) = handles.clips(scan_name(1),:);
%         else
%             handles.clips(scan_name(2:end),:) = [];
%             handles.table1.clips(i,:) = handles.clips(scan_name(1),:);
%         end
%     else
%         handles.clips(size(handles.clips,1)+1,1)= handles.database(patient).day(time_point).parameters(scan(i));
%         handles.clips(size(handles.clips,1),2)= {handles.data_selected(i).image.reco.globalmin};
%         handles.clips(size(handles.clips,1),3)= {handles.data_selected(i).image.reco.globalmax};
%         handles.table1.clips(i,:) = handles.clips(size(handles.clips,1),:);
%     end
%
%
% end
% if handles.data_loaded(1).number_of_scan < 4
%     for i = handles.data_loaded(1).number_of_scan+1:4
%         handles.table1.clips{i,1} = {'NaN'};
%         handles.table1.clips(i,2:3) = {0 0};
%     end
% end

% update MIA_table_pixel_values header
% guidata(hObject, handles);
% col_header(:,1) = {'','','','',''};
% if size(handles.database(patient).day(time_point).parameters(scan),2)>=4
%     col_header(2:5,1) = handles.database(patient).day(time_point).parameters(scan(1:4))';
% else
%     col_header(2:size(handles.database(patient).day(time_point).parameters(scan),2)+1,1) = handles.database(patient).day(time_point).parameters(scan)';
% end
% set(handles.MIA_table_pixel_values, 'ColumnName', col_header);
% set(handles.MIA_table1, 'ColumnName', col_header);

% % reset MIA_plot1
% set(get(handles.MIA_plot1, 'XLabel'), 'String', '');
% set(get(handles.MIA_plot1, 'YLabel'), 'String', '');
% set(get(handles.MIA_plot1, 'ZLabel'), 'String', '');
% if ~isempty(findobj('Tag', 'Colorbar'))
%     cbfreeze('del');
% end
% guidata(hObject, handles);


function handles = MIA_load_axes_PRM(hObject, ~, handles)
% PRM mode i.e. need to open the one parameter (diffusion
% or perfusion or...) for every time point

data_selected = get_data_selected(handles);
if numel(data_selected) ~= 1
    warndlg('In PRM mode you can open only on scan!!', 'Warning');
    return
end
% find indice of the same scan name across each time point to the selected
% patient
data_to_load = find(handles.database.Patient == handles.database.Patient(data_selected) &...
    handles.database.SequenceName == handles.database.SequenceName(data_selected));

if data_to_load <2
    warndlg(strcat({'Need more than one '}, scan_name, ' scan to run the PRM mode') ,'Warning');
    return
end

for i = 1:numel(data_to_load)
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

set(handles.MIA_patient_information_title, 'String', [char(unique(handles.database.Patient(data_to_load))) '_' char(unique(handles.database.Tp(data_selected)))]);
set(handles.MIA_orientation_space_popupmenu, 'String',  char(unique(handles.database.SequenceName(data_to_load),'stable')), 'Value', 1);
handles.data_loaded.number_of_scan = numel(data_to_load);
handles.data_loaded.info_data_loaded = handles.database(data_to_load,:);

for i=1:2 
    stri = num2str(i);
    %  and clear Axes unused
    eval(['cla(handles.MIA_data' stri ');']);
    if i > handles.data_loaded.number_of_scan
        set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off', 'Value', 1);
        set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off', 'Value', 1);
        
    else
        mat_size = handles.data_loaded.Scan(i).V(1).private.dat.dim;
        eval(['set(handles.MIA_data' stri '_title, ''String'', char(handles.database.SequenceName(data_to_load(i))));']);
        
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

% 
% 
% % 
% % scan_selected = 0;
% % scan_name = handles.database(patient).day(time_point).parameters(scan);
% % for ii = 1:numel(handles.database(patient).day)
% %     match = find(strcmp(scan_name, handles.database(patient).day(ii).parameters')==1);
% %     if numel(match) > 1
% %         % the warning message is coded in the MIA_load_axes_Callback
% %         warndlg(['There are more than 1 ' scan_name{:} ' for day ' handles.database(patient).day(ii).date],'Warning');
% %         return
% %     end
% %     if numel(match) == 1
% %         %         tmp_tp = find(match==1);
% %         fid=fopen(fullfile(handles.database(patient).path,handles.database(patient).day(ii).scans_file{match}) ,'r');
% %         if fid>0
% %             fclose(fid);
% %             new = load(fullfile(handles.database(patient).path,handles.database(patient).day(ii).scans_file{match}));
% %             clear tmp;
% %         else
% %             warndlg('No data file in this folder','Warning');
% %             return
% %         end
% %         scan_selected = scan_selected+1;
% %         if scan_selected == 1
% %             handles.data_selected_for_PRM = new.uvascim;
% %             handles.data_selected_for_PRM.image.day = handles.database(patient).day(ii).date;
% %             
% %         else
% %             tmp = new.uvascim.image;
% %             tmp.day = handles.database(patient).day(ii).date;
% %             handles.data_selected_for_PRM.image(scan_selected) = tmp;
% %             
% %             clear tmp
% %         end
% %     end
% % end
% % if scan_selected <2
% %     warndlg(strcat({'Need more than one '}, scan_name, ' scan to run the PRM mode') ,'Warning');
% %     return
% % end
% % handles.data_selected_for_PRM.scan_number = numel(handles.data_selected_for_PRM.image);
% % 
% % % recover information for each scan selected (x, y , echo, slice, expt, z_min, z_max)
% % data_selected_for_PRM_info(1:numel(handles.data_selected_for_PRM.image),1:5) = 1;
% % handles.data_selected_for_PRM.list_scan(1:4) = {''};
% % handles.data_selected_for_PRM.list_scan(1:numel(handles.data_selected_for_PRM.image)) = handles.database(patient).day(time_point).parameters(scan);
% % z_offset = [];
% % 
% % for i=1:numel(handles.data_selected_for_PRM.image)
% %     data_selected_for_PRM_info(i,1:size(size(handles.data_selected_for_PRM.image(i).reco.data),2))=size(handles.data_selected_for_PRM.image(i).reco.data);
% %     tmp_z_min = min(handles.data_selected_for_PRM.image(i).reco.fov_offsets(3,1,:));
% %     tmp_z_max = max(handles.data_selected_for_PRM.image(i).reco.fov_offsets(3,1,:));
% %     data_selected_for_PRM_info(i,6:7) = [tmp_z_min tmp_z_max];
% %     data_selected_for_PRM_info(i,8) = handles.data_selected_for_PRM.image(i).reco.thickness;
% %     data_selected_for_PRM.scan_z_offset(i).values = data_selected_for_PRM_info(i,6):data_selected_for_PRM_info(i,8):data_selected_for_PRM_info(i,7);
% %     data_selected_for_PRM.scan_z_offset(i).values = round(data_selected_for_PRM.scan_z_offset(i).values*100)/100;
% %     if numel(data_selected_for_PRM.scan_z_offset(i).values) ~= size(handles.data_selected_for_PRM.image(i).reco.data,4)
% %         data_selected_for_PRM.scan_z_offset(i).values =  round(squeeze(handles.data_selected_for_PRM.image(i).reco.fov_offsets(3,1,:))'*100)/100;
% %     end
% %     z_offset = [z_offset data_selected_for_PRM.scan_z_offset(i).values]; %#ok<AGROW>
% % end
% % 
% % z_offset = sort(unique(round(z_offset*100)/100));
% % slice_nbr = numel(z_offset);
% % 
% % if slice_nbr > max(data_selected_for_PRM_info(:,4))
% %     data_selected_for_PRM_info(:,4) =  slice_nbr;
% % else
% %     data_selected_for_PRM_info(:,4) =max(data_selected_for_PRM_info(:,4));
% % end
% % 
% % data_selected_for_PRM_info(:,3) = max(data_selected_for_PRM_info(:,3));
% % data_selected_for_PRM_info(:,5) = max(data_selected_for_PRM_info(:,5));
% % 
% % handles.data_selected_for_PRM.scan_info = data_selected_for_PRM_info;
% % handles.data_selected_for_PRM.z_offset = z_offset;
% % handles.data_selected_for_PRM.scan_z_offset = data_selected_for_PRM.scan_z_offset;
% % 
% % % save patient information
% % handles.data_selected_for_PRM.patient_info.name_nbr = patient;
% % handles.data_selected_for_PRM.patient_info.name = handles.database(patient).name;
% % handles.data_selected_for_PRM.patient_info.timepoint_nbr = time_point;
% % handles.data_selected_for_PRM.patient_info.timepoint = handles.database(patient).day(time_point).date;
% % handles.data_selected_for_PRM.patient_info.group =  handles.database(patient).group;
% % set(handles.MIA_patient_information_title, 'String', [handles.data_selected_for_PRM.patient_info.name '_' handles.data_selected_for_PRM.patient_info.timepoint]);
% % 
% % % resize data
% % res_option = get(handles.MIA_resolution_popupmenu, 'Value');
% % handles.resolution_selected = handles.resolution(res_option);
% % % "Original" resolution selected
% % if  handles.resolution_selected == 1
% %     handles.resolution_selected = max(handles.data_selected_for_PRM.scan_info(:,1));
% % else
% %     handles.resolution_selected = handles.resolution(res_option);
% % end
% % handles.data_selected_for_PRM_resized =  handles.data_selected_for_PRM;
% % for ii=1:handles.data_selected_for_PRM.scan_number
% %     imgsize = data_selected_for_PRM_info(ii,:);
% %     factor = [handles.resolution_selected/handles.data_selected_for_PRM.scan_info(ii,1) handles.resolution_selected/handles.data_selected_for_PRM.scan_info(ii,2)];
% %     data = zeros([ceil(factor(1)*imgsize(1)),ceil(factor(2)*imgsize(2)),imgsize(3:5)]);
% %     
% %     for a=1:size(handles.data_selected_for_PRM.image(ii).reco.data,3)
% %         for b=1:size(handles.data_selected_for_PRM.image(ii).reco.data,4)
% %             for c=1:size(handles.data_selected_for_PRM.image(ii).reco.data,5)
% %                 data(:,:,a,b,c) = imresize(handles.data_selected_for_PRM.image(ii).reco.data(:,:,a,b,c),[size(data, 1) size(data,2)],'bilinear');
% %             end
% %         end
% %     end
% %     handles.data_selected_for_PRM_resized.image(ii).reco.data =  data;
% % end
% % 
% % % recover information for each scan selected and resized (x, y , echo, slice, expt, z_min, z_max)
% % data_selected_for_PRM_resized_info(1:numel(handles.data_selected_for_PRM_resized.image),1:7) = 1;
% % for i=1:handles.data_selected_for_PRM.scan_number
% %     data_selected_for_PRM_resized_info(i,1:size(size(handles.data_selected_for_PRM_resized.image(i).reco.data),2))=size(handles.data_selected_for_PRM_resized.image(i).reco.data);
% %     data_selected_for_PRM_resized_info(i,6:8) = data_selected_for_PRM_info(i, 6:8);
% %     
% % end
% % handles.data_selected_for_PRM_resized.scan_info = data_selected_for_PRM_resized_info;
% % 
% % % Align slices
% % if handles.data_selected_for_PRM.scan_number > 1
% %     thickness_sorted = sort(data_selected_for_PRM_resized_info(:,8));
% %     slice_to_remove = [];
% %     for i=1:numel(unique(thickness_sorted))
% %         match_scan = abs(data_selected_for_PRM_resized_info(:,8) - thickness_sorted(i)) < 0.00001;
% %         match_scan = find(match_scan== 1);
% %         for j = 1:numel(match_scan)
% %             scan_nbr = match_scan(j);
% %             matrix_size = size(handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data);
% %             if size(matrix_size,2) == 2
% %                 matrix_size(3) = 1;
% %             end
% %             matrix_size(4) = numel(handles.data_selected_for_PRM_resized.z_offset);
% %             tmp=NaN(matrix_size);
% %             % first data (smaller slice thickness)
% %             if i == 1 && j == 1
% %                 for jj = 1:numel(data_selected_for_PRM.scan_z_offset(scan_nbr).values)
% %                     match =  abs(handles.data_selected_for_PRM_resized.z_offset - data_selected_for_PRM.scan_z_offset(scan_nbr).values(jj)) <  1e-3;
% %                     match_1 = find(match== 1);
% %                     tmp(:,:,:,match_1,:,:) = handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data(:,:,:,jj,:,:); %#ok<FNDSB>
% %                 end
% %             else
% %                 old_match = 0;
% %                 for jj = 1:size(handles.data_selected_for_PRM.image(scan_nbr).reco.data, 4)
% %                     match =  abs(handles.data_selected_for_PRM_resized.z_offset - data_selected_for_PRM.scan_z_offset(scan_nbr).values(jj)) <  1e-3;
% %                     match_1 = find(match== 1);
% %                     if match_1 ~= old_match+1
% %                         for x = (old_match+1):match_1-1
% %                             if (data_selected_for_PRM.scan_z_offset(scan_nbr).values(jj)-data_selected_for_PRM_resized_info(scan_nbr,8)/2 < handles.data_selected_for_PRM_resized.z_offset(x) - min(data_selected_for_PRM_resized_info(:,8))/2) ||...
% %                                     abs((data_selected_for_PRM.scan_z_offset(scan_nbr).values(jj)-data_selected_for_PRM_resized_info(scan_nbr,8)/2)- (handles.data_selected_for_PRM_resized.z_offset(x) - min(data_selected_for_PRM_resized_info(:,8))/2)) <  1e-3
% %                                 tmp(:,:,:,x,:,:) = handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data(:,:,:,jj,:,:);
% %                                 slice_to_remove = [slice_to_remove match_1];
% %                             end
% %                         end
% %                         tmp(:,:,:,match_1,:,:) = handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data(:,:,:,jj,:,:);
% %                     else
% %                         tmp(:,:,:,match_1,:,:) = handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data(:,:,:,jj,:,:);
% %                     end
% %                     for xx = match_1:numel(handles.data_selected_for_PRM_resized.z_offset)
% %                         if (data_selected_for_PRM.scan_z_offset(scan_nbr).values(jj)+data_selected_for_PRM_resized_info(scan_nbr,8)/2 > handles.data_selected_for_PRM_resized.z_offset(xx)+  min(data_selected_for_PRM_resized_info(:,8))/2) ||...
% %                                 abs((data_selected_for_PRM.scan_z_offset(scan_nbr).values(jj)+data_selected_for_PRM_resized_info(scan_nbr,8)/2) - (handles.data_selected_for_PRM_resized.z_offset(xx)+  min(data_selected_for_PRM_resized_info(:,8))/2)) <  1e-3
% %                             tmp(:,:,:,xx,:,:) = handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data(:,:,:,jj,:,:);
% %                         end
% %                     end
% %                     old_match = match_1;
% %                 end
% %             end
% %             handles.data_selected_for_PRM_resized.image(scan_nbr).reco.data = tmp;
% %         end
% %         
% %     end
% %     if ~isempty(slice_to_remove)
% %         slice_to_remove = unique(slice_to_remove);
% %         for i = 1 : numel(handles.data_selected_for_PRM_resized.image)
% %             handles.data_selected_for_PRM_resized.image(i).reco.data(:,:,:,slice_to_remove,:,:) = [];
% %         end
% %         handles.data_selected_for_PRM_resized.z_offset(slice_to_remove) = [];
% %         handles.data_selected_for_PRM_resized.scan_info(:,4) = handles.data_selected_for_PRM_resized.scan_info(1,4) - numel(slice_to_remove);
% %     end
% % end
% % 
% 
% 
% 
% % %set popupmenu(s) (echo and expt for each axes) and clear Axes if needed
% % for i=1:2
% %     stri = num2str(i);
% %     if numel(handles.data_selected_for_PRM.image) > i-1 && size(handles.data_selected_for_PRM.image(i).reco.data, 3) > 1 %numel(scan)
% %         set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max', size(handles.data_selected_for_PRM.image(i).reco.data, 3),...
% %             'SliderStep',[1/(size(handles.data_selected_for_PRM.image(i).reco.data, 3)-1) min(5/(size(handles.data_selected_for_PRM.image(i).reco.data, 3)-1),1)]);
% %     else
% %         set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Visible', 'off');
% %         set(eval(['handles.MIA_data', stri, '_echo_slider']), 'Value', 1);
% %         %  and clear Axes unused
% %         eval(['cla(handles.MIA_data' stri ');']);
% %         eval(['set(handles.MIA_data' stri '_title, ''String'', '''');']);
% %     end
% %     if numel(handles.data_selected_for_PRM.image) > i-1 && size(handles.data_selected_for_PRM.image(i).reco.data, 5) > 1 %numel(scan)
% %         set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'on', 'Value', 1, 'Min', 1, 'Max', size(handles.data_selected_for_PRM.image(i).reco.data, 5),...
% %             'SliderStep',[1/(size(handles.data_selected_for_PRM.image(i).reco.data, 5)-1) min(5/(size(handles.data_selected_for_PRM.image(i).reco.data, 5)-1),1)]);
% %     else
% %         set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Visible', 'off');
% %         set(eval(['handles.MIA_data', stri, '_expt_slider']), 'Value', 1);
% %         
% %     end
% % end
% 
% % resize windows handles.MIA_data1
% set(handles.MIA_data1, 'Position', [0.0188 0.4529 0.2523 0.3140]);
% 
% 
% set(handles.MIA_data2, 'Visible', 'on');
% set(handles.MIA_data2_title, 'Visible', 'on');
% set(handles.MIA_data3, 'Visible', 'off');
% set(handles.MIA_data3_title, 'Visible', 'off');
% set(handles.MIA_data4, 'Visible', 'off');
% set(handles.MIA_data4_title, 'Visible', 'off');
% 
% 
% % set MIA_slider_slice
% % set(handles.MIA_slider_slice, 'Max', data_selected_for_PRM_info(1,4));
% % set(handles.MIA_slider_slice,'Value',1);
% % set(handles.MIA_slider_slice,'Min',1);
% % if data_selected_for_PRM_info(1,4) == 1
% %     set(handles.MIA_slider_slice,'Visible', 'off');
% % else
% %     set(handles.MIA_slider_slice,'Visible', 'on');
% %     set(handles.MIA_slider_slice,'SliderStep',[1/(data_selected_for_PRM_info(1,4)-1) min(5/(data_selected_for_PRM_info(1,4)-1),1)]);
% % end
% 
% % % save clips info of data loaded
% % result = strcmp(handles.clips(:,1), handles.database(patient).day(time_point).parameters(scan(1)));
% % handles.table1.clips{1,1} =  handles.clips(find(result, 1),1);
% % handles.table1.clips(1,2:3) =  handles.clips(find(result, 1),2:3);
% % 
% % % update MIA_table_pixel_values header
% % guidata(hObject, handles);
% % col_header(:,1) = {'','','','',''};
% % if size(handles.database(patient).day(time_point).parameters(scan),2)>=4
% %     col_header(2:5,1) = handles.database(patient).day(time_point).parameters(scan(1:4))';
% % else
% %     col_header(2,1) = strcat(handles.data_selected_for_PRM_resized.image(1).day, '-', handles.database(patient).day(time_point).parameters(scan)');
% %     col_header(3,1) = strcat(handles.data_selected_for_PRM_resized.image(2).day, '-', handles.database(patient).day(time_point).parameters(scan)');
% %     
% % end
% % 
% % set(handles.MIA_table_pixel_values, 'ColumnName', col_header);
% % set(handles.MIA_table1, 'ColumnName', col_header);
% % 
% % % reset MIA_plot1
% % set(get(handles.MIA_plot1, 'XLabel'), 'String', '');
% % set(get(handles.MIA_plot1, 'YLabel'), 'String', '');
% % set(get(handles.MIA_plot1, 'ZLabel'), 'String', '');
% % if ~isempty(findobj('Tag', 'Colorbar'))
% %     cbfreeze('del');
% % end
% 
% guidata(hObject, handles);


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
    end
    
    % update the ROI matrix (new ROI, resized...)
    if isfield(handles.data_loaded, 'ROI')
        handles = MIA_update_VOI_displayed(hObject, eventdata, handles);
    end
    %update MIA_plot1
    if isfield(handles,'ROI_selected_resized')
        handles = MIA_find_VOI_coordonates(hObject,handles);
        if handles.display_option.view_plot == 1
            if ~isempty(handles.data_ploted.coordonates)
                handles =MIA_update_plot1_single(hObject,handles);
            end
        end
    else
        if ~isempty(get(handles.MIA_plot1, 'Children'))
            delete(get(handles.MIA_plot1, 'Children'));
            legend(handles.MIA_plot1,'off');
            hold(handles.MIA_plot1, 'off');
            set(handles.MIA_plot1, 'XTick', []);
            set(handles.MIA_plot1, 'YTick', []);
        end
        if isfield(handles, 'data_ploted')
            handles = rmfield(handles, 'data_ploted');
        end
    end
    % Update the VOI_cluster matrix (new cluster, resized...)
    if isfield(handles, 'ROI_cluster_resized')
        handles = MIA_update_VOI_cluster_displayed(hObject,handles);
    end
end
slice_nbr = get(handles.MIA_slider_slice, 'Value');
% is zommed?
if get(handles.MIA_data1, 'Children')~=0
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


VOI_values = {};
Slice_values= {};

% display every data available (image, ROI, cluster...)

if isfield(handles, 'data_displayed')
    number_of_data_to_displayed = numel(fieldnames(handles.data_displayed));
    axe = zeros(1,size( handles.data_displayed.image,4));
    for i=1:size(handles.data_displayed.image,4)
        
        stri = num2str(i);
        % store current contrast
        current_contrast = get(handles.(sprintf('MIA_data%d', i)), 'Clim');
        if     ~isempty(get(handles.(sprintf('MIA_data%d', i)), 'Children'));
            delete(get(handles.(sprintf('MIA_data%d', i)), 'Children'));
        end
        switch number_of_data_to_displayed
            case 1 % image only   
                % Clip image displayed at the 2 extremum (min, max)
                image_to_display =squeeze(handles.data_displayed.image(:,:,slice_nbr,i,1));
                image_to_display(image_to_display<prctile(image_to_display(:),1)) = prctile(image_to_display(:),1);
                image_to_display(image_to_display>prctile(image_to_display(:),99)) = prctile(image_to_display(:),99);
                image(image_to_display,'CDataMapping','Scaled','Parent', handles.(sprintf('MIA_data%d', i)),'Tag',sprintf('data%d', i));
                
                % apply the colormap selected
                colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
                eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
                
                if isZoomed == true
                    set(handles.(sprintf('MIA_data%d', i)),'XLim', XLim_zoomed);
                    set(handles.(sprintf('MIA_data%d', i)),'YLim', YLim_zoomed);
                    setappdata(handles.(sprintf('MIA_data%d', i)),'matlab_graphics_resetplotview', origInfo);
                else
                    set(handles.(sprintf('MIA_data%d', i)), 'XLim', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 2)]);
                    set(handles.(sprintf('MIA_data%d', i)), 'YLim', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 1)]);
                end
                set(handles.(sprintf('MIA_data%d', i)), 'Visible', 'on', 'XTick' , [], 'YTick', []);
            case 2 % image + ROI
                 % Clip image displayed at the 2 extremum (min, max)
                image_to_display =squeeze(handles.data_displayed.image(:,:,slice_nbr,i,1));
                image_to_display(image_to_display<prctile(image_to_display(:),0.1)) = prctile(image_to_display(:),0.1);
                image_to_display(image_to_display>prctile(image_to_display(:),99.9)) = prctile(image_to_display(:),99.9);
                image(image_to_display,'CDataMapping','Scaled','Parent', handles.(sprintf('MIA_data%d', i)),'Tag',sprintf('data%d', i));
                
                %image(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,1)),'CDataMapping','Scaled','Parent', handles.(sprintf('MIA_data%d', i)),'Tag','data1');

                colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
                eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
                
                if isZoomed == true
                    set(handles.(sprintf('MIA_data%d', i)),'XLim', XLim_zoomed);
                    set(handles.(sprintf('MIA_data%d', i)),'YLim', YLim_zoomed);
                    setappdata(handles.(sprintf('MIA_data%d', i)),'matlab_graphics_resetplotview', origInfo);
                    
                else
                    set(handles.(sprintf('MIA_data%d', i)), 'XLim', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 2)]);
                    set(handles.(sprintf('MIA_data%d', i)), 'YLim', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 1)]);
                end
                hold(handles.(sprintf('MIA_data%d', i)), 'on');
                if strcmp(get(handles.MIA_menu_roi_fill, 'Checked'), 'on')
                    fillroi = true;
                    trans = get(handles.MIA_PRM_slider_trans, 'Value')/100;
                else
                    fillroi = false;
                end
                % ROI if on the slice
                ROI_indices = find(handles.data_loaded.info_data_loaded.Type == 'ROI');
                for x = 1:numel(handles.data_loaded.ROI)
                    if handles.data_displayed.ROI.on_slice(x,slice_nbr) == 1
                        roi_a_appliquer=handles.data_loaded.ROI(x).nii(:,:,slice_nbr);
                        
                        if fillroi
                            roiRGB = repmat(roi_a_appliquer,[1 1 3]) .* permute(repmat(rgb(handles.colors{x}),[size(roi_a_appliquer,1) 1 size(roi_a_appliquer,1)]),[1 3 2]);
                            image(roiRGB,'AlphaData',roi_a_appliquer*trans,'CDataMapping','Scaled','Parent', handles.(sprintf('MIA_data%d', i)),'Tag',sprintf('data%d', i));
                        else
                            contour(handles.(sprintf('MIA_data%d', i)), roi_a_appliquer, 1, 'Color',rgb(handles.colors{x}),...
                                'Visible', 'on',...
                                'tag','ROI_contour');
                        end
                        % calculate the mean value inside the ROI over the
                        % Slice (mean / SD)
                        tmp_slice = squeeze(handles.data_displayed.image(:,:,slice_nbr,i,1));
                        
                        mean_slice = nanmean(tmp_slice(logical(roi_a_appliquer)));
                        sd_slice =  nanstd(tmp_slice(logical(roi_a_appliquer)));
                    else
                        mean_slice = nan;
                        sd_slice =  nan;
                    end
                    % calculate the mean value inside the ROI over the
                    % volume (mean / SD)
                    data_3D = squeeze(handles.data_displayed.image(:,:,:,i)); %(:,:,:,get(handles.(sprintf('MIA_data%d_echo_slider', i)), 'Value'),get(handles.(sprintf('MIA_data%d_expt_slider', i)), 'Value')));
                    mean_VOI = nanmean(data_3D(logical(handles.data_loaded.ROI(x).nii)));
                    sd_VOI =  nanstd(data_3D(logical(handles.data_loaded.ROI(x).nii)));
                    
                   
                    
                    clear data_slice data_3D
                    
                    if i == 1
                        Slice_values(x*2-1,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(x))) '-slice-mean']} ;
                        Slice_values(x*2-1,i+1) = num2cell(mean_slice);
                        Slice_values(x*2,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(x))) '-slice-SD']} ;
                        Slice_values(x*2,i+1) = num2cell(sd_slice);
                        
                        VOI_values(x*2-1,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(x))) '-VOI-mean']} ;
                        VOI_values(x*2-1,i+1) = num2cell(mean_VOI);
                        VOI_values(x*2,1) = {[char(handles.data_loaded.info_data_loaded.SequenceName(ROI_indices(x))) '-VOI-SD']} ;
                        VOI_values(x*2,i+1) = num2cell(sd_VOI);
                    else
                        Slice_values(x*2-1,i+1) = num2cell(mean_slice);
                        Slice_values(x*2,i+1) = num2cell(sd_slice);
                        VOI_values(x*2-1,i+1) = num2cell(mean_VOI);
                        VOI_values(x*2,i+1) = num2cell(sd_VOI);
                    end
                end
                set(handles.(sprintf('MIA_data%d', i)), 'Visible', 'on', 'XTick' , [], 'YTick', []);
                hold(handles.(sprintf('MIA_data%d', i)), 'off');
            case 3 % image + ROI + cluster
                
                %On charge le fichier contenant : "uvascroi" : toutes les
                %donn?es du clustering (classement des pixels),
                %"Informations" : Les informations du clustering (Nombre de
                %classes, patients consid?r?s, cartes utilis?es),
                %"Statistiques" : Les statistiques calcul?es ? partir du
                %classement des pixels, par tranche, volume, etc..
                load(handles.data_selected.roifile{:},'uvascroi','Informations','Statistiques');
                
                
                % image
                if  sum(strcmp(handles.data_selected.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected.image(i).reco.color_map, 'Color')
                    eval(['image(uint8(squeeze(handles.data_selected_resized.image(i).reco.data(:,:,1,slice_nbr,:))),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
                elseif sum(strcmp(handles.data_selected.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected.image(i).reco.color_map, 'Gray')
                    eval(['image(uint8(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:))),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
                    
                else
                    
                    if( max(~isreal(handles.data_displayed.image(:)))==1 )    % complex data
                        eval(['image(squeeze(abs(handles.data_displayed.image(:,:,slice_nbr,i,1))),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
                    else    % real data
                        eval(['image(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,1)),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
                    end
                end
                colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
                eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
                eval(['set(handles.MIA_data' stri ', ''Visible'', ''on'', ''XTick'' , [], ''YTick'', []);']);
                if isZoomed == true
                    eval(['set(handles.MIA_data' stri ',''XLim'', XLim_zoomed);']);
                    eval(['set(handles.MIA_data' stri ',''YLim'', YLim_zoomed);']);
                    eval(['setappdata(handles.MIA_data' stri ',''matlab_graphics_resetplotview'', origInfo);']);
                    
                else
                    eval(['set(handles.MIA_data' stri ', ''XLim'', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 2)]);']);
                    eval(['set(handles.MIA_data' stri ', ''YLim'', [1,size(squeeze(handles.data_displayed.image(:,:,slice_nbr,i,:)), 1)]);']);
                end
                
                % cluster
                eval(['hold(handles.MIA_data' stri ', ''on'');']);
                eval(['image(squeeze(handles.data_displayed.VOI_cluster.data(:,:,slice_nbr,:)), ''CDataMapping'',''Scaled'', ''parent'', handles.MIA_data' stri ', ''AlphaData'',handles.data_displayed.VOI_cluster.trans(:,:,slice_nbr), ''Tag'', ''data' stri '_ROI_cluster'');']);    %Affichage du cluster selectionn??
                eval(['hold(handles.MIA_data' stri ', ''off'');']);
                
                
                
                % Faire le lien entre toutes les images et les images o?
                % apparait le cluster, moins nombreuses
                % PresenceClustToutesTranches est un vecteur de taille
                % ?gale au nombre de tranches. Si la i?me tranche permet d'afficher
                % une ROI, et donc un cluster, alors la i?me composante du
                % vecteur PresenceClustToutesTranches contient un 1. Elle
                % contient un 0 s'il n'y a rien ? afficher
                PresenceClustToutesTranches = zeros(1,size(handles.ROI_cluster_resized.map,3));
                for ccc = 1:size(handles.ROI_cluster_resized.map,3)
                    if sum(sum(handles.data_displayed.VOI_cluster.trans(:,:,ccc))) == 0
                        PresenceClustToutesTranches(ccc) = 0;
                    else
                        PresenceClustToutesTranches(ccc) = 1;
                    end
                end
                
                
                
                %On r?cup?re les bonnes statistiques (ie celles qui
                %correspondent au patient s?lectionn?)
                % En effet, le fichier "Statistiques" charg? pr?cedemment
                % contient les statistiques de tous les patients consid?r?s
                % lors du clustering.
                if exist('Informations', 'var')
                    Sign = Informations.Sign;
                    for j = 1:length(Sign)
                        if handles.data_selected_resized.patient_info.name == Sign(1,j)
                            if handles.data_selected_resized.patient_info.timepoint == Sign(2,j)
                                MoyCartesTranches = Statistiques(j).MoyCartesTranches;
                                ProbTranches = Statistiques(j).ProbTranches;
                                MoyCartesVolume = Statistiques(j).MoyCartesVolume;
                                ProbVolume = Statistiques(j).ProbVolume;
                            end
                        end
                    end
                    
                    %On affiche d?somais les statistiques de volume du patient
                    %selectionn?. On cr?e donc un tableau.
                    
                    %Nom des lignes, avec les couleurs correspondant aux clusters
                    % La fonction colorgen necessite de formuler les
                    % composantes des couleurs en hexad?cimal. On utilise pour
                    % cel? la fonction rgb2hex (r?cup?r?e sur le matlab file
                    % exchange) qui permet de transformer en hexad?cimal une
                    % couleur formul?e en rgb.
                    
                    StatVol = cell(size(MoyCartesVolume,1),size(MoyCartesVolume,2)+2);
                    for ii=1:size(MoyCartesVolume,1)/2
                        colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
                        StatVol{2*ii-1,1} = colergen(rgb2hex(Informations(1).Couleurs(ii,:)), strcat('C_',num2str(ii),'_Volume_Mean'));
                        StatVol{2*ii,1} = colergen(rgb2hex(Informations(1).Couleurs(ii,:)), strcat('C_',num2str(ii),'_Volume_SD'));
                        StatVol{2*ii-1,size(MoyCartesVolume,2)+2} = ProbVolume(ii);
                    end
                    
                    %Affichage des statistiques de volume
                    for iii=1:size(MoyCartesVolume,1)
                        for iiii = 1:size(MoyCartesVolume,2)
                            StatVol{iii,iiii+1} = MoyCartesVolume(iii,iiii);
                        end
                    end
                    
                    %Affichage des noms des colonnes
                    NomsCartes = cell(1,length(Informations.Cartes));
                    for kk = 1:length(Informations.Cartes)
                        NomsCartes{kk+1} = char(Informations.Cartes(kk));
                    end
                    NomsCartes{length(Informations.Cartes)+2} = char('Repartition %');
                    NomsCartes{1} = char('_____________________');
                    
                    %L'image affichee dispose t elle d'un cluster ?
                    % Ce veteur est peut-etre redondant. Peut-etre que le
                    % vecteur PresenceClustToutesTranches suffit. A v?rifier.
                    PresenceCluster = zeros(1,slice_nbr);
                    for ccc = 1:slice_nbr
                        if sum(sum(handles.data_displayed.VOI_cluster.trans(:,:,ccc))) == 0
                            PresenceCluster(ccc) = 0;
                        else
                            PresenceCluster(ccc) = 1;
                        end
                    end
                    
                    %Si pas de cluster, pas de statistiques de la tranche
                    if PresenceCluster(slice_nbr) == 0
                        StatTranche = [];
                    else
                        %Si il y a un cluster a afficher, on affiche les statistiques et les noms des
                        %lignes comme precedemment
                        StatTranche = cell(size(MoyCartesTranches,1),size(MoyCartesTranches,2)+2);
                        for ii=1:size(MoyCartesVolume,1)/2
                            colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
                            StatTranche{2*ii-1,1} = colergen(rgb2hex(Informations(1).Couleurs(ii,:)), strcat('C_',num2str(ii),'_Tranche_Mean'));
                            StatTranche{2*ii,1} = colergen(rgb2hex(Informations(1).Couleurs(ii,:)), strcat('C_',num2str(ii),'_Tranche_SD'));
                            StatTranche{2*ii-1,size(MoyCartesTranches,2)+2} = ProbTranches(ii,sum(PresenceCluster(1:slice_nbr)));
                        end
                        
                        for iii=1:size(MoyCartesTranches,1)
                            for iiii = 1:size(MoyCartesTranches,2)
                                StatTranche{iii,iiii+1} = MoyCartesTranches(iii,iiii,sum(PresenceCluster(1:slice_nbr)));
                            end
                        end
                    end
                    
                    %On laisse une ligne blanche entre les 2 tableaux, pour
                    %plus de clart??
                    Espace = cell(1,size(MoyCartesTranches,2)+2);
                end
                % ROI if on the slice
                eval(['hold(handles.MIA_data' stri ', ''on'');']);
                if strcmp(get(handles.MIA_menu_roi_fill, 'Checked'), 'on')
                    fillroi = true;
                    trans = get(handles.MIA_PRM_slider_trans, 'Value')/100;
                else
                    fillroi = false;
                end
                for x = 1:numel(handles.ROI_selected_resized)
                    if handles.data_displayed.ROI.on_slice(x,slice_nbr) ==1
                        roi_a_appliquer=handles.ROI_selected_resized(x).data.value(:,:,slice_nbr);
                        
                        if fillroi
                            roiRGB = repmat(roi_a_appliquer,[1 1 3]) .* permute(repmat(rgb(handles.colors{x}),[size(roi_a_appliquer,1) 1 size(roi_a_appliquer,1)]),[1 3 2]); %#ok<NASGU>
                            eval(['image(roiRGB,''AlphaData'',roi_a_appliquer*trans,''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri '2'');']);
                        else
                            warning('off')
                            contour(eval(['handles.MIA_data' stri]), roi_a_appliquer, 1, 'Color',rgb(handles.colors{x}),...     % Affichage contour (ROI) et modification du fond de l'image
                                'Visible', 'on',...
                                'tag','ROI_contour');
                            warning('on')
                        end
                    end
                end
                eval(['hold(handles.MIA_data' stri ', ''off'');']);
                
        end
        
        axe(i) = handles.(sprintf('MIA_data%d', i));
        %%%%%%%% activate clic on graph MIA_dataXX_ButtonDownFcn
        set(get(handles.(sprintf('MIA_data%d', i)), 'Children'), 'HitTest', 'off');
        set(handles.(sprintf('MIA_data%d', i)),'ButtonDownFcn', @MIA_clic_on_image);
        set(get(handles.(sprintf('MIA_data%d', i)), 'Children'), 'ButtonDownFcn', @MIA_clic_on_image);
        % update contrast using the values stored (if it is not a new scan
        % loaded)
        if handles.display_option.manual_contrast == 1
            set(handles.(sprintf('MIA_data%d', i)), 'Clim', current_contrast );
        end
        
    end
    
end
if number_of_data_to_displayed == 2 && ~isempty(VOI_values)
    set(handles.MIA_table1, 'Data', [Slice_values' VOI_values']');
    width= get(handles.MIA_table1,'ColumnWidth' );
    width{1} = max(cellfun(@length,Slice_values(:,1)))*6;
    set(handles.MIA_table1,'ColumnWidth', width)
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


%%
scan_of_reference = get(handles.MIA_orientation_space_popupmenu, 'Value');

switch get(hObject, 'Tag')
    case {'MIA_data1_echo_slider', 'MIA_data1_expt_slider'}
        
        data1_echo_nbr = get(handles.MIA_data1_echo_slider, 'Value');
        data1_expt_nbr = get(handles.MIA_data1_expt_slider, 'Value');
        if handles.mode == 1
            handles.data_displayed.image(:,:,:,1) = read_slice(handles.data_loaded.Scan(1).V, handles.data_loaded.Scan(scan_of_reference).V, data1_echo_nbr, data1_expt_nbr);
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
            handles.data_displayed.image(:,:,:,1) = read_slice(handles.data_loaded.Scan(scan_number).V, handles.data_loaded.Scan(scan_of_reference).V, data1_echo_nbr, data1_expt_nbr);    
        end
    case {'MIA_data2_echo_slider', 'MIA_data2_expt_slider'}
        
        data2_echo_nbr = get(handles.MIA_data2_echo_slider, 'Value');
        data2_expt_nbr = get(handles.MIA_data2_expt_slider, 'Value');
        if handles.mode == 1
            handles.data_displayed.image(:,:,:,2) = read_slice(handles.data_loaded.Scan(2).V, handles.data_loaded.Scan(scan_of_reference).V, data2_echo_nbr, data2_expt_nbr);
        else
            if  handles.data_loaded.number_of_scan ==2
                scan_number = 2;
            else
                scan_number = get(handles.MIA_PRM_slider_tp, 'Value');
            end
             handles.data_displayed.image(:,:,:,2) = read_slice(handles.data_loaded.Scan(scan_number).V, handles.data_loaded.Scan(scan_of_reference).V, data2_echo_nbr, data2_expt_nbr);            
        end
    case {'MIA_data3_echo_slider', 'MIA_data3_expt_slider'}
        data3_echo_nbr = get(handles.MIA_data3_echo_slider, 'Value');
        data3_expt_nbr = get(handles.MIA_data3_expt_slider, 'Value');
        handles.data_displayed.image(:,:,:,3) = read_slice(handles.data_loaded.Scan(3).V, handles.data_loaded.Scan(scan_of_reference).V, data3_echo_nbr, data3_expt_nbr);
    case {'MIA_data4_echo_slider', 'MIA_data4_expt_slider'}
        data4_echo_nbr = get(handles.MIA_data4_echo_slider, 'Value');
        data4_expt_nbr = get(handles.MIA_data4_expt_slider, 'Value');
        handles.data_displayed.image(:,:,:,4) = read_slice(handles.data_loaded.Scan(4).V, handles.data_loaded.Scan(scan_of_reference).V, data4_echo_nbr, data4_expt_nbr);
    case {'MIA_new_roi'}
        return
    otherwise
        if isfield(handles, 'data_displayed')
            handles = rmfield(handles, 'data_displayed');
        end
        if handles.mode == 1
            for i=1:handles.data_loaded.number_of_scan
                stri = num2str(i);
                eval(['data' stri '_echo_nbr = get(handles.MIA_data' stri '_echo_slider, ''Value'');']);
                eval(['data' stri '_expt_nbr = get(handles.MIA_data' stri '_expt_slider, ''Value'');']);
                
                eval(['ima' stri '= read_slice(handles.data_loaded.Scan(i).V, handles.data_loaded.Scan(scan_of_reference).V, data' stri '_echo_nbr, data' stri '_expt_nbr);']);
                %      eval(['ima' stri '= squeeze(handles.data_loaded.Scan(i).nii(:,:,:,round(data' stri '_echo_nbr),round(data' stri '_expt_nbr)));']);
                
                %     if sum(strcmp(handles.data_selected.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected.image(i).reco.color_map, 'gray')
                %         eval(['ima' stri '= squeeze(handles.data_selected_resized.image(i).reco.data(:,:,round(data' stri '_echo_nbr),:,:));']);
                %     else
                %         eval(['ima' stri '= squeeze(handles.data_selected_resized.image(i).reco.data(:,:,round(data' stri '_echo_nbr),:,round(data' stri '_expt_nbr)));']);
                
                %     end
                %     if isfield(handles.data_selected_resized.image(i), 'echo_cluster') && ~isempty(handles.data_selected_resized.image(i).reco.echo_cluster)
                %         eval(['title' stri ' = strcat(handles.data_selected.list_scan{i}, ''-'', handles.data_selected_resized.image(i).reco.echo_cluster{data' stri '_echo_nbr});']);
                %     else
                %         eval(['title' stri ' = handles.data_selected.list_scan{i};']);
                %     end
                %     eval(['set(handles.MIA_data' stri '_title, ''String'', {title' stri '});']);
                %eval(['ax(' stri ') = handles.MIA_data' stri ';']);
                
                %set clip
                %     eval(['scan_name = strcmp(handles.clips(:,1)'', title' stri ');']);
                %     tempclip = handles.clips(scan_name == 1,2:3);
                %     tempclip =cell2num(tempclip);  %#ok<NASGU>
                %     eval(['ima' stri '(ima' stri '>tempclip(1,2))=NaN;']);
                %     eval(['ima' stri '(ima' stri '<tempclip(1,1))=NaN;']);
                
                %     if sum(strcmp(handles.data_selected.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected.image(i).reco.color_map, 'Color')
                %         handles.data_displayed.image(:,:,:,i,:) = eval(['ima' num2str(i)]);
                %     else
                handles.data_displayed.image(:,:,:,i) = eval(['ima' num2str(i)]);
                %     end
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
                    if  handles.data_loaded.number_of_scan ==2
                        scan_number = 2;
                    else
                        scan_number = get(handles.MIA_PRM_slider_tp, 'Value');
                    end
                end
                
                
                stri = num2str(i);
                eval(['data' stri '_echo_nbr = get(handles.MIA_data' stri '_echo_slider, ''Value'');']);
                eval(['data' stri '_expt_nbr = get(handles.MIA_data' stri '_expt_slider, ''Value'');']);
                eval(['ima' stri '= read_slice(handles.data_loaded.Scan(scan_number).V, handles.data_loaded.Scan(scan_of_reference).V, data' stri '_echo_nbr, data' stri '_expt_nbr);']);
                
                handles.data_displayed.image(:,:,:,i) = eval(['ima' num2str(i)]);
            end
        end
        if ~strcmp(get(hObject, 'Tag'), 'MIA_PRM_slider_tp')
            if ~size(handles.data_displayed.image, 3) == 1
                set(handles.MIA_slider_slice,'Visible', 'off');
            else
                set(handles.MIA_slider_slice,'Visible', 'on');
                set(handles.MIA_slider_slice,'Value', 1);
                set(handles.MIA_slider_slice, 'max', size(handles.data_displayed.image, 3), 'Value', 1, 'SliderStep',[1/(size(handles.data_displayed.image, 3) -1) min(5/(size(handles.data_displayed.image, 3) -1),1)]);
            end
        end
end


guidata(hObject, handles);




% function MIA_update_axes_PRM(hObject, eventdata, handles)
% handles = guidata(hObject);
% if ~isfield(handles, 'data_selected_for_PRM')
%     return
% end
% slice_nbr = get(handles.MIA_slider_slice, 'Value');
% 
% 
% %save data displayed
% if ~strcmp(get(hObject, 'Tag'), 'MIA_slider_slice')
%     if isfield(handles, 'data_displayed')
%         handles = rmfield(handles, 'data_displayed');
%     end
%     % Update image_displayed matrix
%     if isfield(handles, 'data_selected_for_PRM')
%         handles = MIA_update_image_displayed_PRM(hObject, eventdata, handles);
%     end
%     % update the ROI matrix (new ROI, resized...)
%     if isfield(handles, 'ROI_selected_resized')
%         handles = MIA_update_VOI_displayed_PRM(hObject, eventdata, handles);
%     end
%     %update MIA_plot1
%     if isfield(handles,'ROI_selected_resized')
%         handles = MIA_find_VOI_coordonates(hObject,handles);
%         if ~isempty(handles.data_ploted.coordonates)
%             handles =MIA_update_plot1_PRM(hObject,handles);
%         end
%     else
%         if ~isempty(get(handles.MIA_plot1, 'Children'))
%             delete(get(handles.MIA_plot1, 'Children'));
%             legend(handles.MIA_plot1,'off');
%             hold(handles.MIA_plot1, 'off');
%             set(handles.MIA_plot1, 'XTick', []);
%             set(handles.MIA_plot1, 'YTick', []);
%         end
%         if isfield(handles, 'data_ploted')
%             handles = rmfield(handles, 'data_ploted');
%         end
%     end
%     % Update the PRM matrix (new cluster, resized...)
%     if isfield(handles, 'ROI_PRM_resized')
%         handles = MIA_update_PRM_Overlay_map(hObject,handles);
%     end
% end
% 
% % is zommed?
% if get(handles.MIA_data1, 'Children')~=0
%     origInfo = getappdata(handles.MIA_data1, 'matlab_graphics_resetplotview');
%     if isempty(origInfo)
%         isZoomed = false;
%     elseif isequal(get(handles.MIA_data1,'XLim'), origInfo.XLim) && ...
%             isequal(get(handles.MIA_data1,'YLim'), origInfo.YLim) %&& ...
%         isZoomed = false;
%     else
%         isZoomed = true;
%         XLim_zoomed = get(handles.MIA_data1,'XLim'); %#ok<NASGU>
%         YLim_zoomed = get(handles.MIA_data1,'YLim'); %#ok<NASGU>
%     end
% else
%     isZoomed = false;
% end
% % if the resolution is changed go back unzoomed!
% if strcmp(get(hObject, 'Tag'), 'MIA_resolution_popupmenu')
%     isZoomed = false;
% end
% 
% 
% % display every data available (image, ROI, cluster...)
% if isfield(handles, 'data_displayed')
%     number_of_data_to_displayed = numel(fieldnames(handles.data_displayed));
%     for i=1:2
%         stri = num2str(i);
%         if eval(['~isempty(get(handles.MIA_data' stri ', ''Children''))']);
%             eval(['delete(get(handles.MIA_data' stri ', ''Children''));']);
%         end
%         switch number_of_data_to_displayed
%             case 1 % image only
%                 if( max(~isreal(handles.data_displayed.image(:)))==1 )    % complex data
%                     eval(['image(squeeze(abs(handles.data_displayed.image(:,:,slice_nbr,i))),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
%                 else    % real data
%                     eval(['image(squeeze(handles.data_displayed.image(:,:,slice_nbr,i)),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
%                 end
%                 colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
%                 eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
%                 %                 eval(['colormap(handles.MIA_data' stri ', ''gray'');']);
%                 eval(['set(handles.MIA_data' stri ', ''Visible'', ''on'', ''XTick'' , [], ''YTick'', []);']);
%                 % if image zoomed applay zoom
%                 if isZoomed == true
%                     eval(['set(handles.MIA_data' stri ',''XLim'', XLim_zoomed);']);
%                     eval(['set(handles.MIA_data' stri ',''YLim'', YLim_zoomed);']);
%                 end
%             case 2 % image + ROI -- never used !!
%                 % image
%                 if( max(~isreal(handles.data_displayed.image(:)))==1 )    % complex data
%                     eval(['image(squeeze(abs(handles.data_displayed.image(:,:,slice_nbr,i))),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
%                 else    % real data
%                     eval(['image(squeeze(handles.data_displayed.image(:,:,slice_nbr,i)),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
%                 end
%                 colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
%                 eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
%                 %                 eval(['colormap(handles.MIA_data' stri ', ''gray'');']);
%                 eval(['hold(handles.MIA_data' stri ', ''on'');']);
%                 % ROI if on the slice
%                 for x = 1:numel(handles.ROI_selected_resized)
%                     if findn(handles.ROI_selected_resized(x).data.value(:,:,slice_nbr)==1)
%                         for j = 1:size(handles.data_displayed.ROI.data,1)
%                             roi_contour = handles.data_displayed.ROI.data{j,slice_nbr};
%                             line(roi_contour(1,2:end)', roi_contour(2,2:end)',...
%                                 'parent', eval(['handles.MIA_data' stri]),...
%                                 'Color',rgb(handles.colors{j}),...
%                                 'Visible', 'on',...
%                                 'tag','ROI_contour');
%                         end
%                     end
%                 end
%                 eval(['set(handles.MIA_data' stri ', ''Visible'', ''on'', ''XTick'' , [], ''YTick'', []);']);
%                 eval(['hold(handles.MIA_data' stri ', ''off'');']);
%             case 3 % image + ROI + cluster
%                 % image
%                 if( max(~isreal(handles.data_displayed.image(:)))==1 )    % complex data
%                     eval(['image(squeeze(abs(handles.data_displayed.image(:,:,slice_nbr,i))),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
%                 else    % real data
%                     eval(['image(squeeze(handles.data_displayed.image(:,:,slice_nbr,i)),''CDataMapping'',''Scaled'',''Parent'', handles.MIA_data' stri ',''Tag'',''data' stri ''');']);
%                 end
%                 colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
%                 eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
%                 %                 eval(['colormap(handles.MIA_data' stri ', ''gray'');']);
%                 eval(['set(handles.MIA_data' stri ', ''Visible'', ''on'', ''XTick'' , [], ''YTick'', []);']);
%                 
%                 % if image zoomed applay zoom
%                 if isZoomed == true
%                     eval(['set(handles.MIA_data' stri ',''XLim'', XLim_zoomed);']);
%                     eval(['set(handles.MIA_data' stri ',''YLim'', YLim_zoomed);']);
%                 end
%                 % ROI if on the slice
%                 % Not useful in PRM mode
%                 
%                 % cluster
%                 eval(['hold(handles.MIA_data' stri ', ''on'');']);
%                 eval(['image(squeeze(handles.data_displayed.PRM.data(:,:,slice_nbr,:)), ''CDataMapping'',''Scaled'', ''parent'', handles.MIA_data' stri ', ''AlphaData'',handles.data_displayed.PRM.trans(:,:,slice_nbr), ''Tag'', ''data' stri '_ROI_cluster'');']);
%                 eval(['hold(handles.MIA_data' stri ', ''off'');']);
%         end
%         
%         axe(i) = eval(['handles.MIA_data' stri]);
%         
%         %%%%%%%% activate clic on graph MIA_dataXX_ButtonDownFcn
%         eval(['set(get(handles.MIA_data' stri ', ''Children''), ''HitTest'', ''off'');']);
%         eval(['set(handles.MIA_data' stri ',''ButtonDownFcn'', @MIA_clic_on_image);']);
%         eval(['set(get(handles.MIA_data' stri ', ''Children''), ''ButtonDownFcn'', @MIA_clic_on_image);']);
%     end
% end
% linkaxes(axe, 'xy');
% guidata(hObject, handles);

% function handles = MIA_update_image_displayed_PRM(hObject, eventdata, handles)
% 
% % Display data
% for i=1:2
%     stri = num2str(i);
%     eval(['data' stri '_echo_nbr = get(handles.MIA_data' stri '_echo_slider, ''Value'');']);
%     eval(['data' stri '_expt_nbr = get(handles.MIA_data' stri '_expt_slider, ''Value'');']);
%     if i == 1 % select pre-scan
%         scan_number = get(handles.MIA_PRM_ref_popupmenu, 'Value');
%         if scan_number > 1
%             scan_number = scan_number-1;
%         else %scan number dynamic 1 prior post-time point
%             scan_number = get(handles.MIA_PRM_slider_tp, 'Value') -1;
%             if scan_number == 0  % case PRM_ref = -1 and slider = 1
%                 scan_number = 1;
%             end
%         end
%     else % select post-scan
%         if  handles.data_selected_for_PRM_resized.scan_number ==2
%             scan_number = 2;
%         else
%             scan_number = get(handles.MIA_PRM_slider_tp, 'Value');
%         end
%     end
%     if( max(~isreal(handles.data_selected_for_PRM_resized.image(:)))==1 )    % complex data
%         eval(['ima' stri '= squeeze(abs(handles.data_selected_for_PRM_resized.image(scan_number).reco.data(:,:,data' stri '_echo_nbr,:,data' stri '_expt_nbr)));']);
%     else    % real data
%         eval(['ima' stri '= squeeze(handles.data_selected_for_PRM_resized.image(scan_number).reco.data(:,:,data' stri '_echo_nbr,:,data' stri '_expt_nbr));']);
%     end
%     if isfield(handles.data_selected_for_PRM_resized.image(i), 'echo_cluster') && ~isempty(handles.data_selected_for_PRM_resized.image(i).reco.echo_cluster)
%         eval(['title' stri ' = strcat(handles.data_selected_for_PRM.list_scan{i}, ''-'', handles.data_selected_resized_for_PRM.image(i).reco.echo_cluster{data' stri '_echo_nbr});']);
%     else
%         eval(['title' stri ' = handles.data_selected_for_PRM.list_scan{i};']);
%     end
%     eval(['set(handles.MIA_data' stri '_title, ''String'', {title' stri '});']);
%     eval(['ax(' stri ') = handles.MIA_data' stri ';']);
%     
%     % Update clip
%     eval(['scan_name = strcmp(handles.clips(:,1)'', title' stri ');']);
%     tempclip = handles.clips(scan_name == 1,2:3);
%     tempclip =cell2num(tempclip);  %#ok<NASGU>
%     eval(['ima' stri '(ima' stri '>tempclip(1,2))=NaN;']);
%     eval(['ima' stri '(ima' stri '<tempclip(1,1))=NaN;']);
%     
%     handles.data_displayed.image(:,:,:,i) = eval(['ima' num2str(i)]);
% end
% 
% % linkaxes(ax, 'xy');
% 
% %set title with time point
% tp1 = get(handles.MIA_PRM_ref_popupmenu, 'Value')-1;
% if  handles.data_selected_for_PRM_resized.scan_number ==2
%     tp2 = 2;
% else
%     tp2 = get(handles.MIA_PRM_slider_tp, 'Value');
% end
% if tp1 == 0
%     tp1 = tp2 - 1;
%     if tp1 == 0
%         tp1 = 1;
%     end
% end
% set(handles.MIA_data1_title, 'String', strcat(handles.data_selected_for_PRM_resized.image(tp1).day, '-',get(handles.MIA_data1_title, 'String')));
% set(handles.MIA_data2_title, 'String', strcat(handles.data_selected_for_PRM_resized.image(tp2).day, '-',get(handles.MIA_data2_title, 'String')));



% % Display ROI
% if isfield(handles, 'ROI_selected_resized')
%     handles = MIA_update_VOI_displayed_PRM(hObject, eventdata, handles);
% end
%
% %save data displayed
% if ~strcmp(get(hObject, 'Tag'), 'MIA_slider_slice')
%     if isfield(handles, 'data_displayed')
%         handles = rmfield(handles, 'data_displayed');
%     end
%     handles.data_displayed.image(:,:,:,1) = eval(['ima' num2str(1)]);
%     handles.data_displayed.image(:,:,:,2) = eval(['ima' num2str(2)]);
%
%     %update MIA_plot1 or plot_PRM depending of the mode selected
%     guidata(hObject, handles);
%     if isfield(handles,'ROI_selected_resized')
%         handles = MIA_find_VOI_coordonates(hObject,handles);
%         MIA_update_plot1_PRM(hObject, handles)
%     else % reset MIA_plot1
%         if ~isempty(get(handles.MIA_plot1, 'Children'))
%             delete(get(handles.MIA_plot1, 'Children'));
%             legend(handles.MIA_plot1,'off');
%             hold(handles.MIA_plot1, 'off');
%         end
%         set(handles.MIA_plot1, 'XTick', []);
%         set(handles.MIA_plot1, 'YTick', []);
%         if isfield(handles, 'data_ploted')
%             handles = rmfield(handles, 'data_ploted');
%         end
%         guidata(hObject, handles);
%     end
% end
%
% % Display PRM colormap if mode PRM selected
% if isfield(handles,'ROI_selected_resized')
% MIA_update_PRM_Overlay_map(hObject,handles)
% end

function handles = MIA_update_VOI_displayed(hObject, eventdata, handles)
% slice_nbr = get(handles.MIA_slider_slice, 'Value');
for i = 1:numel(handles.data_loaded.ROI)
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

% function handles = MIA_update_VOI_displayed_PRM(hObject, eventdata, handles)
% 
% for i = 1:numel(handles.ROI_selected_resized)
%     for slice_nbr=1:get(handles.MIA_slider_slice, 'Max');
%         roi_a_appliquer=handles.ROI_selected_resized(i).data.value(:,:,slice_nbr);
%         roi_contour=contourc(double(roi_a_appliquer),1);
%         handles.data_displayed.ROI.data(i,slice_nbr) = {roi_contour};
%     end
% end

% for i = 1:numel(handles.ROI_selected_resized)
%     roi_a_appliquer=handles.ROI_selected_resized(i).data.value(:,:,slice_nbr);
%     roi_contour=contourc(double(roi_a_appliquer),1);
%     if size(roi_contour,2) ~=0
%         for j = 1:2
%             strj = num2str(j);
%             handles.ROI_selected_resized(i).data.hdl_contour= line(roi_contour(1,2:end)', roi_contour(2,2:end)',...
%                 'parent', eval(['handles.MIA_data' strj]),...
%                 'Color',rgb(handles.colors{i}),...
%                 'tag','ROI_contour');
%         end
%     else
%         handles.ROI_selected_resized(i).data.hdl_contour = [];
%     end
% end
% guidata(hObject, handles);

function handles = MIA_find_VOI_coordonates(hObject,handles)

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
        %         if strcmp(get(eval(['handles.MIA_data' stri '_title']), 'String'), handles.histo')
        if handles.mode == 1
            if sum(strcmp(handles.data_selected.list_scan{i}, handles.histo'))>0 && strcmp(handles.data_selected.image(i).reco.color_map, 'Color')
                %             for j = 1:size(handles.data_displayed.image, 3)
                %                 tmp(:,:,j,i) = rgb2gray(uint8(squeeze(handles.data_displayed.image(:, :,j,i,:))));
                %                 index = findn(tmp(:,:,j,i) ~= 0);
                %                 if ~isempty(index)
                %                     for jj = 1:size(index,1)
                %                         if tmp(index(jj,1),index(jj,2),j,i) == 255
                %                             tmp(index(jj,1),index(jj,2),j,i) = 0;
                %                         else
                %                             tmp(index(jj,1),index(jj,2),j,i) = 255 - tmp(index(jj,1),index(jj,2),j,i);
                %                         end
                %                     end
                %                 end
                %             end
                %             tmp(:,:,:,i) = tmp(:,:,:,i).*handles.ROI_selected_resized(ii).data.value;
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
                %             for j = 1:size(handles.data_displayed.image, 3)
                %                 tmp(:,:,j,i) = rgb2gray(uint8(squeeze(handles.data_displayed.image(:, :,j,i,:))));
                %                 index = findn(tmp(:,:,j,i) ~= 0);
                %                 if ~isempty(index)
                %                     for jj = 1:size(index,1)
                %                         if tmp(index(jj,1),index(jj,2),j,i) == 255
                %                             tmp(index(jj,1),index(jj,2),j,i) = 0;
                %                         else
                %                             tmp(index(jj,1),index(jj,2),j,i) = 255 - tmp(index(jj,1),index(jj,2),j,i);
                %                         end
                %                     end
                %                 end
                %             end
                %             tmp(:,:,:,i) = tmp(:,:,:,i).*handles.ROI_selected_resized(ii).data.value;
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
            %handles.data_displayed.image(coordonates(i,1),coordonates(i,2),coordonates(i,3),j);
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
% if handles.data_ploted.coordonates
%return
%end
coordonates = handles.data_ploted.coordonates;

for ii = 1:numel(handles.ROI_selected_resized)
    voi_empty = 0;
    strii=num2str(ii);
    index = findn(coordonates(:,4)== ii);
    if isempty(index)
        voi_empty = 1;
    else
        voi_data =coordonates(index(:,1),:);
        if ~isreal(voi_data) % complex data
            voi_data = abs(voi_data);
        end
        
    end
    if ~isempty(get(handles.MIA_plot1, 'Children'))
        hold(handles.MIA_plot1, 'on');
    end
    if voi_empty == 0
        switch handles.data_loaded(1).number_of_scan
            case 1
                nbin = numel(voi_data(:,5))/(numel(voi_data(:,5))/15);
                if ii > 1
                    hold(handles.MIA_plot1, 'on');
                end
                [f, xi] = histnorm(voi_data(:,5),nbin);
                plot(handles.MIA_plot1,xi,f,...
                    'Color',rgb(handles.colors{ii}),...
                    'Tag', strcat('MIA_plot1_1d', strii));
                clear f xi
            case 2
                scatter(handles.MIA_plot1, voi_data(:,5), voi_data(:,6),...
                    'filled',...
                    'SizeData', 20,...
                    'MarkerFaceColor',rgb(handles.colors{ii}),...
                    'MarkerEdgeColor',rgb(handles.colors{ii}),...
                    'Visible', 'on',...
                    'Tag', strcat('MIA_plot1_2d', strii));
                %                 uistack(findobj('Tag', strcat('MIA_plot1_2d', strii)), 'bottom');
                
            case 3
                scatter3(handles.MIA_plot1, voi_data(:,5), voi_data(:,6), voi_data(:,7),...
                    'filled',...
                    'SizeData', 20,...
                    'MarkerFaceColor',rgb(handles.colors{ii}),...
                    'MarkerEdgeColor',rgb(handles.colors{ii}),...
                    'Tag', strcat('MIA_plot1_3d', strii));
                %                 uistack(findobj('Tag', strcat('MIA_plot1_3d', strii)), 'bottom');
            case 4
                color4d = zeros(size(voi_data(:,8),1),3);
                tmp = jet(256);
                mini=min(voi_data(:,8));
                maxi=max(voi_data(:,8));
                for i = 1:size(voi_data(:,8),1)
                    color4d(i,:) = tmp(round((voi_data(i,8)-mini)/(maxi-mini)*255)+1,:);
                end
                cbfreeze('del');
                scatter3(handles.MIA_plot1, voi_data(:,5), voi_data(:,6), voi_data(:,7), 10, color4d, 'filled',...
                    'SizeData', 30,...
                    'Marker', handles.markers{ii},...
                    'MarkerEdgeColor', 'k',...
                    'Tag', strcat('MIA_plot1_4d', strii));
                %                 uistack(findobj('Tag', strcat('MIA_plot1_4d', strii)), 'bottom');
                if ~strcmp(get(gco, 'Tag'), 'speedy_run_button')
                    colormap(jet);
                    colorbar('peer', handles.MIA_plot1);
                    %% Due to matlab 2015 I had to comment this lines -BL-
                    %                     MIA_plot1_colorbar = cbfreeze( h ,'on');
                    %                     colormap(gray);
                    %                     set(IA_plot1_colorbar, 'YTickLabel', get(IA_plot1_colorbar, 'YTick').*max(voi_data(:,8)));
                    %                     set(get(IA_plot1_colorbar, 'YLabel'), 'String', get(handles.MIA_data4_title, 'String'), 'Rotation', -90, 'VerticalAlignment', 'bottom');
                end
        end
    end
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
legend(handles.MIA_plot1, {handles.ROI_selected_resized.name}', 'Location','NorthEast');
set(get(handles.MIA_plot1, 'XLabel'), 'String', get(handles.MIA_data1_title, 'String'));

set(get(handles.MIA_plot1, 'ZLabel'), 'String', get(handles.MIA_data3_title, 'String'));
hold(handles.MIA_plot1, 'off');



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

% this code worked until 09-2016.... but not anymore, weird?! BL
% GUI_size = get(handles.MIA_GUI, 'Position');
% currPt(1) = currPt(1)/GUI_size(3);
% currPt(2) = currPt(2)/GUI_size(4);
%

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

return
handles = guidata(hObject);

if ~strcmp(get(handles.MIA_GUI,'SelectionType'),'normal')
    %     G=get(handles.MIA_GUI,'userdata');
    G.initpnt=get(gca,'currentpoint');
    G.initClim = get(gca,'Clim');
    set(handles.MIA_GUI,'userdata',G);
    set(handles.MIA_GUI, 'WindowButtonMotionFcn',@MIA_AdjWL);
    
    return
end

if handles.mode == 1
    if size(handles.data_selected_resized.image(1).reco.data,3) == 1 && ...
            size(handles.data_selected_resized.image(1).reco.data, 5) == 1
        return
    end
else
    if size(handles.data_selected_for_PRM_resized.image(1).reco.data,3) == 1 && ...
            size(handles.data_selected_for_PRM_resized.image(1).reco.data, 5) == 1
        return
    end
end
slice_nbre = get(handles.MIA_slider_slice, 'Value');
tag = get(get(hObject, 'Children'), 'Tag');
match_roi = [] ;
if ~ischar(tag)
    for i = 1:numel(tag)
        if strfind(tag{i}, 'data')
            match_data = i;
        elseif strfind(tag{i}, 'ROI_contour')
            match_roi = [match_roi i];
        end
    end
    currPt_on_axe=eval(['get(handles.MIA_' tag{match_data} ',''CurrentPoint'');']);
else
    currPt_on_axe=eval(['get(handles.MIA_' tag ',''CurrentPoint'');']);
end


[pixel_coordinates_2d] = [round(currPt_on_axe(1,1)) round(currPt_on_axe(1,2)) round(currPt_on_axe(1,3))];

voxel = pixel_coordinates_2d(1:2);
if voxel(1) == 0 || voxel(2)==0
    return  %bug somewhere
end

if ~isempty(get(handles.MIA_plot1, 'Children'))
    delete(get(handles.MIA_plot1, 'Children'));
    legend(handles.MIA_plot1,'off');
    hold(handles.MIA_plot1, 'off');
end
legende_txt = [];
if handles.mode == 1
    % single mode
    if size(handles.data_selected_resized.image(1).reco.data,3) > 1 && ...
            size(handles.data_selected_resized.image(1).reco.data, 5) > 1
        
        % select repetition or echoes
    elseif size(handles.data_selected_resized.image(1).reco.data,3) > 1
        % plot echoes
        for i=1:handles.data_loaded(1).number_of_scan
            y_data = []; x_data = [];
            strii = num2str(i);
            if ~isempty(get(handles.MIA_plot1, 'Children'))
                hold(handles.MIA_plot1, 'on');
            end
            if match_roi > 0
                % plot ROI
                for x = 1:numel(match_roi)
                    roi = []; tmp = [];
                    tmp = squeeze(handles.data_selected_resized.image(i).reco.data(:, :, :, slice_nbre, 1));
                    roi = handles.ROI_selected_resized(x).data.value(:,:,slice_nbre);
                    roi(roi == 0) = NaN;
                    for ii = 1:size(tmp,3)
                        y_data(x,ii) = nanmean(nanmean((tmp(:,:,ii).*roi)));
                    end
                end
            else
                y_data(1,:) = squeeze(handles.data_selected_resized.image(i).reco.data(voxel(2), voxel(1), :, slice_nbre, 1));
            end
            y_data(isnan(y_data)) = [];
            if ~isempty(y_data) && numel(y_data) ~= 1
                if ~isfield(handles.data_selected_resized.image(i).reco, 'echotime')
                    if (~isfield(handles.data_selected_resized.image(i).reco, 'Had'))
                        x_data = 1:size(handles.data_selected_resized.image(i).reco.data,3);
                    else
                        continue
                    end
                else
                    if  ((~isempty(regexpi(handles.data_selected_resized.image(1).acq.ppl_name,'\w*dpcasl\w*')))|| (~isempty(regexpi(handles.data_selected_resized.image(1).acq.ppl_name,'\w*Had\w*'))))
                        x_data = 1:size(handles.data_selected_resized.image(i).reco.data,3);
                    else
                        x_data = handles.data_selected_resized.image(i).reco.echotime;
                    end
                end
                %                 user_response = questdlg('Plot by?', 'Bip', 'VOI', 'Voxels', 'Cancel', 'Cancel');
                %                 if strcmp(user_response, 'Cancel') ==1
                %                     return
                %                 end
                user_response = 'VOI';
                
                if strcmp(user_response, 'Voxels') ==1
                    for x = 1:numel(match_roi)
                        roi = []; tmp = [];
                        tmp = squeeze(handles.data_selected_resized.image(i).reco.data(:, :, :, slice_nbre, 1));
                        roi = handles.ROI_selected_resized(x).data.value(:,:,slice_nbre);
                        roi(roi == 0) = NaN;
                        y_data= tmp(:,:,:).*repmat(roi, [1 1 size(tmp,3)]);
                        indice_xy = findn(~isnan(roi) == 1);
                        for y = 1:size(indice_xy,1)
                            if ~isempty(get(handles.MIA_plot1, 'Children'))
                                hold(handles.MIA_plot1, 'on');
                            end
                            plot(handles.MIA_plot1,x_data,squeeze(y_data(indice_xy(y,1),indice_xy(y,2),:))',...
                                'Color',rgb(handles.colors{x}),...
                                'Tag', strcat('MIA_plot1_1d', strii));
                            hold(handles.MIA_plot1, 'on');
                        end
                    end
                else
                    for x = 1:size(y_data,1)
                        if ~isempty(get(handles.MIA_plot1, 'Children'))
                            hold(handles.MIA_plot1, 'on');
                        end
                        plot(handles.MIA_plot1,x_data,y_data(x,:), 'Color',rgb(handles.colors{numel(get(handles.MIA_plot1, 'Children'))+1}),...
                            'Tag', strcat('MIA_plot1_1d', strii));
                        hold(handles.MIA_plot1, 'on');
                        if match_roi > 0
                            legende_txt = [legende_txt strcat(handles.data_selected.list_scan(i), '-', handles.ROI_selected_resized(x).name)];
                        else
                            legende_txt = [legende_txt handles.data_selected.list_scan(i)];
                        end
                        if size(y_data,1) == 1
                            tmp2(i,1:numel(y_data)) =  y_data;
                        end
                    end
                end
                
            end
            hold(handles.MIA_plot1, 'off');
        end
        %plot ratio
        %         if exist('tmp2', 'var') && size(tmp2,1) == 2
        %             figure; plot(x_data, tmp2(2,:)./tmp2(1,:));
        %         end
        %
        if ~isempty(legende_txt)
            legend(handles.MIA_plot1, legende_txt, 'Location','NorthEast');
        end
    elseif size(handles.data_selected_resized.image(1).reco.data, 5) > 1
        % plot repetition
        for i=1:handles.data_loaded(1).number_of_scan
            strii = num2str(i);
            if ~isempty(get(handles.MIA_plot1, 'Children'))
                hold(handles.MIA_plot1, 'on');
            end
            if match_roi > 0
                % plot ROI
                for x = 1:numel(match_roi)
                    roi = []; tmp = [];
                    tmp = squeeze(handles.data_selected_resized.image(i).reco.data(:, :, 1, slice_nbre, :));
                    roi = handles.ROI_selected_resized(x).data.value(:,:,slice_nbre);
                    roi(roi == 0) = NaN;
                    for ii = 1:size(tmp,3)
                        y_data(x,ii) = nanmean(nanmean((tmp(:,:,ii).*roi)));
                    end
                end
            else
                % plot voxel
                clear y_data
                y_data(1,:) = squeeze(handles.data_selected_resized.image(i).reco.data(voxel(2), voxel(1), 1, slice_nbre, 1:eval(['get(handles.MIA_data' strii '_expt_slider, ''max'')'])));
            end
            if ~isempty(y_data) && numel(y_data) ~= 1
                if isfield(handles.data_selected_resized.image(i), 'texte')
                    if strcmp(handles.data_selected.image(i).texte(1:27),'# === DATA DESCRIPTION FILE')
                        timing = strfind(handles.data_selected.image(i).texte,'Scan Duration [sec]');
                        duree_tot =str2double(handles.data_selected.image(i).texte(1,timing+38:timing+41));
                        TR = duree_tot /handles.data_selected_resized.image(i).reco.no_expts;
                        x_data =TR:TR:duree_tot;
                    elseif isfield(handles.data_selected_resized.image(i).reco, 'paramQuantif')
                        method = scan_acqp('CASL scan info',handles.data_selected_resized.image(i).reco.paramQuantif,2);
                        if ~isempty(strfind(method,'PhSw'))
                            pCASL_PhaseStart=scan_acqp('CASL_PhaseStart_Deg=',handles.data_selected_resized.image(i).texte,1);
                            pCASL_PhaseEnd=scan_acqp('CASL_PhaseEnd_Deg=',handles.data_selected_resized.image(i).texte,1);
                            pCASL_PhaseInc=scan_acqp('CASL_PhaseInc_Deg=',handles.data_selected_resized.image(i).texte,1);
                            pCASL_Phase = pCASL_PhaseStart:pCASL_PhaseInc:pCASL_PhaseEnd-1;
                            %x_data = pCASL_Phase(2:end);      %exp 1 is a dummy scan
                            %x_data = mod(x_data,360);
                            x_data = 0:15:360;
                            %y_data=y_data(:,2:end);
                            %y_data=horzcat(y_data(:,26:end),y_data(:,8:26));
                            y_data=horzcat((y_data(:,2:7)+y_data(:,26:end))/2,y_data(:,8:26));
                            %y_data=horzcat(y_data(:,2),y_data(:,27:end),y_data(:,8:26));
                        else
                            timing = sscanf(scan_acqp('##$PVM_ScanTimeStr=',handles.data_selected_resized.image(i).texte,2),'%dh%dm%ds%dms');
                            duree_tot=timing(1)*3600+timing(2)*60+timing(3)+timing(4)/1000; %en seconde
                            TR = duree_tot / handles.data_selected_resized.image(i).reco.no_expts;
                            x_data = TR:TR:duree_tot;
                        end
                    else
                        timing = sscanf(scan_acqp('##$PVM_ScanTimeStr=',handles.data_selected_resized.image(i).texte,2),'%dh%dm%ds%dms');
                        duree_tot=timing(1)*3600+timing(2)*60+timing(3)+timing(4)/1000; %en seconde
                        TR = duree_tot / handles.data_selected_resized.image(i).reco.no_expts;
                        x_data = TR:TR:duree_tot;
                    end
                else
                    x_data = handles.data_selected_resized.image(i).acq.tr:handles.data_selected_resized.image(i).acq.tr:handles.data_selected_resized.image(i).acq.tr*...
                        handles.data_selected_resized.image(i).acq.no_expts;
                end
                for x = 1:size(y_data,1)
                    if ~isempty(get(handles.MIA_plot1, 'Children'))
                        hold(handles.MIA_plot1, 'on');
                    end
                    plot(handles.MIA_plot1,x_data,y_data(x,:), 'Color',rgb(handles.colors{i}),...
                        'Tag', strcat('MIA_plot1_1d', strii));
                    hold(handles.MIA_plot1, 'off');
                end
                legende_txt = [legende_txt strcat(handles.data_selected.list_scan(i), '-x', num2str(voxel(2)), '-y', num2str(voxel(1)), '-z', num2str(slice_nbre))];
            end
            if isfield(handles.data_selected_resized.image(i).reco, 'paramQuantif')
                method = scan_acqp('CASL scan info',handles.data_selected_resized.image(i).reco.paramQuantif,2);
                if ~isempty(strfind(method,'PhSw'))
                    hold(handles.MIA_plot1, 'on');
                    controlPhStart = str2double(scan_acqp('CASL_PhaseCorrC=',handles.data_selected_resized.image(i).texte,0));
                    plot(handles.MIA_plot1,rem(360-controlPhStart,360),(max(y_data(x,:))+min(y_data(x,:)))/2,'+','Color','r','Linewidth',2);
                    hold(handles.MIA_plot1, 'off');
                    figure, plot(x_data,y_data);
                    hold on, plot(rem(360-controlPhStart,360),(max(y_data(x,:))+min(y_data(x,:)))/2,'+','Color','r','Linewidth',2);
                    hold off
                end
            end
            
        end
        if ~isempty(legende_txt)
            legend(handles.MIA_plot1, legende_txt, 'Location','NorthEast');
        end
    end
    
else
    % PRM mode
    warndlg('The click function for the PRM mode is not coded yet', 'Warning');
end






function MIA_table1_add_pixel_value(~,handles,pixel_coordinates)

[pixel_coordinates_2d] = [round(pixel_coordinates(1,1)) round(pixel_coordinates(1,2)) round(pixel_coordinates(1,3))];
% if pixel_coordinates_2d(1) > handles.resolution_selected || pixel_coordinates_2d(2) > handles.resolution_selected ||...
%         pixel_coordinates_2d(1) == 0 || pixel_coordinates_2d(2) == 0
%     return % bug somewhere !!!
% end
if ~isfield(handles, 'data_displayed')
    return
end
values = zeros(1,size(handles.data_displayed.image, 4));
for i = 1:size(handles.data_displayed.image, 4)
    if numel(size(handles.data_displayed.image))>4 && max(max(handles.data_displayed.image(:, :,pixel_coordinates_2d(3),i, 2))) > 0
        histo_data_in_gray = rgb2gray(uint8(squeeze(handles.data_displayed.image(:, :,pixel_coordinates_2d(3),i,:))));
        pixel_value = 255-histo_data_in_gray(pixel_coordinates_2d(2), pixel_coordinates_2d(1));
        if pixel_value == 255
            pixel_value = 0;
        end
        values(i) = pixel_value;
        
    else
        % rgb2ind(uint8(squeeze(handles.data_displayed.image(pixel_coordinates_2d(2), pixel_coordinates_2d(1),pixel_coordinates_2d(3),2, :))), 32)
        values(i) = handles.data_displayed.image(pixel_coordinates_2d(2), pixel_coordinates_2d(1),pixel_coordinates_2d(3),i,1);
        
    end
end
table_data = get(handles.MIA_table_pixel_values, 'data');
table_data(1,2:1+size(handles.data_displayed.image, 4)) = num2cell(values);
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



% --- Executes on button press in MIA_set_clip.
function MIA_set_clip_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_set_clip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'clips')
    return
end

tt = sort(handles.clips(:,1));
clips_sorted = cell(size(handles.clips));
% For remove deleted line
deletedClip = cellfun(@isempty,tt);
if sum(deletedClip) >= 1
    tt = tt(sum(deletedClip)+1:end,1);
    clips_sorted = cell(size(tt,1),size(handles.clips,2));
end
% for i = 1:numel(tt)
%     clips_sorted(i,:) = handles.clips(strcmp(tt(i), handles.clips(:,1)),:);
% end
% handles.clips = clips_sorted;
if isempty(findobj('Tag', 'MIA_clip_table'))
    
    set(handles.MIA_set_clip, 'String', 'Update Clip');
    tmp(1:size(handles.clips,2)) = 1;
    uitable('Tag', 'MIA_clip_table', 'Units', 'normalized','Position', [0.56 0.015 0.44 0.79],...
        'ColumnWidth',{175},...
        'ColumnName',{'Parameters','Clip - Min','Clip - Max'},...
        'ColumnEditable', logical(tmp),...
        'Data',handles.clips);
else
    set(handles.MIA_set_clip, 'String', 'Set Clip');
    %     handles.clips = get(findobj('Tag', 'MIA_clip_table'), 'Data');
    new_data = get(findobj('Tag', 'MIA_clip_table'), 'Data');
    % update names of the parameters
    new_parameters_name = new_data(:,1);
    isnew = find(~strcmp(handles.clips(:,1), new_parameters_name)==1);
    handles.clips = new_data;
    for i = 1:numel(isnew)
        if isempty(new_parameters_name{isnew(i),1})
            old_parameter_name = handles.clips(isnew(i),1);
            handles.clips(isnew(i),:) = [];
        else
            for ii = 1:numel(handles.database)
                for j = 1:numel(handles.database(ii).day)
                    match = find(strcmp(handles.clips(isnew(i),1), handles.database(ii).day(j).parameters') ==1);
                    for jj = 1:numel(match)
                        %update database
                        handles.database(ii).day(j).parameters(match(jj)) = new_parameters_name(isnew(i));
                        %rename corresponding file
                        [PATHSTR,NAME,EXT] = fileparts([handles.database(ii).path   handles.database(ii).day(j).scans_file{match(jj)}]);
                        new_name = [handles.database(ii).name '-' handles.database(ii).day(j).date '-' new_parameters_name{isnew(i)}];
                        new_name = strrep(new_name, '*', 'star');
                        if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
                            warning_text = sprintf('##$ This file no not exist\n##$ %s',...
                                fullfile(PATHSTR,[NAME,EXT]));
                            msgbox(warning_text, 'rename file warning') ;
                        else
                            movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                            handles.database(ii).day(j).scans_file{match(jj)} = [new_name,EXT];
                        end
                    end
                end
            end
            
            %update handles.clips listing
            if find(strcmp(handles.clips(:,1),  new_parameters_name(isnew(i))),1) > 0
                handles.clips(isnew(i),:) = [];
            else
                handles.clips(isnew(i),1) = new_parameters_name(isnew(i));
            end
            % Update scan_liste's name and image title
            if isfield(handles, 'data_displayed')
                if handles.mode == 1
                    for ii=1:handles.data_selected_resized.scan_number
                        if strcmp(handles.data_selected_resized.list_scan{ii}, handles.clips(isnew(i),1))
                            strii = num2str(ii);
                            handles.data_selected_resized.list_scan{ii} = new_parameters_name{isnew(i)};
                            set(eval(['handles.MIA_data' strii '_title']), 'String',new_parameters_name{isnew(i)});
                        end
                    end
                else
                    for ii=1:handles.data_selected_for_PRM_resized.scan_number
                        if strcmp(handles.data_selected_for_PRM_resized.list_scan{ii}, handles.clips(isnew(i),1))
                            strii = num2str(ii);
                            handles.data_selected_for_PRM_resized.list_scan{ii} = new_parameters_name{isnew(i)};
                            set(eval(['handles.MIA_data' strii '_title']), 'String',new_parameters_name{isnew(i)});
                        end
                    end
                end
            end
        end
    end
    %clear table1
    delete(findobj('Tag', 'MIA_clip_table'));
    %update MIA_plot1
    guidata(hObject, handles);
    MIA_update_database_display(hObject, eventdata, handles)
    %update MIA_axes if needed
    if isfield(handles, 'data_selected_resized')
        MIA_update_axes(hObject, eventdata, handles)
    end
    
    % save database
    path_root = pwd;
    handles= guidata(handles.MIA_GUI);
    cd(handles.database(1).databaseinfo.pathname);
    handles.database(1).databaseinfo.VOIs = handles.VOIs;
    handles.database(1).databaseinfo.clips = handles.clips;
    database = handles.database; %#ok<NASGU>
    save(handles.database(1).databaseinfo.filename, 'database');
    cd(path_root);
end





% --- Executes on button press in MIA_table1_add_cluster.
function MIA_table1_add_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_table1_add_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_loaded') ||  ~isfield(handles, 'ROI_selected_resized')
    return
end
if strcmp(get(hObject, 'Tag'), 'MIA_load_axes')
    answer = [];
    for i = 1:numel({handles.ROI_selected_resized.name})
        if strcmp({handles.ROI_selected_resized(i).name}, handles.table1.cluster') == 0
            answer = [answer {handles.ROI_selected_resized(i).name}]; %#ok<AGROW>
        end
    end
    if isempty(answer)
        return
    end
    handles=guidata(handles.MIA_GUI);
else
    answer = inputdlg({'Name of the new cluster?'}, 'Bip', [1 30]);
    if isempty(answer)
        return
    end
end
cluster_name = {};
for i = 1:numel({handles.ROI_selected_resized.name})
    for j = 1:numel(answer)
        cluster_name = [cluster_name strcat(handles.ROI_selected_resized(i).name, '-', answer(j))]; %#ok<AGROW>
    end
end


for i = 1:numel(cluster_name)
    table_data = get(handles.MIA_table1, 'data');
    handles.table1.cluster = [handles.table1.cluster cluster_name(i)];
    handles.table1.cluster_row = [handles.table1.cluster_row size(table_data,1)+1];
    
    table_data(size(table_data,1),1) = {''};
    table_data(size(table_data,1)+1,1) = strcat(...
        '<html><span style="color: #FF0000; font-weight: bold;">', ...
        cluster_name(i), ...
        '</span></html>');
    table_data(size(table_data,1)+1,1) = {'lim Min'};
    table_data(size(table_data,1),2:5) = handles.table1.clips(:,2)';
    table_data(size(table_data,1)+1,1) = {'lim Max'};
    table_data(size(table_data,1),2:5) = handles.table1.clips(:,3)';
    table_data(size(table_data,1)+1,1) = {'Mean'};
    table_data(size(table_data,1)+1,1) = {'Median'};
    table_data(size(table_data,1)+1,1) = {'SD'};
    table_data(size(table_data,1)+1,1) = {'% of ROI'};
    table_data(size(table_data,1)+1,1) = {''};
    set(handles.MIA_table1, 'Data', table_data);
    
    guidata(hObject, handles);
    
end
% upadate data in MIA_table1
MIA_table1_update_data_Callback(hObject,eventdata,handles)
MIA_update_axes(hObject, eventdata, handles)



% --- Executes on button press in MIA_table1_remove_cluster.
function MIA_table1_remove_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_table1_remove_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.table1.cluster)
    return
end

[selection,ok] = listdlg('Name', 'Select the cluster to remove', 'ListString', handles.table1.cluster);
if ok == 0
    return
end
if isfield(handles.data_displayed.image, 'VOI_cluster')
    handles = rmfield(handles, 'ROI_cluster_resized');
    handles.data_displayed = rmfield(handles.data_displayed, 'VOI_cluster');
    %     guidata(hObject, handles)
    %     MIA_update_axes(hObject, eventdata, handles)
end

for i = 1:numel(selection)
    table_data = get(handles.MIA_table1, 'data');
    table_data(handles.table1.cluster_row(selection(i)):handles.table1.cluster_row(selection(i))+7,:) = [];
    for j=1:size(handles.table1.cluster_row,2)
        if handles.table1.cluster_row(j) > handles.table1.cluster_row(selection(i))
            handles.table1.cluster_row(j) = handles.table1.cluster_row(j) - 8;
        end
    end
    set(handles.MIA_table1, 'Data', table_data);
    
    guidata(hObject, handles);
end
handles.table1.cluster(selection) = [];
handles.table1.cluster_row(selection) = [];

if isfield(handles,'ROI_cluster_resized') && isempty(handles.ROI_cluster_resized)
    handles= rmfield(handles, 'ROI_cluster_resized');
end
guidata(hObject, handles);
MIA_table1_update_data_Callback(hObject, eventdata, handles)




% --- Executes on button press in MIA_table1_update_data.
function MIA_table1_update_data_Callback(hObject,eventdata, handles)
%global data_in_table
% hObject    handle to MIA_table1_update_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if isempty(handles.table1.cluster)
    %     if isfield(handles.data_displayed, 'VOI_cluster')
    %         handles.data_displayed= rmfield(handles.data_displayed, 'VOI_cluster');
    %     end
    %     if isfield(handles, 'ROI_cluster_resized')
    %         handles= rmfield(handles, 'ROI_cluster_resized');
    %     end
    %     guidata(hObject, handles);
    %     MIA_update_axes(hObject, eventdata, handles)
    return
end
if ~isfield(handles, 'ROI_selected_resized')
    warndlg('Please load a ROI');
    return
end
table_data = get(handles.MIA_table1, 'data');

clear sub_data

%initialize the ROI_cluster matrix
if size(handles.data_displayed.image,3)==1
    handles.ROI_cluster_resized.map = zeros([size(handles.ROI_selected_resized(1).data.value),1, 3]);
    handles.ROI_cluster_resized.trans = zeros([size(handles.ROI_selected_resized(1).data.value),1]);
else
    handles.ROI_cluster_resized.map = zeros([size(handles.ROI_selected_resized(1).data.value), 3]);
    handles.ROI_cluster_resized.trans = zeros(size(handles.ROI_selected_resized(1).data.value));
end
trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
sep = strfind(handles.table1.cluster, '-');
if ~isempty(cell2mat(sep))
    for i= 1:numel(sep)
        short_cluster_name{i} = handles.table1.cluster{i}(1:sep{i}(end)-1); %#ok<AGROW>
    end
end
for i = 1:numel(handles.table1.cluster)
    no_data = 0;
    exit = 0;
    for j = 1:handles.data_loaded(1).number_of_scan
        if exit == 0
            if ~exist('sub_data', 'var')
                tmp = findn(handles.data_ploted.coordonates(:,4+j) > cell2num(table_data(handles.table1.cluster_row(i)+1,1+j)));
                sub_data = handles.data_ploted.coordonates(tmp(:,1),:);
            else
                tmp = findn(sub_data(:,4+j) > cell2num(table_data(handles.table1.cluster_row(i)+1,1+j)));
                if isempty(tmp)
                    no_data = 0;
                    exit =1;
                    sub_data = tmp;
                else
                    sub_data = sub_data(tmp(:,1),:);
                    no_data = 1;
                end
            end
            if exist('sub_data', 'var') && ~isempty(sub_data)
                tmp = findn(sub_data(:,4+j) < cell2num(table_data(handles.table1.cluster_row(i)+2,1+j)));
                if isempty(tmp)
                    no_data = 0;
                    exit =1;
                    sub_data = tmp;
                else
                    sub_data = sub_data(tmp(:,1),:);
                    no_data = 1;
                end
            end
        end
    end
    % calculate Mean/Median/SD/% of ROI
    for j = 1:handles.data_loaded(1).number_of_scan
        if no_data == 1
            if ~isempty(cell2mat(sep))
                match = strcmp(short_cluster_name{i}, {handles.ROI_selected_resized.name});
                if ~isempty(sub_data)
                    data_in_ROI = findn(sub_data(:,4)== find(match == 1));
                    if ~isempty(data_in_ROI)
                        sub_data = sub_data(data_in_ROI(:,1),:);
                        table_data(handles.table1.cluster_row(i)+3,1+j) = num2cell(mean(sub_data(:,4+j)));
                        table_data(handles.table1.cluster_row(i)+4,1+j) = num2cell(median(sub_data(:,4+j)));
                        table_data(handles.table1.cluster_row(i)+5,1+j) = num2cell(std(sub_data(:,4+j)));
                        table_data(handles.table1.cluster_row(i)+6,2) = num2cell(size(sub_data,1)/sum(handles.data_ploted.coordonates(:,4) == find(match == 1))*100);
                    else
                        table_data(handles.table1.cluster_row(i)+3,1+j) = {'No voxel'};
                        table_data(handles.table1.cluster_row(i)+4,1+j) = {'No voxel'};
                        table_data(handles.table1.cluster_row(i)+5,1+j) = {'No voxel'};
                        table_data(handles.table1.cluster_row(i)+6,2) = {'0'};
                    end
                else
                    table_data(handles.table1.cluster_row(i)+3,1+j) = {'No voxel'};
                    table_data(handles.table1.cluster_row(i)+4,1+j) = {'No voxel'};
                    table_data(handles.table1.cluster_row(i)+5,1+j) = {'No voxel'};
                    table_data(handles.table1.cluster_row(i)+6,2) = {'0'};
                end
            else
                table_data(handles.table1.cluster_row(i)+3,1+j) = num2cell(mean(sub_data(:,4+j)));
                table_data(handles.table1.cluster_row(i)+4,1+j) = num2cell(median(sub_data(:,4+j)));
                table_data(handles.table1.cluster_row(i)+5,1+j) = num2cell(std(sub_data(:,4+j)));
                table_data(handles.table1.cluster_row(i)+6,2) = num2cell(size(sub_data,1)/size(handles.data_ploted.coordonates,1)*100);
            end
        else
            table_data(handles.table1.cluster_row(i)+3,1+j) = {'No voxel'};
            table_data(handles.table1.cluster_row(i)+4,1+j) = {'No voxel'};
            table_data(handles.table1.cluster_row(i)+5,1+j) = {'No voxel'};
            table_data(handles.table1.cluster_row(i)+6,2) = {'0'};
            %             guidata(hObject, handles);
        end
    end
    % drow cluster
    if ~isempty(sub_data)
        if ~isempty(sub_data)
            for ii = 1:numel(sub_data(:,1))
                handles.ROI_cluster_resized.map(sub_data(ii,1), sub_data(ii,2), sub_data(ii,3),:) = handles.colors_rgb(i,:);
                handles.ROI_cluster_resized.trans(sub_data(ii,1), sub_data(ii,2), sub_data(ii,3)) = trans;
            end
        end
    end
    %     guidata(hObject, handles);
    clear sub_data
end
%update MIA_table1
set(handles.MIA_table1, 'Data', table_data);
guidata(hObject, handles);
% MIA_update_axes(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MIA_plot1_select_dots_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_plot1_select_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'ROI_selected_resized')
    warndlg('Please load a ROI');
    return
end

k = waitforbuttonpress; %#ok<NASGU>
point1 = get(gca,'CurrentPoint');
% button down detected
finalRect = rbbox;            %#ok<NASGU> % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
lim_min = min(point1,point2);
lim_max = max(point1,point2);
if lim_min == lim_max
    return
end
if size(handles.data_displayed.image,4) == 4
    data4_min =min(min(get(findobj('Tag', 'data4'), 'CData')));
    data4_max =max(max(get(findobj('Tag', 'data4'), 'CData')));
else
    data4_min = 0;
    data4_max = 0;
end
orientation = num2str(get(handles.MIA_plot1, 'View'));
switch orientation
    case num2str([0 90])
        z = get(handles.MIA_plot1, 'ZTick');
        table1_lim_min = [num2cell(lim_min(1,1)) num2cell(lim_min(1,2)) min(z) data4_min];
        table1_lim_max = [num2cell(lim_max(1,1)) num2cell(lim_max(1,2)) max(z)  data4_max];
    case num2str([0 0])
        y = get(handles.MIA_plot1, 'YTick');
        table1_lim_min = [num2cell(lim_min(1,1))  min(y) num2cell(lim_min(1,3)) data4_min];
        table1_lim_max = [num2cell(lim_max(1,1)) max(y) num2cell(lim_max(1,3)) data4_max];
    case num2str([90 0])
        x = get(handles.MIA_plot1, 'XTick');
        table1_lim_min = [min(x) num2cell(lim_min(1,2)) num2cell(lim_min(1,3)) data4_min];
        table1_lim_max = [max(x) num2cell(lim_max(1,2)) num2cell(lim_max(1,3)) data4_max];
    otherwise
        warndlg('Please swith to a 2D view', 'Warning');
        return
end
if size(handles.data_displayed.image,4) == 1
    table1_lim_min(1,2:4) = {''};
    table1_lim_max(1,2:4) = {''};
end

answer = inputdlg({'Name of the new cluster?'}, 'Bip', [1 30]);
table_data = get(handles.MIA_table1, 'data');
if isempty(answer)
    return
end
cluster_name = {};
for i = 1:numel({handles.ROI_selected_resized.name})
    for j = 1:numel(answer)
        cluster_name = [cluster_name strcat(handles.ROI_selected_resized(i).name, '-', answer(j))]; %#ok<AGROW>
    end
end

for i = 1:numel(cluster_name)
    handles.table1.cluster = [handles.table1.cluster cluster_name(i)];
    handles.table1.cluster_row = [handles.table1.cluster_row size(table_data,1)+1];
    
    table_data(size(table_data,1),1) = {''};
    table_data(size(table_data,1)+1,1) = strcat(...
        '<html><span style="color: #FF0000; font-weight: bold;">', ...
        cluster_name(i), ...
        '</span></html>');
    table_data(size(table_data,1)+1,1) = {'lim Min'};
    table_data(size(table_data,1),2:5) = table1_lim_min;
    table_data(size(table_data,1)+1,1) = {'lim Max'};
    table_data(size(table_data,1),2:5) = table1_lim_max;
    table_data(size(table_data,1)+1,1) = {'Mean'};
    table_data(size(table_data,1)+1,1) = {'Median'};
    table_data(size(table_data,1)+1,1) = {'SD'};
    table_data(size(table_data,1)+1,1) = {'% of ROI'};
    table_data(size(table_data,1)+1,1) = {''};
    set(handles.MIA_table1, 'Data', table_data);
end
guidata(hObject, handles);
% upadate data in MIA_table1
MIA_table1_update_data_Callback(hObject,eventdata,handles)


%         % upadate data in MIA_table1
%         MIA_table1_update_data_Callback(hObject,eventdata,handles)
%
%         %% find voxels
%         dot_selected_with_x_min = findn(handles.data_ploted.coordonates(:,5) > lim_min(1,1));
%         if isempty(dot_selected_with_x_min)
%             return
%         end
%         dot_selected_with_x_min_max = findn(handles.data_ploted.coordonates(dot_selected_with_x_min(:,1),5) < lim_max(1,1));
%          if isempty(dot_selected_with_x_min_max)
%             return
%         end
%         dot_selected_with_x_min_max_y_min  = findn(handles.data_ploted.coordonates(dot_selected_with_x_min(dot_selected_with_x_min_max(:,1),1),6) > lim_min(1,2));
%          if isempty(dot_selected_with_x_min_max_y_min)
%             return
%         end
%         dot_selected  = findn(handles.data_ploted.coordonates(dot_selected_with_x_min(dot_selected_with_x_min_max(dot_selected_with_x_min_max_y_min(:,1),1),1),6) < lim_max(1,2));
%          if isempty(dot_selected)
%             return
%         end
%         tmp = handles.data_ploted.coordonates(dot_selected_with_x_min(dot_selected_with_x_min_max(dot_selected_with_x_min_max_y_min(dot_selected(:,1),1),1),1),:);


% --- Executes on button press in MIA_mode_single_button.
function MIA_mode_single_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_mode_single_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update display mode PRM (multi-time point) --> single time point

handles = MIA_clear_data(hObject, eventdata, handles);

set(handles.MIA_mode_single_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_mode_multi_button, 'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.mode = 1;



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
set(handles.MIA_set_clip, 'Position', [0.4275 0.0364 0.0719 0.025]);

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
% reset every cluster for now
handles.table1.cluster = [];
handles.table1.cluster_row = [];

%uncheck the mask option
set(handles.MIA_menu_define_mask, 'Checked', 'off');

% restore the manual contrast to 0;
handles.display_option.manual_contrast = 0;

if isfield(handles, 'data_selected_resized')
    handles = rmfield(handles, 'data_selected_resized');
end
if isfield(handles, 'data_displayed')
    handles = rmfield(handles, 'data_displayed');
end

set(handles.MIA_menu_image_information,'Checked','off');

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
set(handles.MIA_table1, 'Data', {'', '', '', '', ''});
set(handles.MIA_patient_information_title, 'String', 'No images');

set(handles.MIA_orientation_space_popupmenu, 'String', 'Select orientation')

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
handles = MIA_clear_data(hObject, eventdata, handles);

set(handles.MIA_mode_multi_button, 'backgroundColor', [1,0,0]);
set(handles.MIA_mode_single_button,'backgroundColor', [0.9412, 0.9412, 0.9412]);
handles.mode = 2;


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
set(handles.MIA_set_clip, 'Position', [0.4275 0.330 0.0719 0.025]);
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
set(handles.MIA_menu_image_information,'Checked','off');

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
coordonates = handles.data_ploted.coordonates;
% determine PRM+, PRM0 and PRM-
CI = str2double(get(handles.MIA_PRM_CI, 'String'));
PRMr = coordonates(:,6) > coordonates(:,5)+CI;
PRMb = coordonates(:,6) < coordonates(:,5)-CI;
PRMg = ~or(PRMr, PRMb);
% get(handles.MIA_data1)

scatter(handles.MIA_plot1,coordonates(PRMr,5), coordonates(PRMr,6),...
    'filled',...
    'SizeData', 20,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','r',...
    'Visible', 'on',...
    'Tag', 'MIA_plot1_PRM_r');
hold(handles.MIA_plot1, 'on');
scatter(handles.MIA_plot1,coordonates(PRMb,5), coordonates(PRMb,6),...
    'filled',...
    'SizeData', 20,...
    'MarkerFaceColor','b',...
    'MarkerEdgeColor','b',...
    'Visible', 'on',...
    'Tag', 'MIA_plot1_PRM_b');

scatter(handles.MIA_plot1,coordonates(PRMg,5), coordonates(PRMg,6),...
    'filled',...
    'SizeData', 20,...
    'MarkerFaceColor','g',...
    'MarkerEdgeColor','g',...
    'Visible', 'on',...
    'Tag', 'MIA_plot1_PRM_g');

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
if numel(size(handles.ROI_selected_resized.data.value)) == 2
    handles.ROI_PRM_resized.map=zeros([size(handles.ROI_selected_resized.data.value),1,3]);
else
    handles.ROI_PRM_resized.map=zeros([size(handles.ROI_selected_resized.data.value),3]);
end

handles.ROI_PRM_resized.trans=zeros(size(handles.ROI_selected_resized.data.value));
trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
for i =1:numel(PRMr)
    if PRMr(i) == 1
        handles.ROI_PRM_resized.map(coordonates(i,1), coordonates(i,2), coordonates(i,3),:) = [1 0 0];
    elseif PRMg(i) == 1
        handles.ROI_PRM_resized.map(coordonates(i,1), coordonates(i,2), coordonates(i,3),:) = [0 1 0];
    elseif PRMb(i) ==1
        handles.ROI_PRM_resized.map(coordonates(i,1), coordonates(i,2), coordonates(i,3),:) = [0 0 1];
    end
    handles.ROI_PRM_resized.trans(coordonates(i,1), coordonates(i,2), coordonates(i,3)) = trans;
end

%update table1
table_data = cell(8,3);
table_data(1:8,1) = {'PRMr (%)', 'PRMg (%)', 'PRMb (%)', 'Mean_pre', 'SD_pre', 'Mean_post', 'SD_post', 'Vol (mm3)'};

table_data(1,3) = {sum(PRMr)/numel(coordonates(:,1))*100};
table_data(2,3) = {sum(PRMg)/numel(coordonates(:,1))*100};
table_data(3,3) = {sum(PRMb)/numel(coordonates(:,1))*100};
table_data(4,3) = {mean(coordonates(:,5))};
table_data(5,3) = {std(coordonates(:,5))};
table_data(6,3) = {mean(coordonates(:,6))};
table_data(7,3) = {std(coordonates(:,6))};

voxel_volume = handles.data_selected_for_PRM_resized.image(1).acq.fov(1)/handles.resolution_selected *...
    handles.data_selected_for_PRM_resized.image(1).acq.fov(1)/handles.resolution_selected*...
    handles.data_selected_for_PRM_resized.image(1).acq.thickness;

table_data(8,3) = {numel(coordonates(:,1))*voxel_volume};
set(handles.MIA_table1, 'Data', table_data);

guidata(hObject, handles);

function handles = MIA_update_PRM_Overlay_map(hObject,handles)

slice_nbr = get(handles.MIA_slider_slice, 'Value');
% for i = 1:numel(handles.ROI_selected_resized)
for slice_nbr=1:get(handles.MIA_slider_slice, 'Max');
    handles.data_displayed.PRM.data(:,:,slice_nbr,:)=squeeze(handles.ROI_PRM_resized.map(:,:,slice_nbr,:));
    handles.data_displayed.PRM.trans(:,:,slice_nbr) = handles.ROI_PRM_resized.trans(:,:,slice_nbr);
end
% end



% handles = guidata(hObject);
% if findobj('Tag', 'ROI_contour')
% delete(findobj('Tag', 'ROI_contour'));
% end
% if findobj('Tag', 'data1_PRM')
%     delete(findobj('Tag', 'data1_PRM'));
%     delete(findobj('Tag', 'data2_PRM'));
% end
% %update images with colormap
% slice = get(handles.MIA_slider_slice, 'Value');
%
% hold(handles.MIA_data1, 'on');
% image(squeeze(handles.ROI_PRM_resized.map(:,:,slice,:)), 'parent',...
%     handles.MIA_data1, 'AlphaData',handles.ROI_PRM_resized.trans(:,:,slice), 'Tag', 'data1_PRM');
% hold(handles.MIA_data1, 'off');
%
%
% hold(handles.MIA_data2, 'on');
% image(squeeze(handles.ROI_PRM_resized.map(:,:,slice,:)), 'parent',...
%     handles.MIA_data2, 'AlphaData',handles.ROI_PRM_resized.trans(:,:,slice), 'Tag', 'data2_PRM');
% hold(handles.MIA_data2, 'off');

function handles = MIA_update_VOI_cluster_displayed(hObject,handles)
% slice_nbr = get(handles.MIA_slider_slice, 'Value');
for slice_nbr=1:get(handles.MIA_slider_slice, 'Max');
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


if strcmp(get(handles.MIA_menu_define_mask, 'Checked'), 'on')
    tolerance=get(handles.MIA_PRM_slider_trans, 'Value');
    MIA_mask(hObject, eventdata, handles, tolerance);
elseif strcmp(get(handles.MIA_menu_roi_fill, 'Checked'), 'on')
    MIA_update_axes(hObject, eventdata, handles)
elseif isfield(handles, 'brain_extraction_ROI')
    MIA_Brain_Extraction(hObject, eventdata, handles)
else
    if isfield(handles, 'data_selected_for_PRM') && isfield(handles, 'ROI_PRM_resized')
        trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
        index = findn(handles.ROI_PRM_resized.map ~=0);
        for i = 1: size(index,1)
            handles.ROI_PRM_resized.trans(index(i,1),index(i,2),index(i,3)) = trans;
        end
        guidata(hObject,handles)
        MIA_update_axes(hObject, eventdata, handles)
    elseif isfield(handles, 'data_loaded') && isfield(handles, 'ROI_cluster_resized')
        trans = round(get(handles.MIA_PRM_slider_trans, 'Value'))/100;
        if trans == 0
            trans = 0.01;
        end
        handles.ROI_cluster_resized.trans(handles.ROI_cluster_resized.trans ~=0)=trans;
        guidata(hObject,handles)
        MIA_update_axes(hObject, eventdata, handles)
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


% --- Executes on button press in MIA_speedy_mode.
function MIA_speedy_mode_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_speedy_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    warndlg('No database loaded', 'Warning');
    return
end
speedy_module(hObject, eventdata, handles)

% --- Executes on button press in MIA_table1_store_cluster.
function MIA_table1_store_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_table1_store_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'data_loaded')
    return
end
if isempty(handles.table1.cluster)
    return
end
% define every name of cluster available
tmp= strfind(handles.table1.cluster, '-');
cluster_name=cell(numel(handles.table1.cluster):1);
for i = 1:numel(handles.table1.cluster)
    cluster_name{i} = handles.table1.cluster{i}(tmp{i}+1:end);
end
cluster_name = unique(cluster_name);

[selection,ok] = listdlg('Name', 'Bip', 'ListString', cluster_name,'ListSize', [200 150], 'PromptString', 'Select the cluster to store');
if ok == 0
    return
end
table1_data = get(handles.MIA_table1, 'Data');
parameter_name = unique(get(handles.MIA_table1, 'ColumnName'));
% Remove the fist columnName (useless)
parameter_name = parameter_name(2:end);
for i = 1:numel(selection)
    % First cluster saved
    if ~isfield(handles.database(1).databaseinfo, 'cluster')
        handles.database(1).databaseinfo.cluster(1).cluster_name = cluster_name(i);
        match = strfind(handles.table1.cluster, cluster_name{i});
        done = 0;
        for j=1:numel(match)
            if match{j} >0 && done == 0 % &
                handles.database(1).databaseinfo.cluster(1).parameter_name = [];
                handles.database(1).databaseinfo.cluster(1).min = [];
                handles.database(1).databaseinfo.cluster(1).max = [];
                for jj=1:numel(parameter_name)
                    handles.database(1).databaseinfo.cluster(1).parameter_name =  [handles.database(1).databaseinfo.cluster(1).parameter_name parameter_name(jj)];
                    handles.database(1).databaseinfo.cluster(1).min = [handles.database(1).databaseinfo.cluster(1).min table1_data(handles.table1.cluster_row(j)+1,jj+1)];
                    handles.database(1).databaseinfo.cluster(1).max = [handles.database(1).databaseinfo.cluster(1).max table1_data(handles.table1.cluster_row(j)+2,jj+1)];
                end
                done = 1;
            end
        end
    else
        cluster_nbr = numel(handles.database(1).databaseinfo.cluster)+1;
        handles.database(1).databaseinfo.cluster(cluster_nbr).cluster_name = cluster_name(i);
        match = strfind(handles.table1.cluster, cluster_name{i});
        done = 0;
        for j=1:numel(match)
            if match{j} >0 && done == 0 % &
                for jj=1:numel(parameter_name)
                    handles.database(1).databaseinfo.cluster(cluster_nbr).parameter_name =  [handles.database(1).databaseinfo.cluster(cluster_nbr).parameter_name parameter_name(jj)];
                    handles.database(1).databaseinfo.cluster(cluster_nbr).min = [handles.database(1).databaseinfo.cluster(cluster_nbr).min table1_data(handles.table1.cluster_row(j)+1,jj+1)];
                    handles.database(1).databaseinfo.cluster(cluster_nbr).max = [handles.database(1).databaseinfo.cluster(cluster_nbr).max table1_data(handles.table1.cluster_row(j)+2,jj+1)];
                end
                done = 1;
            end
        end
    end
end

guidata(hObject, handles);



% --- Executes on button press in MIA_table1_load_cluster.
function MIA_table1_load_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_table1_load_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_loaded')
    return
end
if ~isfield(handles.database(1).databaseinfo, 'cluster')
    warndlg('No cluster saved', 'Warning');
    return
end
if ~isfield(handles, 'ROI_selected_resized')
    warndlg('Load a least 1 ROI', 'Warning');
    return
end

cluster_list = [handles.database(1).databaseinfo.cluster.cluster_name];

[selection,ok] = listdlg('Name', 'Bip', 'ListString', cluster_list,'ListSize', [200 150], 'PromptString', 'Select the clusters you want to load');
if ok == 0
    return
end

cluster_name = {};
parameter_loaded_name = unique(get(handles.MIA_table1, 'ColumnName'));
parameter_name = parameter_loaded_name(2:end);
for j = 1:numel(selection)
    if sum(strcmp(handles.database(1).databaseinfo.cluster(selection(j)).parameter_name, parameter_name')) ~=...
            numel(handles.database(1).databaseinfo.cluster(selection(j)).parameter_name)
        parameters_needed = handles.database(1).databaseinfo.cluster(selection(j)).parameter_name;
        message = cell(numel(parameters_needed)+2,1);
        message{1,1} = 'The scans loaded do not correspond to scans used ';
        message{2,1} = 'to defind this cluster. You must load:';
        message{3,1} = '';
        for i = 1:numel(parameters_needed)
            message(3+i,1) = parameters_needed(i);
        end
        warndlg(message, 'Warning');
        return
    end
    for i = 1:numel({handles.ROI_selected_resized.name})
        cluster_name = strcat(handles.ROI_selected_resized(i).name, '-', cluster_list(selection(j)));
        table_data = get(handles.MIA_table1, 'data');
        handles.table1.cluster = [handles.table1.cluster cluster_name];
        handles.table1.cluster_row = [handles.table1.cluster_row size(table_data,1)+1];
        
        table_data(size(table_data,1),1) = {''};
        table_data(size(table_data,1)+1,1) = strcat(...
            '<html><span style="color: #FF0000; font-weight: bold;">', ...
            cluster_name, ...
            '</span></html>');
        table_data(size(table_data,1)+1,1) = {'lim Min'};
        table_data(size(table_data,1),2:2+numel(handles.database(1).databaseinfo.cluster(selection(j)).min)-1) =handles.database(1).databaseinfo.cluster(selection(j)).min;
        table_data(size(table_data,1)+1,1) = {'lim Max'};
        table_data(size(table_data,1),2:2+numel(handles.database(1).databaseinfo.cluster(selection(j)).max)-1) = handles.database(1).databaseinfo.cluster(selection(j)).max;
        table_data(size(table_data,1)+1,1) = {'Mean'};
        table_data(size(table_data,1)+1,1) = {'Median'};
        table_data(size(table_data,1)+1,1) = {'SD'};
        table_data(size(table_data,1)+1,1) = {'% of ROI'};
        table_data(size(table_data,1)+1,1) = {''};
        set(handles.MIA_table1, 'Data', table_data);
        
        guidata(hObject, handles);
    end
end
% upadate data in MIA_table1
MIA_table1_update_data_Callback(hObject,eventdata,handles)


% --------------------------------------------------------------------
function MIA_menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MIA_3D_rendering_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_3D_rendering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles.data_loaded, 'ROI')
    return
end

% Display the surface
figure;
transparency = 0.2;
if handles.mode == 1
    if isfield(handles, 'ROI_cluster_resized')
        D = handles.ROI_cluster_resized.map;
        D(isnan(D)) = 0;
        aa = reshape(D, [size(D,1)*size(D,2)*size(D,3) 3]);
        
        classes = unique(aa,'rows');
        for i = 1:size(classes,1)
            if sum(classes(i,:) == [0 0 0]) ~= 3
                % resize
                sub_class = D(:,:,:,1) == classes(i,1) &  D(:,:,:,2) == classes(i,2) & D(:,:,:,3) == classes(i,3);
                
                sub_class = padarray(sub_class,[1 1 1],'both');
                
                % Create an isosurface
                Ds = smooth3(sub_class, 'gaussian', 5);
                
                surface = isosurface(Ds,0.3333);
                if ~isempty(surface.vertices)
                    
                    surface.vertices(:,1) = surface.vertices(:,1)*handles.data_selected_resized.image(1).acq.fov(1)/handles.data_selected_resized.image(1).acq.no_samples;
                    surface.vertices(:,2) = surface.vertices(:,2)*handles.data_selected_resized.image(1).acq.fov(2)/handles.data_selected_resized.image(1).acq.no_views;
                    surface.vertices(:,3) = surface.vertices(:,3)*handles.data_selected_resized.image(1).acq.thickness;
                    
                    % subplot(1,2,1);
                    hiso = patch('Vertices',surface.vertices,...
                        'Faces',surface.faces,...
                        'FaceColor',classes(i,:),...
                        'EdgeColor','none');
                    alpha(0.7)
                end
            end
        end
    else
        for i=1:numel(handles.data_loaded.ROI)
            D =  size(handles.data_loaded.ROI(i).nii);
            
            % resize
            D = padarray(D,[1 1 1],'both');
            
            % Create an isosurface
            Ds = smooth3(D, 'gaussian', 5);
            % Ds = D;
            surface = isosurface(Ds,0.3333);
            surface.vertices(:,1) = surface.vertices(:,1)*handles.data_selected_resized.image(1).acq.fov(1)/handles.data_selected_resized.image(1).acq.no_samples;
            surface.vertices(:,2) = surface.vertices(:,2)*handles.data_selected_resized.image(1).acq.fov(2)/handles.data_selected_resized.image(1).acq.no_views;
            surface.vertices(:,3) = surface.vertices(:,3)*handles.data_selected_resized.image(1).acq.thickness;
            %             surface.vertices(:,1) = surface.vertices(:,1)*abs(handles.data_loaded.ROI(i).V.mat(1,1);
            %             surface.vertices(:,2) = surface.vertices(:,2)*abs(handles.data_loaded.ROI(i).V.mat(2,2));
            %             surface.vertices(:,3) = surface.vertices(:,3)*abs(handles.data_loaded.ROI(i).V.mat(3,3));
            
            % subplot(1,2,1);
            hiso = patch('Vertices',surface.vertices,...
                'Faces',surface.faces,...
                'FaceColor',handles.colors_rgb(i,:),...
                'EdgeColor','none');
            alpha_value = transparency+i*0.1;
            alpha(alpha_value)
        end
        
    end
else
    rgb_color  = [1 0 0; 0 1 0; 0 0 1];
    
    test = reshape(handles.ROI_PRM_resized.map, [size(handles.ROI_PRM_resized.map, 1)*size(handles.ROI_PRM_resized.map, 2)*size(handles.ROI_PRM_resized.map, 3) 3]);
    for i= [1 3]
        %
        %         if PRMr(i) == 1
        %             handles.ROI_PRM_resized.map(coordonates(i,1), coordonates(i,2), coordonates(i,3),:) = [1 0 0];
        %         elseif PRMg(i) == 1
        %             handles.ROI_PRM_resized.map(coordonates(i,1), coordonates(i,2), coordonates(i,3),:) = [0 1 0];
        %         elseif PRMb(i) ==1
        %             handles.ROI_PRM_resized.map(coordonates(i,1), coordonates(i,2), coordonates(i,3),:) = [0 0 1];
        %         end
        %         findn(test == 1 0 0);
        D = zeros([size(test,1) 1]);
        %red
        if i == 1
            D(test(:,1) == 1 & test(:,2) == 0 & test(:,3) == 0) = 1;
            %green
        elseif i == 2
            D(test(:,1) == 0 & test(:,2) == 1 & test(:,3) == 0) = 1;
            % blue
        elseif i == 3
            D(test(:,1) == 0 & test(:,2) == 0 & test(:,3) == 1) = 1;
        end
        
        D = reshape(D, [size(handles.ROI_PRM_resized.map, 1), size(handles.ROI_PRM_resized.map,2 ), size(handles.ROI_PRM_resized.map, 3)]);
        
        % resize
        D = padarray(D,[1 1 1],'both');
        
        % remove isolated modification
        D = bwareaopen(D, 5);
        % Create an isosurface
        Ds = smooth3(D, 'gaussian', 5);
        % Ds = D;
        surface = isosurface(Ds,0.3333);
        if ~isempty(surface.vertices)
            surface.vertices(:,1) = surface.vertices(:,1)*handles.data_selected_for_PRM_resized.image(1).acq.fov(1)/handles.data_selected_for_PRM_resized.image(1).acq.no_samples;
            surface.vertices(:,2) = surface.vertices(:,2)*handles.data_selected_for_PRM_resized.image(1).acq.fov(2)/handles.data_selected_for_PRM_resized.image(1).acq.no_views;
            surface.vertices(:,3) = surface.vertices(:,3)*handles.data_selected_for_PRM_resized.image(1).acq.thickness;
            
            
            %         handles.colors ={'b', 'g', 'y', 'm', 'c', 'r', 'k', 'w'};
            % handles.colors_rgb = [0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 0 0; 1 1 1];
            
            % subplot(1,2,1);
            hiso = patch('Vertices',surface.vertices,...
                'Faces',surface.faces,...
                'FaceColor',rgb_color(i,:),...
                'EdgeColor','none');
            
            alpha(0.5);
        end
    end
    
    
end

view(45,30)
axis tight
daspect([1,1,.4])
lightangle(45,30);
% set(gcf,'Renderer','zbuffer');
lighting gouraud %phong
% isonormals(Ds,hiso)

set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)
set(gca, 'Visible', 'on', 'XTick' , [], 'YTick', [], 'ZTick', []);


% % Reconstruct the volume and display it as montage
% OV = surface2volume(surface,[],1);
% nDims = size(OV);
% subplot(1,2,2);
% figure;montage(reshape(OV,nDims(1),nDims(2),1,nDims(3)),[0 1]);


% --- Executes on button press in MIA_table1_save_cluster_as_VOI.
function MIA_table1_save_cluster_as_VOI_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_table1_save_cluster_as_VOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.table1.cluster)
    return
end

[selection,ok] = listdlg('Name', 'Select the cluster to save', 'ListString', handles.table1.cluster);
if ok == 0
    return
end
[VOI_number, ok] = listdlg('PromptString','Select the new scan name:',...
    'Name', 'Select a Name',...
    'SelectionMode','single',...
    'ListSize', [400 300],...
    'ListString',handles.VOIs');
if ok == 0
    return
end

if strcmp('Other', handles.VOIs(VOI_number)) == 1
    ROI_name=inputdlg('Name of the new VOI','VOI''s name',1,{'New VOI'});
    handles.VOIs = [handles.VOIs(1:end-1), ROI_name, 'Other'];
else
    ROI_name = handles.VOIs(VOI_number);
end
% determine patient id and time point
id = handles.data_selected.patient_info.name_nbr;

time_point = handles.data_selected.patient_info.timepoint_nbr;

for i=1:get(handles.MIA_slider_slice, 'Max')
    rgb = squeeze(handles.data_displayed.VOI_cluster(selection).data(:,:,i,1))+squeeze(handles.data_displayed.VOI_cluster(selection).data(:,:,i,2))+squeeze(handles.data_displayed.VOI_cluster(selection).data(:,:,i,3));
    rgb(rgb>0)=1;
    uvascroi(i).value=rgb;
    uvascroi(i).area_in_pixel=sum(sum(uvascroi(i).value));
    uvascroi(i).position=[0,0; 0, handles.resolution_selected; handles.resolution_selected, handles.resolution_selected; handles.resolution_selected, 0];
    uvascroi(i).name=ROI_name{:};
    uvascroi(i).colorindex=1;
    %%% case roi drown on histo data
    if isfield(handles.data_selected.image(1).reco, 'color_map')
        uvascroi(i).no_samples = size(handles.data_selected.image.reco.data,1);
        uvascroi(i).no_views=size(handles.data_selected.image.reco.data,2);
        uvascroi(i).fov_offsets=handles.data_selected.image(1).reco.fov_offsets;
    else
        uvascroi(i).displayedecho=handles.data_selected.image(1).reco.displayedecho;
        uvascroi(i).displayedslice=handles.data_selected.image(1).reco.displayedslice;
        uvascroi(i).displayedexpt=handles.data_selected.image(1).reco.displayedexpt;
        uvascroi(i).fov_orientation=handles.data_selected.image(1).reco.fov_orientation(1:9,handles.data_selected.image(1).reco.displayedecho,i,handles.data_selected.image(1).reco.displayedexpt);
        uvascroi(i).no_samples=handles.data_selected.image(1).reco.no_samples;
        if isfield(handles.data_selected.image(1).reco, 'no_views')
            uvascroi(i).no_views=handles.data_selected.image(1).reco.no_views;
        end
        uvascroi(i).fov=handles.data_selected.image(1).acq.fov;
        if isfield(handles.data_selected.image(1).reco, 'no_views')
            uvascroi(i).area_in_mmxmm=(uvascroi(i).fov(1)/uvascroi(i).no_samples)*(uvascroi(i).fov(2)/uvascroi(i).no_views)*uvascroi(i).area_in_pixel;
        end
        uvascroi(i).reco_number=handles.data_selected.image(1).reco_number;
        uvascroi(i).scan_number=handles.data_selected.image(1).scan_number;
        uvascroi(i).filename=handles.data_selected.image(1).filename;
        uvascroi(i).filenames=sprintf('%s\n',handles.data_selected.image(1).filename);
        uvascroi(i).fov_offsets=handles.data_selected.image(1).reco.fov_offsets(1:3,handles.data_selected.image(1).reco.displayedecho,i,handles.data_selected.image(1).reco.displayedexpt);
        
    end
    uvascroi(i).thickness = handles.data_selected.image(1).reco.thickness;
    uvascroi(i).date=date;
    uvascroi(i).affiche=0;
end
if ~isempty(uvascroi)
    pathname = [handles.database(id).path handles.database(id).name handles.database(id).day(time_point).date '-VOI-' ROI_name{:} '.mat'];
    save(pathname,'uvascroi');
    handles.database(id).day(time_point).VOIs = [handles.database(id).day(time_point).VOIs ROI_name];
    voi_file = [handles.database(id).name handles.database(id).day(time_point).date '-VOI-' ROI_name{:} '.mat'];
    handles.database(id).day(time_point).VOIs_file = [handles.database(id).day(time_point).VOIs_file {voi_file}];
end

guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles)


% --------------------------------------------------------------------
function MIA_plot1_plot_all_voxels_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_plot1_plot_all_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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




% --------------------------------------------------------------------
function MIA_menu_save_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_save_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles.database.Properties.UserData , 'db_filename')
    path_root = pwd;
    [filename,pathname] = uiputfile('.mat','Please enter the name of the database', path_root );%pathstr(1:file_sep(end)));
    if pathname == 0
        return
    end
    handles.database.Properties.UserData.db_filename = filename;
    handles.database.Properties.UserData.MIA_root_path = pathname;
end

handles.database.Properties.UserData.VOIs = handles.VOIs;
handles.database.Properties.UserData.histo = handles.histo;
database = handles.database; %#ok<NASGU>
save(fullfile(handles.database.Properties.UserData.MIA_root_path, handles.database.Properties.UserData.db_filename), 'database');

% save handles
guidata(hObject, handles)

msgbox('Done', 'Information') ;

% --- Executes on button press in MIA_new_roi.
function MIA_new_roi_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_new_roi (see GCBO)
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
    [x, y, ~] = size(handles.data_selected_for_PRM_resized.image(which_image).reco.data);
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
    if sum(strcmp(newVOI_name, ROI_listing) > 0)
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
                handles.data_loaded.ROI(end).nii = read_volume(handles.data_loaded.ROI(end).V , handles.data_loaded.Scan(which_image).V(1));
            else
                handles.data_loaded.ROI.V = spm_vol(fullfilename(handles, ROI_idex, '.nii'));
                handles.data_loaded.ROI.nii = read_volume(handles.data_loaded.ROI.V , handles.data_loaded.Scan(which_image).V(1));
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
            list_of_choices = {'Union', 'Instersection', 'Delete'};
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
            set(handles.MIA_menu_view_voxel_on_map, 'Check', 'off');
            handles.display_option.view_pixel_on_map = 0;
            guidata(hObject, handles);
            
            eval(['hroi=impoly(handles.MIA_data' image_number ',[]);']);
            
            set(handles.MIA_menu_view_voxel_on_map, 'Check', 'on');
            handles.display_option.view_pixel_on_map = 1;
            guidata(hObject, handles);
        else
            eval(['hroi=impoly(handles.MIA_data' image_number ',[]);']);
            
        end
    case 'Delete'
        
        %         handles.ROI_selected(ROI_selected_nbr).data(match_slice) = [];
        %         % delete the roi file and the roi from the database if empty
        %         if isempty(handles.ROI_selected(ROI_selected_nbr).data)
        %             handles.ROI_selected(ROI_selected_nbr) = [];
        %             if isfield(handles, 'ROI_selected_resized')
        %                 handles.ROI_selected_resized(ROI_selected_nbr) = [];
        %             end
        %             % remove from database
        %             set(handles.MIA_scans_list, 'Value', VOI_number);
        %             MIA_update_database_display(hObject, eventdata, handles);
        %             MIA_remove_scan_Callback(hObject, eventdata, handles)
        %             handles = guidata(hObject);
        %             % delete file
        %             pathname = [handles.database(patient).path VOI_file_name{:}];
        %             delete(pathname);
        %             if isempty(handles.ROI_selected)
        %                 handles = rmfield(handles, 'ROI_selected');
        %                 if isfield(handles, 'ROI_selected_resized')
        %                     handles = rmfield(handles, 'ROI_selected_resized');
        %                 end
        %                 if isfield(handles, 'ROI_PRM_resized')
        %                     handles = rmfield(handles, 'ROI_PRM_resized');
        %                 end
        %                 guidata(hObject, handles)
        %                 MIA_update_axes(hObject, eventdata, handles)
        %                 return
        %             end
        %             guidata(hObject, handles)
        %             MIA_update_axes(hObject, eventdata, handles)
        %             return
        %         end
        %         uvascroi = handles.ROI_selected(ROI_selected_nbr).data;
end




position = getPosition(hroi);
Scan_of_reference_selected = get(handles.MIA_orientation_space_popupmenu, 'Value');
% Scan_of_reference_listing = get(handles.MIA_orientation_space_popupmenu, 'String');
% Scan_of_reference = Scan_of_reference_listing(Scan_of_reference_selected,:);

if isempty(position) || size(position, 1) ==1
    % case of delete?
    return
end
switch ROI_case
    case 'Delete'
        
    case 'New slice'  % new ROI slice added to an existing ROI
        handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = createMask(hroi);
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(handles.data_loaded.ROI(ROI_loaded_idex).nii, handles.data_loaded.Scan(Scan_of_reference_selected).V);
        
    case 'New ROI'
        % create the VOI matrix
        ROI_matrix=zeros([x y z]);
        ROI_matrix(:,:,slice_nbr)=createMask(hroi);
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(ROI_matrix, handles.data_loaded.Scan(Scan_of_reference_selected).V);
        
    case 'Union'
        handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) + createMask(hroi);
        %         handles.data_loaded.ROI(ROI_loaded_idex).nii(handles.data_loaded.ROI(ROI_loaded_idex).nii > 0 ) = 1;
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(handles.data_loaded.ROI(ROI_loaded_idex).nii, handles.data_loaded.Scan(Scan_of_reference_selected).V);
    case 'Instersection'
        handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr) = or(handles.data_loaded.ROI(ROI_loaded_idex).nii(:,:,slice_nbr), createMask(hroi)) - createMask(hroi);
        % transform the ROI_matrix in order to match to the nii hearder
        % (rotation/translation)
        ROI_matrix = write_volume(handles.data_loaded.ROI(ROI_loaded_idex).nii, handles.data_loaded.Scan(Scan_of_reference_selected).V);
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
if exist([handles.database.Properties.UserData.MIA_data_path, 'ROI_data', filesep], 'dir') ~= 7
    status = mkdir([handles.database.Properties.UserData.MIA_data_path, 'ROI_data', filesep]);
    if status == 0
        warndlg('You do not have the right to write in the folder!', 'Warning');
        return
    end
end
V_ROI.fname =  [handles.database.Properties.UserData.MIA_data_path, 'ROI_data', filesep, file_name, '.nii'] ;

V_ROI = rmfield(V_ROI,'private'); % Delete old nifti header. Will be recreated to match new image properties
V_ROI.dt(1) = spm_type(outputDatatype); % save images in specified format
if spm_type(outputDatatype,'intt')
    V_ROI = rmfield(V_ROI,'pinfo'); % integer datatype : let spm_write_vol decide on scaling
else
    V_ROI.pinfo(1:2) = [1;0];       % do not apply any scaling when saving as float data
end
% save the ROI in nii file (could be a new ROI or and old but updated)
spm_write_vol(V_ROI,ROI_matrix);
% spm_jsonwrite(fullfilename(handles, roi_index, '.json'), hroi_data);
% load the new ROI
if isfield(handles.data_loaded, 'ROI')
    if sum(ROI_loaded_idex) ~= 0 % update the structre for an updated ROI
        handles.data_loaded.ROI(ROI_loaded_idex).V = spm_vol(V_ROI.fname);
        handles.data_loaded.ROI(ROI_loaded_idex).nii = read_volume(handles.data_loaded.ROI(ROI_loaded_idex).V , handles.data_loaded.Scan(Scan_of_reference_selected).V);
    else %% add new ROI to the data_loaded_ROI structure (another ROI is already loaded)
        handles.data_loaded.ROI(numel(handles.data_loaded.ROI)+1).V = spm_vol(V_ROI.fname);
        handles.data_loaded.ROI(end).nii = read_volume(handles.data_loaded.ROI(end).V , handles.data_loaded.Scan(Scan_of_reference_selected).V);
    end
else %% add the new ROI as load ROI
    handles.data_loaded.ROI.V = spm_vol(V_ROI.fname);
    handles.data_loaded.ROI.nii = read_volume(handles.data_loaded.ROI.V , handles.data_loaded.Scan(Scan_of_reference_selected).V);
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






% % save handles
% guidata(hObject, handles)
%
% % find every VOI loaded already
% if isfield(handles, 'ROI_selected_resized')
%     voi_to_load = [{handles.ROI_selected_resized.name} newVOI_name];
% else
%     voi_to_load = newVOI_name;
% end
%
% % select patient/timepoint/Voi mode and select the good one
% % patient
% set(handles.MIA_name_list, 'Value', patient);
% MIA_update_database_display(hObject, eventdata, handles);
%
% % time point
% set(handles.MIA_time_points_list, 'Value', timepoint);
% MIA_update_database_display(hObject, eventdata, handles);
%
% % VOI
% set(handles.MIA_scans_list, 'Value', 1);
% set(handles.MIA_scan_VOIs_button, 'Value', 1)
% MIA_update_database_display(hObject, eventdata, handles);
%
% % Select the good VOI
% match_VOI_to_open = [];
% for i = 1:numel(voi_to_load)
%     if strcmp(voi_to_load(i), newVOI_name) ==1 && new ==1
%     else
%         match_VOI_to_open = [match_VOI_to_open find(strcmp(get(handles.MIA_scans_list, 'String'), voi_to_load(i)) ==1)];
%     end
% end
% match_VOI_to_open = sort(unique(match_VOI_to_open));
% set(handles.MIA_scans_list, 'Value', match_VOI_to_open);
% MIA_update_database_display(hObject, eventdata, handles);
%
% % Load the VOIs
% MIA_load_axes_Callback(hObject, eventdata, handles)
%
% % update the handles
% handles = guidata(hObject);


% if new == 0
%     for i = 1:numel(handles.ROI_selected)
%         % find the good VOI
%         if strcmp(handles.ROI_selected(i).data(1).name, newVOI_name) ==1
%             match_ROI_selected = [match_ROI_selected 1];
%             % check if the VOI to drow is already on this slice
%             z_offset = [handles.ROI_selected(i).data.fov_offsets];
%             if handles.mode == 1
%                 match_slice = find(abs(z_offset(3,:) - handles.data_selected_resized.z_offset(slice_nbr)) < 1e-3,1);
%             else
%                 match_slice = find(abs(z_offset(3,:) - handles.data_selected_for_PRM_resized.z_offset(slice_nbr))< 1e-3,1);
%
%             end
%
%             if ~isempty(match_slice)
%                 list_of_choices = {'Modify', 'Union', 'Instersection', 'Delete'};
%                 [user_response,ok] = listdlg('Name', 'Bip', 'ListString', list_of_choices,'ListSize', [200 150], 'PromptString', 'What do you want to do?', 'SelectionMode', 'single');
%                 if isempty(user_response) || ok == 0
%                     return
%                 end
%                 user_response = list_of_choices{user_response};
%             end
%         else
%             match_ROI_selected = [match_ROI_selected 0];
%         end
%     end
% end

% drow voi// warning: need to uncheck the view voxel on map before drawing
% the new VOI
% if handles.mode == 1
%     scan_nbr = handles.data_loaded(1).number_of_scan;
% else
%
%     scan_nbr = 2;
% end
% ROI_selected_nbr = find(match_ROI_selected ==1);




% create the VOI matrix (using the same resolution as the old roi
%         [x_roi, ~, ~] = size(handles.ROI_selected(ROI_selected_nbr).data(1).value);
%         % need to interpolate the new roi and the position in case of differences between old
%         factor = x_roi/x;
%         %         bw = imresize(bw,factor,'bilinear');
%         for iii = 1:numel(handles.ROI_selected(ROI_selected_nbr).data)
%             handles.ROI_selected(ROI_selected_nbr).data(iii).position =  handles.ROI_selected(ROI_selected_nbr).data(iii).position.*factor;
%         end

%         done = 0;
%         switch user_response
%             case 'Modify'
%
%             case 'Union'
%                 %                 z_offset = [handles.ROI_selected(ROI_selected_nbr).data.fov_offsets];
%                 %                 nbr = find(z_offset(3,:) == handles.data_selected_resized.image(1).reco.fov_offsets(3,1,slice_nbr,1),1);
%                 %                 old_mask = handles.ROI_selected_resized(ROI_selected_nbr).data(match_slice).value;
%                 old_mask = handles.ROI_selected_resized(ROI_selected_nbr).data.value;
%                 bw = bw+logical(old_mask(:,:,slice_nbr));
%                 bw(bw>1)=1;
%             case 'Instersection'
%                 %                 z_offset = [handles.ROI_selected(ROI_selected_nbr).data.fov_offsets];
%                 %                 nbr = find(z_offset(3,:) == handles.data_selected_resized.image(1).reco.fov_offsets(3,1,slice_nbr,1),1);
%                 %                 old_mask = handles.ROI_selected(ROI_selected_nbr).data(match_slice).value;
%                 %                 bw = or(logical(old_mask), bw)-bw;
%                 old_mask = handles.ROI_selected_resized(ROI_selected_nbr).data.value;
%                 bw = or(logical(old_mask(:,:,slice_nbr)), bw)-bw;
%                 bw(bw<0)=0;
%                 if sum(sum(bw)) == 0
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice) = [];
%                     if isfield(handles, 'ROI_selected_resized')
%                         handles.ROI_selected_resized(ROI_selected_nbr) = [];
%                     end
%                     % delete the roi file and the roi from the database if empty
%                     if isempty(handles.ROI_selected(ROI_selected_nbr).data)
%                         handles.ROI_selected(ROI_selected_nbr) = [];
%                         % delete file
%                         pathname = [handles.database(patient).path VOI_file_name{:}];
%                         delete(pathname);
%                         if isempty(handles.ROI_selected)
%                             handles = rmfield(handles, 'ROI_selected');
%                             if isfield(handles, 'ROI_selected_resized')
%                                 handles = rmfield(handles, 'ROI_selected_resized');
%                             end
%                             % remove frome database
%                             set(handles.MIA_scans_list, 'Value', VOI_number);
%                             MIA_update_database_display(hObject, eventdata, handles);
%                             MIA_remove_scan_Callback(hObject, eventdata, handles)
%                             guidata(hObject, handles)
%                             MIA_update_axes(hObject, eventdata, handles)
%                             return
%                         end
%                         % remove frome database
%                         set(handles.MIA_scans_list, 'Value', VOI_number);
%                         MIA_update_database_display(hObject, eventdata, handles);
%                         MIA_remove_scan_Callback(hObject, eventdata, handles)
%                         guidata(hObject, handles)
%                         MIA_update_axes(hObject, eventdata, handles)
%                         return
%                     end
%                     uvascroi = handles.ROI_selected(ROI_selected_nbr).data;
%                     done = 1;
%                 end
%             case 'New'
%                 match_slice = size(handles.ROI_selected(ROI_selected_nbr).data, 2)+1;
%                 handles.ROI_selected(ROI_selected_nbr).data(match_slice) =  handles.ROI_selected(ROI_selected_nbr).data(match_slice-1);
%         end
%         if done == 0
%
%             % update the resolusion of the ROI (cf. handles.resolution)
%             for iii = 1:numel(handles.ROI_selected(ROI_selected_nbr).data)
%                 handles.ROI_selected(ROI_selected_nbr).data(iii).no_samples=x;
%                 handles.ROI_selected(ROI_selected_nbr).data(iii).no_views=y;
%                 if ~strcmp(user_response, 'New')
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).value =  handles.ROI_selected_resized(ROI_selected_nbr).data.value(:,:,match_slice);
%                 end
%             end
%             handles.ROI_selected(ROI_selected_nbr).data(match_slice).position = position;
%             handles.ROI_selected(ROI_selected_nbr).data(match_slice).value = bw;
%             handles.ROI_selected(ROI_selected_nbr).data(match_slice).area_in_pixel = sum(bw(:));
%             handles.ROI_selected(ROI_selected_nbr).data(match_slice).displayedslice = slice_nbr;
%
%             if handles.mode == 1
%                 if numel(size(handles.data_selected_resized.image(which_image).reco.fov_offsets)) == 4
%                     match = abs((round(squeeze(handles.data_selected_resized.image(which_image).reco.fov_offsets(3,1,:,1))*100)/100) - handles.data_selected_resized.z_offset(slice_nbr)) < 1e-3;
%                 else
%                     match = abs((round(squeeze(handles.data_selected_resized.image(which_image).reco.fov_offsets(3,1,:,1))*100)/100) - handles.data_selected_resized.z_offset(slice_nbr)) < 1e-3;
%                 end
%                 %match = find(match==1,1,'first');
%                 if numel(size(handles.data_selected_resized.image(which_image).reco.fov_offsets)) == 2
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_offsets(:,1) =handles.data_selected_resized.image(which_image).reco.fov_offsets(:,match);
%                 elseif handles.data_selected.image(which_image).reco.no_expts == 1 && handles.data_selected.image(which_image).reco.no_echoes == 1
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_offsets =handles.data_selected_resized.image(which_image).reco.fov_offsets(:,1,match);
%                 elseif handles.data_selected.image(which_image).reco.no_expts > 1 && handles.data_selected.image(which_image).reco.no_echoes == 1
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_offsets = handles.data_selected.image(which_image).reco.fov_offsets(1:3,1,slice_nbr,1);
%                 end
%                 if numel(size(handles.data_selected_resized.image(which_image).reco.fov_offsets)) == 4
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_offsets = handles.data_selected_resized.image(which_image).reco.fov_offsets(:,1,match,1);
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_orientation = handles.data_selected.image(which_image).reco.fov_orientation(1:9,1,match,1);
%                 else
%                     if isfield(handles.data_selected.image(which_image).reco, 'fov_orientation');
%                         handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_orientation = handles.data_selected.image(which_image).reco.fov_orientation(1:9,1,match);
%                     end
%                 end
%             else
%                 if numel(size(handles.data_selected_for_PRM_resized.image(which_image).reco.fov_offsets)) == 4
%                     match = abs(squeeze(handles.data_selected_for_PRM_resized.image(which_image).reco.fov_offsets(3,1,:,1)) - handles.data_selected_for_PRM_resized.z_offset(slice_nbr)) < 1e-3;
%                 else
%                     match = abs(squeeze(handles.data_selected_for_PRM_resized.image(which_image).reco.fov_offsets(3,1,:)) - handles.data_selected_for_PRM_resized.z_offset(slice_nbr)) < 1e-3;
%                 end
%                 handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_offsets(:,1) =handles.data_selected_for_PRM_resized.image(which_image).reco.fov_offsets(:,match);
%                 if numel(size(handles.data_selected_for_PRM_resized.image(which_image).reco.fov_offsets)) == 4
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_orientation = handles.data_selected_for_PRM_resized.image(1).reco.fov_orientation(1:9,1,match,1);
%                 else
%                     if isfield(handles.data_selected_for_PRM_resized.image(which_image).reco, 'fov_orientation')
%                         handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov_orientation = handles.data_selected_for_PRM_resized.image(which_image).reco.fov_orientation(1:9,1,match);
%                     end
%                 end
%             end
%             if isfield(handles.ROI_selected(ROI_selected_nbr).data(match_slice), 'fov')
%                 handles.ROI_selected(ROI_selected_nbr).data(match_slice).area_in_mmxmm = ...
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov(1)/handles.ROI_selected(ROI_selected_nbr).data(match_slice).no_samples*...
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).fov(2)/handles.ROI_selected(ROI_selected_nbr).data(match_slice).no_views*...
%                     handles.ROI_selected(ROI_selected_nbr).data(match_slice).area_in_pixel;
%             end
%             uvascroi = handles.ROI_selected(ROI_selected_nbr).data;
%         end

% if new == 1

% save handles
%     guidata(hObject, handles)

%     MIA_update_database_display(hObject, eventdata, handles);
% end



%
% % Select the VOI included the new one
% match_VOI_to_open = [];
% for i = 1:numel(voi_to_load)
%     match_VOI_to_open = [match_VOI_to_open find(strcmp(get(handles.MIA_scans_list, 'String'), voi_to_load(i)) ==1)];
% end
% match_VOI_to_open = sort(unique(match_VOI_to_open));
% set(handles.MIA_scans_list, 'Value', match_VOI_to_open);
% MIA_update_database_display(hObject, eventdata, handles);
% save handles

% % Load the VOIs
% MIA_load_axes_Callback(hObject, eventdata, handles)
%
% % update database display
% MIA_update_database_display(hObject, eventdata, handles);

% function  handles =IA_save_automatic_ROI_form_Atlas(hObject, eventdata, handles, ROI_name)
function  handles = MIA_save_automatic_ROI_form_Map(hObject, eventdata, handles, Map_name, patient_selected, day_selected, threshold, ROI_name)
% create and save new ROI (ROI_name) using a scan (MAP_name) and a
% threshold [value_min value_max]on this map for a specific patient (patient_selected)
% and day (day_selected)

Image_data = load(char(fullfile(handles.database(patient_selected).path,...
    handles.database(patient_selected).day(day_selected).scans_file(ismember(handles.database(patient_selected).day(day_selected).parameters, Map_name) ==1 ))));

indx = findn(Image_data.uvascim.image.reco.data > threshold(1) & Image_data.uvascim.image.reco.data < threshold(2));

ROI_data = zeros(size(Image_data.uvascim.image.reco.data));

for ii = 1:numel(indx(:,1))
    ROI_data(indx(ii,1), indx(ii,2), 1, indx(ii,4),:) = 1;
end
if strcmp(ROI_name, 'Atlas_Fit')
    ROI_data = imfill(ROI_data);
    ROI_data = imdilate(ROI_data,strel('disk',5));
end
% remove CSF ROI
if strcmp(ROI_name, 'Atlas_GM') || strcmp(ROI_name, 'Atlas_Striat-R') || strcmp(ROI_name, 'Atlas_Striat-L') || strcmp(ROI_name, 'Atlas_WM')
    Atlas_CSF_ROI = load(char(fullfile(handles.database(patient_selected).path,...
        handles.database(patient_selected).day(day_selected).VOIs_file(ismember(handles.database(patient_selected).day(day_selected).VOIs, 'Atlas_CSF') ==1 ))));
end
% remove WM ROI
if  strcmp(ROI_name, 'Atlas_Striat-R') || strcmp(ROI_name, 'Atlas_Striat-L') || strcmp(ROI_name, 'Atlas_Cortex-L')  || strcmp(ROI_name, 'Atlas_Cortex-L')
    Atlas_WM_ROI = load(char(fullfile(handles.database(patient_selected).path,...
        handles.database(patient_selected).day(day_selected).VOIs_file(ismember(handles.database(patient_selected).day(day_selected).VOIs, 'Atlas_WM') ==1 ))));
end
uvascroi(size(ROI_data,4)).value =zeros([size(ROI_data,1) size(ROI_data,2)]);

for i=1:size(ROI_data,4)
    data2D = squeeze(ROI_data(:,:,i));
    se = strel('disk',1);
    uvascroi(i).value=imerode(data2D,se);
    if strcmp(ROI_name, 'Atlas_CSF') || strcmp(ROI_name, 'Atlas_WM')
        uvascroi(i).value = imdilate(data2D,se);
    end
    % remove CSF ROI
    if strcmp(ROI_name, 'Atlas_GM') || strcmp(ROI_name, 'Atlas_Striat-R')  || strcmp(ROI_name, 'Atlas_Striat-L') || strcmp(ROI_name, 'Atlas_WM')
        uvascroi(i).value = uvascroi(i).value .* ~Atlas_CSF_ROI.uvascroi(i).value;
    end
    
    % remove WM ROI
    if strcmp(ROI_name, 'Atlas_Striat-R') || strcmp(ROI_name, 'Atlas_Striat-L')  || strcmp(ROI_name, 'Atlas_Cortex-L') || strcmp(ROI_name, 'Atlas_Cortex-L')
        uvascroi(i).value = uvascroi(i).value .* ~Atlas_WM_ROI.uvascroi(i).value;
    end
    
    uvascroi(i).area_in_pixel=sum(sum(uvascroi(i).value));
    uvascroi(i).position=[0,0; 0, size(ROI_data,1); size(ROI_data,1), size(ROI_data,1); size(ROI_data,1), 0];
    uvascroi(i).name=ROI_name;
    uvascroi(i).colorindex=1;
    uvascroi(i).displayedecho=Image_data.uvascim.image.reco.displayedecho;
    uvascroi(i).displayedslice=Image_data.uvascim.image.reco.displayedslice;
    uvascroi(i).displayedexpt=Image_data.uvascim.image.reco.displayedexpt;
    uvascroi(i).fov_orientation=Image_data.uvascim.image.reco.fov_orientation(1:9,Image_data.uvascim.image.reco.displayedecho,i,Image_data.uvascim.image.reco.displayedexpt);
    uvascroi(i).no_samples=Image_data.uvascim.image.reco.no_samples;
    if isfield(Image_data.uvascim.image.reco, 'no_views')
        uvascroi(i).no_views=Image_data.uvascim.image.reco.no_views;
    end
    uvascroi(i).fov=Image_data.uvascim.image.acq.fov;
    if isfield(Image_data.uvascim.image.reco, 'no_views')
        uvascroi(i).area_in_mmxmm=(uvascroi(i).fov(1)/uvascroi(i).no_samples)*(uvascroi(i).fov(2)/uvascroi(i).no_views)*uvascroi(i).area_in_pixel;
    end
    uvascroi(i).reco_number=Image_data.uvascim.image.reco_number;
    uvascroi(i).scan_number=Image_data.uvascim.image.scan_number;
    uvascroi(i).filename=Image_data.uvascim.image.filename;
    uvascroi(i).filenames=sprintf('%s\n',Image_data.uvascim.image.filename);
    uvascroi(i).fov_offsets=Image_data.uvascim.image.reco.fov_offsets(1:3,Image_data.uvascim.image.reco.displayedecho,i,Image_data.uvascim.image.reco.displayedexpt);
    
end
uvascroi(i).thickness = Image_data.uvascim.image.reco.thickness;
uvascroi(i).date=date;
uvascroi(i).affiche=0;

if ~isempty(uvascroi)
    pathname = [handles.database(patient_selected).path handles.database(patient_selected).name handles.database(patient_selected).day(day_selected).date '-VOI-' ROI_name '.mat'];
    save(pathname,'uvascroi');
    handles.database(patient_selected).day(day_selected).VOIs = [handles.database(patient_selected).day(day_selected).VOIs {ROI_name}];
    voi_file = [handles.database(patient_selected).name handles.database(patient_selected).day(day_selected).date '-VOI-' ROI_name '.mat'];
    handles.database(patient_selected).day(day_selected).VOIs_file = [handles.database(patient_selected).day(day_selected).VOIs_file {voi_file}];
    % Check if this new VOI exist already
    % if not update the handles.VOIs field
    if sum(strcmp(handles.VOIs', {ROI_name})) == 0
        handles.VOIs =  [handles.VOIs(1:end-1) {ROI_name} handles.VOIs(end-1:end)];
    end
end


guidata(findobj('Tag', 'MIA_GUI'), handles);
MIA_update_database_display(hObject, eventdata, handles)




% --------------------------------------------------------------------
function MIA_compute_parametric_maps_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_compute_parametric_maps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% define all the scan available

parameters_list = char(unique(handles.database.SequenceName))';

compute_parametric_maps_menu(handles, parameters_list)


function handles = MIA_add_new_scan_to_database(handles, data, patient, time_point, new_parameter_name, new_parameter_file)

%add to the database
handles.database(patient).day(time_point).parameters = [handles.database(patient).day(time_point).parameters {new_parameter_name}];
handles.database(patient).day(time_point).scans_file = [handles.database(patient).day(time_point).scans_file new_parameter_file];
% update handles.clips
if sum(strcmp(handles.clips(:,1), new_parameter_name)) ==0
    handles.clips(numel(handles.clips(:,1))+1,1) = {new_parameter_name};
    if isfield(data, 'clip') && ~isempty(data.clip)
        handles.clips(numel(handles.clips(:,1)),2) = num2cell(data.clip(1));
        handles.clips(numel(handles.clips(:,1)),3) = num2cell(data.clip(2));
    else
        handles.clips(numel(handles.clips(:,1)),2) = num2cell(data.reco.globalmin);
        handles.clips(numel(handles.clips(:,1)),3) = num2cell(data.reco.globalmax);
    end
end



% --------------------------------------------------------------------
function MIA_cloneScanVoi_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_cloneScanVoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

patient = get(handles.MIA_name_list, 'Value');
time_point= get(handles.MIA_time_points_list, 'Value');
scan = get(handles.MIA_scans_list, 'Value');
Scan_VOI_option = get(handles.MIA_scan_VOIs_button, 'Value');
% Clone scan else clone VOI
if Scan_VOI_option == 0
    file_name = handles.database(patient).day(time_point).scans_file{scan};
    new_file_name = [spm_str_manip(file_name,'r'),'-cloned.mat'];
    parameter_name =  handles.database(patient).day(time_point).parameters(scan);
    new_parameter_name = strcat(parameter_name, '-cloned');
    % update databas
    handles.database(patient).day(time_point).scans_file = [handles.database(patient).day(time_point).scans_file {new_file_name}];
    handles.database(patient).day(time_point).parameters = [handles.database(patient).day(time_point).parameters new_parameter_name];
    
    % Check if this new parameter exist already
    % if not update the handles.clip  fields
    if sum(strcmp(handles.clips(:,1), new_parameter_name)) == 0
        % clone the clips value
        match_clip = find(strcmp(handles.clips(:,1)', parameter_name) == 1);
        handles.clips(numel(handles.clips(:,1))+1,1) = new_parameter_name;
        handles.clips(numel(handles.clips(:,1)),2:3) =  handles.clips(match_clip,2:3); %#ok<FNDSB>
    end
else
    file_name = handles.database(patient).day(time_point).VOIs_file{scan};
    new_file_name = [spm_str_manip(file_name,'r'),'-cloned.mat'];
    parameter_name =  handles.database(patient).day(time_point).VOIs(scan);
    new_parameter_name = strcat(parameter_name, '-cloned');
    % update databas
    handles.database(patient).day(time_point).VOIs_file = [handles.database(patient).day(time_point).VOIs_file {new_file_name}];
    handles.database(patient).day(time_point).VOIs = [handles.database(patient).day(time_point).VOIs new_parameter_name];
    
    % Check if this new VOI exist already
    % if not update the handles.VOIs field
    if sum(strcmp(handles.VOIs', new_parameter_name)) == 0
        handles.VOIs =  [handles.VOIs(1:end-1) new_parameter_name handles.VOIs(end-1:end)];
    end
end
% create the new file
copyfile([handles.database(patient).path file_name], [handles.database(patient).path new_file_name], 'f');
% rename the ROI name in the uvascroi sturcture
if Scan_VOI_option == 1
    load([handles.database(patient).path new_file_name])
    for i=1:numel(uvascroi) %#ok<NODEF>
        uvascroi(i).name = char(new_parameter_name); %#ok<AGROW>
    end
    save([handles.database(patient).path new_file_name],'uvascroi');
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

MIA_menu_save_database_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function MIA_remove_omit_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_remove_omit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end

patient = get(handles.MIA_name_list, 'Value');

omit_obj = findobj(handles.MIA_name_right_click, 'Label', 'Omit');
omit = get(omit_obj, 'Checked');

if strcmp('off', omit)
    for i = 1:numel(patient)
        handles.database(patient(i)).omit = 1;
    end
else
    for i = 1:numel(patient)
        handles.database(patient(i)).omit = 0;
    end
end
guidata(hObject, handles);

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
function MIA_menu_image_information_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_image_information (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_displayed')
    return
end
if handles.mode == 1
    time_point =  handles.data_selected.patient_info.timepoint_nbr;
    patient =  handles.data_selected.patient_info.name_nbr;
else
    time_point =  handles.data_selected_for_PRM.patient_info.timepoint_nbr;
    patient =  handles.data_selected_for_PRM.patient_info.name_nbr;
end

show_anything = strcmp(get(handles.MIA_menu_image_information,'Checked'),'on');
if show_anything == 0
    set(handles.MIA_menu_image_information,'Checked','on');
else
    set(handles.MIA_menu_image_information,'Checked','off');
    if handles.mode == 1
        for i = 1:handles.data_loaded(1).number_of_scan
            strii = num2str(i);
            delete(findobj(allchild(0),'Tag',eval(['''image_infos_data' strii ''''])));
            guidata(hObject, handles);
        end
        return
    else
        for i = 1:2
            strii = num2str(i);
            delete(findobj(allchild(0),'Tag',eval(['''image_infos_data' strii ''''])));
            guidata(hObject, handles);
        end
        return
    end
end

% Display or hide image information
if handles.mode == 1
    for i = 1:handles.data_loaded(1).number_of_scan
        strii = num2str(i);
        if isfield(handles.data_selected.image(i).reco, 'paramQuantif')
            %             info_image = handles.data_selected.image(i).reco.paramQuantif;
            %             info_image = [info_image 'slice offset :', num2str(handles.data_selected.scan_z_offset(i).values,'%g ')];
            %             info_image = [info_image 'slice thickness :',num2str(handles.data_selected.image(i).reco.thickness)];
            %             info_image = [info_image 'Matrix :',[num2str(handles.data_selected.image(i).reco.no_samples) 'x' num2str(handles.data_selected.image(i).reco.no_views)]];
            %             nbr_of_line = numel(strfind(info_image, '##$'));
            
            info_image = handles.data_selected.image(i).reco.paramQuantif;
            if size(info_image,1) ==2;
                info_image = char(info_image);
                info_image(size(info_image,1)+1, 1:length(['slice offset : ' num2str(handles.data_selected.scan_z_offset(i).values,'%g ')])) =...
                    ['slice offset : ' num2str(handles.data_selected.scan_z_offset(i).values,'%g ')];
                
                info_image(size(info_image,1)+1, 1:length(['slice thickness : ' num2str(handles.data_selected.image(i).reco.thickness)])) = ['slice thickness : ' num2str(handles.data_selected.image(i).reco.thickness)];
                nbr_of_line = size(info_image, 1);
            else
                info_image = [info_image 'slice offset :', num2str(handles.data_selected.scan_z_offset(i).values,'%g ')];
                info_image = [info_image 'slice thickness :',num2str(handles.data_selected.image(i).reco.thickness)];
                nbr_of_line = numel(strfind(info_image, '##$'));
            end
            
            %             if length(info_image) > 400
            title = eval(['get(handles.MIA_data' strii '_title, ''String'');']);
            msgbox(info_image, title{:}) ;
            %             elseif nbr_of_line > 0
            %                 text(20,5+9*nbr_of_line,info_image,...
            %                     'color','red',...
            %                     'Parent', eval(strcat('handles.MIA_data', strii)),...
            %                     'Tag', strcat('image_infos_data', strii),...
            %                     'Units','Pixels',...
            %                     'FontWeight','Normal',...
            %                     'FontSize',8,...
            %                     'Interpreter','None',...
            %                     'ButtonDownFcn','');
            %             end
        elseif isfield(handles.data_selected.image(i).reco, 'iminfos') && ~isfield(handles.data_selected.image(i).reco, 'paramQuantif')
            info_image = handles.data_selected.image(i).reco.iminfos;
            info_image = [info_image 'slice offset :', num2str(handles.data_selected.scan_z_offset(i).values,'%g ')];
            info_image = [info_image 'slice thickness :',num2str(handles.data_selected.image(i).reco.thickness)];
            info_image = [info_image 'Matrix :',[num2str(handles.data_selected.image(i).reco.no_samples) 'x' num2str(handles.data_selected.image(i).reco.no_views)]];
            posx=20;
            n_info_lines = size(info_image,2);
            posy=5+15*n_info_lines;
            hdl_info_text = zeros(1,n_info_lines);
            for this_line=1:n_info_lines
                hdl_info_text(this_line)=text(posx,posy,info_image(:,this_line),...
                    'color','white',...
                    'Parent', eval(strcat('handles.MIA_data', strii)),...
                    'Tag', strcat('image_infos_data', strii),...
                    'Units','Pixels',...
                    'FontWeight','Normal',...
                    'FontSize',8,...
                    'Interpreter','None',...
                    'ButtonDownFcn','');
                posy=posy-15;
            end
        end
    end
else
    for i = 1:2
        strii = num2str(i);
        if i==1
            tp = get(handles.MIA_PRM_ref_popupmenu, 'Value');
            if tp > 1
                tp = tp-1;
            else
                tp = get(handles.MIA_PRM_slider_tp, 'Value') -1;
                if tp == 0  % case PRM_ref = -1 and slider = 1
                    tp = 1;
                end
            end
        elseif i==2
            tp = get(handles.MIA_PRM_slider_tp, 'Value');
        end
        if isfield(handles.data_selected_for_PRM_resized.image(i).reco, 'paramQuantif')
            info_image = handles.data_selected_for_PRM.image(tp).reco.paramQuantif;
            if size(info_image,1) ==2;
                info_image = char(info_image);
                info_image(size(info_image,1)+1, 1:length(['slice offset : ' num2str(handles.data_selected_for_PRM.scan_z_offset(i).values,'%g ')])) =...
                    ['slice offset : ' num2str(handles.data_selected_for_PRM.scan_z_offset(i).values,'%g ')];
                
                info_image(size(info_image,1)+1, 1:length(['slice thickness : ' num2str(handles.data_selected_for_PRM.image(i).reco.thickness)])) = ['slice thickness : ' num2str(handles.data_selected_for_PRM.image(i).reco.thickness)];
                %                 nbr_of_line = size(info_image, 1);
            else
                info_image = [info_image 'slice offset :', num2str(handles.data_selected_for_PRM.scan_z_offset(i).values,'%g ')];
                info_image = [info_image 'slice thickness :',num2str(handles.data_selected_for_PRM.image(i).reco.thickness)];
                %                 nbr_of_line = numel(strfind(info_image, '##$'));
            end
            
            
            %             if length(info_image) > 400
            title = eval(['get(handles.MIA_data' strii '_title, ''String'');']);
            msgbox(info_image, title{:}) ;
            
            %             elseif nbr_of_line > 0
            %                 text(20,5+9*nbr_of_line,info_image,...
            %                     'color','red',...
            %                     'Parent', eval(strcat('handles.MIA_data', strii)),...
            %                     'Tag', strcat('image_infos_data', strii),...
            %                     'Units','Pixels',...
            %                     'FontWeight','Normal',...
            %                     'FontSize',8,...
            %                     'Interpreter','None',...
            %                     'ButtonDownFcn','');
            %             end
        elseif isfield(handles.data_selected_for_PRM.image(tp).reco, 'iminfos') && ~isfield(handles.data_selected_for_PRM.image(tp).reco, 'paramQuantif')
            info_image = handles.data_selected_for_PRM.image(tp).reco.iminfos;
            info_image = [info_image 'slice offset :', num2str(handles.data_selected_for_PRM.scan_z_offset(i).values,'%g ')];
            info_image = [info_image 'slice thickness :',num2str(handles.data_selected_for_PRM.image(i).reco.thickness)];
            posx=20;
            n_info_lines = size(info_image,2);
            posy=5+15*n_info_lines;
            hdl_info_text = zeros(1,n_info_lines);
            for this_line=1:n_info_lines
                hdl_info_text(this_line)=text(posx,posy,info_image(:,this_line),...
                    'color','white',...
                    'Parent', eval(strcat('handles.MIA_data', strii)),...
                    'Tag', strcat('image_infos_data', strii),...
                    'Units','Pixels',...
                    'FontWeight','Normal',...
                    'FontSize',8,...
                    'Interpreter','None',...
                    'ButtonDownFcn','');
                posy=posy-15;
            end
        end
    end
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function MIA_menu_export_to_tiff_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_export_to_tiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_displayed')
    return
end
export_mode_list = {'All slices', 'Slice displayed'};

[mode_selected,ok] = listdlg('Name', 'Export Mode?', 'ListString',  export_mode_list,'ListSize', [200 150],...
    'PromptString', 'Please select the export mode');
if ok == 0
    return
end
switch mode_selected
    case 1
        slice_nbr = 1:size(handles.data_displayed.image,3);
        end_of_the_tiff_name = 'all_slices';
    case 2
        slice_nbr = get(handles.MIA_slider_slice, 'Value');
        end_of_the_tiff_name = strcat('slice-', num2str(slice_nbr));
end
if handles.mode == 1
    patient_name = handles.data_selected.patient_info.name;
    tp = handles.data_selected.patient_info.timepoint;
    tiff_path = handles.database(handles.data_selected.patient_info.name_nbr).path;
else
    patient_name = handles.data_selected_for_PRM.patient_info.name;
    tp = handles.data_selected_for_PRM.patient_info.timepoint;
    tiff_path = handles.database(handles.data_selected_for_PRM.patient_info.name_nbr).path;
end


tagstruct.ImageLength = size(handles.data_displayed.image,1);
tagstruct.ImageWidth = size(handles.data_displayed.image,2);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
% tagstruct.RowsPerStrip    = 16; %%%%
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

for i = 1:size(handles.data_displayed.image,4)
    strii = num2str(i);
    tiff_name = strcat(tiff_path, patient_name, '-',tp , '-', eval(['get(handles.MIA_data' strii '_title, ''String'')']), '-', end_of_the_tiff_name, '.tif');
    tiff_name = strrep(tiff_name, '*', 'star');
    for j = 1:numel(slice_nbr)
        if j == 1
            t = Tiff(tiff_name{:},'w');
        else
            t = Tiff(tiff_name{:},'a');
        end
        t.setTag(tagstruct);
        %                 t.write(int32(handles.data_displayed.image(:,:,slice_nbr(j),i)));
        t.write(single(handles.data_displayed.image(:,:,slice_nbr(j),i)));
        t.close();
    end
end





% --------------------------------------------------------------------
function MIA_menu_export_montage_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_export_montage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Ancien code (Celui de Benjamin)
% patient_list =get(handles.MIA_name_list, 'String');
% parameters_list = [];
% VOIs_list = [];
% time_point_list  = [];
% for i=1:numel(get(handles.MIA_name_list, 'String'))
%     for j = 1:numel(handles.database(i).day)
%         parameters_list = [parameters_list handles.database(i).day(j).parameters];
%         VOIs_list = [VOIs_list handles.database(i).day(j).VOIs];
%         time_point_list  = [time_point_list {handles.database(i).day(j).date}];
%     end
% end
% parameters_list = unique(parameters_list);
% VOIs_list = unique(VOIs_list);
% time_point_list  = unique(time_point_list);
%
%
% [patient_selected, ok1] = listdlg('PromptString','Select the patient:',...
%     'Name', 'Selecting...',...
%     'ListSize', [400 300],...
%     'ListString',patient_list');
% if ok1 == 0
%     return
% end
% if numel(patient_selected) == 1
%     [time_point_selected, ok1] = listdlg('PromptString','Select the time point:',...
%         'Name', 'Selecting...',...
%         'ListSize', [400 300],...
%         'ListString', {handles.database(patient_selected).day.date}');
% else
%     [time_point_selected, ok1] = listdlg('PromptString','Select the time point:',...
%         'Name', 'Selecting...',...
%         'ListSize', [400 300],...
%         'ListString', time_point_list');
%
% end
% if ok1 == 0
%     return
% end
% if numel(patient_selected) == 1 && numel(time_point_selected) == 1
%     [parameters_selected, ok1] = listdlg('PromptString','Select the parameters:',...
%         'Name', 'Selecting...',...
%         'ListSize', [400 300],...
%         'ListString',parameters_list');
% else
%     [parameters_selected, ok1] = listdlg('PromptString','Select the parameters:',...
%         'Name', 'Selecting...',...
%         'ListSize', [400 300],...
%         'ListString',parameters_list');
% end
% if ok1 == 0
%     return
% end
%
%
%
% if numel(patient_selected) == 1 && numel(time_point_selected) && numel(parameters_selected) == 1
%     [VOI_selected, ok1] = listdlg('PromptString','Select the VOI:',...
%         'Name', 'Selecting...',...
%         'ListSize', [400 300],...
%         'ListString',VOIs_list');
% else
%     [VOI_selected, ok1] = listdlg('PromptString','Select the VOI:',...
%         'Name', 'Selecting...',...
%         'ListSize', [400 300],...
%         'ListString',VOIs_list');
% end
% if ok1 == 0
%     return
% end
%
%
%
%
%
% [slice_selected, ok1] = listdlg('PromptString','Select the slice:',...
%     'Name', 'Selecting...',...
%     'ListSize', [400 300],...
%     'ListString',[{'All'}, {'Define by User'}]');
% if ok1 == 0
%     return
% end
% if slice_selected == 2
%     slice_selected = inputdlg('Enter the slice number separted by a space!!', 'bip', 1);
%     slice_selected = str2double(slice_selected{:});
% else
%     slice_selected = inf;
% end
%
% if strcmp(handles.database(1).path(end), filesep) == 0
%     path = strcat(handles.database(1).path, filesep);
% else
%     path = handles.database(1).path;
% end
% resolution = [256 256];
% image_for_montage = [];
% titles_for_images={};
% for i = 1:numel(patient_selected)
%     for j=1:numel(time_point_selected)
%         for x = 1:numel(parameters_selected)
%             if time_point_selected(j) > numel(handles.database(patient_selected(i)).day)
%                 continue
%             end
%             % load the the scan
%             parameter_nbr = strcmp(parameters_list(parameters_selected(x)), handles.database(patient_selected(i)).day(time_point_selected(j)).parameters) == 1;
%             scan_filename = strcat(path, handles.database(patient_selected(i)).day(time_point_selected(j)).scans_file(parameter_nbr));
%             if ~isempty(scan_filename) && sum(parameter_nbr) == 1
%                 fid=fopen(scan_filename{:} ,'r');
%                 if fid>0
%                     fclose(fid);
%                     load(scan_filename{:});
%                 else
%                     warning_text = sprintf('##$ Can not open the following map\n##$ %s\n##$',...
%                         scan_filename{:});
% %                    logbook{numel(logbook)+1} =warning_text;
%                     continue
%                 end
%                 if slice_selected == inf
%                     warndlg('Mutli slice not coded yet!!', 'Warning');
%                     return
%                 else
%                     if size(uvascim.image.reco.data, 4) >=max(slice_selected)
%                         for ii = 1:numel(slice_selected)
%                             if isempty(image_for_montage)
%                                 image_for_montage(:,:,1) = imresize((uvascim.image.reco.data(:,:,1,slice_selected(ii),1)- cell2num(handles.clips(strcmp(handles.clips(:,1), parameters_list(parameters_selected(x))) == 1,2)))/...
%                                     cell2num(handles.clips(strcmp(handles.clips(:,1), parameters_list(parameters_selected(x))) == 1,3)),  resolution, 'nearest');
%                                 titles_for_images = {[handles.database(patient_selected(i)).name '-'...
%                                     handles.database(patient_selected(i)).day(time_point_selected(j)).date '-' parameters_list{parameters_selected(x)}...
%                                     '(' num2str(cell2num(handles.clips(strcmp(parameters_list(parameters_selected(x)), handles.clips(:,1)),2)))...
%                                     '-' num2str(cell2num(handles.clips(strcmp(parameters_list(parameters_selected(x)), handles.clips(:,1)),3))) ')']};
%                                 %                                 titles_for_images = {[ handles.database(patient_selected(i)).name '-'...
%                                 %                                    parameters_list{parameters_selected(x)} ]};
%                             else
%                                 image_for_montage(:,:,size(image_for_montage,3)+1) =  imresize((uvascim.image.reco.data(:,:,1,slice_selected(ii),1)- cell2num(handles.clips(strcmp(handles.clips(:,1), parameters_list(parameters_selected(x))) == 1,2)))/...
%                                     cell2num(handles.clips(strcmp(handles.clips(:,1), parameters_list(parameters_selected(x))) == 1,3)),  resolution, 'nearest');
%                                 titles_for_images(size(titles_for_images,2)+1) ={[handles.database(patient_selected(i)).name '-'...
%                                     handles.database(patient_selected(i)).day(time_point_selected(j)).date '-' parameters_list{parameters_selected(x)}...
%                                     '(' num2str(cell2num(handles.clips(strcmp(parameters_list(parameters_selected(x)), handles.clips(:,1)),2)))...
%                                     '-' num2str(cell2num(handles.clips(strcmp(parameters_list(parameters_selected(x)), handles.clips(:,1)),3))) ')']};
%                                 %                                 titles_for_images(size(titles_for_images,2)+1) ={[handles.database(patient_selected(i)).name '-'...
%                                 %                                   parameters_list{parameters_selected(x)}]};
%                                 %
%                             end
%                         end
%                     else
%                         warning_text = sprintf('##$ No slice numbre %d in the following map\n##$ %s \n##$',...
%                             max(slice_selected), scan_filename{:});
%                         %logbook{numel(logbook)+1} =warning_text;
%                         continue
%                     end
%                 end
%             end
%
%         end
%     end
% end
% Nrows = ceil(sqrt(size(image_for_montage,3)));
% Ncols = ceil(size(image_for_montage,3)/Nrows);
% f=figure;
% h=montage(permute(image_for_montage, [1 2 4 3]),...
%     'DisplayRange', [0 1], 'Size', [Nrows Ncols]);
% x=get(h,'XData');
% y=get(h,'YData');
% titles_to_display = cell(1,Nrows*Ncols);
% if size(titles_to_display,2) > size(image_for_montage,3)
%     titles_to_display(size(image_for_montage,3)+1:end)={''};
% end
% titles_to_display(1,1:size(titles_for_images,2)) = titles_for_images;
% titles_to_display=reshape(titles_to_display, [Ncols Nrows])';
% for i = 1:Nrows
%     for j=1:Ncols
%         if isempty(isempty(titles_to_display{i,j}))
%         end
%         text((j*x(2)/Ncols)-x(2)/Ncols, (i*y(2)/Nrows)-y(2)/Nrows,titles_to_display{i,j}, 'BackgroundColor', [1 1 1])
%     end
% end



%% Nouveau code (Celui de Clement)

% Creation des listes proposees a l'utilisateur
Patients_list = get(handles.MIA_name_list, 'String');
List = [];
Param = [];
Pat = [];
Tim = [];
TpInd = [];
PatInd = [];
for i=1:length(Patients_list)
    for j=1:length(handles.database(i).day)
        List = [List strcat(Patients_list(i),'_',handles.database(i).day(j).date)];
        Pat = [Pat Patients_list(i)];
        Tim = [Tim str2num(handles.database(i).day(j).date)];
        TpInd = [TpInd j]; % Indice du timepoint relatif aux autres timepoints du meme patient
        PatInd = [PatInd i]; % Indice du patient
    end
end

% On propose a l'utilisateur de choisir un ou plusieurs couples Patient/Timepoint
[patients_selected, ok1] = listdlg('PromptString','Select the patient(s):',...
    'Name', 'Selecting...',...
    'ListSize', [400 300],...
    'ListString',List');
if ok1 == 0
    return
end


% On recupere les parametres disponibles pour les patients selectionnes. Si
% un parametre n'est disponible que pour une partie des patients
% selectionnes, il est quand meme affiche.
for i=1:length(patients_selected)
    Param = [Param handles.database(PatInd(i)).day(TpInd(i)).parameters];
end
Param = unique(Param);

% On propose a l'utilisateur de choisir un ou plusieurs parametres
[parameters_selected, ok1] = listdlg('PromptString','Select the parameter(s):',...
    'Name', 'Selecting...',...
    'ListSize', [400 300],...
    'ListString',Param');
if ok1 == 0
    return
end


% On laisse la possibilite a l'utilisateur d'afficher une ROI (ROI ou
% cluster) sur les cartes.
Clust = questdlg('Would you like to print a ROI ?',...
    'Clusters',...
    'Yes','No','Cancel','Cancel');

switch Clust
    case 'No'
        %Quelles tranches afficher ? L'utilisateur peut choisir d'afficher
        %toutes les tranches de tous les patients selectionnes (plusieurs
        %figures seront creees) ou uniquement la tranche centrale de chaque
        %patient (une seule figure est creee avec un subplot).
        [slice_selected, ok1] = listdlg('PromptString','Select the slice:',...
            'Name', 'Selecting...',...
            'ListSize', [400 300],...
            'ListString',[{'All'}, {'The central slice'}]');
        if ok1 == 0
            return
        end
        
    case 'Yes'
        %Quelles tranches afficher ? L'utilisateur peut choisir d'afficher
        %toutes les tranches de tous les patients selectionnes (plusieurs
        %figures seront creees), uniquement la tranche centrale de chaque
        %patient (une seule figure est creee avec un subplot), ou
        %uniquement la tranche ou le cluster est le plus grand pour chaque
        %patient (une seule figure est creee avec un subplot).
        [slice_selected, ok1] = listdlg('PromptString','Select the slice:',...
            'Name', 'Selecting...',...
            'ListSize', [400 300],...
            'ListString',[{'All'}, {'The central slice'},{'Biggest ROI Slice'}]');
        if ok1 == 0
            return
        end
        
    case 'Cancel'
        return
end




if slice_selected == 1 % ie si l'utilisateur a decide d'afficher toutes les tranches
    
    if strcmp(Clust,'Yes') % ie si l'utilisateur a decide d'afficher une ROI (ou cluster)
        ClusterList = [];
        %On recupere la liste des clusters disponibles a l'affichage.
        for i=1:length(patients_selected)
            ClusterList = [ClusterList handles.database(PatInd(patients_selected(i))).day.VOIs];
        end
        ClusterList = unique(ClusterList);
        
        % L'utilisateur choisit la ROI ou le cluster a afficher.
        [cluster_selected, ok1] = listdlg('PromptString','Select the cluster:',...
            'Name', 'Selecting...',...
            'ListSize', [400 300],...
            'ListString',ClusterList');
        if ok1 == 0
            return
        end
    end
    
    
    
    
    
    
    
    for i = 1:length(patients_selected)
        for j=1:length(parameters_selected)
            A = List(patients_selected(i));
            B = Param(parameters_selected(j));
            h = figure('Name',strcat(A{1},'_',B{1}),'units','normalized','outerposition',[0 0 1 1]);
            % On recupere l'indice du fichier representant le bon parametre
            [Max,Ind] = max(strcmp(Param(parameters_selected(j)),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters));
            if Max == 0 % ie si le parametre n'est pas disponible pour ce patient
                axis off
            else
                fichier = strcat(handles.database(1).databaseinfo.pathname,'MIA_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind));
                load(fichier{1});
                NbSlices = size(uvascim.image.reco.data,4);
                % En fonction du nombre de tranches a afficher, on definit
                % la taille du subplot. Il est possible d'augmenter le
                % nombre de cas ci dessous.
                if NbSlices < 13
                    lignes = 3;
                    colonnes = 4;
                elseif NbSlices < 26
                    lignes = 5;
                    colonnes = 5;
                elseif NbSlices < 50
                    lignes = 7;
                    colonnes = 7;
                elseif NbSlices < 73
                    lignes = 8;
                    colonnes = 9;
                else
                    warndlg('You cannot display more than 72 slices')
                end
                
                if strcmp(Clust,'Yes') % Si l'utilisateur a choisi d'afficher une ROI (ou cluster) (on est toujours dans le cas ou on affiche toutes les tranches du patient)
                    % On recupere l'indice du fichier representant la bonne
                    % ROI
                    [Max2,Ind2] = max(strcmp(ClusterList(cluster_selected),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).VOIs));
                    if Max2 == 0 %ie si la ROI n'est pas disponible pour ce patient
                        %Si la ROi n'est pas dicponible, on affiche les
                        %cartes sans ROI.
                        for k = 1:NbSlices
                            subplot(lignes,colonnes,k)
                            imshow(uvascim.image.reco.data(:,:,1,k,1),[min(min(uvascim.image.reco.data(:,:,1,k,1))) max(max(uvascim.image.reco.data(:,:,1,k,1)))])
                            title(strcat('Tranche : ',num2str(k)));
                        end
                    else
                        fichier2 = strcat(handles.database(1).databaseinfo.pathname,'MIA_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).VOIs_file(Ind2));
                        load(fichier2{1})
                        Offsets = {uvascroi.fov_offsets}; % On recupere les offsets des ROI, afin de les afficher sur les bonnes tranches
                        for k = 1:NbSlices
                            subplot(lignes,colonnes,k)
                            PrintFlag = 0;
                            for l=1:length(Offsets)
                                if uvascim.image.reco.fov_offsets(3,1,k) == Offsets{l}(3) % Si la tranche admet une ROI
                                    FacteurZoom = size(uvascim.image.reco.data,1)/size(uvascroi(l).value,1); % On mets les deux images (Tranche et ROI) a la meme taille
                                    %image(imresize(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),1/FacteurZoom),'CDataMapping','Scaled')
                                    
                                    %Affichage de la carte :
                                    imshow(imresize(uvascim.image.reco.data(:,:,1,k,1),1/FacteurZoom),[min(min(uvascim.image.reco.data(:,:,1,k,1))) max(max(uvascim.image.reco.data(:,:,1,k,1)))])
                                    title(strcat('Tranche : ',num2str(k)));
                                    hold on
                                    axis off
                                    if length(size(uvascroi(l).value)) == 2 % Si c'est une ROI (pas un cluster)
                                        contour(uvascroi(l).value,1,'LineColor','r') % On l'affiche
                                        axis off
                                        hold off
                                        PrintFlag = 1;
                                        
                                    elseif length(size(uvascroi(l).value)) == 3 % Si c'est un cluster (Pas une ROI)
                                        %On definit une matrice A qui
                                        %permettra de rendre transparent
                                        %tout ce qui n'est pas cluster,
                                        %afin de superposer la carte et le
                                        %cluster.
                                        A = (uvascroi(l).value == 0);
                                        A = sum(A,3);
                                        B =(A==0);
                                        A = (B==0);
                                        %image(imresize(uvascroi(k).value,4),'CDataMapping','Scaled','AlphaData',imresize(A,4));
                                        image(uvascroi(l).value,'CDataMapping','Scaled','AlphaData',A); % On l'affiche
                                        hold off
                                        PrintFlag = 1;
                                        axis off
                                    else
                                        warndlg('The data in the selected ROI has not the good size')
                                    end
                                elseif PrintFlag == 0 && l == length(Offsets)
                                    % Si on a passe toutes les tranches de
                                    % ROI sans pouvoir en afficher une,
                                    % alors que le patient admet une ROI,
                                    % alors on affiche la carte sans ROI.
                                    % Cette situation est due ? une
                                    % mauvaise correspondance entre les
                                    % offsets des tranches de la ROI et
                                    % ceux des tranches de la carte.
                                    imshow(uvascim.image.reco.data(:,:,1,k,1),[min(min(uvascim.image.reco.data(:,:,1,k,1))) max(max(uvascim.image.reco.data(:,:,1,k,1)))])
                                    title(strcat('Tranche : ',num2str(k)));
                                    axis off
                                end
                            end
                        end
                    end
                elseif strcmp(Clust,'No') % Si l'utilisateur a choisi de ne pas afficher de ROI (on est toujours dans le cas ou on affiche toutes les tranches du patient)
                    for k = 1:NbSlices
                        subplot(lignes,colonnes,k)
                        imshow(uvascim.image.reco.data(:,:,1,k,1),[min(min(uvascim.image.reco.data(:,:,1,k,1))) max(max(uvascim.image.reco.data(:,:,1,k,1)))])
                        title(strcat('Tranche : ',num2str(k)));
                    end
                end
                
            end
            
        end
    end
    
elseif slice_selected == 2 && strcmp(Clust,'No') % ie l'utilisateur a choisi d'afficher la tranche centrale et de ne pas afficher de ROI
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    
    
    for i=1:length(patients_selected)
        for j=1:length(parameters_selected)
            % On cree une figure en subplot avec une ligne par patient et
            % une colonne par parametre.
            subplot(length(patients_selected),length(parameters_selected),(i-1)*length(parameters_selected)+j)
            annotation('textbox',[0 (1-i/(length(patients_selected)+1)) 0.1 0.1],'String',List(patients_selected(i)),'Interpreter','none');
            annotation('textbox',[j/(length(parameters_selected)+1) 0.9 0.1 0.1],'String',Param(parameters_selected(j)),'Interpreter','none');
            % On recupere l'indice du fichier representant le bon parametre
            [Max,Ind] = max(strcmp(Param(parameters_selected(j)),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters));
            if Max == 0 % ie si le parametre n'est pas disponible
                axis off
            else
                fichier = strcat(handles.database(1).databaseinfo.pathname,'MIA_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind));
                load(fichier{1});
                TrancheCentrale = floor(size(uvascim.image.reco.data,4)/2);
                %image(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),'CDataMapping','Scaled')
                imshow(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),[min(min(uvascim.image.reco.data(:,:,1,TrancheCentrale,1))) max(max(uvascim.image.reco.data(:,:,1,TrancheCentrale,1)))])
                title(strcat('Tranche : ',num2str(TrancheCentrale)));
                %image(handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters,'CDataMapping','Scaled');
            end
        end
    end
    
elseif slice_selected ==2 && strcmp(Clust,'Yes') % ie l'utilisateur a choisi d'afficher la tranche centrale et d'afficher une ROI
    ClusterList = [];
    % On recupere la liste des ROI disponibles
    for i=1:length(patients_selected)
        ClusterList = [ClusterList handles.database(PatInd(patients_selected(i))).day.VOIs];
    end
    ClusterList = unique(ClusterList);
    
    %L'utilisateur choisit la ROI qu'il veut afficher
    [cluster_selected, ok1] = listdlg('PromptString','Select the cluster:',...
        'Name', 'Selecting...',...
        'ListSize', [400 300],...
        'ListString',ClusterList');
    if ok1 == 0
        return
    end
    
    
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    
    
    for i=1:length(patients_selected)
        for j=1:length(parameters_selected)
            subplot(length(patients_selected),length(parameters_selected),(i-1)*length(parameters_selected)+j)
            annotation('textbox',[0 (1-i/(length(patients_selected)+1)) 0.1 0.1],'String',List(patients_selected(i)),'Interpreter','none');%Patients_list(patients_selected(i)));
            annotation('textbox',[j/(length(parameters_selected)+1) 0.9 0.1 0.1],'String',Param(parameters_selected(j)),'Interpreter','none');
            % On recupere l'indice du fichier representant le bon parametre
            [Max,Ind] = max(strcmp(Param(parameters_selected(j)),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters));
            if Max == 0
                axis off
            else
                fichier = strcat(handles.database(1).databaseinfo.pathname,'MIA_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind));
                load(fichier{1});
                TrancheCentrale = floor(size(uvascim.image.reco.data,4)/2);
                % On recupere l'indice du fichier representant la bonne ROI
                [Max2,Ind2] = max(strcmp(ClusterList(cluster_selected),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).VOIs));
                if Max2 == 0 %Si la ROI n'est pas disponible pour ce patient
                    %image(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),'CDataMapping','Scaled')
                    imshow(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),[min(min(uvascim.image.reco.data(:,:,1,TrancheCentrale,1))) max(max(uvascim.image.reco.data(:,:,1,TrancheCentrale,1)))])
                    title(strcat('Tranche : ',num2str(TrancheCentrale)));
                else
                    fichier2 = strcat(handles.database(1).databaseinfo.pathname,'MIA_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).VOIs_file(Ind2));
                    load(fichier2{1})
                    Offsets = {uvascroi.fov_offsets}; % On recupere les offsets des ROI, afin de les afficher sur les bonnes tranches
                    PrintFlag = 0;
                    for k=1:length(Offsets)
                        if uvascim.image.reco.fov_offsets(3,1,TrancheCentrale) == Offsets{k}(3)
                            FacteurZoom = size(uvascim.image.reco.data,1)/size(uvascroi(TrancheCentrale).value,1); % On mets les deux images (Tranche et ROI) a la meme taille
                            %image(imresize(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),1/FacteurZoom),'CDataMapping','Scaled')
                            
                            %On affiche la carte
                            imshow(imresize(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),1/FacteurZoom),[min(min(uvascim.image.reco.data(:,:,1,TrancheCentrale,1))) max(max(uvascim.image.reco.data(:,:,1,TrancheCentrale,1)))])
                            title(strcat('Tranche : ',num2str(TrancheCentrale)));
                            axis off
                            %title(strcat(handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind),'Tranche : ',num2str(TrancheCentrale)));
                            hold on
                            if length(size(uvascroi(k).value)) == 2 % Si c'est une ROI (pas un cluster)
                                contour(uvascroi(k).value,1,'LineColor','r') %On l'affiche
                                axis off
                                hold off
                                PrintFlag = 1;
                                
                            elseif length(size(uvascroi(k).value)) == 3 % Si c'est un cluster (pas une ROI)
                                %On definit une matrice A qui permettra de
                                %rendre transparent tout ce qui n'est pas
                                %cluster, afin de superposer la carte et le
                                %cluster.
                                A = (uvascroi(k).value == 0);
                                A = sum(A,3);
                                B =(A==0);
                                A = (B==0);
                                %image(imresize(uvascroi(k).value,4),'CDataMapping','Scaled','AlphaData',imresize(A,4));
                                image(uvascroi(k).value,'CDataMapping','Scaled','AlphaData',A);
                                hold off
                                PrintFlag = 1;
                                axis off
                            else
                                warndlg('The data in the selected ROI has not the good size')
                            end
                        elseif PrintFlag == 0 && k == length(Offsets)
                            % Si on a passe toutes les tranches de
                            % ROI sans pouvoir en afficher une,
                            % alors que le patient admet une ROI,
                            % alors on affiche la carte sans ROI.
                            % Cette situation est due ? une
                            % mauvaise correspondance entre les
                            % offsets des tranches de la ROI et
                            % ceux des tranches de la carte.
                            %image(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),'CDataMapping','Scaled')
                            imshow(uvascim.image.reco.data(:,:,1,TrancheCentrale,1),[min(min(uvascim.image.reco.data(:,:,1,TrancheCentrale,1))) max(max(uvascim.image.reco.data(:,:,1,TrancheCentrale,1)))])
                            title(strcat('Tranche : ',num2str(TrancheCentrale)));
                            %title(strcat(handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind),'Tranche : ',num2str(TrancheCentrale)));
                        end
                    end
                end
            end
        end
    end
    
    
elseif slice_selected == 3 % ie l'utilisateur a choisi d'afficher la tranche ou le cluster est le plus gros (sous entendu il a egalement choisi d'afficher une ROI/cluster)
    
    ClusterList = [];
    %On recupere la liste des ROI disponibles a l'affichage
    for i=1:length(patients_selected)
        ClusterList = [ClusterList handles.database(PatInd(patients_selected(i))).day.VOIs];
    end
    ClusterList = unique(ClusterList);
    
    %On propose a l'utilisateur de choisir la ROI a afficher
    [cluster_selected, ok1] = listdlg('PromptString','Select the cluster:',...
        'Name', 'Selecting...',...
        'ListSize', [400 300],...
        'ListString',ClusterList');
    if ok1 == 0
        return
    end
    
    
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    
    
    
    for i=1:length(patients_selected)
        for j=1:length(parameters_selected)
            subplot(length(patients_selected),length(parameters_selected),(i-1)*length(parameters_selected)+j)
            annotation('textbox',[0 (1-i/(length(patients_selected)+1)) 0.1 0.1],'String',List(patients_selected(i)),'Interpreter','none');%Patients_list(patients_selected(i)));
            annotation('textbox',[j/(length(parameters_selected)+1) 0.9 0.1 0.1],'String',Param(parameters_selected(j)),'Interpreter','none');
            % On recupere l'indice du fichier representant la bonne ROI
            [Max2,Ind2] = max(strcmp(ClusterList(cluster_selected),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).VOIs));
            if Max2 ==0 % Si la ROI n'est pas disponible pour ce patient, on va juste afficher la carte
                % On recupere l'indice du fichier representant le bon
                % parametre
                [Max3,Ind3] = max(strcmp(Param(parameters_selected(j)),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters));
                if Max3 == 0 %Si le parametre n'est pas disponible pour ce patient
                    axis off
                else
                    fichier = strcat(handles.database(1).databaseinfo.pathname,'MIA_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind3));
                    load(fichier{1});
                    Tranche = floor(size(uvascim.image.reco.data,4)/2);
                    %image(uvascim.image.reco.data(:,:,1,Tranche,1),'CDataMapping','Scaled')
                    imshow(uvascim.image.reco.data(:,:,1,Tranche,1),[min(min(uvascim.image.reco.data(:,:,1,Tranche,1))) max(max(uvascim.image.reco.data(:,:,1,Tranche,1)))])
                    title(strcat('Tranche :',num2str(Tranche)))
                    axis off
                end
            else % Si la ROI est disponible
                fichier2 = strcat(handles.database(1).databaseinfo.pathname,'Image_Analyses_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).VOIs_file(Ind2));
                load(fichier2{1})
                area = zeros(length(uvascroi),1);
                for l=1:length(uvascroi)
                    area(l) = uvascroi(l).area_in_pixel;
                end
                
                % On recupere l'indice du fichier representant le bon
                % parametre
                [Max,Ind] = max(strcmp(Param(parameters_selected(j)),handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters));
                if Max == 0 %Si le parametre n'est pas disponible pour ce patient
                    axis off
                else
                    fichier = strcat(handles.database(1).databaseinfo.pathname,'Image_Analyses_data/',handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind));
                    load(fichier{1});
                    
                    PrintFlag =0 ;
                    %On recupere l'indice de la tranche (de la ROI) qui
                    %admet la plus grande ROI.
                    [Max,Ind3] = max(area);
                    offset = uvascroi(Ind3).fov_offsets(3); % On recupere les offsets des ROI, afin de les afficher sur les bonnes tranches
                    for ll=1:length(uvascim.image.reco.fov_offsets(3,1,:,1))
                        Tranche = 0;
                        if uvascim.image.reco.fov_offsets(3,1,ll) == offset
                            FacteurZoom = size(uvascim.image.reco.data,1)/size(uvascroi(Ind3).value,1);  % On mets les deux images (Tranche et ROI) a la meme taille
                            %image(imresize(uvascim.image.reco.data(:,:,1,ll,1),1/FacteurZoom),'CDataMapping','Scaled')
                            %On affiche la carte
                            imshow(imresize(uvascim.image.reco.data(:,:,1,ll,1),1/FacteurZoom),[min(min(uvascim.image.reco.data(:,:,1,ll,1))) max(max(uvascim.image.reco.data(:,:,1,ll,1)))])
                            %title(strcat(handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind),'Tranche : ',num2str(ll)));
                            title(strcat('Tranche :',num2str(ll)))
                            axis off
                            hold on
                            if length(size(uvascroi(Ind3).value)) == 2    %C'est une ROI (Pas un cluster)
                                contour(uvascroi(Ind3).value,1,'LineColor','r') % On l'affiche
                                axis off
                                hold off
                                PrintFlag = 1;
                                Tranche = ll;
                            elseif length(size(uvascroi(Ind3).value)) == 3   %C'est un Cluster (Pas une ROI)
                                %On definit une matrice A qui permettra de
                                %rendre transparent tout ce qui n'est pas
                                %cluster, afin de superposer la carte et le
                                %cluster.
                                A = (uvascroi(Ind3).value == 0);
                                A = sum(A,3);
                                B =(A==0);
                                A = (B==0);
                                %image(imresize(uvascroi(k).value,4),'CDataMapping','Scaled','AlphaData',imresize(A,4));
                                image(uvascroi(Ind3).value,'CDataMapping','Scaled','AlphaData',A);
                                hold off
                                PrintFlag = 1;
                                Tranche = ll;
                            else
                                warndlg('The data in the selected ROI has not the good size')
                            end
                            %image(handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).parameters,'CDataMapping','Scaled');
                        elseif PrintFlag == 0 && ll == length(uvascim.image.reco.fov_offsets(3,1,:))
                            % Si on a passe toutes les tranches de
                            % ROI sans pouvoir en afficher une,
                            % alors que le patient admet une ROI,
                            % alors on affiche la carte sans ROI.
                            % Cette situation est due ? une
                            % mauvaise correspondance entre les
                            % offsets des tranches de la ROI et
                            % ceux des tranches de la carte.
                            if Tranche == 0
                                Tranche = floor(size(uvascim.image.reco.data,4)/2);
                            end
                            %image(uvascim.image.reco.data(:,:,1,Tranche,1),'CDataMapping','Scaled')
                            imshow(uvascim.image.reco.data(:,:,1,Tranche,1),[min(min(uvascim.image.reco.data(:,:,1,Tranche,1))) max(max(uvascim.image.reco.data(:,:,1,Tranche,1)))])
                            title(strcat('Tranche :',num2str(Tranche)))
                            axis off
                            %title(strcat(handles.database(PatInd(patients_selected(i))).day(TpInd(patients_selected(i))).scans_file(Ind),'Tranche : ',num2str(Tranche)));
                        end
                        
                    end
                end
            end
        end
    end
    
end








% --------------------------------------------------------------------
function MIA_math_module_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_math_module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'database')
    return
end

% define all the scan available
% parameters_list = [];
% VOIs_list = [];
% for i=1:numel(get(handles.MIA_name_list, 'String'))
%     for j = 1:numel(handles.database(i).day)
%         parameters_list = [parameters_list handles.database(i).day(j).parameters];
%         VOIs_list = [VOIs_list handles.database(i).day(j).VOIs];
%     end
% end

rows = handles.database.Type == 'Scan';
vars = {'SequenceName'};
new_parameters_list = handles.database{rows, vars}; % Ou alors handles.database(rows, :) ?

new_VOIs_list = handles.database{handles.database.Type == 'ROI', {'SequenceName'}}; %Ou alors handles.database(handles.database.Type == 'ROI', :) ?


%parameters_list = unique(parameters_list);
%VOIs_list = unique(VOIs_list);
%Math_module(handles, parameters_list, VOIs_list)
Math_module(handles, char(new_parameters_list), char(new_VOIs_list))


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
switch selection,
    case 'Parameters'
        listing = [];
        for i=1:numel(get(handles.MIA_name_list, 'String'))
            for j = 1:numel(handles.database(i).day)
                listing = [listing handles.database(i).day(j).parameters];
            end
        end
        listing = unique(listing);
    case 'VOIs'
        listing = [];
        for i=1:numel(get(handles.MIA_name_list, 'String'))
            for j = 1:numel(handles.database(i).day)
                listing = [listing handles.database(i).day(j).VOIs]; %#ok<AGROW>
            end
        end
        listing = unique(listing);
        if isempty(listing)
            listing = {''};
        end
    case 'Cancel'
        return
end

[list_to_remove,ok] = listdlg('Name', 'Question?', 'ListString', listing','ListSize', [300 400],...
    'PromptString', 'Select the parameters/VOIs you want to remove');
if ok == 0 || isempty(list_to_remove)
    return
end

user_response = questdlg('Do you want to delete the corresponding file(s) too (.mat)??', 'Warning', 'Yes', 'No', 'Cancel', 'Cancel');
if strcmp(user_response, 'Cancel')
    return
end
if strcmp(user_response, 'Yes')
    user_response = questdlg('Do you REALLY want to delete the corresponding file(s) (.mat)??', 'Warning', 'No', 'Yes', 'Cancel', 'Cancel');
end
% define all the scan available
for i = 1:numel(handles.database)
    for j = 1:numel(handles.database(i).day)
        for ii = 1:numel(list_to_remove)
            switch selection,
                case 'Parameters'
                    if strcmp(user_response, 'Yes')
                        if sum(strcmp(handles.database(i).day(j).parameters, listing{list_to_remove(ii)}))>0
                            old_path = pwd;
                            cd(handles.database(i).path);
                            delete(handles.database(i).day(j).scans_file{strcmp(handles.database(i).day(j).parameters, listing{list_to_remove(ii)})});
                            cd(old_path);
                        end
                    end
                    handles.database(i).day(j).scans_file(strcmp(handles.database(i).day(j).parameters, listing{list_to_remove(ii)})) = [];
                    handles.database(i).day(j).parameters(strcmp(handles.database(i).day(j).parameters, listing{list_to_remove(ii)})) = [];
                case 'VOIs'
                    if strcmp(user_response, 'Yes')
                        if sum(strcmp(handles.database(i).day(j).VOIs_file, listing{list_to_remove(ii)}))>0
                            old_path = pwd;
                            cd(handles.database(i).path);
                            delete(handles.database(i).day(j).VOIs_file{strcmp(handles.database(i).day(j).parameters, listing{list_to_remove(ii)})});
                            cd(old_path);
                        end
                    end
                    handles.database(i).day(j).VOIs_file(strcmp(handles.database(i).day(j).VOIs, listing{list_to_remove(ii)})) = [];
                    handles.database(i).day(j).VOIs(strcmp(handles.database(i).day(j).VOIs, listing{list_to_remove(ii)})) = [];
            end
        end
    end
end
for ii = 1:numel(list_to_remove)
    switch selection
        case 'Parameters'
            handles.clips(strcmp(handles.clips(:,1),listing{list_to_remove(ii)}),:) = [];
        case 'VOIs'
            handles.VOIs(strcmp(handles.VOIs,listing{list_to_remove(ii)})) = [];
    end
end

guidata(hObject, handles)
MIA_update_database_display(hObject, eventdata, handles)
msgbox('Done', 'Message') ;


% --------------------------------------------------------------------
function MIA_menu_define_mask_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_define_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_displayed')
    msgbox('Please load 1 scan', 'Message') ;
    return
end
if handles.data_loaded(1).number_of_scan ~=1
    msgbox('Please load 1 scan only', 'Message') ;
    return
end
if strcmp(get(handles.MIA_menu_define_mask, 'Checked'), 'off')
    set(handles.MIA_menu_define_mask, 'Checked', 'on')
    handles.disk_strel = strel('disk',5);
    handles.image_resize_original = handles.data_selected_resized.image.reco.data;
    guidata(hObject, handles)
else
    
    listing = handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).parameters;
    [scan_list,ok] = listdlg('Name', 'Question?', 'ListString', listing','ListSize', [300 400],...
        'PromptString', 'Select select the scans to mask');
    if ok == 0 || isempty(scan_list)
        handles.data_selected_resized.image.reco.data = handles.image_resize_original;
        if isfield(handles, 'color_mask3d')
            handles= rmfield(handles, 'color_mask3d');
            handles= rmfield(handles, 'disk_strel');
            handles= rmfield(handles, 'image_resize_original');
        end
        set(handles.MIA_menu_define_mask, 'Checked', 'off')
        guidata(hObject, handles)
        MIA_update_axes(hObject, eventdata, handles)
        return
    end
    
    selection_overwrite = questdlg('Do you want to overwrite the files (.mat) with the scan masked?',...
        'Warning',...
        'Yes','No', 'Cancel', 'Cancel');
    if strcmp(selection_overwrite, 'Yes')
        selection_overwrite = questdlg('Do you REALLY want to overwrite the files(.mat) with the scan masked?',...
            'Warning',...
            'No','Yes', 'Cancel', 'Cancel');
    end
    mask = double(squeeze(handles.color_mask3d(:,:,1,:,1)));
    mask(mask==0) = NaN;
    logbook = {};
    for i = 1:numel(scan_list)
        load(fullfile(handles.database(handles.data_selected.patient_info.name_nbr).path,...
            handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{scan_list(i)}));
        [x, y, echo, slice, repet]=size(uvascim.image.reco.data);
        if slice ~= size(mask, 3)
            logbook{numel(logbook)+1} =[handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{scan_list(i)}...
                '-mask -Skip (not the same slice nbr)'];
            continue
        end
        mask3d_resized=  imresize(mask, [x, y],'bilinear');
        mask5d = repmat(mask3d_resized,[1 1 1 echo repet]);
        mask5d = permute(mask5d, [1 2 4 3 5]);
        
        uvascim.image.reco.data=uvascim.image.reco.data.* mask5d;
        switch selection_overwrite
            case 'Yes'
                % save file
                [path, name, ext] = fileparts(fullfile(handles.database(handles.data_selected.patient_info.name_nbr).path,...
                    handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{scan_list(i)}));
                scan_name = [name ext];
                save(fullfile([path, scan_name]),'uvascim');
                MIA_update_database_display(hObject, eventdata, handles);
                clear uvascim;
            case 'No'
                % save file
                [path, name, ext] = fileparts(fullfile(handles.database(handles.data_selected.patient_info.name_nbr).path,...
                    handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{scan_list(i)}));
                scan_name = [name '-masked' ext];
                save(fullfile(path, scan_name),'uvascim');
                % add to the database and update handles.clips
                handles =  MIA_add_new_scan_to_database(handles, uvascim.image, handles.data_selected.patient_info.name_nbr, handles.data_selected.patient_info.timepoint_nbr,...
                    [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).parameters{scan_list(i)} '-masked'], [name '-masked' ext]);
                guidata(hObject, handles)
                MIA_update_database_display(hObject, eventdata, handles);
                clear uvascim;
            case 'Cancel'
                handles= rmfield(handles, 'color_mask3d');
                set(handles.MIA_menu_define_mask, 'Checked', 'off')
                handles= rmfield(handles, 'disk_strel');
                handles.data_selected_resized.image.reco.data = handles.image_resize_original;
                handles= rmfield(handles, 'image_resize_original');
                guidata(hObject, handles)
                MIA_update_axes(hObject, eventdata, handles)
                return
        end
        logbook{numel(logbook)+1} =[handles.database(handles.data_selected.patient_info.name_nbr).name...
            handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).date...
            [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{scan_list(i)}...
            '-mask'] '-Done'];
    end
    
    handles= rmfield(handles, 'image_resize_original');
    handles= rmfield(handles, 'color_mask3d');
    set(handles.MIA_menu_define_mask, 'Checked', 'off')
    handles= rmfield(handles, 'disk_strel');
    guidata(hObject, handles)
    if ~isempty(logbook)
        listdlg('ListString', logbook', 'ListSize',[400 500], 'Name', 'logbook');
    else
        msgbox('Done', 'logbook') ;
    end
    return
end

% retrieve value of clip for min and max of slider
ListClips = handles.clips;
SelectedScan = handles.data_selected.list_scan;
IndxScanClip = strcmp(SelectedScan,ListClips);
[min,maxval] = ListClips{IndxScanClip,2:3};

% Ask for value of min and max
prompt = {'Min','Max'};
dlg_title = 'Please initate the min and max values?';
num_lines = 1;
def = {sprintf('%d',min),sprintf('%d',maxval)};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer = inputdlg(prompt,dlg_title,num_lines,def, options);
if isempty(answer)
    return
end

min = str2double(answer{1});
maxval =str2double(answer{2});

set(handles.MIA_PRM_slider_trans,'Value',(maxval-min)/2);
set(handles.MIA_PRM_slider_trans,'max',maxval);
set(handles.MIA_PRM_slider_trans,'min',min);
MIA_mask(hObject, eventdata, handles, (maxval-min)/2)

function MIA_mask(hObject, eventdata, handles, tolerance)

if isfield(handles, 'color_mask3d')
    handles= rmfield(handles, 'color_mask3d');
end
[x, y, echo, slice, repet]=size(handles.data_selected_resized.image.reco.data);

n_echoes = get(handles.MIA_data1_echo_slider, 'Value');
n_slices = get(handles.MIA_slider_slice, 'Value');
n_expts = get(handles.MIA_data1_expt_slider, 'Value');
color_mask3d = false([x y 1 n_slices 1]);
for j = 1:slice,
    im = handles.image_resize_original(:,:,n_echoes,j,n_expts);
    % coordinates of a list of reference point
    %     list_ind = 1;
    %Create the binary mask
    g = double(im);
    %     max_ref = max(g(list_ind));
    max_ref =  min(g(:));
    color_mask = (g-max_ref) <= tolerance;
    % Connected component labelling and brain selection
    color_mask = ~color_mask;
    label = bwlabel(color_mask, 4);
    STATS = regionprops(label, 'Area'); %#ok<MRPBW>
    
    [taille_du_plus_grans_label,plus_grand_label] = max([STATS.Area]); %#ok<ASGLU>
    
    if ~isempty(plus_grand_label)
        color_mask (label ~= plus_grand_label) = 0;
    end;
    color_mask = imerode(color_mask,handles.disk_strel);
    color_mask = imdilate(color_mask,handles.disk_strel);
    warning('off');
    color_mask = imclose(color_mask,handles.disk_strel);
    warning('on'); %#ok<WNON>
    handles.color_mask3d(:,:,1,j,1) = color_mask;
end
if length(size(handles.image_resize_original)) == 4
    color_mask5d = repmat(handles.color_mask3d,[1 1 echo 1]);
else
    color_mask5d = repmat(handles.color_mask3d,[1 1 echo 1 repet]);
end

handles.data_selected_resized.image.reco.data=handles.image_resize_original .* color_mask5d;
guidata(hObject, handles)
MIA_update_axes(hObject, eventdata, handles)

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

handles = guidata(hObject);
if ~isfield(handles, 'database')
    return
end
namelist = {handles.database(:).name};
namelist_up = sort(namelist);
for i = 1:numel(handles.database)
    match = strmatch(namelist_up(i), {handles.database(:).name}', 'exact');
    database_tmp(i) = handles.database(match);
end
handles.database = database_tmp;

guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);

% --------------------------------------------------------------------
function MIA_sort_name_down_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_name_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

if ~isfield(handles, 'database')
    return
end
namelist = {handles.database(:).name};
namelist_up = sort(namelist);
patient_number = numel(handles.database)+1;
for i = 1:numel(handles.database)
    match = strmatch(namelist_up(i), {handles.database(:).name}', 'exact');
    database_tmp(patient_number-i) = handles.database(match);
end
handles.database = database_tmp;

guidata(hObject, handles);
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

handles = guidata(hObject);
if ~isfield(handles, 'database')
    return
end
patient = get(handles.MIA_name_list, 'Value');
time_point_list = {handles.database(patient).day.date};
time_point_list_up = sort(time_point_list);
for i = 1:numel(handles.database(patient).day)
    match = strmatch(time_point_list_up(i), {handles.database(patient).day.date}', 'exact');
    database_tmp(i) = handles.database(patient).day(match);
end
handles.database(patient).day = database_tmp;

guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);

% --------------------------------------------------------------------
function MIA_sort_time_point_down_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_time_point_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

if ~isfield(handles, 'database')
    return
end
patient = get(handles.MIA_name_list, 'Value');
time_point_list = {handles.database(patient).day.date};
time_point_list_up = sort(time_point_list);
time_point_number = numel(handles.database(patient).day)+1;
for i = 1:numel(handles.database(patient).day)
    match = strmatch(time_point_list_up(i), {handles.database(patient).day.date}', 'exact');
    database_tmp(time_point_number-i) = handles.database(patient).day(match);
end
handles.database(patient).day = database_tmp;

guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);

% --------------------------------------------------------------------
function MIA_sort_scan_up_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_sort_scan_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = guidata(hObject);
if ~isfield(handles, 'database')
    return
end
patient = get(handles.MIA_name_list, 'Value');
time_point = get(handles.MIA_time_points_list, 'Value');
if get(handles.MIA_scan_VOIs_button, 'Value') == 0
    parameters_list = handles.database(patient).day(time_point).parameters;
else
    parameters_list = handles.database(patient).day(time_point).VOIs;
end
parameters_list_up = sort(parameters_list);


if get(handles.MIA_scan_VOIs_button, 'Value') == 0
    for i = 1:numel(handles.database(patient).day(time_point).parameters)
        match = strmatch(parameters_list_up(i), handles.database(patient).day(time_point).parameters', 'exact');
        if numel(match) > 1
            MIA_warning_duplicate_scan_fcn(handles, parameters_list_up(i),patient,time_point)
            return
        end
        database_tmp_parameters(i) = handles.database(patient).day(time_point).parameters(match);
        database_tmp_scans_file(i)= handles.database(patient).day(time_point).scans_file(match);
    end
else
    for i = 1:numel(handles.database(patient).day(time_point).VOIs)
        match = strmatch(parameters_list_up(i), handles.database(patient).day(time_point).VOIs', 'exact');
        if numel(match) > 1
            MIA_warning_duplicate_VOI_fcn(handles, parameters_list_up(i),patient,time_point)
            return
        end
        database_tmp_VOIs(i) = handles.database(patient).day(time_point).VOIs(match);
        database_tmp_VOIs_file(i)= handles.database(patient).day(time_point).VOIs_file(match);
    end
end

if get(handles.MIA_scan_VOIs_button, 'Value') == 0
    handles.database(patient).day(time_point).parameters = database_tmp_parameters;
    handles.database(patient).day(time_point).scans_file = database_tmp_scans_file;
else
    handles.database(patient).day(time_point).VOIs =database_tmp_VOIs;
    handles.database(patient).day(time_point).VOIs_file = database_tmp_VOIs_file;
end

guidata(hObject, handles);
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
patient = get(handles.MIA_name_list, 'Value');
time_point = get(handles.MIA_time_points_list, 'Value');
if get(handles.MIA_scan_VOIs_button, 'Value') == 0
    parameters_list = handles.database(patient).day(time_point).parameters;
else
    parameters_list = handles.database(patient).day(time_point).VOIs;
end
parameters_list_up = sort(parameters_list);

time_point_number = length(parameters_list_up)+1;
if get(handles.MIA_scan_VOIs_button, 'Value') == 0
    for i = 1:numel(handles.database(patient).day(time_point).parameters)
        match = strmatch(parameters_list_up(i), handles.database(patient).day(time_point).parameters', 'exact');
        database_tmp_parameters(time_point_number-i) = handles.database(patient).day(time_point).parameters(match);
        database_tmp_scans_file(time_point_number-i)= handles.database(patient).day(time_point).scans_file(match);
    end
else
    for i = 1:numel(handles.database(patient).day(time_point).VOIs)
        match = strmatch(parameters_list_up(i), handles.database(patient).day(time_point).VOIs', 'exact');
        database_tmp_VOIs(time_point_number-i) = handles.database(patient).day(time_point).VOIs(match);
        database_tmp_VOIs_file(time_point_number-i)= handles.database(patient).day(time_point).VOIs_file(match);
    end
end

if get(handles.MIA_scan_VOIs_button, 'Value') == 0
    handles.database(patient).day(time_point).parameters = database_tmp_parameters;
    handles.database(patient).day(time_point).scans_file = database_tmp_scans_file;
else
    handles.database(patient).day(time_point).VOIs =database_tmp_VOIs;
    handles.database(patient).day(time_point).VOIs_file = database_tmp_VOIs_file;
end

guidata(hObject, handles);
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
switch selection,
    case 'Parameters'
        listing = [];
        for i=1:numel(get(handles.MIA_name_list, 'String'))
            for j = 1:numel(handles.database(i).day)
                listing = [listing handles.database(i).day(j).parameters];
            end
        end
        listing = unique(listing);
    case 'VOIs'
        listing = [];
        for i=1:numel(get(handles.MIA_name_list, 'String'))
            for j = 1:numel(handles.database(i).day)
                listing = [listing handles.database(i).day(j).VOIs]; %#ok<AGROW>
            end
        end
        listing = unique(listing);
        if isempty(listing)
            listing = {''};
        end
    case 'Cancel'
        return
end

[list_to_rename,ok] = listdlg('Name', 'Question?', 'ListString', listing','SelectionMode','single','ListSize', [500 500],...
    'PromptString', 'Select the parameters/VOIs you want to rename');
if ok == 0 || isempty(list_to_rename)
    return
end

name_option = [listing 'Other'];
old_scan_name = listing(list_to_rename);

[new_scan_name, ok1] = listdlg('PromptString',sprintf('Select the new name for replace %s',old_scan_name{:}),...
    'Name', 'Select a Name',...
    'SelectionMode','single',...
    'ListSize', [500 500],...
    'ListString',name_option);

if ok1 == 0
    return
end

switch selection
    case 'VOIs'
        if strcmp('Other', name_option(new_scan_name)) == 1
            newparameter = inputdlg('Name of the new VOI ', 'Question?', 1, {''});
        else
            newparameter = name_option(new_scan_name);
        end
    case 'Parameters'
        if strcmp('Other',name_option(new_scan_name)) == 1
            newparameter = inputdlg('Name of the new Scan ', 'Question?', 1, {''});
            if isempty(newparameter)
                return
            end
            if sum(strcmp(newparameter, handles.clips(:,1))) == 0
                handles.clips(size(handles.clips,1)+1,1) = newparameter;
                handles.clips(size(handles.clips,1),2) = handles.clips(strcmp(old_scan_name, handles.clips(:,1)),2);
                handles.clips(size(handles.clips,1),3) = handles.clips(strcmp(old_scan_name, handles.clips(:,1)),3);
            end
        else
            newparameter = name_option(new_scan_name);
        end
end

for patient = 1:numel(handles.database)
    for timepoint = 1:numel(handles.database(patient).day)
        for ii = 1:numel(list_to_rename)
            % rename the scan file
            switch selection
                case 'VOIs'
                    id_param = strcmp(handles.database(patient).day(timepoint).VOIs, listing{list_to_rename(ii)});
                    if sum(id_param)>0
                        [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).VOIs_file{id_param}]);
                        new_name = [handles.database(patient).name '-' handles.database(patient).day(timepoint).date '-' newparameter{:}];
                        if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
                            warning_text = sprintf('##$ This file no not exist\n##$ %s',...
                                fullfile(PATHSTR,[NAME,EXT]));
                            msgbox(warning_text, 'rename file warning') ;
                        elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
                            user_response = questdlg('The new file name for this ROI exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
                            if strcmp(user_response, 'No')
                                return
                            end
                            user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
                            % overwrite the file if requested
                            if strcmp(user_response, 'Yes')
                                if ~strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]))
                                    movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                                    handles.database(patient).day(timepoint).VOIs_file{id_param} = [new_name,EXT];
                                    handles.database(patient).day(timepoint).VOIs{id_param} = newparameter{1};
                                end
                                %% Add for modify name in ROI struct
                                load([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{id_param}])
                                for itVoi = 1 : numel(uvascroi)
                                    uvascroi(itVoi).name = newparameter{1};
                                end
                                save([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{id_param}],'uvascroi')
                            end
                            
                        else
                            movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                            handles.database(patient).day(timepoint).VOIs_file{id_param} = [new_name,EXT];
                            handles.database(patient).day(timepoint).VOIs{id_param} = newparameter{1};
                            %% Add for modify name in ROI
                            load([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{id_param}])
                            for itRoi = 1 : numel(uvascroi)
                                uvascroi(itRoi).name = newparameter{1};
                            end
                            save([PATHSTR filesep handles.database(patient).day(timepoint).VOIs_file{id_param}],'uvascroi')
                        end
                    end
                case 'Parameters'
                    id_param = strcmp(handles.database(patient).day(timepoint).parameters, listing{list_to_rename(ii)});
                    if sum(id_param)>0
                        [PATHSTR,NAME,EXT] = fileparts([handles.database(patient).path   handles.database(patient).day(timepoint).scans_file{id_param}]);
                        new_name = [handles.database(patient).name '-' handles.database(patient).day(timepoint).date '-' newparameter{:}];
                        new_name = strrep(new_name, '*', 'star');
                        if  exist(fullfile(PATHSTR,[NAME,EXT]), 'file') == 0
                            warning_text = sprintf('##$ This file no not exist\n##$ %s',...
                                fullfile(PATHSTR,[NAME,EXT]));
                            msgbox(warning_text, 'rename file warning') ;
                        elseif  exist(fullfile(PATHSTR,[new_name,EXT]), 'file') == 2
                            user_response = questdlg('The new file name exist already for this patient/time point. Do you want to overwrite this file or no?','Warning', 'Yes', 'No', 'No');
                            if strcmp(user_response, 'No')
                                return
                            end
                            user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
                            % overwrite the file if requested
                            if strcmp(user_response, 'Yes')
                                if ~strcmp(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]))
                                    movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                                    handles.database(patient).day(timepoint).scans_file{id_param} = [new_name,EXT];
                                    handles.database(patient).day(timepoint).parameters{id_param} = newparameter{1};
                                end
                            end
                        else
                            movefile(fullfile(PATHSTR,[NAME,EXT]), fullfile(PATHSTR,[new_name,EXT]), 'f');
                            handles.database(patient).day(timepoint).scans_file{id_param} = [new_name,EXT];
                            handles.database(patient).day(timepoint).parameters{id_param} = newparameter{1};
                        end
                    end
            end
        end
    end
end

guidata(hObject, handles);

% Add new name to handles.VOIs
if( strcmp('VOIs',selection) && strcmp('Other', handles.VOIs(new_scan_name)) )
    handles.VOIs=[handles.VOIs(1:end-1), newparameter, handles.VOIs(end-1:end)];
end

% check if the old VOIs exist somewhere else, if not remove it
vois_list = [];
for i=1:numel(handles.database)
    for j = 1:numel(handles.database(i).day)
        vois_list = [vois_list handles.database(i).day(j).VOIs]; %#ok<AGROW>
        
    end
end
vois_list= unique(vois_list);

if sum(strcmp(old_scan_name, vois_list')) == 0
    match = find(strcmp(old_scan_name,handles.VOIs), 1);
    handles.VOIs(match) = [];
end

% check if the old parameter exist somewhere else, if not remove it
parameters_list = [];
for i=1:numel(handles.database)
    for j = 1:numel(handles.database(i).day)
        parameters_list = [parameters_list handles.database(i).day(j).parameters]; %#ok<AGROW>
        
    end
end
parameters_list= unique(parameters_list);

if sum(strcmp(old_scan_name, parameters_list')) == 0
    match = find(strcmp(old_scan_name,handles.clips(:,1)), 1);
    handles.clips(match,:) = [];
end

guidata(hObject, handles)

%%% update graph and display
MIA_update_database_display(hObject, eventdata, handles)
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
    end;
catch
    
end;



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function MIA_GUI_WindowButtonUpFcn(hObject, eventdata, handles)

% hObject    handle to MIA_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if ~strcmp(get(handles.MIA_GUI,'SelectionType'),'normal')
    set(handles.MIA_GUI,'WindowButtonMotionFcn',   @(hObject,eventdata)MIA2('MIA_GUI_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject)));
end
% save contrast
handles.display_option.manual_contrast = 1;
guidata(hObject, handles)


% --------------------------------------------------------------------
function MIA_copy_ScanVoi_to_other_tp_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_copy_ScanVoi_to_other_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
patient_selected = get(handles.MIA_name_list, 'Value');
time_point_selected = get(handles.MIA_time_points_list, 'Value');
scan_voi_selected = get(handles.MIA_scans_list, 'Value');
[time_point,ok] = listdlg('PromptString', 'Select 1 or several time point',...
    'Name', 'Question?',...
    'ListSize', [300 500],...
    'ListString', {handles.database(patient_selected).day.date}');
if ok == 0
    return
end
logbook = {};
for i=1:numel(time_point)
    for j = 1:numel(scan_voi_selected)
        % copie scan
        if get(handles.MIA_scan_VOIs_button, 'Value') == 0
            scan_name = handles.database(patient_selected).day(time_point_selected).parameters(scan_voi_selected(j));
            file_name_to_copy = handles.database(patient_selected).day(time_point_selected).scans_file(scan_voi_selected(j));
            listing = handles.database(patient_selected).day(time_point(i)).parameters';
            % copie VOI
        else
            scan_name = handles.database(patient_selected).day(time_point_selected).VOIs(scan_voi_selected(j));
            file_name_to_copy = handles.database(patient_selected).day(time_point_selected).VOIs_file(scan_voi_selected(j));
            listing = handles.database(patient_selected).day(time_point(i)).VOIs';
        end
        PATHSTR = handles.database(patient_selected).path;
        if sum(strcmp(scan_name, listing)) == 0
            filename = [handles.database(patient_selected).name, '-', ...
                handles.database(patient_selected).day(time_point(i)).date, '-',...
                scan_name{:}, '-ROI.mat'];
            filename = regexprep(filename,'[: ]','-');
            filename = regexprep(filename,'*','star');
            filename = regexprep(filename,'(','');
            filename = regexprep(filename,')','');
            filename = regexprep(filename,',','');
            filename = regexprep(filename,'<','');
            filename = regexprep(filename,'>','');
            filename = regexprep(filename,'/','-');
            if  exist(fullfile(PATHSTR,file_name_to_copy{:}), 'file') == 0
                warning_text = sprintf('Can not %s because it do not exist\n',...
                    scan_name{:});
                logbook{numel(logbook)+1} =warning_text;
                continue
            elseif  exist(fullfile(PATHSTR,filename), 'file') == 2
                
                user_response = questdlg(['The file: ' filename ' already exist. Do you want to overwrite this file or no?'],'Warning', 'Yes', 'No', 'No');
                if strcmp(user_response, 'No')
                    return
                end
                user_response = questdlg('Do you REALLY want to overwrite this file or no?','Warning', 'No', 'Yes', 'No');
                % overwrite the file if requested
                if strcmp(user_response, 'Yes')
                    copyfile(fullfile(PATHSTR,file_name_to_copy{:}), fullfile(PATHSTR,filename), 'f');
                    if get(handles.MIA_scan_VOIs_button, 'Value') == 0
                        if isempty(handles.database(patient_selected).day(time_point(i)).parameters)
                            handles.database(patient_selected).day(time_point(i)).parameters = scan_name;
                            handles.database(patient_selected).day(time_point(i)).scans_file =  {filename};
                        else
                            handles.database(patient_selected).day(time_point(i)).parameters(end+1) = scan_name;
                            handles.database(patient_selected).day(time_point(i)).scans_file(end+1) =  {filename};
                        end
                    else
                        if isempty(handles.database(patient_selected).day(time_point(i)).VOIs_file)
                            handles.database(patient_selected).day(time_point(i)).VOIs = scan_name;
                            handles.database(patient_selected).day(time_point(i)).VOIs_file =  {filename};
                        else
                            handles.database(patient_selected).day(time_point(i)).VOIs(end+1) = scan_name;
                            handles.database(patient_selected).day(time_point(i)).VOIs_file(end+1) =  {filename};
                        end
                    end
                end
            else
                copyfile(fullfile(PATHSTR,file_name_to_copy{:}), fullfile(PATHSTR,filename), 'f');
                if get(handles.MIA_scan_VOIs_button, 'Value') == 0
                    if isempty(handles.database(patient_selected).day(time_point(i)).parameters)
                        handles.database(patient_selected).day(time_point(i)).parameters = scan_name;
                        handles.database(patient_selected).day(time_point(i)).scans_file =  {filename};
                    else
                        handles.database(patient_selected).day(time_point(i)).parameters(end+1) = scan_name;
                        handles.database(patient_selected).day(time_point(i)).scans_file(end+1) =  {filename};
                    end
                else
                    if isempty(handles.database(patient_selected).day(time_point(i)).VOIs_file)
                        handles.database(patient_selected).day(time_point(i)).VOIs = scan_name;
                        handles.database(patient_selected).day(time_point(i)).VOIs_file =  {filename};
                    else
                        handles.database(patient_selected).day(time_point(i)).VOIs(end+1) = scan_name;
                        handles.database(patient_selected).day(time_point(i)).VOIs_file(end+1) =  {filename};
                    end
                end
            end
            warning_text = sprintf('%s was copied to %s\n',...
                scan_name{:},...
                handles.database(patient_selected).day(time_point(i)).date);
            logbook{numel(logbook)+1} =warning_text;
        else
            warning_text = sprintf('%s already exist for %s\n',...
                scan_name{:},...
                handles.database(patient_selected).day(time_point(i)).date);
            logbook{numel(logbook)+1} =warning_text;
        end
    end
end

%update handes
guidata(hObject, handles)

if ~isempty(logbook)
    listdlg('ListString', logbook', 'ListSize',[250 350], 'Name', 'logbook');
else
    msgbox('Done', 'logbook') ;
end

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
%
%
% MRIManager.MRIFileManagerFM
if isfield(handles, 'database')
    MIA_root_path = handles.database.Properties.UserData.MIA_root_path;
    %     MIA_data_path = handles.database.Properties.UserData.MIA_data_path ;
else
    MIA_root_path = uigetdir(pwd, 'Select the directory to save the your new projet');
    if sum(MIA_root_path) == 0
        return
    end
    
    MIA_root_path = [MIA_root_path filesep];
    %% create the output folder ('MIA_data')
    MIA_data_path = [MIA_root_path, 'MIA_data', filesep];
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
    % create the database structure
    handles.database =  table;%
    handles.database = cell2table(cell(0,8));
    handles.database.Properties.VariableNames = {'Group','Patient', 'Tp', 'Path', 'Filename', 'Type', 'IsRaw', 'SequenceName'};
    handles.database.Properties.UserData.MIA_root_path = MIA_root_path;
    handles.database.Properties.UserData.MIA_data_path = MIA_data_path;
    
    
end
%% create the tmp folder
MIA_tmp_folder = [MIA_root_path, 'tmp'];
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
MRIFileManager.FileManagerFrame.main({MIA_tmp_folder_for_java, namExport})
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

handles = guidata(handles.MIA_GUI);
MIA_tmp_folder = [handles.database.Properties.UserData.MIA_root_path, 'tmp'];
log_file =struct2table(jsondecode(char(data_loaded)));

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
                    if ~isempty(handles.database)
                        
                        %% check if a scan with the same SequenceName exist for this patient at this time point. If so, add suffix to the SequenceName (ie. SequenceName(X)
                        if sum(handles.database.Patient ==  char(name_selected) & handles.database.Tp ==  char(tp_selected) &  handles.database.SequenceName == char(clean_variable_name(char(json_data.ProtocolName), ''))) == 1
                            %            idx = find(ismember(handles.database.SequenceName, char(clean_variable_name(char(json_data.ProtocolName), ''))));
                            idx =  strfind(cellstr(handles.database.SequenceName), char(clean_variable_name(char(json_data.ProtocolName), '')));
                            nbr_of_seq = sum([idx{:}]);
                            %             seq_name = handles.database.SequenceName(handles.database.Patient ==  name_selected & handles.database.Tp ==  answer(2) &  handles.database.SequenceName == char(clean_variable_name(char(json_data.ProtocolName), '')));
                            %             seq_name = handles.database.SequenceName((handles.database.Patient) ==  name_selected & (handles.database.Tp ==  answer(2)) &  startsWith(char(handles.database.SequenceName),char(json_data.ProtocolName)));
                            seq_name = [char(json_data.ProtocolName) '(' num2str(nbr_of_seq+1) ')'];
                            file_name = strcat(name_selected, '-', tp_selected,'-',seq_name,'_',datestr(now,'yyyymmdd-HHMMSSFFF'));
                            %             new_data = {'Patient','Tp','nii','json','type', 'SequenceName', 'Group';
                            %                 name_selected,answer(2),fullfile([char(file_name), '.nii']),fullfile([char(file_name), '.json']), 'Scan', seq_name, 'Undefined'};
                        else
                            seq_name = clean_variable_name(char(json_data.ProtocolName), '');
                            file_name = strcat(name_selected, '-', tp_selected,'-',seq_name,'_',datestr(now,'yyyymmdd-HHMMSSFFF'));
                            %             new_data = {'Patient','Tp','nii','json','type', 'SequenceName', 'Group';
                            %                 name_selected,answer(2),fullfile([char(file_name), '.nii']),fullfile([char(file_name), '.json']), 'Scan', seq_name, 'Undefined'};
                        end
                    else
                        seq_name = clean_variable_name(char(json_data.ProtocolName), '');
                        file_name = strcat(name_selected , '-', tp_selected,'-',seq_name,'_',datestr(now,'yyyymmdd-HHMMSSFFF'));
                        
                    end
                    new_data = table(categorical(cellstr('Undefined')), categorical(cellstr(name_selected)), categorical(cellstr(tp_selected)), categorical(cellstr(handles.database.Properties.UserData.MIA_data_path)), categorical(cellstr(file_name)),...
                        categorical(cellstr('Scan')), categorical(1), categorical(cellstr(seq_name)),...
                        'VariableNames', {'Group','Patient', 'Tp', 'Path', 'Filename', 'Type', 'IsRaw', 'SequenceName'});
                    
                    %% add data to the database
                    handles.database = [handles.database; new_data];
                    
                    %%save the files (nii + json)
                    movefile(fullfile(MIA_tmp_folder, [NAME, '.nii']), [char(handles.database.Path(end)), char(handles.database.Filename(end)), '.nii']);
                    movefile(fullfile(MIA_tmp_folder, [NAME, '.json']), [char(handles.database.Path(end)), char(handles.database.Filename(end)), '.json']);
                    
                    %             % update the json file
                    %             log_file_to_update(index_data_to_import(m),:) = [];
                    %             spm_jsonwrite(fullfile(MIA_tmp_folder, log_files), table2struct(log_file_to_update));
                    guidata(hObject, handles);
                end
            end
        end
    end
end

% delete temps files and folder
rmdir(MIA_tmp_folder, 's')

%   MIA_add_name_Callback(hObject, eventdata, handles)

guidata(hObject, handles);

MIA_update_database_display(hObject, eventdata, handles);
msgbox('Import Done!', 'Message') ;


% --------------------------------------------------------------------
function MIA_menu_Clean_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Clean_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%% create a listing of all file present in the database

warning_result = questdlg('This function will move every file present in your hard drive and not present in your database, do you want to continue?',...
    'Warning', 'Yes', 'No', 'No');
if isempty(warning_result) || strcmp(warning_result, 'No')
    return
end

scan_list = [];
VOIs_list = [];
for i=1:size(handles.database,2)
    for j = 1:numel(handles.database(i).day)
        scan_list = [scan_list handles.database(i).day(j).scans_file]; %#ok<AGROW>
        VOIs_list = [VOIs_list handles.database(i).day(j).VOIs_file]; %#ok<AGROW>
    end
end
database_file = cell([scan_list VOIs_list])';

%%% create a listing of all matlab file in the Image_analyses_data folder
cellf = @(fun, arr) cellfun(fun, num2cell(arr), 'uniformoutput',0);
hard_drive_files=cellf(@(f) char(f.tochar()), java.io.File(handles.database(1).path).list());
hard_drive_files = cell2struct(hard_drive_files, 'name' , size(hard_drive_files,1));

if exist(fullfile(handles.database(1).path, 'to_trash'), 'dir') ~= 7
    status = mkdir(fullfile(handles.database(1).path, 'to_trash'));
    if status == 0
        warndlg('You do not the right to write in the folder!', 'Warning');
        return
    end
end

%% create a listing of all file how need to be deleted
file_moved = 0;
for i=1:size(hard_drive_files,1)
    if strcmp(hard_drive_files(i).name, 'to_trash') || strcmp(hard_drive_files(i).name, 'desktop.ini')
        continue
    end
    if sum(strcmp(database_file, hard_drive_files(i).name)) == 0
        file_moved = file_moved+1;
        FileRename([handles.database(1).path, hard_drive_files(i).name],...
            [handles.database(1).path, 'to_trash', filesep,hard_drive_files(i).name])
    end
end

msgbox([num2str(file_moved), ' Files have been moved to ', handles.database(1).path, 'to_trash'], 'Message') ;




% --------------------------------------------------------------------
function MIA_menu_export_PRM_to_tiff_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_export_PRM_to_tiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% need to be in PRM mode
if handles.mode == 1
    warndlg('This function works only in PRM mode!', 'Warning');
    return
end
if ~isfield(handles, 'data_displayed')
    warndlg('No data loaded!', 'Warning');
    return
end
if ~isfield(handles.data_displayed, 'PRM')
    warndlg('No VOI loaded!', 'Warning');
    return
end

slice_nbr = get(handles.MIA_slider_slice, 'Value');

[image_selected,ok] = listdlg('Name', 'Bip',...
    'ListString', [get(handles.MIA_data1_title, 'String') get(handles.MIA_data2_title, 'String') 'Other']','ListSize', [200 150], 'PromptString', 'Select the background image');
if ok == 0
    return
end
if image_selected == 3
    %     [image_selected,ok] = listdlg('Name', 'Bip',...
    %     'ListString', handles.database(handles.data_selected_for_PRM.patient_info.name_nbr).day(handles.data_selected_for_PRM.patient_info.timepoint_nbr).pae        ,'ListSize', [200 150], 'PromptString', 'Select the background image');
    
end
%%% display image
figure; h = imagesc(handles.data_displayed.image(:,:,slice_nbr,image_selected));
colormap_selected = handles.colormap(get(handles.MIA_colormap_popupmenu,'Value'));
eval(['colormap(handles.MIA_data' stri ', ''' colormap_selected{:} ''');']);
colormap( 'gray');
hold on

%%% display the PRM
test = squeeze(handles.data_displayed.PRM.data(:,:,slice_nbr,:));
h = imshow(test);
hold off
%%% create a transparent mask for the PRM
index = findn(sum(test,3)==0);
transpartant = ones([size(handles.data_displayed.PRM.data, 1) size(handles.data_displayed.PRM.data, 2)]);
for ii = 1: size(index,1)
    transpartant(index(ii,1),index(ii,2)) = 0;
end
%%% apply the transpartent mask
set(h, 'AlphaData',transpartant);




% --------------------------------------------------------------------
function MIA_menu_Importing_form_another_database_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_menu_Importing_form_another_database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname]=uigetfile('*.mat','Load database','MultiSelect','off');

if pathname == 0
    return
else
    if ~strcmp(class(filename),'double') %#ok<STISA>
        
        tmp_database = load(fullfile(pathname,filename));
        data_listing = {};
        patient_id_timepoint = [];
        for i = 1:numel(tmp_database.database)
            for j=1:numel(tmp_database.database(i).day)
                data_listing(size(data_listing, 1)+1,:) = {[tmp_database.database(i).name, char(tmp_database.database(i).day(j).date)]};
                patient_id_timepoint(size(patient_id_timepoint, 1)+1,:) = [i, j];
            end
        end
        
        [data_to_import, ok] = listdlg('PromptString','Select a the data you want to import',...
            'Name', 'Select a data',...
            'ListSize', [400 300],...
            'ListString',data_listing);
        
        if ok == 0
            return
        end
        
        % for old database
        for i = 1:numel(data_to_import)
            
            % upadate the database
            if sum(strcmp({handles.database.name}', tmp_database.database(patient_id_timepoint(data_to_import(i),1)).name)) == 0
                
                handles.database(numel(handles.database)+1).name =  tmp_database.database(patient_id_timepoint(data_to_import(i),1)).name;
                handles.database(numel(handles.database)).path =  handles.database(1).path;
                handles.database(numel(handles.database)).group =  tmp_database.database(patient_id_timepoint(data_to_import(i),1)).group;
                handles.database(numel(handles.database)).omit =  tmp_database.database(patient_id_timepoint(data_to_import(i),1)).omit;
                handles.database(numel(handles.database)).day = ...
                    tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2));
                
            else
                potistion_in_database = strcmp({handles.database.name}', tmp_database.database(patient_id_timepoint(data_to_import(i),1)).name);
                handles.database(potistion_in_database).day(numel(handles.database(potistion_in_database).day)+1) = ...
                    tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2));
            end
            % copy the ROI files
            for j=1:numel(tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).VOIs_file)
                if exist(char(fullfile(tmp_database.database(patient_id_timepoint(data_to_import(i),1)).path, tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).VOIs_file(j))) , 'file') == 2
                    copyfile(char(fullfile(tmp_database.database(patient_id_timepoint(data_to_import(i),1)).path, tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).VOIs_file(j))),...
                        char(fullfile(handles.database(1).path, tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).VOIs_file(j))), 'f');
                end
            end
            % copy the parameter files
            for j=1:numel(tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).scans_file)
                if exist(char(fullfile(tmp_database.database(patient_id_timepoint(data_to_import(i),1)).path, tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).scans_file(j))) , 'file') == 2
                    copyfile(char(fullfile(tmp_database.database(patient_id_timepoint(data_to_import(i),1)).path, tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).scans_file(j))),...
                        char(fullfile(handles.database(1).path, tmp_database.database(patient_id_timepoint(data_to_import(i),1)).day(patient_id_timepoint(data_to_import(i),2)).scans_file(j))), 'f');
                end
            end
            guidata(hObject, handles);
            MIA_update_database_display(hObject, eventdata, handles);
        end
    end
end
num_timepoint = 0;
for i=1:numel(handles.database)
    num_timepoint = num_timepoint+ numel({handles.database(i).day.date});
end
title = ['database name: ',filename,'; ', num2str(numel(handles.database)), ' patients and ',  num2str(num_timepoint), ' ','time points'];
set(handles.MIA_GUI, 'Name', title);
set(handles.MIA_name_list, 'Value', 1);


guidata(hObject, handles);
MIA_update_database_display(hObject, eventdata, handles);
msgbox('Done', 'Message') ;



% --- Executes on button press in MIA_Brain_Extraction.
function MIA_Brain_Extraction_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_Brain_Extraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'data_loaded')
    warndlg('Please open 1 scan (in Single mode)!', 'Warning');
    return
end
if handles.data_loaded(1).number_of_scan > 1
    warndlg('Please open 1 scan only!', 'Warning');
    return
end
if sum(strcmp([handles.data_selected.list_scan{:} '-masked'], handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).parameters)) == 1
    warndlg('Brain extraction already performed on this scan!', 'Warning');
    return
end

I = squeeze(handles.data_selected_resized.image.reco.data(:,:,1,:,1));
%%radius of the structural element (in pixel)
radius = round(size(handles.data_selected_resized.image.reco.data,1)*size(handles.data_selected_resized.image.reco.data,2) / 4.3691e+03); %16/(256*256) =2.4414e-04      16 = un nombre optimal pour les images en 256x256
voxeldim = [handles.data_selected_resized.image.reco.fov(1)/size(handles.data_selected_resized.image.reco.data,1) handles.data_selected_resized.image.reco.fov(2)/size(handles.data_selected_resized.image.reco.data,2) handles.data_selected_resized.image.reco.no_slices/size(handles.data_selected_resized.image.reco.data,4)];

selection = questdlg('Mice or Rat?',...
    'Warning',...
    'Mice','Rat', 'Human', 'Rat');
if isempty(selection)
    return
end
switch selection,
    case 'Mice'
        BrSize = [100 550];
    case 'Rat'
        BrSize = [1200 4400];
    case 'Human'
        BrSize = [120000 440000];
end


maxiter = 200;
OptMask = 0;


% this code come from
% https://sites.google.com/site/chuanglab/software/3d-pcnn
% email : k.chuang@uq.edu.au
% Reference & Citation
% N Chou, J Wu, J Bai, A Qiu, KH Chuang. ???Robust automatic rodent brain extraction using 3D Pulse-coupled NeuralNetworks (PCNN).??? IEEE Trans Imag Proc 20(9):2554-64, 2011.
% [I_border, G_I, optG] = PCNN3D(I, p, voxeldim, BrSize, maxiter, OptMask)
%
% I: 3D matrix of the brain image
% p: radius of the structural element (in pixel)
% voxeldim: [x, y, z] voxel dimensions of the data set (in mm)
% BrSize: [min max] brain volume sizes (in mm3) for estimating optimal iteration; default values [100 550] is for mouse brain. For adult rat brain, one may try [1200 4400]
% maxiter: maximum iteration. Default is 200
% OptMask: output only the mask based on the optimal guess. Note: the guess may not always give the best result.
%
% The outputs:
%    I_border: A cell array that stores the binary mask at each iteration.
%    G_I: the plot of volume of mask (in mm3) vs. iterations, like the following.
%    optG: the suggested optimal iteration.

[handles.brain_extraction_ROI, handles.brain_extraction_G_I, handles.brain_extraction_Best_It] = PCNN3D(I, radius, voxeldim, BrSize, maxiter, OptMask);
close(gcf)

handles.brain_extraction_Original_data = handles.data_selected_resized.image.reco.data;
set(handles.MIA_PRM_slider_trans,'Value',handles.brain_extraction_Best_It);
set(handles.MIA_PRM_slider_trans,'max',find(handles.brain_extraction_G_I ~= 0, 1, 'last' )-1);
set(handles.MIA_PRM_slider_trans,'min',find(handles.brain_extraction_G_I ~= 0, 1 ));

handles.StopButton = uicontrol('Style','Pushbutton','Position',[250 70 80 40],'Tag', 'Brain_Extraction_Stop', 'String','Done','Callback',@MIA_Brain_Extraction_StopCallback);
handles.CancelButton = uicontrol('Style','Pushbutton','Position',[450 70 80 40],'Tag', 'Brain_Extraction_Cancel', 'String','Cancel','Callback',@MIA_Brain_Extraction_CancelCallback);

guidata(hObject, handles)
% update current display
MIA_PRM_slider_trans_Callback(hObject, eventdata, handles)




function MIA_Brain_Extraction(hObject, eventdata, handles)

Iteration_number = round(get(handles.MIA_PRM_slider_trans,'Value'));
Mask = reshape(full(cell2mat(handles.brain_extraction_ROI{Iteration_number})), [size(handles.data_selected_resized.image.reco.data,1) size(handles.data_selected_resized.image.reco.data,2) handles.data_selected_resized.image.reco.no_slices]);
handles.data_selected_resized.image.reco.data=handles.brain_extraction_Original_data  .* repmat(permute(Mask, [1 2 4 3 5]), [1 1 size(handles.brain_extraction_Original_data, 3) 1 size(handles.brain_extraction_Original_data, 5)]);
% test replace all values = 0 to NaN. Maybe better for co-registration
handles.data_selected_resized.image.reco.data(handles.data_selected_resized.image.reco.data == 0) = NaN;

guidata(hObject, handles)

MIA_update_axes(hObject, eventdata, handles)

function MIA_Brain_Extraction_StopCallback(hObject, eventdata, ~)
%// Retrieve elements from handles structure.

handles = guidata(hObject);
logbook = {};
Iteration_number = round(get(handles.MIA_PRM_slider_trans,'Value'));
Mask = reshape(full(cell2mat(handles.brain_extraction_ROI{Iteration_number})), [size(handles.data_selected_resized.image.reco.data,1) size(handles.data_selected_resized.image.reco.data,2) handles.data_selected_resized.image.reco.no_slices]);

% save file

load(fullfile(handles.database(handles.data_selected.patient_info.name_nbr).path,...
    handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{strcmp(handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).parameters, handles.data_selected.list_scan{1})}));
[x, y, echo, ~, repet]=size(uvascim.image.reco.data); %#ok<NODEF>

mask3d_resized=  imresize(Mask, [x, y],'bilinear');
mask5d = repmat(mask3d_resized,[1 1 1 echo repet]);
mask5d = permute(mask5d, [1 2 4 3 5]);

uvascim.image.reco.data=uvascim.image.reco.data.* mask5d;

parametre_nbr = strcmp(handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).parameters, handles.data_selected.list_scan{1});
[path, name, ext] = fileparts(fullfile(handles.database(handles.data_selected.patient_info.name_nbr).path,...
    handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).scans_file{parametre_nbr}));
scan_name = [name '-masked' ext];
save(fullfile(path, scan_name),'uvascim');

% add to the database and update handles.clips
handles =  MIA_add_new_scan_to_database(handles, uvascim.image, handles.data_selected.patient_info.name_nbr, handles.data_selected.patient_info.timepoint_nbr,...
    [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).parameters{parametre_nbr} '-masked'], [name '-masked' ext]);
guidata(hObject, handles)
MIA_update_database_display(hObject, eventdata, handles);
% handles = guidata(hObject);
clear uvascim;

% save the ROI
new_ROI_name = 'Mask';
for i=1:get(handles.MIA_slider_slice, 'Max')
    uvascroi(i).value = squeeze(mask3d_resized(:,:,i));
    uvascroi(i).area_in_pixel=sum(sum(uvascroi(i).value));
    uvascroi(i).position=[0,0; 0, handles.resolution_selected; handles.resolution_selected, handles.resolution_selected; handles.resolution_selected, 0];
    uvascroi(i).name={new_ROI_name};
    uvascroi(i).colorindex=1;
    %%% case roi drown on histo data
    uvascroi(i).displayedecho=handles.data_selected.image(1).reco.displayedecho;
    uvascroi(i).displayedslice=handles.data_selected.image(1).reco.displayedslice;
    uvascroi(i).displayedexpt=handles.data_selected.image(1).reco.displayedexpt;
    uvascroi(i).fov_orientation=handles.data_selected.image(1).reco.fov_orientation(1:9,1,i,1);
    uvascroi(i).no_samples=handles.data_selected.image(1).reco.no_samples;
    if isfield(handles.data_selected.image(1).reco, 'no_views')
        uvascroi(i).no_views=handles.data_selected.image(1).reco.no_views;
    end
    uvascroi(i).fov=handles.data_selected.image(1).acq.fov;
    if isfield(handles.data_selected.image(1).reco, 'no_views')
        uvascroi(i).area_in_mmxmm=(uvascroi(i).fov(1)/uvascroi(i).no_samples)*(uvascroi(i).fov(2)/uvascroi(i).no_views)*uvascroi(i).area_in_pixel;
    end
    if isfield(handles.data_selected.image(1), 'reco_number')
        uvascroi(i).reco_number=handles.data_selected.image(1).reco_number;
    end
    if isfield(handles.data_selected.image(1), 'scan_number')
        uvascroi(i).scan_number=handles.data_selected.image(1).scan_number;
    end
    uvascroi(i).filename=handles.data_selected.image(1).filename;
    uvascroi(i).filenames=sprintf('%s\n',handles.data_selected.image(1).filename);
    uvascroi(i).fov_offsets=handles.data_selected.image(1).reco.fov_offsets(1:3,1,i,1);
    
    %     end
    uvascroi(i).thickness = handles.data_selected.image(1).reco.thickness;
    uvascroi(i).date=date;
    uvascroi(i).affiche=0;
end
if ~isempty(uvascroi)
    pathname = [handles.database(handles.data_selected.patient_info.name_nbr).path handles.database(handles.data_selected.patient_info.name_nbr).name handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).date '-VOI-' new_ROI_name '.mat'];
    save(pathname,'uvascroi');
    handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs = [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs {new_ROI_name}];
    voi_file = [handles.database(handles.data_selected.patient_info.name_nbr).name handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).date '-VOI-' new_ROI_name '.mat'];
    handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs_file = [handles.database(handles.data_selected.patient_info.name_nbr).day(handles.data_selected.patient_info.timepoint_nbr).VOIs_file {voi_file}];
end
% add cluster to the database
handles.VOIs = handles.VOIs(1:end-1);
handles.VOIs(end+1) = {new_ROI_name};
handles.VOIs = unique(handles.VOIs);
handles.VOIs(end+1) = {'Other'};

% delete variables used during the brain extraction function
handles = rmfield(handles, 'brain_extraction_Original_data');
handles = rmfield(handles, 'brain_extraction_ROI');
handles = rmfield(handles, 'brain_extraction_G_I');
handles = rmfield(handles, 'brain_extraction_Best_It');

% delete buttons
% delete buttons
if ~isempty(findobj('Tag', 'Brain_Extraction_Stop'))
    delete(findobj('Tag', 'Brain_Extraction_Stop'))
end
if ~isempty(findobj('Tag', 'Brain_Extraction_Cancel'))
    delete(findobj('Tag', 'Brain_Extraction_Cancel'))
end

guidata(handles.MIA_GUI, handles)




function MIA_Brain_Extraction_CancelCallback(hObject, eventdata, ~)

%// Retrieve elements from handles structure.
handles = guidata(hObject);
if isfield(handles, 'brain_extraction_Original_data')
    handles.data_selected_resized.image.reco.data=handles.brain_extraction_Original_data;
    handles = rmfield(handles, 'brain_extraction_Original_data');
end

if isfield(handles, 'brain_extraction_ROI')
    handles = rmfield(handles, 'brain_extraction_ROI');
end
if isfield(handles, 'brain_extraction_G_I')
    
    handles = rmfield(handles, 'brain_extraction_G_I');
end
if isfield(handles, 'brain_extraction_Best_It')
    handles = rmfield(handles, 'brain_extraction_Best_It');
end
if isfield(handles, 'StopButton')
    handles = rmfield(handles, 'StopButton');
end
if isfield(handles, 'CancelButton')
    handles = rmfield(handles, 'CancelButton');
end

guidata(hObject, handles)

MIA_update_axes(hObject, eventdata, handles)
% delete buttons
if ~isempty(findobj('Tag', 'Brain_Extraction_Stop'))
    delete(findobj('Tag', 'Brain_Extraction_Stop'))
end
if ~isempty(findobj('Tag', 'Brain_Extraction_Cancel'))
    delete(findobj('Tag', 'Brain_Extraction_Cancel'))
end

% --- Executes on button press in MIA_test_button.
function MIA_test_button_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_test_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'database')
    return
end

handles.database
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
% set(handles.MIA_GUI,'WindowButtonUpFcn',   @(hObject,eventdata)MIA2('MIA_slider_release_click',hObject,eventdata,guidata(hObject)));
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
if isfield(handles,'data_selected') == 0
    [filename,pathname] = uigetfile('.mat',...
        'Select a cluster to print its statistics',...
        strcat(pwd,'/ClustersGMM'));
    load(strcat(pathname,filename),'Informations','Statistiques');
else
    load(handles.data_selected.roifile{:},'uvascroi','Informations', 'Statistiques');
end
Signatures(Informations, Statistiques)


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
%      handles.data_loaded.Scan(i).nii = read_volume(handles.data_loaded.Scan(i).V, handles.data_loaded.Scan(scan_of_reference).V);
% end
if isfield(handles.data_loaded, 'ROI')
    for i=1:numel(handles.data_loaded.ROI)
        handles.data_loaded.ROI(i).nii = read_volume(handles.data_loaded.ROI(i).V, handles.data_loaded.Scan(scan_of_reference).V);
    end
end
% % update slider_slice
% set(handles.MIA_slider_slice, 'max', size(handles.data_loaded.Scan(1).nii, 3), 'Value', 1, 'SliderStep',[1/(size(handles.data_loaded.Scan(1).nii, 3) -1) min(5/(size(handles.data_loaded.Scan(1).nii, 3) -1),1)]);
%
% %update handes
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


% --- Executes on button press in MIA_pipeline_creator.
function MIA_pipeline_creator_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_pipeline_creator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MIA_pipeline_creator
if ~isfield(handles, 'database')
    warndlg('No database loaded', 'Warning');
    return
end
MIA_pipeline(hObject, eventdata, handles)

function nii_json_fullfilename = fullfilename(handles, nii_index, ext)

nii_json_fullfilename = [char(handles.database.Path(nii_index)) char(handles.database.Filename(nii_index)) ext];

function MIA_remove_scan(hObject, eventdata, handles, nii_index)
% delete the file (nii/json) from the hard drive
for i=1:numel(nii_index)
    % delete the nii file
    delete(fullfilename(handles, nii_index(i), '.nii'))
    % delete the json file
    delete(fullfilename(handles, nii_index(i), '.json'))
end


% remove the scan from the database
handles.database(nii_index,:) = [];
switch get(hObject, 'Tag')
    case 'MIA_remove_scan'
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
