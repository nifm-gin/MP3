function varargout = speedy_module(varargin)
% SPEEDY_MODULE MATLAB code for speedy_module.fig
%      SPEEDY_MODULE, by itself, creates a new SPEEDY_MODULE or raises the existing
%      singleton*.
%
%      H = SPEEDY_MODULE returns the handle to a new SPEEDY_MODULE or the handle to
%      the existing singleton*.
%
%      SPEEDY_MODULE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPEEDY_MODULE.M with the given input arguments.
%
%      SPEEDY_MODULE('Property','Value',...) creates a new SPEEDY_MODULE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speedy_module_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speedy_module_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speedy_module

% Last Modified by GUIDE v2.5 23-Feb-2017 16:12:40


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @speedy_module_OpeningFcn, ...
    'gui_OutputFcn',  @speedy_module_OutputFcn, ...
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


% --- Executes just before speedy_module is made visible.
function speedy_module_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to speedy_module (see VARARGIN)

% Choose default command line output for speedy_module
handles.output = hObject;

% Save the Image_Analysis_module handles in the handles of the speedy module
handles.MIA_data = varargin{3};

handles.mean_by_mode = 'ROI';
handles.color_RGB = ...
    [ 1 0 0; ...                  Red
    0 1 0; ...                  Lime
    0 0 1; ...                  Blue
    1 1 0;  ...                 Yellow
    0 1 1; ...                  Aqua
    1 0 1; ...                  Fuchsia
    0 0.5 0.5; ...              Teal
    1 0.6471 0; ...             Orange
    0.5 0 0.5; ...              Purple
    0 0.5 0; ...                Green
    0.5 0 0; ...                Maroon
    0 0 0.5; ...                Navy
    0.5 0.5 0; ...              Olive
    0.7529 0.7529 0.7529; ...   Silver
    0.5 0.5 0.5; ...            Gray
    0 0 0; ...                 Black
    1 1 1];                     %White

handles.color_shortName = {'k','r','g','b','y','c','m'};
set(findobj('Tag', 'speedy_histogram_based_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_id'), 'Value', 0)

handles.speedy_graph1_tooltip = uicontrol(handles.speedy_GUI, 'Style', 'text',...
    'Visible', 'off', 'Units', 'characters');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes speedy_module wait for user response (see UIRESUME)
% uiwait(handles.speedy_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = speedy_module_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in speedy_patient_name_listbox.
function speedy_patient_name_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_patient_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_patient_name_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_patient_name_listbox

IA_update_database_display(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function speedy_patient_name_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_patient_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_time_points_listbox.
function speedy_time_points_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_time_points_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_time_points_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_time_points_listbox

IA_update_database_display(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function speedy_time_points_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_time_points_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_ROIs_listbox.
function speedy_ROIs_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_ROIs_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_ROIs_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_ROIs_listbox

IA_update_database_display(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function speedy_ROIs_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_ROIs_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_clusters_listbox.
function speedy_clusters_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_clusters_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_clusters_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_clusters_listbox


% --- Executes during object creation, after setting all properties.
function speedy_clusters_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_clusters_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_parameters_listbox.
function speedy_parameters_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_parameters_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_parameters_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_parameters_listbox

IA_update_database_display(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.template_dataset.
function speedy_parameters_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_parameters_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_group_name_listbox.
function speedy_group_name_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_group_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_group_name_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_group_name_listbox

% get group info

IA_update_database_display(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function speedy_group_name_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_group_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in speedy_getData.
function speedy_getData_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_getData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global
global voxels_database

resolution = get(handles.speedy_resolution_popupmenu, 'Value');
if resolution ~= 1
    switch resolution
        case 2
            factor = [64 64];
            resolution_selected = 64;
        case 3
            factor = [128 128];
            resolution_selected = 128;
        case 4
            factor = [256 256];
            resolution_selected = 256;
        case 5
            factor = [512 512];
            resolution_selected = 512;
    end
end
logbook = {};

list_sorted =  sort(handles.MIA_data.clips(:,1));
[parameters_selected,ok] = listdlg('Name', 'Which parameter?', 'ListString',  list_sorted,'ListSize', [300 300],...
    'PromptString', 'Please select the parameter you want to analyse?');
if ok == 0
    return
end


[VOI_selected,ok] = listdlg('Name', 'Which VOI?', 'ListString',  handles.MIA_data.VOIs(1:end-2)','ListSize', [300 300],...
    'PromptString', 'Please select the VOI you want to analyse?');
if ok == 0
    return
end

[~,NAME,EXT] = fileparts([handles.MIA_data.database(1).databaseinfo.pathname handles.MIA_data.database(1).databaseinfo.filename]);

n=0;
for patient = 1:numel(handles.MIA_data.database)
    for timepoint = 1:numel(handles.MIA_data.database(patient).day)
        for parameter = 1:numel(parameters_selected)
            for roi = 1:numel(handles.MIA_data.database(patient).day(timepoint).VOIs)
                n=n+1;
            end
        end
    end
end
voxels_database = cell(n,2);
n=0;
for patient = 1:numel(handles.MIA_data.database)
%for patient = 22:37
    for timepoint = 1:numel(handles.MIA_data.database(patient).day)
        for parameter = 1:numel(parameters_selected)
            set(handles.speedy_info_text, 'String', ['Current Prossess : ', handles.MIA_data.database(patient).name, '-',...
                handles.MIA_data.database(patient).day(timepoint).date, '-', list_sorted{parameters_selected(parameter)}]);
            drawnow;
            
            match_parameter = strcmp(handles.MIA_data.database(patient).day(timepoint).parameters', list_sorted(parameters_selected(parameter),1));
            if sum(match_parameter) ~= 1
                continue
            end
            fid=fopen(fullfile(handles.MIA_data.database(patient).path, handles.MIA_data.database(patient).day(timepoint).scans_file{match_parameter}),'r');
            if fid>0
                fclose(fid);
                data_loaded = load(fullfile(handles.MIA_data.database(patient).path, handles.MIA_data.database(patient).day(timepoint).scans_file{match_parameter}));
                if ~or(isfield(data_loaded, 'uvascim'), isfield(data_loaded, 'so2struct'))
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient).day(timepoint).scans_file{match_parameter} '-is not a correct file'];
                    clear data_loaded
                    continue
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient).day(timepoint).scans_file{match_parameter} '-is not a correct file'];
                clear data_loaded
                continue
            end
            
            if resolution > 1 &&...
                    (size(data_loaded.uvascim.image.reco.data, 1) ~= factor(1) ||...
                    size(data_loaded.uvascim.image.reco.data, 2) ~= factor(2))
                tmp = zeros([factor size(data_loaded.uvascim.image.reco.data, 3) size(data_loaded.uvascim.image.reco.data, 4) size(data_loaded.uvascim.image.reco.data, 5)]);
                for echo = 1:size(data_loaded.uvascim.image.reco.data, 3)
                    for expt = 1:size(data_loaded.uvascim.image.reco.data, 5)
                        if numel(size(data_loaded.uvascim.image.reco.data(:,:,echo,:,expt))) == 3
                            tmp(:,:,echo,:,expt) = imresize3d(squeeze(data_loaded.uvascim.image.reco.data(:,:,echo,:,expt)), [],...
                                [factor size(data_loaded.uvascim.image.reco.data, 4)],'cubic');
                        else
                            tmp(:,:,echo,:,expt) = imresize(squeeze(data_loaded.uvascim.image.reco.data(:,:,echo,:,expt)),...
                                factor,'cubic');
                        end
                    end
                end
                data_loaded.uvascim.image.reco.data = tmp;
                clear tmp
                % resize data
            end
            parameter_size = size(data_loaded.uvascim.image.reco.data);
            %             for roi = 1:numel(handles.MIA_data.database(patient).day(timepoint).VOIs)
            for roi = 1:numel(VOI_selected)
                % load ROI
                %                 fid=fopen(fullfile(handles.MIA_data.database(patient).path, handles.MIA_data.database(patient).day(timepoint).VOIs_file{roi}),'r');
                if sum(strcmp(handles.MIA_data.database(patient).day(timepoint).VOIs, handles.MIA_data.VOIs(VOI_selected(roi))))  == 1
                    fid=fopen(fullfile(handles.MIA_data.database(patient).path,...
                        handles.MIA_data.database(patient).day(timepoint).VOIs_file{ strcmp(handles.MIA_data.database(patient).day(timepoint).VOIs, handles.MIA_data.VOIs(VOI_selected(roi)))}),'r');
                    
                    
                    
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient).name, '-', handles.MIA_data.database(patient).day(timepoint).date, '-No ', handles.MIA_data.VOIs{VOI_selected(roi)}, '-VOI'];
                    clear ROI_loaded
                    continue
                end
                %                 fid=fopen(fullfile(handles.MIA_data.database(patient).path, handles.MIA_data.database(patient).day(timepoint).VOIs_file{roi}),'r');
                
                if fid>0
                    fclose(fid);
                    %                     ROI_loaded = load(fullfile(handles.MIA_data.database(patient).path, handles.MIA_data.database(patient).day(timepoint).VOIs_file{roi}));
                    ROI_loaded = load( fullfile(handles.MIA_data.database(patient).path,...
                        handles.MIA_data.database(patient).day(timepoint).VOIs_file{ strcmp(handles.MIA_data.database(patient).day(timepoint).VOIs, handles.MIA_data.VOIs(VOI_selected(roi)))}));
                    
                    if ~isfield(ROI_loaded, 'uvascroi')
                        logbook{numel(logbook)+1} = [handles.MIA_data.database(patient).name, '-', handles.MIA_data.database(patient).day(timepoint).date, '-', handles.MIA_data.VOIs{VOI_selected(roi)}, '-is not a correct file'];
                        %                         [handles.MIA_data.database(patient).day(timepoint).VOIs_file{roi} '-is not a correct file'];
                        clear ROI_loaded
                        continue
                    end
                else
                    logbook{numel(logbook)+1} = [handles.MIA_data.database(patient).name, '-', handles.MIA_data.database(patient).day(timepoint).date, '-', handles.MIA_data.VOIs{VOI_selected(roi)}, '-is not a correct file'];
                    %                     logbook{numel(logbook)+1} =[handles.MIA_data.database(patient).day(timepoint).VOIs_file{roi} '-is not a correct file'];
                    clear ROI_loaded
                    continue
                end
                
                for roi_nbr = 1:numel(ROI_loaded.uvascroi)
                    %resize ROI in order to match to the parameter size
                    match_slice = abs(squeeze(data_loaded.uvascim.image.reco.fov_offsets(3,1,:,1))-ROI_loaded.uvascroi(roi_nbr).fov_offsets(3)) <  1e-5;
                    if match_slice == 0
                        match_slice = abs(squeeze(data_loaded.uvascim.image.reco.fov_offsets(3,1,:,1))-ROI_loaded.uvascroi(roi_nbr).fov_offsets(3)) <  0.5;
                    end
                    if  sum(match_slice) == 1 %&& data_loaded.uvascim.image.reco.thickness == ROI_loaded.uvascroi(roi_nbr).thickness
                        %                         roi_resize = imresize(ROI_loaded.uvascroi(roi_nbr).value(:,:),[parameter_size(1) parameter_size(2)],'bilinear');
                        roi_resize = logical(imresize(ROI_loaded.uvascroi(roi_nbr).value(:,:),[parameter_size(1) parameter_size(2)],'bilinear'));
                        
                        voxel_position = findn(roi_resize == 1);
                        if isempty(voxel_position)
                            continue
                        end
                        for echo = 1:size(data_loaded.uvascim.image.reco.data, 3)
                            for expt = 1:size(data_loaded.uvascim.image.reco.data, 5)
                                data_2d = squeeze(data_loaded.uvascim.image.reco.data(:,:,echo, match_slice, expt));
                                if ~exist('voxels_database_tmp_cell','var')
                                    % store the cell's data
                                    voxels_database_tmp_cell = [repmat({handles.MIA_data.database(patient).group}, [size(voxel_position,1), 1])...
                                        repmat({handles.MIA_data.database(patient).name}, [size(voxel_position,1), 1])...
                                        repmat({handles.MIA_data.database(patient).day(timepoint).date}, [size(voxel_position,1), 1])...
                                        repmat(handles.MIA_data.database(patient).day(timepoint).parameters(match_parameter), [size(voxel_position,1), 1])...
                                        repmat(handles.MIA_data.VOIs(VOI_selected(roi)), [size(voxel_position,1), 1])];
                                    % store the values data
                                    voxels_database_tmp_double = [voxel_position(:,1) voxel_position(:,2) repmat(echo, [size(voxel_position,1), 1])...
                                        repmat(ROI_loaded.uvascroi(roi_nbr).fov_offsets(3), [size(voxel_position,1), 1])...
                                        repmat(expt, [size(voxel_position,1), 1]) data_2d(voxel_position(:,1)+(voxel_position(:,2)-1)*size(data_2d,1))];
                                else
                                    % store the cell's data
                                    voxels_database_tmp_cell=  cat(1,voxels_database_tmp_cell,[repmat({handles.MIA_data.database(patient).group}, [size(voxel_position,1), 1])...
                                        repmat({handles.MIA_data.database(patient).name}, [size(voxel_position,1), 1])...
                                        repmat({handles.MIA_data.database(patient).day(timepoint).date}, [size(voxel_position,1), 1])...
                                        repmat(handles.MIA_data.database(patient).day(timepoint).parameters(match_parameter), [size(voxel_position,1), 1])...
                                        repmat(handles.MIA_data.VOIs(VOI_selected(roi)), [size(voxel_position,1), 1])]);
                                    % store the values data
                                    voxels_database_tmp_double =  cat(1,voxels_database_tmp_double,[voxel_position(:,1) voxel_position(:,2) repmat(echo, [size(voxel_position,1), 1])...
                                        repmat(ROI_loaded.uvascroi(roi_nbr).fov_offsets(3), [size(voxel_position,1), 1])...
                                        repmat(expt, [size(voxel_position,1), 1]) data_2d(voxel_position(:,1)+(voxel_position(:,2)-1)*size(data_2d,1))]);
                                end
                                clear data_2d
                                
                            end
                        end
                        clear roi_resize voxel_position
                    elseif sum(match_slice) == 1 && data_loaded.uvascim.image.reco.thickness ~= ROI_loaded.uvascroi(roi_nbr).thickness
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient).name '-' handles.MIA_data.database(patient).day(timepoint).date...
                            '- slice thickness missmatch between ' handles.MIA_data.database(patient).day(timepoint).parameters{match_parameter} ' and ' ROI_loaded.uvascroi(roi_nbr).name{:}];
                    end
                end
            end
            if exist('voxels_database_tmp_cell', 'var')
                n=n+1;
                voxels_database_tmp{n,1} = voxels_database_tmp_cell;
                voxels_database_tmp{n,2}=voxels_database_tmp_double;
                clear voxels_database_tmp_cell voxels_database_tmp_double
            end
        end
    end
end
if ~exist('voxels_database_tmp', 'var')
    % display info
    if ~isempty(logbook)
        listdlg('ListString', logbook', 'ListSize',[250 350], 'Name', 'logbook')
    end
    return
end
if isempty(voxels_database_tmp(:,1))
    logbook{numel(logbook)+1} ='No data collected (No parameters/VOIs';
end
if n < size(voxels_database_tmp, 1)
    voxels_database_tmp(n+1:end,:) = [];
end

handles.MIA_data.database(1).databaseinfo.voxels_database_need_to_update = 0;
handles.MIA_data.database(1).databaseinfo.voxels_database_voxels_to_update = [];
guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data);


if size(voxels_database_tmp,1) >= 2
    set(handles.speedy_info_text, 'String', 'Current Prossess : Merging data');
    drawnow
    cell_tmp = cat(1,voxels_database_tmp{:,1});
    double_tmp = cat(1,voxels_database_tmp{:,2});
    clear voxels_database_tmp
    voxels_database =  dataset({cell_tmp(:,1), 'Group'},{cell_tmp(:,2), 'Patient'}, {cell_tmp(:,3), 'Timepoint'},...
        {cell_tmp(:,4), 'Parameter'}, {cell_tmp(:,5), 'VOI'}...
        , {double_tmp(:,1), 'Coord_X'}, {double_tmp(:,2), 'Coord_Y'}, {double_tmp(:,3), 'Echo'}, {double_tmp(:,4), 'Coord_Z'},...
        {double_tmp(:,5), 'Expt'}, {double_tmp(:,6), 'Value'});
    clear cell_tmp double_tmp
    
    handles.MIA_data.database(1).databaseinfo.voxels_database_filename = [NAME, '-voxels_database', EXT];
    
    guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data);
    % update title
    set(handles.speedy_GUI, 'Name',  handles.MIA_data.database(1).databaseinfo.voxels_database_filename);
    
    % % transforme data to nominal data
    voxels_database.Group = nominal(voxels_database.Group);
    voxels_database.Patient = nominal(voxels_database.Patient);
    voxels_database.Timepoint = nominal(voxels_database.Timepoint);
    voxels_database.Parameter = nominal(voxels_database.Parameter);
    voxels_database.VOI = nominal(voxels_database.VOI);
    %% test add resolution to the dataset
    if exist('resolution_selected', 'var')
        voxels_database.Properties.UserData = resolution_selected;
    end
    % save voxels_database
    save(fullfile(handles.MIA_data.database(1).databaseinfo.pathname,  handles.MIA_data.database(1).databaseinfo.voxels_database_filename),'voxels_database', '-v7.3');
end
% end
handles.voxels_database_shortlist = unique(voxels_database(:,1:5));

% update the handles

parameters_list = cellstr(unique(voxels_database.Parameter));
for i =1:size(parameters_list,1)
    min(i) = handles.MIA_data.clips(strcmp(handles.MIA_data.clips(:,1), parameters_list(i)),2);
    max(i) =  handles.MIA_data.clips(strcmp(handles.MIA_data.clips(:,1), parameters_list(i)),3);
end
handles.clips =  dataset({nominal(parameters_list), 'Parameter'},...
    {cell2num(min)', 'Clip_Min'},...
    {cell2num(max)', 'Clip_Max'});
guidata(hObject, handles);

set(handles.speedy_info_text, 'String', 'Current Prossess : Get/Update data done!');
%update display
IA_update_database_display(hObject, eventdata, handles)

% update clipped by another parameter popupmenu menu
set(handles.speedy_clipped_by_another_parameter_popupmenu, 'String', ['No', cellstr(handles.clips.Parameter)'], 'Value',  1);

% display info
if ~isempty(logbook)
    listdlg('ListString', logbook', 'ListSize',[250 350], 'Name', 'logbook')
end

% --- Executes on selection change in speedy_measure_listbox.
function speedy_measure_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_measure_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_measure_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_measure_listbox


% --- Executes during object creation, after setting all properties.
function speedy_measure_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_measure_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_resolution_popupmenu.
function speedy_resolution_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_resolution_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_resolution_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_resolution_popupmenu


% --- Executes during object creation, after setting all properties.
function speedy_resolution_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_resolution_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function IA_update_database_display(hObject, eventdata, handles)

% get patients info
group_contents = cellstr(get(handles.speedy_group_name_listbox,'String'));
group_name = group_contents(get(handles.speedy_group_name_listbox,'Value'));
if ~strcmp(group_name, 'all')
    for i=1:numel(group_name)
        if i == 1
            tmp = handles.voxels_database_shortlist(handles.voxels_database_shortlist.Group == group_name{i},:);
        else
            tmp =  union(tmp, handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Group == group_name{i},:));
        end
    end
    handles.voxels_database_shortlist_subgroup = tmp;
else
    handles.voxels_database_shortlist_subgroup = handles.voxels_database_shortlist;
    if strcmp(get(hObject, 'Tag'), 'speedy_group_name_listbox')
        set(handles.speedy_group_name_listbox, 'Value',1);
        set(handles.speedy_patient_name_listbox, 'Value', 1);
        set(handles.speedy_time_points_listbox, 'Value', 1);
        set(handles.speedy_parameters_listbox, 'Value', 1);
        set(handles.speedy_ROIs_listbox, 'Value',1);
    end
end

% get patients info
patient_contents = cellstr(get(handles.speedy_patient_name_listbox,'String'));
patient_name = patient_contents(get(handles.speedy_patient_name_listbox,'Value'));
if ~strcmp(patient_name, 'all')
    for i=1:numel(patient_name)
        if i == 1
            tmp =  handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Patient == patient_name{i},:);
        else
            tmp =  union(tmp, handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Patient == patient_name{i},:));
        end
    end
    handles.voxels_database_shortlist_subgroup =tmp;
end
clear tmp
% get timepoint info
timepoint_contents = cellstr(get(handles.speedy_time_points_listbox,'String'));
timepoint_name = timepoint_contents(get(handles.speedy_time_points_listbox,'Value'));
if ~strcmp(timepoint_name, 'all')
    for i=1:numel(timepoint_name)
        if i == 1
            tmp =  handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Timepoint == timepoint_name{i},:);
        else
            tmp =  union(tmp,handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Timepoint == timepoint_name{i},:));
        end
    end
    handles.voxels_database_shortlist_subgroup =tmp;
    
end
clear tmp

% get parameter info
parameter_contents = cellstr(get(handles.speedy_parameters_listbox,'String'));
parameter_name = parameter_contents(get(handles.speedy_parameters_listbox,'Value'));
if ~strcmp(parameter_name, 'all')
    for i=1:numel(parameter_name)
        if i == 1
            tmp = handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Parameter == parameter_name{i},:);
        else
            tmp = union(tmp, handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.Parameter == parameter_name{i},:));
        end
    end
    handles.voxels_database_shortlist_subgroup =tmp;
end
clear tmp

% get ROI info
ROI_contents = cellstr(get(handles.speedy_ROIs_listbox,'String'));
ROI_name = ROI_contents(get(handles.speedy_ROIs_listbox,'Value'));
if ~strcmp(ROI_name, 'all')
    for i=1:numel(ROI_name)
        if i == 1
            tmp = handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.VOI == ROI_name{i},:);
        else
            tmp = union(tmp, handles.voxels_database_shortlist_subgroup(handles.voxels_database_shortlist_subgroup.VOI == ROI_name{i},:));
        end
    end
    handles.voxels_database_shortlist_subgroup = tmp;
end


% update display
group_contents =  ['all' cellstr(unique(handles.voxels_database_shortlist_subgroup.Group)')];
patient_contents = ['all' cellstr(unique(handles.voxels_database_shortlist_subgroup.Patient)')];
timepoint_contents = ['all' cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)')];
parameter_contents =  ['all' cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter)')];
ROI_contents = ['all' cellstr(unique(handles.voxels_database_shortlist_subgroup.VOI)')];
PRM_ref_contents = cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)');

match_value = [];
for i = 1:numel(group_name)
    match_value = [match_value find(strcmp(group_contents, group_name{i}),1)];
end
set(handles.speedy_group_name_listbox, 'String', group_contents, 'Value',  match_value);

match_value = [];
for i = 1:numel(patient_name)
    match_value = [match_value find(strcmp(patient_contents, patient_name{i}),1)];
end
set(handles.speedy_patient_name_listbox, 'String', patient_contents, 'Value',  match_value);

match_value = [];
for i = 1:numel(timepoint_name)
    match_value = [match_value find(strcmp(timepoint_contents, timepoint_name{i}),1)];
end
set(handles.speedy_time_points_listbox, 'String', timepoint_contents, 'Value',  match_value);

match_value = [];
for i = 1:numel(parameter_name)
    match_value = [match_value find(strcmp(parameter_contents, parameter_name{i}),1)];
end
set(handles.speedy_parameters_listbox, 'String', parameter_contents, 'Value',  match_value);


match_value = [];
for i = 1:numel(ROI_name)
    match_value = [match_value find(strcmp(ROI_contents, ROI_name{i}),1)];
end
set(handles.speedy_ROIs_listbox, 'String', ROI_contents, 'Value',  match_value);


if isfield(handles.MIA_data.database(1).databaseinfo, 'cluster')
    clusters_list = ['all' 'No cluster' handles.MIA_data.database(1).databaseinfo.cluster.cluster_name];
else
    clusters_list = [{'all'} {'No cluster'}];
end
set(handles.speedy_clusters_listbox, 'String', clusters_list');

% match_value = [];
% for i = 1:numel(timepoint_name)
%     match_value = [match_value find(strcmp(PRM_ref_contents, timepoint_name{i}),1)];
% end
set(handles.speedy_PRM_ref_popupmenu, 'String', PRM_ref_contents');


% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in speedy_go_button.
function speedy_go_button_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global voxels_database

measure_mode_listing = get(handles.speedy_measure_listbox, 'String');
measure_mode = measure_mode_listing{get(handles.speedy_measure_listbox, 'Value'),1};

if ~isempty(get(handles.speedy_graph1, 'Children'))
    delete(get(handles.speedy_graph1, 'Children'));
    legend(handles.speedy_graph1,'off');
    hold(handles.speedy_graph1, 'off');
end
if get(handles.speedy_clip_to_use_multiple, 'Value')
    warndlg('Multiple clip not coded yet, please select Unique clip','Warning');
    return
end

mean_by_mode = find([get(handles.speedy_mean_by_group, 'Value') get(handles.speedy_mean_by_subject, 'Value')...
    get(handles.speedy_mean_by_slice , 'Value')  get(handles.speedy_histogram_based_by_group, 'Value') get(handles.speedy_histogram_based_by_ROI, 'Value') ...
    get(handles.speedy_4Ddata_by_echo_id, 'Value')  get(handles.speedy_4Ddata_by_echo_group, 'Value')] == 1);

time_point_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint));
table_data = dataset;
legend_txt = {};
set(handles.speedy_info_text, 'String', 'Current Prossess : Plotting....');
switch measure_mode
    case 'Mean+/-SD'
        [table_data, legend_txt, x_label, x_tick, x_lim, XTickLabel,y_label, y_tick, y_lim, yTickLabel] = speedy_measure_by_mean_mode(handles, mean_by_mode, time_point_listing, table_data,legend_txt);
    case 'Correlation'
        if sum(mean_by_mode == [1 3])
            warndlg('This Correlation mode is not coded yet','Warning');
            return
        end
        parameter_nbr = get(handles.speedy_parameters_listbox, 'Value');
        
        if numel(parameter_nbr) ~= 2 || sum(parameter_nbr == 1) ==1
            warndlg('Please selected 2 and only parameters in order to use the "Correlation" mode','Warning');
            return
        end
        [table_data, legend_txt, x_label, x_tick, x_lim, XTickLabel,y_label, y_tick, y_lim, yTickLabel] = speedy_measure_correlation_mode(handles, mean_by_mode, time_point_listing, parameter_nbr, table_data,legend_txt);
    case 'Volume+/-SD'
       [table_data, legend_txt, x_label, x_tick, x_lim, y_tick, y_lim, y_label, XTickLabel, yTickLabel] = speedy_measure_by_volume_mode(handles, mean_by_mode, time_point_listing, table_data,legend_txt);
    case 'PRM'
        parameter_nbr = get(handles.speedy_parameters_listbox, 'Value');
        if numel(parameter_nbr) ~= 1
            warndlg('Please selected only 1 parameters in order to use the "PRM" mode','Warning');
            return
        end
        ROI_nbr = get(handles.speedy_ROIs_listbox, 'Value');
        if numel(ROI_nbr) ~= 1
            warndlg('Please selected only 1 ROI in order to use the "PRM" mode','Warning');
            return
        end
        [table_data, legend_txt, x_label, x_tick, x_lim, XTickLabel,y_label, y_tick, y_lim, yTickLabel] = speedy_measure_by_PRM_mode(handles, mean_by_mode);
    case 'Clustering: GMM'
        speedy_Clustering_GMM_mode(hObject, eventdata, handles);
        return
    case 'Clustering: MMSD'
        speedy_Clustering_MMSD_mode(hObject, eventdata, handles);
        return
end

hold(handles.speedy_graph1, 'off');
% set axes and legend
set(get(handles.speedy_graph1, 'XLabel'), 'String', x_label, 'FontUnits', 'Normalized', 'FontSize', 0.04);
set(handles.speedy_graph1, 'XTick', x_tick, 'FontUnits', 'Normalized', 'FontSize', 0.03)%, 'XLim', x_lim);
set(handles.speedy_graph1, 'XTickLabel', XTickLabel, 'FontUnits', 'Normalized', 'FontSize', 0.03);

set(get(handles.speedy_graph1, 'YLabel'), 'String', y_label, 'FontUnits', 'Normalized', 'FontSize', 0.03);
set(handles.speedy_graph1, 'YTick', y_tick, 'FontUnits', 'Normalized', 'FontSize', 0.03)%, 'YLim', y_lim);
set(handles.speedy_graph1, 'YTickLabel', yTickLabel, 'FontUnits', 'Normalized', 'FontSize', 0.03);
legend(handles.speedy_graph1, legend_txt, 'Location','NorthEast');

% update table
switch measure_mode
    case 'Mean+/-SD'
        if mean_by_mode == 1
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:4)),  num2cell(double(table_data(:,5:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames),...
                'FontSize', 0.05);
        else
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:5)),  num2cell(double(table_data(:,6:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames),...
                'FontSize', 0.05);
        end
    case 'Volume+/-SD'
        if mean_by_mode == 2
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:5)),  num2cell(double(table_data(:,6:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames),...
                'FontSize', 0.05);
        end
    case 'Correlation'
        if mean_by_mode == 5
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:6)),  num2cell(double(table_data(:,7:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames),...
                'FontSize', 0.05);
        elseif mean_by_mode == 2
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:4)),  num2cell(double(table_data(:,5:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames),...
                'FontSize', 0.05);
        else
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:5)),  num2cell(double(table_data(:,6:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames),...
                'FontSize', 0.05);
        end
    case 'PRM'
        if mean_by_mode == 2
            set(handles.speedy_table1, 'data', [cellstr(table_data(:,1:4)),  num2cell(double(table_data(:,5:end)))], 'ColumnName', cellstr(table_data.Properties.VarNames));
        end
end
handles.data_ploted = table_data;
% Update handles structure

guidata(hObject, handles); 


function speedy_Clustering_MMSD_mode(~,~,~)

       
        %Delete temporary directories when function is finished or interrupted
        deleteTemp = onCleanup(@() cleanMeUp());
        
        %Select steps of the pipeline by toggling checkboxes on
        h = figure('name','Select steps of the pipeline');
        uicontrol(h, 'Style', 'Text',...
            'String','Please select the steps of the pipeline to go through:',...
            'Position', [30 260 400 50]);
            
        %Learning mode checkbox
        a = uicontrol(h, 'Style', 'Checkbox',...
            'String','Learning',...
            'Position', [30 250 130 20]);
        
        %Prediction mode checkbox
        b = uicontrol(h, 'Style', 'Checkbox',...
            'String','Prediction',...
            'Position', [30 200 130 20]);
        
        %Comparison mode checkbox
        c = uicontrol(h, 'Style', 'Checkbox',...
            'String','Comparison',...
            'Position', [30 150 130 20]);
        
        %Testing mode checkbox
%         d = uicontrol(h, 'Style', 'Checkbox',...
%             'String','Test',...
%             'Position', [300 250 130 20]);
        
        %Testing/Comparison mode checkbox
%         e = uicontrol(h, 'Style', 'Checkbox',...
%             'String','Test/Comparison',...
%             'Position', [300 200 130 20]);
            
        %OK button
        f = uicontrol(h, 'Style','togglebutton',...
            'Position',[30 30 130 20],...
            'String','OK','CallBack',@ok);
            
        %Cancel button
        g = uicontrol(h, 'Style','togglebutton',...
            'Position',[300 30 130 20],...
            'String','Cancel','CallBack',@cancel);

        
        %Wait for user's response
        uiwait;
        
            %If user closes the window
                if ishandle(h)==0
                    return
                end
        
            %If user presses OK button
                if f.Value == 1
                    learn = a.Value;
                    pred = b.Value;
                    comp = c.Value;
                    %test= d.Value;
                    %testcomp = e.Value;
                    close
                    if learn==0 && pred==0 && comp==0 %&& test==0 && testcomp==0
                        msgbox('No step selected.', 'Warning')
                        return
                    end
            %If user presses 'Cancel' button
                elseif g.Value == 1
                    close
                    return
                end
                
                %Demander à l'utilisateur de régler le nombre de classes
                %(pas fonctionnel)
                
%         prompt = 'How many clusters?';
%         dlg_title = 'Input number of clusters';
%         num_lines = 1;
%         K = inputdlg(prompt,dlg_title,num_lines);
%         fid = fopen('tmp.txt','w');
%         %formatSpec = '%d';
%         fprintf(fid,'%d',K);
%         fclose(fid);
%         !/data/home/sandra/Documents/Stage_Sofiane/Code_R_local/script_local.sh
        
        
        %Create temporary folder
        mkdir('/data/home/sandra/Documents/Stage_Sofiane/temp/');
            
        %Select data for learning part
        if learn == 1
            [filename1,pathname1] = uigetfile('*_ROIContro.csv',...
                'Select files for learning part',...
                '/data/home/sandra/Documents/Stage_Sofiane/Database/Projet_Gliome',...
                'MultiSelect', 'on');
            %If no file is selected
            if isequal(filename1,0)
                quest = questdlg('No file selected for learning part. Continue anyway?',...
                'No file selected.', ...
                'Yes', 'No', 'Cancel', 'Cancel');
                switch quest
                    case 'No'
                        return
                    case 'Cancel'
                        return
                end
            %If data is selected 
            else
                %Create temporary folder
                mkdir('/data/home/sandra/Documents/Stage_Sofiane/temp/learn/');
                %Case 1: one file is selected
                if iscell(filename1) == 0
                    path_to_file = fullfile(pathname1,filename1);
                    copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
                %Case 2: multiple files are selected
                else
                    numfiles = numel(filename1);
                    for i = 1:numfiles
                        path_to_file = fullfile(pathname1, filename1{i});
                        copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
                    end
                end
            end
        end
        
        %Select data for prediction part
        if pred == 1
            [filename2,pathname2] = uigetfile('*_ROIBrain.csv',...
                'Select files for prediction part',...
                '/data/home/sandra/Documents/Stage_Sofiane/Database/Projet_Gliome',...
                'MultiSelect', 'on');
            %If no file is selected
            if isequal(filename1,0)
                quest = questdlg('No file selected for prediction part. Continue anyway?',...
                'No file selected.', ...
                'Yes', 'No', 'Cancel', 'Cancel');
                switch quest
                    case 'No'
                        return
                    case 'Cancel'
                        return
                end
            %If data is selected 
            else
                %Create temporary folder
                mkdir('/data/home/sandra/Documents/Stage_Sofiane/temp/pred/');
                %Case 1: one file is selected
                if iscell(filename2) == 0
                    path_to_file = fullfile(pathname2,filename2);
                    copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
                %Case 2: multiple files are selected
                else
                    numfiles = numel(filename2);
                    for i = 1:numfiles
                        path_to_file = fullfile(pathname2, filename2{i});
                        copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
                    end
                end
            end
        end
        
        %Select data for comparison part
        if comp == 1
            [filename3,pathname3] = uigetfile('*_ROITumorflair.csv',...
                'Select Data Files for comparison part',...
                '/data/home/sandra/Documents/Stage_Sofiane/Database/Projet_Gliome',...
                'MultiSelect', 'on');
            %If no file is selected
            if isequal(filename3,0)
                quest = questdlg('No file selected for comparison part. Continue anyway?',...
                'No file selected.', ...
                'Yes', 'No', 'Cancel', 'Cancel');
                switch quest
                    case 'No'
                        return
                    case 'Cancel'
                        return
                end
            %If data is selected 
            else
                %Create temporary folder
                mkdir('/data/home/sandra/Documents/Stage_Sofiane/temp/comp/');
                %Case 1: one file is selected
                if iscell(filename3) == 0
                    path_to_file = fullfile(pathname3,filename3);
                    copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
                %Case 2: multiple files are selected
                else
                    numfiles = numel(filename3);
                    for i = 1:numfiles
                        path_to_file = fullfile(pathname3, filename3{i});
                        copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
                    end
                end
            end
        end
        

%         %Select data for testing part
%         if test == 1
%             [filename4,pathname4] = uigetfile('*_ROIBrain.csv',...
%                 'Select Data Files for testing part',...
%                 '/data/home/sandra/Documents/Stage_Sofiane/Data',...
%                 'MultiSelect', 'on');
%             %If no file is selected
%             if isequal(filename4,0)
%                 quest = questdlg('No file selected for testing part. Continue anyway?',...
%                 'No file selected.', ...
%                 'Yes', 'No', 'Cancel', 'Cancel');
%                 switch quest
%                     case 'No'
%                         return
%                     case 'Cancel'
%                         return
%                 end
%             %If data is selected 
%             else
%                 %Create temporary folder
%                 mkdir('/data/home/sandra/Documents/Stage_Sofiane/temp/testing/');
%                 %Case 1: one file is selected
%                 if iscell(filename4) == 0
%                     path_to_file = fullfile(pathname4,filename4);
%                     copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
%                 %Case 2: multiple files are selected
%                 else
%                     numfiles = numel(filename4);
%                     for i = 1:numfiles
%                         path_to_file = fullfile(pathname4, filename4{i});
%                         copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
%                     end
%                 end
%             end
%         end
%         
%         %Select data for test/comparison part
%         if testcomp == 1
%             [filename5,pathname5] = uigetfile('*_ROITumorflair.csv',...
%                 'Select Data Files for test/comparison part',...
%                 '/data/home/sandra/Documents/Stage_Sofiane/Data',...
%                 'MultiSelect', 'on');
%             %If no file is selected
%             if isequal(filename5,0)
%                 quest = questdlg('No file selected for test/comparison part. Continue anyway?',...
%                 'No file selected.', ...
%                 'Yes', 'No', 'Cancel', 'Cancel');
%                 switch quest
%                     case 'No'
%                         return
%                     case 'Cancel'
%                         return
%                 end
%             %If data is selected 
%             else
%                 %Create temporary folder
%                 mkdir('/data/home/sandra/Documents/Stage_Sofiane/temp/testcomp/');
%                 %Case 1: one file is selected
%                 if iscell(filename5) == 0
%                     path_to_file = fullfile(pathname5,filename5);
%                     copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
%                 %Case 2: multiple files are selected
%                 else
%                     numfiles = numel(filename5);
%                     for i = 1:numfiles
%                         path_to_file = fullfile(pathname5, filename5{i});
%                         copyfile(path_to_file, '/data/home/sandra/Documents/Stage_Sofiane/temp');
%                     end
%                 end
%             end
%         end
            
        
        %Choice between local and cluster computing
        choice = questdlg('How would you like to perform the clustering?',...
            'Choose Clustering Mode', ...
            'Local', 'IN2P3', 'Cancel', 'Cancel');
        %Handle user response
        switch choice
            case 'Local'
                %Run local shell script
%                 setenv('LD_LIBRARY_PATH',getenv('PATH'));
%                 system('cp -r /data/home/sandra/Documents/Stage_Sofiane/temp /data/home/sandra/Documents/Stage_Sofiane/Code_R_local/Data');
%                 cd('/data/home/sandra/Documents/Stage_Sofiane/Code_R_local');
%                 system('./script_local.sh');
                
                  setenv('LD_LIBRARY_PATH',getenv('PATH'));
                  system('cp -r /data/home/sandra/Documents/Stage_Sofiane/temp /data/home/sandra/Documents/Stage_Sofiane/CodeAlexis/Data');
                  cd('/data/home/sandra/Documents/Stage_Sofiane/CodeAlexis');
                  system('./Script_v2016_02_08.sh');
                  %system('./script_local.sh');
            case 'IN2P3'
                %Run cluster shell script
                system('cp -r /data/home/sandra/Documents/Stage_Sofiane/temp /data/home/sandra/Documents/Stage_Sofiane/Code_R_cluster/Data');
                cd('~/Documents/Stage_Sofiane/Code_R_cluster');
                system('scp -r /data/home/sandra/Documents/Stage_Sofiane/temp skerbal@ccage7.in2p3.fr:/sps/inter/gin/skerbal/Code_R/Data || ./delete_temp.sh');
                system('./submit_job.sh || ./delete_temp.sh');
                
            case 'Cancel'
        end
        

        
        
function speedy_Clustering_GMM_mode(hObject,eventdata,handles)
global voxels_database
%global data_in_table

if isempty(voxels_database.Properties.UserData)
    warning('You chose the original resolution for each image. If those resolutions are differents, it will impact your results.')
end

Normalisation_mode = questdlg('How would you normalize your data?', 'Warning', 'Patient-by-Patient', 'All database', 'None', 'None');
if isempty(Normalisation_mode)
    return
end
Clipping_mode = questdlg('Do you want to clip your data?', 'Warning', 'Yes', 'No', 'No');
if isempty(Clipping_mode)
    return
end
Heuristique_mode = questdlg('Do you want to automaticaly find the optimal number of clusters? (Be careful, it takes time!)', 'Warning', 'Yes', 'No', 'No');
if isempty(Heuristique_mode)
    return
end
if strcmp(Heuristique_mode, 'Yes')
    question = strcat('Slope heuristic : How many clusters do you want to test?');
    NbClassesTestees  = inputdlg(question,'Clustering''s option',1,{'9'});
    if isempty(NbClassesTestees) || ~isnumeric(str2double(NbClassesTestees{:}))
        return
    end
    NbClassesTestees = str2double(NbClassesTestees{:});
end

if strcmp(Heuristique_mode, 'No')
    question = strcat('How many clusters do you want?');
    NbClust = inputdlg(question,'Clustering''s option',1,{'5'});
    if isempty(NbClust) || ~isnumeric(str2double(NbClust{:}))
        return
    end
    NbClust = str2double(NbClust{:});
end


% create the new dataset format
VOI_list = unique(handles.voxels_database_shortlist_subgroup.VOI);
Id_list = unique(handles.voxels_database_shortlist_subgroup.Patient);
Tp_list = unique(handles.voxels_database_shortlist_subgroup.Timepoint);
All_Pameter_list = unique(handles.voxels_database_shortlist_subgroup.Parameter);
data_in_table = dataset;
template_dataset= voxels_database(1,:);
template_dataset.Value = [];
template_dataset.Echo = [];
template_dataset.Expt = [];
template_dataset.Parameter = [];
for ii = 1:numel(All_Pameter_list)
    para_name = char(All_Pameter_list(ii));
    para_name = clean_variable_name(para_name);
    eval(['template_dataset.' para_name '=nan;']);
end
template_dataset(1,:) = [];

% complet the new dataset format
warning off
for i=1:numel(Id_list)
    for j=1:numel(Tp_list)
        for x=1:numel(VOI_list)
            tmp_data = template_dataset;
            tmp_database = voxels_database(voxels_database.Patient == Id_list(i) &...
                voxels_database.Timepoint == Tp_list(j) &...
                voxels_database.VOI == VOI_list(x),:);
            if ~isempty(tmp_database)
                for ii = 1:numel(All_Pameter_list)
                    para_name = char(All_Pameter_list(ii));
                    para_name = clean_variable_name(para_name);
                    if sum(tmp_database.Parameter == All_Pameter_list(ii)) ~= 0
                        tmp = tmp_database(tmp_database.Parameter == All_Pameter_list(ii),:);
                        if isempty(tmp_data)
                            tmp_data.Group = tmp.Group;
                            tmp_data.Patient = tmp.Patient;
                            tmp_data.Timepoint = tmp.Timepoint;
                            tmp_data.VOI = tmp.VOI;
                            tmp_data.Coord_X = tmp.Coord_X;
                            tmp_data.Coord_Y = tmp.Coord_Y;
                            tmp_data.Coord_Z =tmp.Coord_Z;
                            eval(['tmp_data.' para_name '=tmp.Value;']);
                        else
                            if size(tmp_data,1) == size(tmp,1)
                                eval(['tmp_data.' para_name '=tmp.Value;']);
                            end
                        end
                    else
                        if isempty(tmp_data)
                            tmp = tmp_database(tmp_database.Parameter == tmp_database.Parameter(1),:);
                            tmp_data.Group = tmp.Group;
                            tmp_data.Patient = tmp.Patient;
                            tmp_data.Timepoint = tmp.Timepoint;
                            tmp_data.VOI = tmp.VOI;
                            tmp_data.Coord_X = tmp.Coord_X;
                            tmp_data.Coord_Y = tmp.Coord_Y;
                            tmp_data.Coord_Z =tmp.Coord_Z;
                            eval(['tmp_data.' para_name '=nan([size(tmp,1),1]);']);
                        else
                            eval(['tmp_data.' para_name '=nan([size(tmp_data,1),1]);']);  %join
                        end
                    end
                    clear tmp
                end
                %% Clip data
                if strcmp(Clipping_mode, 'Yes')
                    for iii=1:numel(All_Pameter_list)
                        para_name = char(All_Pameter_list(iii));
                        para_name = clean_variable_name(para_name);
                        % tmp_data.ADC(tmp_data.ADC< handles.clips.Clip_Min(handles.clips.Parameter == char(All_Pameter_list(iii))) = NaN;
                        eval(['tmp_data.' para_name '(tmp_data.' para_name '< handles.clips.Clip_Min(handles.clips.Parameter == char(All_Pameter_list(iii)))) = NaN;']);
                        eval(['tmp_data.' para_name '(tmp_data.' para_name '> handles.clips.Clip_Max(handles.clips.Parameter == char(All_Pameter_list(iii)))) = NaN;']);
                    end
                end
                if strcmp(Normalisation_mode, 'Patient-by-Patient')
                    for iii=1:numel(All_Pameter_list)
                        para_name = char(All_Pameter_list(iii));
                        para_name = clean_variable_name(para_name);
                        % Normalize data
                        % data_normalize = data - mean /SD
                        eval(['tmp_mean = nanmean(tmp_data.' para_name ');']);
                        eval(['tmp_std = nanstd(tmp_data.' para_name ');']);
                        eval(['tmp_data.' para_name '=(tmp_data.' para_name ' - tmp_mean )/ tmp_std;']);
                    end
                end
                if isempty(data_in_table)
                    data_in_table = tmp_data;
                else
                    data_in_table = union(data_in_table, tmp_data);
                end
                clear tmp_data
            end
        end
    end
end

% normalize data (data pooled)
% data_normalize = data - mean /SD
if strcmp(Normalisation_mode, 'All database')
    for i=1:numel(All_Pameter_list)
        para_name = char(All_Pameter_list(i));
        para_name = clean_variable_name(para_name);
        eval(['tmp_mean = nanmean(data_in_table.' para_name ');']);
        eval(['tmp_std = nanstd(data_in_table.' para_name ');']);
        eval(['data_in_table.' para_name '=(data_in_table.' para_name ' - tmp_mean )/ tmp_std;']);
    end
end
warning on %#ok<WNON>

% release memory
clear tmp_database
%remove voxel which contains at list 1 NaN value
data_in_table(sum(isnan(double(data_in_table(:,8:size(data_in_table,2)))),2)>0,:)=[];
options = statset ( 'maxiter', 1000);

%% test with GMM
if strcmp(Heuristique_mode, 'Yes')
    
    % L'heuristique de pente utilise le coefficient directeur de le
    % regression lin�aire de la vraisemblance en fonction du
    % nombre de classes du modele. Afin d'avoir une regression
    % significative, on effectue cette regression sur au minimum 5 points.
    % Il faut donc calculer la vraisemblance sur 5 classes de plus que la
    % derni�re � tester.
    ptsheurist = NbClassesTestees + 5;
    
    
    %Vecteur pour stocker la logvraisemblance
    loglike = zeros(1,ptsheurist);
    
    %On stocke les modeles calcules pour ne pas avoir a les recalculer une
    %fois le nombre de classes optimal trouve.
    modeles = cell(1,ptsheurist);
    
    parfor kk=1:ptsheurist
        
        %L'option "Replicate,10" signifie que l'on va calculer 10 fois le
        %modele en modifiant l'initialisation. Le modele renvoye est celui
        %de plus grande vraisemblance.
        modeles{kk} = fitgmdist(double(data_in_table(:,8:size(data_in_table,2))),kk,'Options',options, 'Regularize', 1e-5,'Replicates',10);
        
        loglike(kk) = -modeles{kk}.NegativeLogLikelihood;
        
        %La ligne suivante permet uniquement de suivre l'avancement du
        %calcul des modeles
        disp(strcat('Modele_', num2str(kk)))
    end
    NbCartes = size(data_in_table,2)-7;
    
    %Le vecteur alpha contient les coefficients directeurs des regresions
    %lineaires de la logvraisemblance en fonction du nombre de classes du
    %modele. Sa ieme composante contient le coefficient directeur de la
    %regression lineaire de la log vraisemblance en fonction du nombre de
    %classes du modele en ne prenant pas en compte les i-1 premiers points.
    alpha = zeros(NbClassesTestees,2);
    
    %Le vecteur eqbic contient pour chaque valeur alpha l'equivalent BIC
    %applique a chaque valeur de la log vraisemblance. On obtient donc une
    %matrice ou chaque ligne correspond a l'equivalent BIC applique en
    %chaque valeur de la log vraisemblance pour une valeur de alpha. On
    %passe ainsi d'une ligne a l'autre en modifiant alpha. Dans l'optique
    %de tracer les courbes uniquement a partir du point i, la matrice est initialisee a la valeur NaN.
    eqbic = NaN(NbClassesTestees,length(loglike));
    
    %Le vecteur eqbic2 est similaire au vecteur eqbic mais avec un autre
    %critere.
    eqbic2 = NaN(NbClassesTestees,length(loglike));
    
    for j = 1:NbClassesTestees
        %La regression lineaire
        alpha(j,:) = polyfit(j:ptsheurist,loglike(j:end),1);
        for i=j:length(loglike)
            %eqbic2(j,i) = 2*alpha(j,1)*(i-1+NbCartes*i+(1+NbCartes)*NbCartes/2*i)-loglike(i);
            eqbic(j,i) = 2*alpha(j,1)*i-loglike(i);
        end
    end
    %figure
    %plot(eqbic2.')
    figure
    plot(eqbic.')
    [M,I] = nanmin(eqbic,2);  %Pour chacune des courbes de l'eqbic, l'indice pour lequel le minimum est atteint est considere comme etant le nombre optimal de clusters
    figure
    plot(0:10,0:10,'r')
    hold on
    plot(I,'b')
    k = 0;
    % On vient de tracer le nombre de clusters optimal pour chaque courbe,
    % donc pour chaque coefficient directeur, donc pour chaque point i
    % debut de la regression lineaire. Le nombre de classes optimal global
    % est le point k minimum pour lequel f(k) = k, soit l'intersection de
    % la courbe tracee et de le bissectrice du plan.
    for i=1:length(I)
        if I(i) == i && k == 0
            k = i;
        end
    end
    
    if k == 0
        warndlg('Cannot find an optimal number of cluster, try again and test higher numbers of clusters','Cannot find a number of clusters');
        return
    end
    gmfit = modeles{k};

end


if strcmp(Heuristique_mode, 'No')
    gmfit = fitgmdist(double(data_in_table(:,8:size(data_in_table,2))),NbClust,'Options',options, 'Regularize', 1e-5,'Replicates',1);
    k = NbClust;
end


%On classe tous les pixels en fonction des gaussiennes trouvees par le
%modele
data_in_table.cluster = cluster(gmfit,double(data_in_table(:,8:size(data_in_table,2))));

% save(char(inputdlg('Sous quel nom voulez vous sauvegarder le modèle ?','Sauvegarde Modèle',1,{strcat(num2str(length(unique(data_in_table.Patient))),'Patients','_',num2str(max(data_in_table.cluster(:))),'Clusters',num2str(size(data_in_table,2)-8),'Cartes','Modele')})),'gmfit')
%assignin('base','data_in_table', data_in_table);
%% test to apply a model to new data
% load('model_6maps_k6.mat')
% data_in_table.cluster = cluster(gmfit,double(data_in_table(:,8:size(data_in_table,2))));
% data_in_table( data_in_table.cluster == 1,:) = [];

% test with GMM
%gmfit = fitgmdist(double(data_in_table(:,8:size(data_in_table,2))),k,'Options',options);%, 'Regularize', 1e-5);
%data_in_table.cluster = cluster(gmfit ,double(data_in_table(:,8:size(data_in_table,2))));

%  
% % %% GMM + BIC analyis
% figure;
% for k = 1:9
%     gmfit = fitgmdist(double(data_in_table(:,8:size(data_in_table,2))),k,'Options',options);%, 'Regularize', 1e-5);
% 
%     result(k).model = gmfit;
% end
% plot(cell2mat( arrayfun(@(c) result(c).model.BIC, (1:length(result)).', 'Uniform', 0) ))
% gmfit = result(6).model;

%% save the model
% save('model_6maps_k5.mat', 'gmfit');


%% replace by a white cluster which is the last one in the handles.color_RGB
data_in_table.cluster(isnan(data_in_table.cluster)) = size(handles.color_RGB,1);

Couleurs = handles.color_RGB;
Cartes = nominal(data_in_table.Properties.VarNames(8:end-1));

% A partir du classement des pixels, on calcule les statistiques du
% clustering
[MoyCartesTranches, ProbTranches, MoyCartesVolume, ProbVolume, Ecart_Type_Global, Sign, MoyGlobal] = AnalyseClusterGMM(data_in_table);

ROI = char(unique(data_in_table.VOI));

%On cree 2 strutures que l'on va sauvegarder avec chaque uvascroi. Ces
%structures contiennent les informations et les statistiques du clustering.
Informations = struct('Couleurs',Couleurs,'Cartes', Cartes , 'Modele', gmfit, 'Sign', Sign,'ROI',ROI);
Statistiques = struct('MoyCartesTranches', MoyCartesTranches , 'ProbTranches', ProbTranches , 'MoyCartesVolume', MoyCartesVolume , 'ProbVolume', ProbVolume, 'Ecart_Type_Global', Ecart_Type_Global,'MoyGlobal', MoyGlobal);
% IA_patient_list = {handles.MIA_data.database.name};
NomDossier = [];
for i = 1:length(Informations.Cartes)
    NomDossier = [NomDossier '_' char(Informations.Cartes(i))];
end
NomDossier2 = [num2str(length(Informations.Cartes)) 'Cartes' filesep NomDossier];
logbook = {};
answer = inputdlg('Comment voulez vous nommer ce clustering ?', 'Choix du nom du clustering', 1,{strcat(num2str(k),'C_',NomDossier)});
if ~exist(strcat(handles.MIA_data.database(1).path,NomDossier2), 'dir')
    mkdir(strcat(handles.MIA_data.database(1).path,NomDossier2));
end
save([strcat( handles.MIA_data.database(1).path,NomDossier2), filesep,answer{1} '-clusterGMM.mat'],'Informations', 'Statistiques');

for i = 1:numel(Id_list)
    for j = 1:numel(Tp_list)
        handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));
        tmp_database = data_in_table(data_in_table.Patient == Id_list(i) &...
            data_in_table.Timepoint == Tp_list(j),:);
        if ~isempty(tmp_database)
            
            % load ROI and get info
            patient_nbr = strcmp(cellstr(unique(tmp_database.Patient)), {handles.MIA_data.database.name});
            time_point_nbr =  strcmp(cellstr(unique(tmp_database.Timepoint)),{handles.MIA_data.database(patient_nbr).day.date});
            VOI_nbr  = strcmp(cellstr(unique(tmp_database.VOI)),handles.MIA_data.database(patient_nbr).day(time_point_nbr).VOIs);
            ROI_filename = fullfile(handles.MIA_data.database(1).path, handles.MIA_data.database(patient_nbr).day(time_point_nbr).VOIs_file(VOI_nbr));
            load(ROI_filename{:})
%             [row, colomn] =size(uvascroi(1).value);
            row = voxels_database.Properties.UserData;
             colomn = row;
            z_offset =[uvascroi.fov_offsets];
            z_offset =  z_offset(3,:);
            output_data = nan([row,colomn,numel(uvascroi) 3]);
            for z=1:size(tmp_database,1)
                output_data(tmp_database.Coord_X(z),tmp_database.Coord_Y(z),abs(z_offset-tmp_database.Coord_Z(z))<  1e-5,:) = handles.color_RGB(tmp_database.cluster(z),:);
                
            end
            
            %                 test(:,:,k,:) = squeeze(output_data(:,:,3,:));
            
            for z = 1:numel(z_offset)
                  uvascroi(z).value =squeeze(output_data(:,:,z,:));
            end
            
            
            % save ROI
             [pathstr,name,ext] = fileparts(ROI_filename{:});
            %save([pathstr, filesep, name '-clusterGMM' ext], 'uvascroi','Informations', 'Statistiques');
             save([pathstr, filesep, name '-' answer{1} '-clusterGMM' ext],  'uvascroi','Informations', 'Statistiques');
            % update database
            handles.MIA_data.database(patient_nbr).day(time_point_nbr).VOIs_file = [handles.MIA_data.database(patient_nbr).day(time_point_nbr).VOIs_file {[name '-' answer{1} '-clusterGMM' ext]}];
            handles.MIA_data.database(patient_nbr).day(time_point_nbr).VOIs = [handles.MIA_data.database(patient_nbr).day(time_point_nbr).VOIs strcat(cellstr(unique(tmp_database.VOI)), '-', answer{1},'-clusterGMM')];
            % Check if this new VOI exist already
            % if not update the handles.VOIs field
            if sum(strcmp(handles.MIA_data.VOIs', strcat(cellstr(unique(tmp_database.VOI)), '-clusterGMM'))) == 0
                handles.MIA_data.VOIs =  [handles.MIA_data.VOIs(1:end-2) strcat(cellstr(unique(tmp_database.VOI)), '-', answer{1}, '-clusterGMM') handles.MIA_data.VOIs(end-1:end)];
            end
            
            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
             MIA('MIA_update_database_display', hObject, eventdata, handles.MIA_data);
%             MIA('MIA_update_database_display', hObject, eventdata, findobj('Tag', 'MIA_GUI'));
            tmp =cellstr(unique(tmp_database.VOI));
            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient_nbr).name handles.MIA_data.database(patient_nbr).day(time_point_nbr).date '-' tmp{:}, '-',answer{1}, '-Done'];
        end
    end
end
% end
% for i = 1:size(output_data,3)
%                 figure;imshow3D(test);
% end
if ~isempty(logbook)
    listdlg('ListString', logbook', 'ListSize',[350 350], 'Name', 'logbook');
else
    msgbox('Done', 'logbook') ;
end


%save(strcat(char(data_in_table.Group(1)),'_',char(data_in_table.Patient(1)),'_',char(data_in_table.Timepoint(1)),'_',num2str(max(data_in_table.cluster(:))),'Clusters',num2str(size(data_in_table,2)-8),'Cartes','Signature'),'ProbVolume')
%export(data_in_table,'file','data_in_table.txt');
%save(char(inputdlg('Sous quel nom voulez vous sauvegarder le fichier contenant les données des patients ?','Sauvegarde',1,{strcat(num2str(length(unique(data_in_table.Patient))),'Patients','_',num2str(max(data_in_table.cluster(:))),'Clusters',num2str(size(data_in_table,2)-8),'Cartes')})),'data_in_table')
%save(char(inputdlg('Sous quel nom voulez vous sauvegarder le fichier contenant les informations necessaires à l''affichage et l''analyse des résultats ?','Sauvegarde',1,{strcat(num2str(length(unique(data_in_table.Patient))),'Patients','_',num2str(max(data_in_table.cluster(:))),'Clusters',num2str(size(data_in_table,2)-8),'Cartes')})),'uvascroi','gmfit','Couleurs','Cartes','Patients','Groupes','ProbVolume','MoyCartesVolume','ProbTranches','MoyCartesTranches');


return






function [table_data, legend_txt, x_label, x_tick, x_lim, XTickLabel,y_label, y_tick, y_lim, yTickLabel] = speedy_measure_by_PRM_mode(handles, mean_by_mode)
global voxels_database

PRM_CI = str2double(get(handles.edit3, 'String'));
if PRM_CI == 0
    warndlg('Please define a correct CI','Warning')
    return
end
PRM_ref = get(handles.speedy_PRM_ref_popupmenu, 'Value');
if PRM_ref ==  numel(get(handles.speedy_PRM_ref_popupmenu, 'String'))
    warndlg('The PRM ref = the last day!','Warning')
    return
end
number_of_tp = PRM_ref+1 : numel(get(handles.speedy_PRM_ref_popupmenu, 'String'));
table_data = dataset;

switch mean_by_mode
    case {1, 2}
        group_listing = nominal(cellstr(unique(handles.voxels_database_shortlist_subgroup.Group)));
        patient_listing = nominal(cellstr(unique(handles.voxels_database_shortlist_subgroup.Patient)));
        Timepoint_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint));
        parameter_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        VOI_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.VOI));
        
        for g = 1:numel(group_listing)
            mean_by_group_data_r = nan([numel(patient_listing) numel(Timepoint_listing)-1]);
            mean_by_group_data_g = nan([numel(patient_listing) numel(Timepoint_listing)-1]);
            mean_by_group_data_b = nan([numel(patient_listing) numel(Timepoint_listing)-1]);
            for i = 1:numel(patient_listing)
                clear tmp_data
                %create sub_voxels_database
                sub_voxels_database = voxels_database(voxels_database.Patient == patient_listing(i) & voxels_database.Parameter == parameter_listing(1),:);
                
                VOI_list = unique(handles.voxels_database_shortlist_subgroup.VOI);
                Id_list = unique(handles.voxels_database_shortlist_subgroup.Patient);
                data_for_1_patient= sub_voxels_database(1,1:end-1);
                data_for_1_patient.Timepoint = [];
                data_for_1_patient.Parameter = [];
                data_for_1_patient.Echo = [];
                data_for_1_patient.Expt = [];
                data_for_1_patient(1,:) = [];
                
                % complet the new dataset format
                warning off
                
                for x=1:numel(VOI_listing)
                    for ii = 1:numel(Timepoint_listing)
                        tmp_database = sub_voxels_database(sub_voxels_database.Patient == Id_list(i) &...
                            sub_voxels_database.Timepoint == Timepoint_listing(ii) &...
                            sub_voxels_database.VOI == VOI_list(x),:);
                        
                        
                        
                        if ~isempty(tmp_database)
                            para_name = char(Timepoint_listing(ii));
                            para_name = clean_variable_name(para_name);
                            
                            if isempty(data_for_1_patient)
                                data_for_1_patient.Group = tmp_database.Group;
                                data_for_1_patient.Patient = tmp_database.Patient;
                                data_for_1_patient.VOI = tmp_database.VOI;
                                data_for_1_patient.Coord_X = tmp_database.Coord_X;
                                data_for_1_patient.Coord_Y = tmp_database.Coord_Y;
                                data_for_1_patient.Coord_Z =tmp_database.Coord_Z;
                                eval(['data_for_1_patient.' para_name '=tmp_database.Value;']);
                                
                            else
                                if size(data_for_1_patient,1) == size(tmp_database,1)
                                    eval(['data_for_1_patient.' para_name '=tmp_database.Value;']);
                                else
                                    %% problem when the number of voxel pre and post differ!!
                                    %% this code works but is really slow!!
                                    %                                 if size(data_for_1_patient,1) < size(tmp_database,1)
                                    %                                     for xx =  1:size(data_for_1_patient,1)
                                    %                                         tmp = tmp_database.Value(tmp_database.Coord_X == data_for_1_patient.Coord_X(xx) & tmp_database.Coord_Y == data_for_1_patient.Coord_Y(xx) & tmp_database.Coord_Z == data_for_1_patient.Coord_Z(xx));
                                    %                                         if ~isempty(tmp)
                                    %                                             eval(['data_for_1_patient.' para_name '(xx) = tmp_database.Value(tmp_database.Coord_X == data_for_1_patient.Coord_X(xx) & tmp_database.Coord_Y == data_for_1_patient.Coord_Y(xx) & tmp_database.Coord_Z == data_for_1_patient.Coord_Z(xx));']);
                                    %                                         end
                                    %                                     end
                                    %                                 else
                                    %                                     for xx =  1:size(tmp_database,1)
                                    %                                        tmp = data_for_1_patient(data_for_1_patient.Coord_X(:) == tmp_database.Coord_X(xx) & data_for_1_patient.Coord_Y(:) == tmp_database.Coord_Y(xx) & ...
                                    %                                             data_for_1_patient.Coord_Z(:) == tmp_database.Coord_Z(xx),:);
                                    %                                         if ~isempty(tmp)
                                    %                                              eval(['data_for_1_patient.' para_name '(data_for_1_patient.Coord_X(:) == tmp_database.Coord_X(xx) & data_for_1_patient.Coord_Y(:) == tmp_database.Coord_Y(xx) & data_for_1_patient.Coord_Z(:) == tmp_database.Coord_Z(xx)) = tmp_database.Value(xx);']);
                                    %
                                    %                                         end
                                    %                                     end
                                    %
                                    %                                 end
                                end
                            end
                        end
                    end
                end
                Ref_name = char(Timepoint_listing(PRM_ref));
                Ref_name = clean_variable_name(Ref_name);
                clear data_for_PRM
                data_for_PRM = dataset;
                data_to_plot = nan([3, numel(number_of_tp)]);
                for ii = 1:numel(number_of_tp) %#ok<FXSET>
                    % ref = TP1
                    eval(['data_for_PRM.TP1 = data_for_1_patient.' Ref_name ';']);
                    % Post = TP2
                    para_name = char(Timepoint_listing(number_of_tp(ii)));
                    para_name = clean_variable_name(para_name);
                    if sum(strcmp(para_name, data_for_1_patient.Properties.VarNames')) ~= 0
                        eval(['data_for_PRM.TP2 = data_for_1_patient.' para_name ';']);
                    else
                        if ~isempty(table_data)
                            table_data.Group(size(table_data,1)+1,:) = {char(unique(data_for_1_patient.Group))};
                            table_data.Patient(size(table_data,1),:) = {char(unique(data_for_1_patient.Patient))};
                            table_data.VOI(size(table_data,1),:) = {char(unique(data_for_1_patient.VOI))};
                            table_data.Parameter(size(table_data,1),:) = {strcat(char(Timepoint_listing(PRM_ref)), 'vs.', char(Timepoint_listing(number_of_tp(ii))))};
                            table_data.PRMred(size(table_data,1),:) = NaN;
                            table_data.PRMgreen(size(table_data,1),:) = NaN;
                            table_data.PRMblue(size(table_data,1),:) =NaN;
                        else
                            table_data.Group = {char(unique(data_for_1_patient.Group))};
                            table_data.Patient = {char(unique(data_for_1_patient.Patient))};
                            table_data.VOI = {char(unique(data_for_1_patient.VOI))};
                            table_data.Parameter = {strcat(char(Timepoint_listing(PRM_ref)), 'vs.', char(Timepoint_listing(number_of_tp(ii))))};
                            table_data.PRMred = NaN;
                            table_data.PRMgreen = NaN;
                            table_data.PRMblue = NaN;
                        end
                        continue
                        
                    end
                    % remove NaNand clip data
                    data_for_PRM(isnan(data_for_PRM.TP1) | isnan(data_for_PRM.TP2), :) = [];
                    data_for_PRM(data_for_PRM.TP1 < handles.clips.Clip_Min(handles.clips.Parameter == parameter_listing ) | data_for_PRM.TP2 < handles.clips.Clip_Min(handles.clips.Parameter == parameter_listing ),:)= [];
                    data_for_PRM(data_for_PRM.TP2 > handles.clips.Clip_Max(handles.clips.Parameter == parameter_listing ) | data_for_PRM.TP2 > handles.clips.Clip_Max(handles.clips.Parameter == parameter_listing ),:)= [];
                    data_for_PRM.diff = data_for_PRM.TP2 - data_for_PRM.TP1;
                    % calculate PRM
                    data_for_PRM.blue = data_for_PRM.diff < -PRM_CI;
                    data_for_PRM.green =  abs(data_for_PRM.diff) < PRM_CI;
                    data_for_PRM.red = data_for_PRM.diff > PRM_CI;
                    
                    if ~isempty(table_data)
                        table_data.Group(size(table_data,1)+1,:) = {char(unique(data_for_1_patient.Group))};
                        table_data.Patient(size(table_data,1),:) = {char(unique(data_for_1_patient.Patient))};
                        table_data.VOI(size(table_data,1),:) = {char(unique(data_for_1_patient.VOI))};
                        table_data.Parameter(size(table_data,1),:) = {strcat(char(Timepoint_listing(PRM_ref)), 'vs.', char(Timepoint_listing(number_of_tp(ii))))};
                        table_data.PRMred(size(table_data,1),:) = sum(data_for_PRM.red)/size(data_for_PRM,1)*100;
                        table_data.PRMgreen(size(table_data,1),:) = sum(data_for_PRM.green)/size(data_for_PRM,1)*100;
                        table_data.PRMblue(size(table_data,1),:) = sum(data_for_PRM.blue)/size(data_for_PRM,1)*100;
                        
                    else
                        table_data.Group = {char(unique(data_for_1_patient.Group))};
                        table_data.Patient = {char(unique(data_for_1_patient.Patient))};
                        table_data.VOI = {char(unique(data_for_1_patient.VOI))};
                        table_data.Parameter = {strcat(char(Timepoint_listing(PRM_ref)), 'vs.', char(Timepoint_listing(number_of_tp(ii))))};
                        table_data.PRMred = sum(data_for_PRM.red)/size(data_for_PRM,1)*100;
                        table_data.PRMgreen = sum(data_for_PRM.green)/size(data_for_PRM,1)*100;
                        table_data.PRMblue = sum(data_for_PRM.blue)/size(data_for_PRM,1)*100;
                    end
                    data_to_plot(1,ii) = sum(data_for_PRM.red)/size(data_for_PRM,1)*100;
                    data_to_plot(2,ii) = sum(data_for_PRM.green)/size(data_for_PRM,1)*100;
                    data_to_plot(3,ii) = sum(data_for_PRM.blue)/size(data_for_PRM,1)*100;
                end
                
                if mean_by_mode == 2
                    if isempty(get(handles.speedy_graph1, 'Children'))
                        hold(handles.speedy_graph1, 'off');
                        plot(handles.speedy_graph1,   number_of_tp, data_to_plot(1,:), 'color', 'r');
                        hold(handles.speedy_graph1, 'on');
                    else
                        plot(handles.speedy_graph1,   number_of_tp, data_to_plot(1,:), 'color', 'r');
                    end
                    plot(handles.speedy_graph1,   number_of_tp, data_to_plot(3,:), 'color', 'b');
                elseif mean_by_mode == 1
                    mean_by_group_data_r(size(mean_by_group_data_r,1)+1,:) = data_to_plot(1,:);
                    mean_by_group_data_g(size(mean_by_group_data_g,1)+1,:) = data_to_plot(2,:);
                    mean_by_group_data_b(size(mean_by_group_data_b,1)+1,:) = data_to_plot(3,:);
                end
            end
            if mean_by_mode == 1
                if isempty(get(handles.speedy_graph1, 'Children'))
                    hold(handles.speedy_graph1, 'off');
                    errorbar(handles.speedy_graph1,number_of_tp, nanmean(mean_by_group_data_r), nanstd(mean_by_group_data_r), 'Color', 'r');
                    hold(handles.speedy_graph1, 'on');
                else
                    errorbar(handles.speedy_graph1,number_of_tp, nanmean(mean_by_group_data_r), nanstd(mean_by_group_data_r), 'Color', 'r');
                end
                errorbar(handles.speedy_graph1,number_of_tp, nanmean(mean_by_group_data_b), nanstd(mean_by_group_data_b), 'Color', 'b');
            end
        end
        
        x_label = 'Time Point';
        x_tick =     number_of_tp(1): number_of_tp(end);
        x_lim =   number_of_tp(1)-1: number_of_tp(end)+1;
        XTickLabel = x_tick';
        
        y_label = ['PRM-' parameter_listing{1}];
        set(get(handles.speedy_graph1, 'YLabel'), 'String', y_label, 'FontUnits', 'Normalized', 'FontSize', 0.03);
        y_tick = get(handles.speedy_graph1, 'YTick');
        y_lim= get(handles.speedy_graph1, 'YLim');
        yTickLabel = get(handles.speedy_graph1, 'yTickLabel');
        tmp= dataset2cell(unique(table_data(:,1:3)));
        tmp = tmp(2:end,:);
        if mean_by_mode == 2
            legend_txt = cell([size(tmp,1)*2 1]);
            for y = 1:size(tmp,1)
                legend_txt(y*2-1,:) = {strcat(tmp{y,:}, 'PRM+')};
                legend_txt(y*2,:) = {strcat(tmp{y,:}, 'PRM-')};
            end
        else
            legend_txt = cell([numel(group_listing)*2 1]);
            for y=1: numel(group_listing)
                legend_txt(y*2-1,:) = strcat(cellstr(group_listing(y,:)), 'PRM+');
                legend_txt(y*2,:) =  strcat(cellstr(group_listing(y,:)), 'PRM-');
            end
        end
        
end




function   [table_data, legend_txt, x_label, x_tick, x_lim, XTickLabel,y_label, y_tick, y_lim, yTickLabel] = speedy_measure_correlation_mode(handles, mean_by_mode, time_point_listing,parameter_nbr, table_data,legend_txt)
global voxels_database
switch mean_by_mode
    
    
    case {4, 2}
        % case 4 Correlation of 2 parameters for the pooled voxels inside the same group/ROI
        % case 2 % Correlation of 2 parameters averaged subject by subject
        group_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Group));
        parameter_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        Timepoint_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint));
        patient_listing = nominal(cellstr(unique(handles.voxels_database_shortlist_subgroup.Patient)));
        VOI_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.VOI));
        for i = 1:numel(group_listing)
            for j = 1:numel(time_point_listing)
                clear tmp_data
                if get(handles.speedy_patient_name_listbox, 'Value') == 1
                    tmp_data = voxels_database;
                else
                    for ii = 1:numel(patient_listing)
                        if ii == 1
                            tmp_data =voxels_database(voxels_database.Patient == patient_listing(ii),:);
                        else
                            tmp_data = union(tmp_data, voxels_database(voxels_database.Patient == patient_listing(ii),:));
                        end
                    end
                end
                for jj =  1:numel(VOI_listing)
                    if i*j*jj > 1
                        hold(handles.speedy_graph1, 'on');
                    end
                    x_data = tmp_data(tmp_data.Group==group_listing{i} &...
                        tmp_data.Timepoint==time_point_listing{j}&...
                        tmp_data.Group==group_listing{i}&...
                        tmp_data.Parameter == parameter_listing{1} &...
                        tmp_data.VOI ==  VOI_listing{jj},:);
                    y_data = tmp_data(tmp_data.Group==group_listing{i} &...
                        tmp_data.Timepoint==time_point_listing{j}&...
                        tmp_data.Group==group_listing{i}&...
                        tmp_data.Parameter == parameter_listing{2} &...
                        tmp_data.VOI ==  VOI_listing{jj}, :);
                    
                    if isempty(y_data)
                        continue
                    end
                    % refine dataset if the ROI is not present on the same
                    % number of slices
                    if size(x_data,1) < size(y_data,1)
                        y_data= y_data(ismember(y_data(:,[1:3, 5:10]), x_data(:,[1:3, 5:10])),:);
                        x_data= x_data(ismember(x_data(:,[1:3, 5:10]), y_data(:,[1:3, 5:10])),:);
                    elseif size(x_data,1) > size(y_data,1)
                        x_data= x_data(ismember(x_data(:,[1:3, 5:10]), y_data(:,[1:3, 5:10])),:);
                        y_data= y_data(ismember(y_data(:,[1:3, 5:10]), x_data(:,[1:3, 5:10])),:);
                    end
                    
                    if isempty(x_data) || isempty(y_data)
                        warndlg('No commun voxels has been found! Did both parameters have the same resolution?','Warning');
                        continue
                    end
                    % clip unique
                    if get(handles.speedy_clip_to_use_unique, 'Value')
                        match_voxel =...
                            x_data.Value >= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),2)) &...
                            x_data.Value <= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),3)) & ...
                            y_data.Value >= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),2)) &...
                            y_data.Value <= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),3));
                    else
                        
                    end
                    
                    if mean_by_mode == 4
                        color_map = rand(1,3);
                        scatter(handles.speedy_graph1, x_data.Value(match_voxel), y_data.Value(match_voxel),...
                            'filled',...
                            'SizeData', 20,...
                            'MarkerFaceColor',color_map,...
                            'MarkerEdgeColor',color_map,...
                            'Visible', 'on');
                        hold(handles.speedy_graph1, 'on');
                        mdl =LinearModel.fit(x_data.Value(match_voxel), y_data.Value(match_voxel),'linear');
                        b0=3; b1=4;
                        x= linspace(min(x_data.Value(match_voxel)),max(x_data.Value(match_voxel)), 5); % Adapt n for resolution of graph
                        y= mdl.Coefficients.Estimate(2)*x+mdl.Coefficients.Estimate(1);
                        
                        plot(handles.speedy_graph1,x,y,'Color',color_map)
                        
                        %%%%%%%% activate clic on graph
                        set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
                        set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
                        set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
                        %                         mdl.Rsquared
                        
                        % save data in a table
                        if i*j*jj == 1
                            table_data.Group = nominal(group_listing{i});
                            table_data.Timepoint = nominal(time_point_listing{j});
                            table_data.VOI = nominal(VOI_listing{jj});
                            table_data.Parameter1 = nominal(parameter_listing(1));
                            table_data.Parameter2 = nominal(parameter_listing(2));
                            table_data.Voxels_nbr =  numel(x_data.Value);
                            table_data.Voxels_used=  sum(match_voxel== 1);
                            table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
                            table_data.R2 =  mdl.Rsquared.Adjusted;
                            table_data.Pvalue = mdl.Coefficients.pValue(1);
                        else
                            tmp = dataset;
                            tmp.Group = nominal(group_listing{i});
                            tmp.Timepoint = nominal(time_point_listing{j});
                            tmp.VOI = nominal(VOI_listing{jj});
                            tmp.Parameter1 = nominal(parameter_listing(1));
                            tmp.Parameter2 = nominal(parameter_listing(2));
                            tmp.Voxels_nbr =  numel(x_data.Value);
                            tmp.Voxels_used=  sum(match_voxel== 1);
                            tmp.Pc_voxels_excluded =  ((tmp.Voxels_nbr - tmp.Voxels_used) / tmp.Voxels_nbr)*100;
                            tmp.R2 =  mdl.Rsquared.Adjusted;
                            tmp.Pvalue = mdl.Coefficients.pValue(1);
                            table_data = union(table_data,tmp);
                        end
                        tmp = [group_listing{i} '-' time_point_listing{j} '-' VOI_listing{jj}];
                        legend_txt(size(legend_txt,1)+1,1) ={tmp};
                        tmp = [group_listing{i} '-' time_point_listing{j} '-' VOI_listing{jj} '-fit'];
                        legend_txt(size(legend_txt,1)+1,1) ={tmp};
                        
                        
                    else
                        tmp = unique(x_data);
                        x_table_data = unique(tmp(:,[1:3, 5]));
                        y_table_data = unique(tmp(:,[1:3, 5]));
                        warning('off')
                        for x = 1:size(x_table_data,1)
                            x_table_data.Value(x) = mean(x_data.Value(ismember(x_data(:,[1:3, 5]), x_table_data(x,1:4))&match_voxel,:));
                            x_table_data.SD(x)= std(x_data.Value(ismember(x_data(:,[1:3, 5]), x_table_data(x,1:4))&match_voxel,:));
                            x_table_data.Voxels_nbr(x) =  numel(x_data.Value(ismember(x_data(:,[1:3, 5]), x_table_data(x,1:4)),:));
                            x_table_data.Voxels_used(x)= numel(x_data.Value(ismember(x_data(:,[1:3, 5]), x_table_data(x,1:4))&match_voxel,:));
                            x_table_data.Pc_voxels_excluded(x) =  ((x_table_data.Voxels_nbr(x) - x_table_data.Voxels_used(x)) / x_table_data.Voxels_nbr(x))*100;
                            
                            y_table_data.Value(x)= mean(y_data.Value(ismember(y_data(:,[1:3, 5]), y_table_data(x,1:4))&match_voxel,:));
                            y_table_data.SD(x)= std(y_data.Value(ismember(y_data(:,[1:3, 5]), y_table_data(x,1:4))&match_voxel,:));
                            y_table_data.Voxels_nbr(x) =  numel(y_data.Value(ismember(x_data(:,[1:3, 5]), y_table_data(x,1:4)),:));
                            y_table_data.Voxels_used(x)= numel(y_data.Value(ismember(x_data(:,[1:3, 5]), y_table_data(x,1:4))&match_voxel,:));
                            y_table_data.Pc_voxels_excluded(x) =  ((y_table_data.Voxels_nbr(x) - y_table_data.Voxels_used(x)) / y_table_data.Voxels_nbr(x))*100;
                        end
                        warning('stats:dataset:subsasgn:DefaultValuesAddedVariable', 'on')
                        parameter_listing_for_table = clean_variable_name(parameter_listing);
                        if i*j*jj == 1
                            
                            table_data=  x_table_data(:,1:4);
                            eval(['table_data.', regexprep(parameter_listing_for_table{1},'-','_'), '=x_table_data.Value;']);
                            eval(['table_data.', [regexprep(parameter_listing_for_table{1},'-','_') '_SD'], '=x_table_data.SD;']);
                            eval(['table_data.', regexprep(parameter_listing_for_table{2},'-','_'), '=y_table_data.Value;']);
                            eval(['table_data.', [regexprep(parameter_listing_for_table{2},'-','_') '_SD'], '=y_table_data.SD;'])
                            table_data.Vox_nbr=x_table_data.Voxels_nbr;
                            table_data.Vox_used=x_table_data.Voxels_nbr;
                            table_data.pc_Vox_ex=x_table_data.Pc_voxels_excluded;
                        else
                            tmp =  x_table_data(:,1:4);
                            eval(['tmp.',regexprep(parameter_listing_for_table{1},'-','_'), '=x_table_data.Value;']);
                            eval(['tmp.', [regexprep(parameter_listing_for_table{1},'-','_') '_SD'], '=x_table_data.SD;']);
                            eval(['tmp.', regexprep(parameter_listing_for_table{2},'-','_'), '=y_table_data.Value;']);
                            eval(['tmp.', [regexprep(parameter_listing_for_table{2},'-','_') '_SD'], '=y_table_data.SD;'])
                            tmp.Vox_nbr=x_table_data.Voxels_nbr;
                            tmp.Vox_used=x_table_data.Voxels_nbr;
                            tmp.pc_Vox_ex=x_table_data.Pc_voxels_excluded;
                            
                            table_data = union(table_data,tmp);
                        end
                    end
                    
                end
            end
        end
        if mean_by_mode == 2
            color_map = rand(1,3);
            parameter_listing_for_table = clean_variable_name(parameter_listing);
            scatter(handles.speedy_graph1,   eval(['table_data.', regexprep(parameter_listing_for_table{1},'-','_')]), eval(['table_data.', regexprep(parameter_listing_for_table{2},'-','_')]),...
                'filled',...
                'SizeData', 20,...
                'MarkerFaceColor',color_map,...
                'MarkerEdgeColor',color_map,...
                'Visible', 'on');
            hold(handles.speedy_graph1, 'on');
            mdl =LinearModel.fit(eval(['table_data.', regexprep(parameter_listing_for_table{1},'-','_')]), eval(['table_data.', regexprep(parameter_listing_for_table{2},'-','_')]),'linear');
            b0=3; b1=4;
            x= linspace(min(eval(['table_data.', regexprep(parameter_listing_for_table{1},'-','_')])),max(eval(['table_data.',regexprep(parameter_listing_for_table{1},'-','_')])), 5); % Adapt n for resolution of graph
            y= mdl.Coefficients.Estimate(2)*x+mdl.Coefficients.Estimate(1);
            plot(handles.speedy_graph1,x,y,'Color',color_map)
            table_data.R2 =  repmat(mdl.Rsquared.Adjusted, [size( table_data,1),1]);
            table_data.Pvalue = repmat(mdl.Coefficients.pValue(1), [size( table_data,1),1]);
            
        end
        x_label = parameter_listing(1);
        x_tick =    double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),2)): double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),3))/5 : double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),3));
        x_lim = [double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),2)) double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),3))];
        XTickLabel =  x_tick';
        
        y_label = parameter_listing(2);
        y_tick =    double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),2)): double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),3))/5 : double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),3));
        y_lim = [double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),2)) double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),3))];
        yTickLabel =  y_tick';
        
    case 5 % Correlation of 2 parameters inside the same ROI
        group_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Group));
        patient_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Patient));
        parameter_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        VOI_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.VOI));
        for i = 1:numel(group_listing)
            for j = 1:numel(time_point_listing)
                for ii = 1:numel(patient_listing)
                    for jj =  1:numel(VOI_listing)
                        if i*j*ii*jj > 1
                            hold(handles.speedy_graph1, 'on');
                        end
                        x_data = voxels_database(voxels_database.Group==group_listing{i} &...
                            voxels_database.Patient==patient_listing{ii}&...
                            voxels_database.Timepoint==time_point_listing{j} &...
                            voxels_database.Parameter == parameter_listing{1} &...
                            voxels_database.VOI ==  VOI_listing{jj}, :);
                        y_data = voxels_database(voxels_database.Group==group_listing{i} &...
                            voxels_database.Patient==patient_listing{ii}&...
                            voxels_database.Timepoint==time_point_listing{j} &...
                            voxels_database.Parameter == parameter_listing{2} &...
                            voxels_database.VOI ==  VOI_listing{jj}, :);
                        if isempty(y_data)
                            continue
                        end
                        % refine dataset if the ROI is not present on the same
                        % number of slices
                        if size(x_data,1) < size(y_data,1)
                            y_data= x_data(ismember(y_data(:,[1:3, 5:10]), x_data(:,[1:3, 5:10])),:);
                        elseif size(x_data,1) > size(y_data,1)
                            x_data= x_data(ismember(x_data(:,[1:3, 5:10]), y_data(:,[1:3, 5:10])),:);
                        end
                        
                        % clip unique
                        if get(handles.speedy_clip_to_use_unique, 'Value')
                            match_voxel =x_data.Value >= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),2)) &...
                                x_data.Value <= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(1)),3)) & ...
                                y_data.Value >= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),2)) &...
                                y_data.Value <= double(handles.clips(handles.clips.Parameter ==  nominal(parameter_listing(2)),3));
                        else
                            
                        end
                        color_map = rand(1,3);
                        
                        scatter(handles.speedy_graph1, x_data.Value(match_voxel), y_data.Value(match_voxel),...
                            'filled',...
                            'SizeData', 20,...
                            'MarkerFaceColor',color_map,...
                            'MarkerEdgeColor',color_map,...
                            'Visible', 'on');
                        hold(handles.speedy_graph1, 'on');
                        mdl =LinearModel.fit(x_data.Value(match_voxel), y_data.Value(match_voxel),'linear');
                        b0=3; b1=4;
                        x= linspace(min(x_data.Value(match_voxel)),max(x_data.Value(match_voxel)), 5); % Adapt n for resolution of graph
                        y= mdl.Coefficients.Estimate(2)*x+mdl.Coefficients.Estimate(1);
                        plot(handles.speedy_graph1,x,y,'Color',color_map)
                        
                        %%%%%%%% activate clic on graph
                        set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
                        set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
                        set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
                        %                         mdl.Rsquared
                        
                        % save data in a table
                        if isempty(table_data)
                            table_data.Group = nominal(group_listing{i});
                            table_data.Patient  = nominal(patient_listing{ii});
                            table_data.Timepoint = nominal(time_point_listing{j});
                            table_data.VOI = nominal(VOI_listing{jj});
                            table_data.Parameter1 = nominal(parameter_listing(1));
                            table_data.Parameter2 = nominal(parameter_listing(2));
                            table_data.Voxels_nbr =  numel(x_data.Value);
                            table_data.Voxels_used=  sum(match_voxel== 1);
                            table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
                            table_data.R2 =  mdl.Rsquared.Adjusted;
                            table_data.Pvalue = mdl.Coefficients.pValue(1);
                        else
                            tmp = dataset;
                            tmp.Group = nominal(group_listing{i});
                            tmp.Patient  = nominal(patient_listing{ii});
                            tmp.Timepoint = nominal(time_point_listing{j});
                            tmp.VOI = nominal(VOI_listing{jj});
                            tmp.Parameter1 = nominal(parameter_listing(1));
                            tmp.Parameter2 = nominal(parameter_listing(2));
                            tmp.Voxels_nbr =  numel(x_data.Value);
                            tmp.Voxels_used=  sum(match_voxel== 1);
                            tmp.Pc_voxels_excluded =  ((tmp.Voxels_nbr - tmp.Voxels_used) / tmp.Voxels_nbr)*100;
                            tmp.R2 =  mdl.Rsquared.Adjusted;
                            tmp.Pvalue = mdl.Coefficients.pValue(1);
                            table_data = union(table_data,tmp);
                        end
                        tmp = [group_listing{i} '-' patient_listing{i} time_point_listing{j} '-' VOI_listing{jj}];
                        legend_txt(size(legend_txt,1)+1,1) ={tmp};
                        tmp = [group_listing{i} '-' patient_listing{i} time_point_listing{j} '-' VOI_listing{jj} '-fit'];
                        legend_txt(size(legend_txt,1)+1,1) ={tmp};
                    end
                end
            end
        end
        x_label = parameter_listing(1);
        x_tick = get(handles.speedy_graph1, 'XTick');
        x_lim =  get(handles.speedy_graph1, 'XLim');
        XTickLabel =  get(handles.speedy_graph1, 'XTickLabel');
        y_label = parameter_listing(2);
        
    otherwise
        return
end

function[table_data, legend_txt, x_label, x_tick, x_lim, y_tick, y_lim, y_label, XTickLabel, yTickLabel] = speedy_measure_by_volume_mode(handles, mean_by_mode, time_point_listing, table_data,legend_txt)
global voxels_database
y_tick = [];
y_lim = '';
yTickLabel = '';
y_label = '';
voxel_size = inputdlg({'Voxel size in µm3? (default 234x234x800'}, 'Bip', [1 30], {'0.0438'});
if isempty(voxel_size)
    return
end
voxel_size = str2double(voxel_size);
% The volume function does not need to compute the ROI volume on each
% Parameter, so we keep only on iteration of each conditions
% This code does no delete the original shorlist_subgroup (no guidata in this function)
handles.voxels_database_shortlist_subgroup.Parameter(1:end) = handles.voxels_database_shortlist_subgroup.Parameter(1);
handles.voxels_database_shortlist_subgroup = unique(handles.voxels_database_shortlist_subgroup);

switch mean_by_mode
    case 1 % mean by Group
        % 1st mean by ROI
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            tmp_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            % clip unique
            if get(handles.speedy_clip_to_use_unique, 'Value')
                y_data = size(tmp_data, 1)*voxel_size;
                sd_data =0;
            else
                
            end
            x_data=find(strcmp(time_point_listing, cellstr(handles.voxels_database_shortlist_subgroup.Timepoint(i))) == 1)+1;
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            if i == 1
                tmp = unique(tmp_data);
                table_data = unique(tmp(:,1:5));
                table_data.Value = y_data;
                table_data.SD = sd_data;
                table_data.Voxels_nbr =  numel(tmp_data.Value);
                table_data.Voxels_used=  sum(double(tmp_data(tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                    tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
            else
                tmp = unique(tmp_data);
                tmp2 = unique(tmp(:,1:5));
                tmp2.Value = y_data;
                tmp2.SD = sd_data;
                tmp2.Voxels_nbr =  numel(tmp_data.Value);
                tmp2.Voxels_used=  sum(double(tmp_data(tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                    tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                tmp2.Pc_voxels_excluded =  ((tmp2.Voxels_nbr - tmp2.Voxels_used) / tmp2.Voxels_nbr)*100;
                table_data = union(table_data,tmp2);
            end
            clear x_data y_data sd_data
        end
        % 2nd mean by group
        group_listing = cellstr(unique(table_data.Group));
        curve_listing = cellstr(unique(table_data(:,[1 4 5])));
        for i=1:size(curve_listing,1)
            for j = 1:numel(time_point_listing)
                if sum(table_data.Group == curve_listing{i,1} & table_data.Parameter == curve_listing{i,2} & table_data.VOI == curve_listing{i,3} &...
                        table_data.Timepoint == time_point_listing{j}) > 1
                    
                    % collect y_data
                    tmp_data=  table_data(table_data.Group == curve_listing{i,1} & table_data.Parameter == curve_listing{i,2} & table_data.VOI == curve_listing{i,3} &...
                        table_data.Timepoint == time_point_listing{j},:);
                    y_data(j) = nanmean(double(tmp_data.Value));
                    sd_data(j) = nanstd(double(tmp_data.Value));
                    Animal_nbr = numel(double(unique(tmp_data.Patient)));
                    Mean_voxels_nbr = nanmean(double(unique(tmp_data.Voxels_nbr)));
                    Mean_Voxels_used =nanmean(double(unique(tmp_data.Voxels_used)));
                    Mean_pc_voxels_excluded = nanmean(double(unique(tmp_data.Pc_voxels_excluded)));
                else
                    y_data(j) = NaN;
                    sd_data(j) = NaN;
                    Animal_nbr = NaN;
                    Mean_voxels_nbr = NaN;
                    Mean_Voxels_used = NaN;
                    Mean_pc_voxels_excluded = NaN;
                end
                x_data(j)= j+1;
                if i == 1 && j == 1
                    table_data_tmp = dataset;
                    table_data_tmp.Group = nominal(curve_listing{i,1});
                    table_data_tmp.Timepoint = nominal(time_point_listing{j});
                    table_data_tmp.Parameter = nominal(curve_listing{i,2});
                    table_data_tmp.VOI = nominal(curve_listing{i,3});
                    table_data_tmp.Animal_nbr = numel(double(unique(tmp_data.Patient)));
                    table_data_tmp.Value = y_data(j);
                    table_data_tmp.SD =  sd_data(j);
                    table_data_tmp.Mean_voxels_nbr = Mean_voxels_nbr;
                    table_data_tmp.Mean_Voxels_used = Mean_Voxels_used;
                    table_data_tmp.Mean_pc_voxels_excluded = Mean_pc_voxels_excluded;
                else
                    tmp = dataset;
                    tmp.Group = nominal(curve_listing{i,1});
                    tmp.Timepoint = nominal(time_point_listing{j});
                    tmp.Parameter = nominal(curve_listing{i,2});
                    tmp.VOI = nominal(curve_listing{i,3});
                    tmp.Animal_nbr = Animal_nbr;
                    tmp.Value = y_data(j);
                    tmp.SD=  sd_data(j);
                    tmp.Mean_voxels_nbr =Mean_voxels_nbr;
                    tmp.Mean_Voxels_used =Mean_Voxels_used;
                    tmp.Mean_pc_voxels_excluded = Mean_pc_voxels_excluded ;
                    
                    table_data_tmp = union(table_data_tmp, tmp);
                end
                clear tmp_data
            end
            errorbar(handles.speedy_graph1,x_data, y_data,sd_data, 'Color', rand(1,3));
            legend_txt(size(legend_txt,1)+1,1) ={[curve_listing{i,:}]};
            %%%%%%%% activate clic on graph
            set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
            set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
            set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
        end
        table_data = table_data_tmp;
        x_label = {'Timepoint'};
        x_tick = 1:numel(cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)))+2;
        x_lim = [1 numel(time_point_listing)+2];
        XTickLabel = [' ' time_point_listing' ' ']';
        y_label ={'ROI size in µm3'};
        
    case 2 % mean by ROI
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            tmp_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            
            % clip unique
            if get(handles.speedy_clip_to_use_unique, 'Value')
                y_data = size(tmp_data,1)*voxel_size;
                sd_data = 0;
            else
                
            end
            
            x_data=find(strcmp(time_point_listing, cellstr(handles.voxels_database_shortlist_subgroup.Timepoint(i))) == 1)+1;
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            if i == 1
                tmp = unique(tmp_data);
                table_data = unique(tmp(:,1:5));
                table_data.Value = y_data;
                table_data.SD = sd_data;
                table_data.Voxels_nbr =  numel(tmp_data.Value);
                table_data.Voxels_used=  sum(double(tmp_data(tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                    tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
            else
                tmp = unique(tmp_data);
                tmp2 = unique(tmp(:,1:5));
                tmp2.Value = y_data;
                tmp2.SD = sd_data;
                tmp2.Voxels_nbr =  numel(tmp_data.Value);
                tmp2.Voxels_used=  sum(double(tmp_data(tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                    tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                tmp2.Pc_voxels_excluded =  ((tmp2.Voxels_nbr - tmp2.Voxels_used) / tmp2.Voxels_nbr)*100;
                
                table_data = union(table_data,tmp2);
            end
            clear x_data y_data sd_data
        end
        
        % plot ROI-by-ROI
        % loop on curve
        curve_listing = cellstr(unique(table_data(:,[1 2 4 5])));
        for i = 1:size(curve_listing)
            % loop on time point
            for j = 1:numel(time_point_listing)
                if sum(table_data.Group == curve_listing{i,1} & table_data.Patient == curve_listing{i,2} &...
                        table_data.Parameter == curve_listing{i,3} & table_data.VOI == curve_listing{i,4} &...
                        table_data.Timepoint == time_point_listing{j}) == 1
                    
                    % collect y_data
                    y_data(j) = table_data.Value(table_data.Group == curve_listing{i,1} & table_data.Patient == curve_listing{i,2} &...
                        table_data.Parameter == curve_listing{i,3} & table_data.VOI == curve_listing{i,4} &...
                        table_data.Timepoint == time_point_listing{j});
                    sd_data(j) = table_data.SD(table_data.Group == curve_listing{i,1} & table_data.Patient == curve_listing{i,2} &...
                        table_data.Parameter == curve_listing{i,3} & table_data.VOI == curve_listing{i,4} &...
                        table_data.Timepoint == time_point_listing{j});
                else
                    y_data(j) = NaN;
                    sd_data(j) =  NaN;
                end
                x_data(j)= j+1;
            end
            errorbar(handles.speedy_graph1,x_data, y_data,sd_data, 'Color', rand(1,3));
            legend_txt(size(legend_txt,1)+1,1) ={[curve_listing{i,:}]};
            %%%%%%%% activate clic on graph
            set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
            set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
            set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
        end
        
        x_label = {'Timepoint'};
        x_tick = 1:numel(cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)))+2;
        x_lim = [1 numel(time_point_listing)+2];
        XTickLabel = [' ' time_point_listing' ' ']';
        y_label ={'ROI size in µm3'};
        
    case 3 % mean by Slice
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            tmp_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            z_offset = unique(tmp_data.Coord_Z);
            for j = 1:numel(z_offset)
                
                % clip unique
                if get(handles.speedy_clip_to_use_unique, 'Value')
                    y_data = size(tmp_data,1)*voxel_size;
                    sd_data = 0;
                else
                    
                end
                
                x_data=find(strcmp(time_point_listing, cellstr(handles.voxels_database_shortlist_subgroup.Timepoint(i))) == 1)+1;
                if ~isempty(table_data)
                    hold(handles.speedy_graph1, 'on');
                end
                errorbar(handles.speedy_graph1,x_data, y_data,sd_data, 'Color', rand(1,3))
                tmp = cellstr(handles.voxels_database_shortlist_subgroup(i,:));
                legend_txt(size(legend_txt,1)+1,1) ={[tmp{:} 'z' num2str(z_offset(j))]};
                
                %%%%%%%% activate clic on graph
                set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
                set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
                set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
                
                if isempty(table_data)
                    tmp = unique(tmp_data(tmp_data.Coord_Z == z_offset(j), :));
                    table_data = unique(tmp(:,1:5));
                    table_data.Coord_Z =  unique(tmp.Coord_Z);
                    table_data.Value = y_data;
                    table_data.SD = sd_data;
                    table_data.Voxels_nbr =  numel(tmp_data.Value(tmp_data.Coord_Z == z_offset(j)));
                    table_data.Voxels_used =  numel(double(tmp_data(tmp_data.Coord_Z == z_offset(j) & tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                        tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                    table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
                else
                    tmp = unique(tmp_data(tmp_data.Coord_Z == z_offset(j), :));
                    tmp2 = unique(tmp(:,1:5));
                    tmp2.Coord_Z =  unique(tmp.Coord_Z);
                    tmp2.Value = y_data;
                    tmp2.SD = sd_data;
                    tmp2.Voxels_nbr =  numel(tmp_data.Value(tmp_data.Coord_Z == z_offset(j)));
                    tmp2.Voxels_used =  numel(double(tmp_data(tmp_data.Coord_Z == z_offset(j) & tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                        tmp2.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                    tmp2.Pc_voxels_excluded =  ((tmp2.Voxels_nbr - tmp2.Voxels_used) / tmp2.Voxels_nbr)*100;
                    table_data = union(table_data,tmp2);
                end
            end
        end
        x_label = {'Timepoint'};
        x_tick = 1:numel(cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)))+2;
        x_lim = [1 numel(time_point_listing)+2];
        XTickLabel = [' ' time_point_listing' ' ']';
        y_label ={'ROI size in µm3'};
        
    case 4 % mean by Voxel and by Group
        % Non available
        
    case 5 % mean by Voxel and by ROI
        % Non available
end

function   [table_data, legend_txt, x_label, x_tick, x_lim, XTickLabel,y_label, y_tick, y_lim, yTickLabel] = speedy_measure_by_mean_mode(handles, mean_by_mode, time_point_listing, table_data,legend_txt)
global voxels_database
clipped_by_another_parameter_status = get(handles.speedy_clipped_by_another_parameter_popupmenu, 'String');
clipped_by_another_parameter_status = clipped_by_another_parameter_status(get(handles.speedy_clipped_by_another_parameter_popupmenu, 'Value'));

if ~strcmp(clipped_by_another_parameter_status, 'No')
    
    % create the new dataset format
    VOI_list = unique(handles.voxels_database_shortlist_subgroup.VOI);
    Id_list = unique(handles.voxels_database_shortlist_subgroup.Patient);
    Tp_list = unique(handles.voxels_database_shortlist_subgroup.Timepoint);
    All_Pameter_list = [clipped_by_another_parameter_status unique(handles.voxels_database_shortlist_subgroup.Parameter)];
    data_in_table = dataset;
    template_dataset= voxels_database(1,:);
    template_dataset.Value = [];
    template_dataset.Echo = [];
    template_dataset.Expt = [];
    template_dataset.Parameter = [];
    for ii = 1:numel(All_Pameter_list)
        para_name = char(All_Pameter_list(ii));
        para_name = clean_variable_name(para_name);
        eval(['template_dataset.' para_name '=nan;']);
    end
    template_dataset(1,:) = [];
    
    % complet the new dataset format
    warning off
    for i=1:numel(Id_list)
        
        for j=1:numel(Tp_list)
            for x=1:numel(VOI_list)
                tmp_data = template_dataset;
                tmp_database = voxels_database(voxels_database.Patient == Id_list(i) &...
                    voxels_database.Timepoint == Tp_list(j) &...
                    voxels_database.VOI == VOI_list(x),:);
                if ~isempty(tmp_database)
                    for ii = 1:numel(All_Pameter_list)
                        para_name = char(All_Pameter_list(ii));
                        para_name = clean_variable_name(para_name);
                        if sum(tmp_database.Parameter == All_Pameter_list(ii)) ~= 0
                            tmp = tmp_database(tmp_database.Parameter == All_Pameter_list(ii),:);
                            if isempty(tmp_data)
                                tmp_data.Group = tmp.Group;
                                tmp_data.Patient = tmp.Patient;
                                tmp_data.Timepoint = tmp.Timepoint;
                                tmp_data.VOI = tmp.VOI;
                                tmp_data.Coord_X = tmp.Coord_X;
                                tmp_data.Coord_Y = tmp.Coord_Y;
                                tmp_data.Coord_Z =tmp.Coord_Z;
                                eval(['tmp_data.' para_name '=tmp.Value;']);
                            else
                                if size(tmp_data,1) == size(tmp,1)
                                    eval(['tmp_data.' para_name '=tmp.Value;']);
                                end
                            end
                        else
                            if isempty(tmp_data)
                                tmp = tmp_database(tmp_database.Parameter == tmp_database.Parameter(1),:);
                                tmp_data.Group = tmp.Group;
                                tmp_data.Patient = tmp.Patient;
                                tmp_data.Timepoint = tmp.Timepoint;
                                tmp_data.VOI = tmp.VOI;
                                tmp_data.Coord_X = tmp.Coord_X;
                                tmp_data.Coord_Y = tmp.Coord_Y;
                                tmp_data.Coord_Z =tmp.Coord_Z;
                                eval(['tmp_data.' para_name '=nan([size(tmp,1),1]);']);
                            else
                                eval(['tmp_data.' para_name '=nan([size(tmp_data,1),1]);']);  %join
                            end
                        end
                        clear tmp
                    end
                    if isempty(data_in_table)
                        data_in_table = tmp_data;
                    else
                        data_in_table = union(data_in_table, tmp_data);
                    end
                    clear tmp_data
                end
            end
        end
    end
    
    clipped_var_name = char(clipped_by_another_parameter_status);
    clipped_var_name = clean_variable_name(clipped_var_name);
    
    clipped_by_another_parameter_min_value = str2double(get(handles.speedy_clipped_by_another_parameter_min_value, 'String'));
    clipped_by_another_parameter_max_value = str2double(get(handles.speedy_clipped_by_another_parameter_max_value, 'String'));
end



switch mean_by_mode
    case 1 % mean by Group
        % 1st mean by ROI
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            
            if ~strcmp(clipped_by_another_parameter_status, 'No')
                
                tmp_data = data_in_table(data_in_table.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                    data_in_table.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                    data_in_table.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                    data_in_table.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            else
                
                tmp_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                    voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                    voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                    voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                    voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
                
            end
            % clip unique
            if get(handles.speedy_clip_to_use_unique, 'Value')
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    para_name = char(handles.voxels_database_shortlist_subgroup.Parameter(i));
                    para_name = clean_variable_name(para_name); %#ok<NASGU>
                    
                    match_voxel = eval(['tmp_data.' char(para_name) '>= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))']) &...
                        eval(['tmp_data.' char(para_name) '<= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3))']) &...
                        eval(['tmp_data.' char(clipped_var_name) '>= clipped_by_another_parameter_min_value']) &...
                        eval(['tmp_data.' char(clipped_var_name) '<= clipped_by_another_parameter_max_value']);
                    
                    y_data = nanmean(double(eval(['tmp_data.' char(para_name) '(match_voxel)'])));
                    sd_data = nanstd(double(eval(['tmp_data.' char(para_name) '(match_voxel)'])));
                else
                    match_voxel =tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2)) &...
                        tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3));
                    
                    y_data = nanmean(double(tmp_data.Value(match_voxel)));
                    sd_data = nanstd(double(tmp_data.Value(match_voxel)));
                    
                end
            else
                
            end
            x_data=find(strcmp(time_point_listing, cellstr(handles.voxels_database_shortlist_subgroup.Timepoint(i))) == 1)+1;
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            if i == 1
                
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    table_data = unique(tmp_data(:,1:4));
                    table_data.Parameter = handles.voxels_database_shortlist_subgroup.Parameter(i);
                else
                    tmp = unique(tmp_data);
                    table_data = unique(tmp(:,1:5));
                end
                table_data.Value = y_data;
                table_data.SD = sd_data;
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    table_data.Voxels_nbr =size(tmp_data, 1);
                else
                    table_data.Voxels_nbr =  sum(~isnan(tmp_data.Value)); % ; %numel(tmp_data.Value);
                end
                
                table_data.Voxels_used= sum(match_voxel);
                table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
            else
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    tmp2 = unique(tmp_data(:,1:4));
                    tmp2.Parameter = handles.voxels_database_shortlist_subgroup.Parameter(i);
                else
                    tmp = unique(tmp_data);
                    tmp2 = unique(tmp(:,1:5));
                end
                tmp2.Value = y_data;
                tmp2.SD = sd_data;
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    tmp2.Voxels_nbr =size(tmp_data, 1);
                else
                    tmp2.Voxels_nbr =  sum(~isnan(tmp_data.Value)); % size(tmp_data, 1) ; %numel(tmp_data.Value);
                end
                
                tmp2.Voxels_used= sum(match_voxel);
                
                tmp2.Pc_voxels_excluded =  ((tmp2.Voxels_nbr - tmp2.Voxels_used) / tmp2.Voxels_nbr)*100;
                %                 table_data = union(table_data,tmp2);
                table_data(size(table_data,1)+1,:) = tmp2;
                
            end
            clear x_data y_data sd_data
        end
        % 2nd mean by group
        group_listing = cellstr(unique(table_data.Group));
        if ~strcmp(clipped_by_another_parameter_status, 'No')
            curve_listing = cellstr(unique(table_data(:,[1 5 4])));
        else
            curve_listing = cellstr(unique(table_data(:,[1 4 5])));
        end
        for i=1:size(curve_listing,1)
            for j = 1:numel(time_point_listing)
                if sum(table_data.Group == curve_listing{i,1} & table_data.Parameter == curve_listing{i,2} & table_data.VOI == curve_listing{i,3} &...
                        table_data.Timepoint == time_point_listing{j}) > 1
                    
                    % collect y_data
                    tmp_data=  table_data(table_data.Group == curve_listing{i,1} & table_data.Parameter == curve_listing{i,2} & table_data.VOI == curve_listing{i,3} &...
                        table_data.Timepoint == time_point_listing{j},:);
                    y_data(j) = nanmean(double(tmp_data.Value));
                    sd_data(j) = nanstd(double(tmp_data.Value));
                    Animal_nbr = numel(double(tmp_data.Patient));
                    Mean_voxels_nbr = nanmean(double(tmp_data.Voxels_nbr));
                    Mean_Voxels_used =nanmean(double(tmp_data.Voxels_used));
                    Mean_pc_voxels_excluded = nanmean(double(tmp_data.Pc_voxels_excluded));
                else
                    y_data(j) = NaN;
                    sd_data(j) = NaN;
                    Animal_nbr = NaN;
                    Mean_voxels_nbr = NaN;
                    Mean_Voxels_used = NaN;
                    Mean_pc_voxels_excluded = NaN;
                end
                x_data(j)= j+1;
                if i == 1 && j == 1
                    table_data_tmp = dataset;
                    table_data_tmp.Group = nominal(curve_listing{i,1});
                    table_data_tmp.Timepoint = nominal(time_point_listing{j});
                    table_data_tmp.Parameter = nominal(curve_listing{i,2});
                    table_data_tmp.VOI = nominal(curve_listing{i,3});
                    table_data_tmp.Animal_nbr = numel(double(unique(tmp_data.Patient)));
                    table_data_tmp.Value = y_data(j);
                    table_data_tmp.SD =  sd_data(j);
                    table_data_tmp.Mean_voxels_nbr = Mean_voxels_nbr;
                    table_data_tmp.Mean_Voxels_used = Mean_Voxels_used;
                    table_data_tmp.Mean_pc_voxels_excluded = Mean_pc_voxels_excluded;
                else
                    tmp = dataset;
                    tmp.Group = nominal(curve_listing{i,1});
                    tmp.Timepoint = nominal(time_point_listing{j});
                    tmp.Parameter = nominal(curve_listing{i,2});
                    tmp.VOI = nominal(curve_listing{i,3});
                    tmp.Animal_nbr = Animal_nbr;
                    tmp.Value = y_data(j);
                    tmp.SD=  sd_data(j);
                    tmp.Mean_voxels_nbr =Mean_voxels_nbr;
                    tmp.Mean_Voxels_used =Mean_Voxels_used;
                    tmp.Mean_pc_voxels_excluded = Mean_pc_voxels_excluded ;
                    
                    %                     table_data_tmp = union(table_data_tmp, tmp);
                    table_data_tmp(size(table_data_tmp,1)+1,:) = tmp;
                    
                end
                clear tmp_data
            end
            if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                errorbar(handles.speedy_graph1,abs(x_data), abs(y_data),sd_data,...
                    'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:))
            else
                errorbar(handles.speedy_graph1,abs(x_data), abs(y_data),sd_data, 'Color', rand(1,3))
            end
            legend_txt(size(legend_txt,1)+1,1) ={[curve_listing{i,:}]};
            %%%%%%%% activate clic on graph
            set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
            set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
            set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
        end
        hold(handles.speedy_graph1, 'off');
        
        
        table_data = table_data_tmp;
        x_label = {'Timepoint'};
        x_tick = 1:numel(cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)))+2;
        x_lim = [1 numel(time_point_listing)+2];
        XTickLabel = [' ' time_point_listing' ' ']';
        y_label =cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        y_temp = min(table_data.Value)-2*max(table_data.SD):((max(table_data.Value)+2*max(table_data.SD))-(min(table_data.Value)-2*max(table_data.SD)))/5:max(table_data.Value)+2*max(table_data.SD);
        if y_temp(1) < 1
            y_tick = y_temp;
        else
            y_tick = unique(round(y_temp));
        end
        y_lim = [min(y_tick) max(y_tick)];
        yTickLabel = num2str(y_tick');
        
    case 2 % mean by ROI
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            if ~strcmp(clipped_by_another_parameter_status, 'No')
                
                tmp_data = data_in_table(data_in_table.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                    data_in_table.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                    data_in_table.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                    data_in_table.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            else
                
                tmp_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                    voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                    voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                    voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                    voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i),:);
            end
            
            % clip unique
            if get(handles.speedy_clip_to_use_unique, 'Value')
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    para_name = char(handles.voxels_database_shortlist_subgroup.Parameter(i));
                    para_name = clean_variable_name(para_name); %#ok<NASGU>
                    
                    match_voxel = eval(['tmp_data.' char(para_name) '>= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))']) &...
                        eval(['tmp_data.' char(para_name) '<= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3))']) &...
                        eval(['tmp_data.' char(clipped_var_name) '>= clipped_by_another_parameter_min_value']) &...
                        eval(['tmp_data.' char(clipped_var_name) '<= clipped_by_another_parameter_max_value']);
                    
                    y_data = nanmean(double(eval(['tmp_data.' char(para_name) '(match_voxel)'])));
                    sd_data = nanstd(double(eval(['tmp_data.' char(para_name) '(match_voxel)'])));
                else
                    match_voxel =tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2)) &...
                        tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3));
                    
                    y_data = nanmean(double(tmp_data.Value(match_voxel)));
                    sd_data = nanstd(double(tmp_data.Value(match_voxel)));
                    
                end
            else
                
            end
            
            %             % clip unique
            %             if get(handles.speedy_clip_to_use_unique, 'Value')
            %                 match_voxel =tmp_data >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2)) &...
            %                     tmp_data <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3));
            %
            %                 y_data = nanmean(double(tmp_data(match_voxel)));
            %                 sd_data = nanstd(double(tmp_data(match_voxel)));
            %             else
            %
            %             end
            
            x_data=find(strcmp(time_point_listing, cellstr(handles.voxels_database_shortlist_subgroup.Timepoint(i))) == 1)+1;
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            if i == 1
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    table_data = unique(tmp_data(:,1:4));
                    table_data.Parameter = handles.voxels_database_shortlist_subgroup.Parameter(i);
                else
                    table_data =  handles.voxels_database_shortlist_subgroup(i,1:5);
                end
                table_data.Value = y_data;
                table_data.SD = sd_data;
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    table_data.Voxels_nbr =size(tmp_data,1);
                else
                    table_data.Voxels_nbr =  sum(~isnan(tmp_data.Value)); %size(tmp_data,1);
                end
                
                table_data.Voxels_used=  sum(match_voxel);
                table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
            else
                
                
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    tmp2 = unique(tmp_data(:,1:4));
                    tmp2.Parameter = handles.voxels_database_shortlist_subgroup.Parameter(i);
                else
                    tmp2 =  handles.voxels_database_shortlist_subgroup(i,1:5);
                end
                tmp2.Value = y_data;
                tmp2.SD = sd_data;
                if ~strcmp(clipped_by_another_parameter_status, 'No')
                    tmp2.Voxels_nbr =  size(tmp_data,1);
                else
                    tmp2.Voxels_nbr =  sum(~isnan(tmp_data.Value)); %  size(tmp_data,1);
                end
                
                
                tmp2.Voxels_used=  sum(match_voxel);
                
                tmp2.Pc_voxels_excluded =  ((tmp2.Voxels_nbr - tmp2.Voxels_used) / tmp2.Voxels_nbr)*100;
                
                %table_data = union(table_data,tmp2);
                table_data(size(table_data,1)+1,:) = tmp2;
            end
            clear x_data y_data sd_data
        end
        
        % plot ROI-by-ROI
        % loop on curve
        if ~strcmp(clipped_by_another_parameter_status, 'No')
            curve_listing = cellstr(unique(table_data(:,[1 2 5 4])));
        else
            curve_listing = cellstr(unique(table_data(:,[1 2 4 5])));
        end
        for i = 1:size(curve_listing)
            % loop on time point
            for j = 1:numel(time_point_listing)
                if sum(table_data.Group == curve_listing{i,1} & table_data.Patient == curve_listing{i,2} &...
                        table_data.Parameter == curve_listing{i,3} & table_data.VOI == curve_listing{i,4} &...
                        table_data.Timepoint == time_point_listing{j}) == 1
                    
                    % collect y_data
                    y_data(j) = table_data.Value(table_data.Group == curve_listing{i,1} & table_data.Patient == curve_listing{i,2} &...
                        table_data.Parameter == curve_listing{i,3} & table_data.VOI == curve_listing{i,4} &...
                        table_data.Timepoint == time_point_listing{j});
                    sd_data(j) = table_data.SD(table_data.Group == curve_listing{i,1} & table_data.Patient == curve_listing{i,2} &...
                        table_data.Parameter == curve_listing{i,3} & table_data.VOI == curve_listing{i,4} &...
                        table_data.Timepoint == time_point_listing{j});
                else
                    y_data(j) = NaN;
                    sd_data(j) =  NaN;
                end
                x_data(j)= j+1;
            end
            
            if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                errorbar(handles.speedy_graph1,x_data, abs(y_data),abs(sd_data),...
                    'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:))
            else
                errorbar(handles.speedy_graph1,abs(x_data), abs(y_data),sd_data, 'Color', rand(1,3))
            end
            
            legend_txt(size(legend_txt,1)+1,1) ={[curve_listing{i,:}]};
            %%%%%%%% activate clic on graph
            set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
            set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
            set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
        end
        
        x_label = {'Timepoint'};
        x_tick = 1:numel(cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)))+2;
        x_lim = [1 numel(time_point_listing)+2];
        XTickLabel = [' ' time_point_listing' ' ']';
        y_label =cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        y_temp = min(table_data.Value)-2*max(table_data.SD):((max(table_data.Value)+2*max(table_data.SD))-(min(table_data.Value)-2*max(table_data.SD)))/5:max(table_data.Value)+2*max(table_data.SD);
        if y_temp(1) < 1
            y_tick = y_temp;
        else
            y_tick = round(y_temp);
        end
        y_lim = [min(y_tick) max(y_tick)];
        yTickLabel = num2str(y_tick');
    case 3 % mean by Slice
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            tmp_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            z_offset = unique(tmp_data.Coord_Z);
            for j = 1:numel(z_offset)
                
                % clip unique
                if get(handles.speedy_clip_to_use_unique, 'Value')
                    y_data = nanmean(double(tmp_data(tmp_data.Coord_Z == z_offset(j) & tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                        tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                    sd_data = nanstd(double(tmp_data(tmp_data.Coord_Z == z_offset(j) & tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                        tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                else
                    
                end
                
                x_data=find(strcmp(time_point_listing, cellstr(handles.voxels_database_shortlist_subgroup.Timepoint(i))) == 1)+1;
                x_data = repmat(x_data, size(y_data));
                if ~isempty(table_data)
                    hold(handles.speedy_graph1, 'on');
                end
                if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                    errorbar(handles.speedy_graph1,x_data, y_data,sd_data,...
                        'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:))
                else
                    errorbar(handles.speedy_graph1,x_data, y_data,sd_data, 'Color', rand(1,3))
                end
                tmp = cellstr(handles.voxels_database_shortlist_subgroup(i,:));
                legend_txt(size(legend_txt,1)+1,1) ={[tmp{:} 'z' num2str(z_offset(j))]};
                
                %%%%%%%% activate clic on graph
                set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
                set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
                set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
                
                if isempty(table_data)
                    tmp = unique(tmp_data(tmp_data.Coord_Z == z_offset(j), :));
                    table_data = unique(tmp(:,1:5));
                    table_data.Coord_Z =  unique(tmp.Coord_Z);
                    table_data.Value = y_data;
                    table_data.SD = sd_data;
                    table_data.Voxels_nbr =  numel(tmp_data.Value(tmp_data.Coord_Z == z_offset(j)));
                    table_data.Voxels_used =  numel(double(tmp_data(tmp_data.Coord_Z == z_offset(j) & tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                        tmp_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                    table_data.Pc_voxels_excluded =  ((table_data.Voxels_nbr - table_data.Voxels_used) / table_data.Voxels_nbr)*100;
                else
                    tmp = unique(tmp_data(tmp_data.Coord_Z == z_offset(j), :));
                    tmp2 = unique(tmp(:,1:5));
                    tmp2.Coord_Z =  unique(tmp.Coord_Z);
                    tmp2.Value = y_data;
                    tmp2.SD = sd_data;
                    tmp2.Voxels_nbr =  numel(tmp_data.Value(tmp_data.Coord_Z == z_offset(j)));
                    tmp2.Voxels_used =  numel(double(tmp_data(tmp_data.Coord_Z == z_offset(j) & tmp_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2))&...
                        tmp2.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3)), 'Value')));
                    tmp2.Pc_voxels_excluded =  ((tmp2.Voxels_nbr - tmp2.Voxels_used) / tmp2.Voxels_nbr)*100;
                    table_data = union(table_data,tmp2);
                end
            end
        end
        x_label = {'Timepoint'};
        x_tick = 1:numel(cellstr(unique(handles.voxels_database_shortlist_subgroup.Timepoint)))+2;
        x_lim = [1 numel(time_point_listing)+2];
        XTickLabel = [' ' time_point_listing' ' ']';
        y_label =cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        y_temp = min(table_data.Value)-2*max(table_data.SD):((max(table_data.Value)+2*max(table_data.SD))-(min(table_data.Value)-2*max(table_data.SD)))/5:max(table_data.Value)+2*max(table_data.SD);
        if y_temp(1) < 1
            y_tick = y_temp;
        else
            y_tick = round(y_temp);
        end
        y_lim = [min(y_tick) max(y_tick)];
        yTickLabel = num2str(y_tick');
        
    case 4 % histogram-based mean by group
        group_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Group));
        parameter_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.Parameter));
        patient_listing = nominal(cellstr(unique(handles.voxels_database_shortlist_subgroup.Patient)));
        VOI_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.VOI));
        y_max = 0;
        for i = 1:numel(group_listing)
            for j = 1:numel(time_point_listing)
                for ii = 1:numel(parameter_listing)
                    clear tmp_data
                    if get(handles.speedy_patient_name_listbox, 'Value') == 1
                        tmp_data = voxels_database;
                    else
                        for iii = 1:numel(patient_listing)
                            if iii == 1
                                tmp_data =voxels_database(voxels_database.Patient == patient_listing(iii),:);
                            else
                                tmp_data = union(tmp_data, voxels_database(voxels_database.Patient == patient_listing(iii),:));
                            end
                        end
                    end
                    for jj =  1:numel(VOI_listing)
                        if i*j*ii*jj > 1
                            hold(handles.speedy_graph1, 'on');
                        end
                        y_data = tmp_data(tmp_data.Group==group_listing{i} &...
                            tmp_data.Timepoint==time_point_listing{j} &...
                            tmp_data.Parameter == parameter_listing{ii} &...
                            tmp_data.VOI ==  VOI_listing{jj}, :);
                        
                        if isempty(y_data)
                            continue
                        end
                        % clip unique
                        if get(handles.speedy_clip_to_use_unique, 'Value')
                            match_voxel =y_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(ii),2)) &...
                                y_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(ii),3));
                        else
                            
                        end
                        if sum(match_voxel) > 1000
                            [f, xi] = histnorm(y_data.Value(match_voxel),round(sum(match_voxel)/50));
                        else
                            [f, xi] = histnorm(y_data.Value(match_voxel),round(sum(match_voxel)/10));
                        end
                        if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                            plot(handles.speedy_graph1,xi,f,...
                                'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:));
                        else
                            plot(handles.speedy_graph1,xi,f,...
                                'Color',rand(1,3));
                        end
                        y_max = max(y_max, max(f));
                        
                        %             hist(handles.speedy_graph1,y_data.Value(match_voxel),round(sum(match_voxel)/5), 'FaceColor',rand(1,3))
                        
                        %%%%%%%% activate clic on graph
                        set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
                        set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
                        set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
                        
                        % save data in a table
                        if i == 1
                            table_data = y_data(match_voxel,:);
                        else
                            table_data = union(table_data,y_data(match_voxel,:));
                        end
                        tmp = [group_listing{i} '-' time_point_listing{j} '-' parameter_listing{ii} '-' VOI_listing{jj}];
                        legend_txt(size(legend_txt,1)+1,1) ={tmp};
                    end
                end
            end
        end
        x_label = parameter_listing';
        x_tick = get(handles.speedy_graph1, 'XTick');
        x_lim =  get(handles.speedy_graph1, 'XLim');
        XTickLabel =  get(handles.speedy_graph1, 'XTickLabel');
        y_label = {'% of voxel'};
        y_tick = 0:y_max/5:y_max;
        y_lim = [0 y_max];
        yTickLabel = num2str((round(y_tick*1000)/10)');
        
    case 5 % histogram-based mean by ROI
        y_max = 0;
        for i = 1:size(handles.voxels_database_shortlist_subgroup,1)
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            y_data = voxels_database(voxels_database.Group==handles.voxels_database_shortlist_subgroup.Group(i) &...
                voxels_database.Patient == handles.voxels_database_shortlist_subgroup.Patient(i) &...
                voxels_database.Timepoint==handles.voxels_database_shortlist_subgroup.Timepoint(i) &...
                voxels_database.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i) &...
                voxels_database.VOI == handles.voxels_database_shortlist_subgroup.VOI(i), :);
            
            % clip unique
            if get(handles.speedy_clip_to_use_unique, 'Value')
                match_voxel =y_data.Value >= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),2)) &...
                    y_data.Value <= double(handles.clips(handles.clips.Parameter == handles.voxels_database_shortlist_subgroup.Parameter(i),3));
            else
                
            end
            if round(sum(match_voxel)/5) > 1
                [f, xi] = histnorm(y_data.Value(match_voxel),round(sum(match_voxel)/5));
                if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                    plot(handles.speedy_graph1,xi,f,...
                        'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:));
                else
                    plot(handles.speedy_graph1,xi,f,...
                        'Color',rand(1,3));
                end
                tmp = cellstr(handles.voxels_database_shortlist_subgroup(i,:));
                legend_txt(size(legend_txt,1)+1,1) ={[tmp{:}]};
            end
            y_max = max(y_max, max(f));
            %             hist(handles.speedy_graph1,y_data.Value(match_voxel),round(sum(match_voxel)/5), 'FaceColor',rand(1,3))
            
            %%%%%%%% activate clic on graph
            set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
            set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
            set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
            
            % save data in a table
            if i == 1
                table_data = y_data(match_voxel,:);
            else
                table_data = union(table_data,y_data(match_voxel,:));
            end
        end
        x_label = {'Voxels'};
        x_tick = get(handles.speedy_graph1, 'XTick');
        x_lim =  get(handles.speedy_graph1, 'XLim');
        XTickLabel =  get(handles.speedy_graph1, 'XTickLabel');
        y_label = {'% of voxel'};
        y_tick = 0:y_max/5:y_max;
        y_lim = [0 y_max];
        yTickLabel = num2str((round(y_tick*1000)/10)');
    case 6 % 4D data, plot echo curves for each patient
        y_max = 0;
        group_list_agument = short_list_agument_funtion(handles, 'Group');
        patient_listing_agument = short_list_agument_funtion(handles,'Patient');
        time_point_listing_agument = short_list_agument_funtion(handles,'Timepoint');
        parameter_listing_agument = short_list_agument_funtion(handles,'Parameter');
        VOI_listing_agument = short_list_agument_funtion(handles,'VOI');
        
        
        eval(['Voxels_to_plot = voxels_database(', group_list_agument, '&' patient_listing_agument, '&'  time_point_listing_agument,'&', parameter_listing_agument,'&', VOI_listing_agument ',:);']);
        Voxels_to_plot = Voxels_to_plot(:, [1:5 8 6:7 9:11]); %#ok<NODEF>
        
        [dot_to_plot,~,pos] = unique(Voxels_to_plot(:,1:6),{'Group', 'Patient', 'Timepoint', 'Parameter', 'Echo', 'VOI'});
        dot_to_plot.Value =accumarray(pos, Voxels_to_plot.Value, [], @mean);
        curve_nbr = unique(dot_to_plot(:,1:5));
        echo_times = unique(dot_to_plot.Echo);
        %         test = cellstr(curve_nbr);
        for i = 1:size(curve_nbr,1)
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            curve_to_plot = dot_to_plot.Value(dot_to_plot.Group == curve_nbr.Group(i) &...
                dot_to_plot.Patient == curve_nbr.Patient(i) &...
                dot_to_plot.Timepoint == curve_nbr.Timepoint(i) &...
                dot_to_plot.Parameter == curve_nbr.Parameter(i)&...
                dot_to_plot.VOI == curve_nbr.VOI(i));
            if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                
                plot(handles.speedy_graph1,echo_times,curve_to_plot,...
                    'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:));
            else
                plot(handles.speedy_graph1,echo_times',curve_to_plot,...
                    'Color',rand(1,3));
            end
            tmptext = cellstr(curve_nbr(i,1:5));
            legend_txt(size(legend_txt,1)+1,1) = {strcat(tmptext{:})};
            drawnow
            
            
            %             test2(i,:) =  curve_to_plot';
        end
        
        y_max = max(y_max, max(dot_to_plot.Value));
        
        
        %%%%%%%% activate clic on graph
        set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
        set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
        set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
        
        % save data in a table
        table_data = Voxels_to_plot(:, [1:5 7:9 6 10:11]);
        
        x_label = {'Echo Number'};
        x_tick = get(handles.speedy_graph1, 'XTick');
        x_lim =  get(handles.speedy_graph1, 'XLim');
        XTickLabel =  get(handles.speedy_graph1, 'XTickLabel');
        y_label = {'SI'};
        y_tick = 0:y_max/5:y_max;
        y_lim = [0 y_max];
        yTickLabel = num2str((round(y_tick*1000)/10)');
    case 7 % 4D data, plot echo curves mean by group / ROI
        
        y_max = 0;
        group_list_agument = short_list_agument_funtion(handles, 'Group');
        patient_listing_agument = short_list_agument_funtion(handles,'Patient');
        time_point_listing_agument = short_list_agument_funtion(handles,'Timepoint');
        parameter_listing_agument = short_list_agument_funtion(handles,'Parameter');
        VOI_listing_agument = short_list_agument_funtion(handles,'VOI');
        
        
        eval(['Voxels_to_plot = voxels_database(', group_list_agument, '&' patient_listing_agument, '&'  time_point_listing_agument,'&', parameter_listing_agument,'&', VOI_listing_agument ',:);']);
        Voxels_to_plot = Voxels_to_plot(:, [1:5 8 6:7 9:11]); %#ok<NODEF>
        
        [dot_to_plot,~,pos] = unique(Voxels_to_plot(:,1:6),{'Group', 'Patient', 'Timepoint', 'Parameter', 'Echo', 'VOI'});
        dot_to_plot.Value =accumarray(pos, Voxels_to_plot.Value, [], @mean);
        
        echo_times = unique(dot_to_plot.Echo);
        
        [dot_to_plot2,~,pos2] = unique(dot_to_plot(:,1:6),{'Group', 'Timepoint', 'Parameter', 'VOI', 'Echo'});
        dot_to_plot2.Mean =accumarray(pos2, dot_to_plot.Value, [], @mean);
        dot_to_plot2.SD =accumarray(pos2, dot_to_plot.Value, [], @std);
        curve_nbr = unique(dot_to_plot2(:,1:5));
        for i = 1:size(curve_nbr,1)
            if i > 1
                hold(handles.speedy_graph1, 'on');
            end
            if numel(get(handles.speedy_graph1, 'Children'))+1 < 8
                errorbar(handles.speedy_graph1,...
                    echo_times,dot_to_plot2.Mean(dot_to_plot2.Group == curve_nbr.Group(i) &...
                    dot_to_plot2.Patient == curve_nbr.Patient(i) &...
                    dot_to_plot2.Timepoint == curve_nbr.Timepoint(i) &...
                    dot_to_plot2.Parameter == curve_nbr.Parameter(i)&...
                    dot_to_plot2.VOI == curve_nbr.VOI(i)),...
                    dot_to_plot2.SD(dot_to_plot2.Group == curve_nbr.Group(i) &...
                    dot_to_plot2.Patient == curve_nbr.Patient(i) &...
                    dot_to_plot2.Timepoint == curve_nbr.Timepoint(i) &...
                    dot_to_plot2.Parameter == curve_nbr.Parameter(i)&...
                    dot_to_plot2.VOI == curve_nbr.VOI(i)),...
                    'Color',handles.color_RGB(numel(get(handles.speedy_graph1, 'Children'))+1,:));
            else
                errorbar(handles.speedy_graph1,...
                    echo_times,dot_to_plot2.Mean(dot_to_plot2.Group == curve_nbr.Group(i) &...
                    dot_to_plot2.Patient == curve_nbr.Patient(i) &...
                    dot_to_plot2.Timepoint == curve_nbr.Timepoint(i) &...
                    dot_to_plot2.Parameter == curve_nbr.Parameter(i)&...
                    dot_to_plot2.VOI == curve_nbr.VOI(i)),...
                    dot_to_plot2.SD(dot_to_plot2.Group == curve_nbr.Group(i) &...
                    dot_to_plot2.Patient == curve_nbr.Patient(i) &...
                    dot_to_plot2.Timepoint == curve_nbr.Timepoint(i) &...
                    dot_to_plot2.Parameter == curve_nbr.Parameter(i)&...
                    dot_to_plot2.VOI == curve_nbr.VOI(i)),...
                    'Color',rand(1,3));
            end
            tmptext = cellstr(curve_nbr(i,1:5));
            legend_txt(size(legend_txt,1)+1,1) = {strcat(tmptext{:})};
            drawnow
        end
        
        
        y_max = max(y_max, max(dot_to_plot2.Mean));
        %%%%%%%% activate clic on graph
        set(get(handles.speedy_graph1, 'Children'), 'HitTest', 'off');
        set(handles.speedy_graph1,'ButtonDownFcn', @speedy_graph1_clic);
        set(get(handles.speedy_graph1, 'Children'), 'ButtonDownFcn', @speedy_graph1_clic);
        
        % save data in a table
        table_data = Voxels_to_plot(:, [1:4 6:9 5 10:11]);
        
        x_label = {'Echo Number'};
        x_tick = get(handles.speedy_graph1, 'XTick');
        x_lim =  get(handles.speedy_graph1, 'XLim');
        XTickLabel =  get(handles.speedy_graph1, 'XTickLabel');
        y_label = {'SI'};
        y_tick = 0:y_max/5:y_max;
        y_lim = [0 y_max];
        yTickLabel = num2str((round(y_tick*1000)/10)');
        
end


function short_list_agument = short_list_agument_funtion(handles, row_name) %#ok<INUSL>

eval(['row_listing = cellstr(unique(handles.voxels_database_shortlist_subgroup.' row_name '));']);
for g = 1:numel(row_listing) %#ok<NODEF>
    if g == 1
        short_list_agument = char(strcat('(voxels_database.', row_name, '== ', '''', row_listing(g,:),  '''', ')'));
    else
        short_list_agument = [short_list_agument(1:size(short_list_agument,2)-1) ' | ' char(strcat('voxels_database.', row_name, '== ', '''', row_listing(g,:),'''', ')'))];
    end
end

% --- Executes on mouse press over axes background.
function speedy_graph1_clic(hObject, eventdata, handles)
% hObject    handle to patient_graph1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
if ~strcmp(get(gcbf, 'SelectionType'), 'normal') ||...
        get(handles.speedy_histogram_based_by_group, 'Value') == 1 ||...
        get(handles.speedy_mean_by_group, 'Value') == 1
    return
end

% position = get(handles.speedy_graph1, 'CurrentPoint')

tooltip= get(handles.speedy_graph1_tooltip, 'String');


if get(handles.speedy_histogram_based_by_ROI, 'Value') == 1
    tmp = handles.voxels_database_shortlist_subgroup;
    curve_listing = [char(tmp.Group)  char(tmp.Patient) char(tmp.Timepoint)  char(tmp.Parameter)  char(tmp.VOI)];
    match_data = handles.voxels_database_shortlist_subgroup(strcmp(curve_listing, tooltip), :);
else
    match_space = strfind(tooltip{:}, ' ');
    value = str2double(tooltip{:}(1:match_space(1)));
    match_data = handles.data_ploted(abs(handles.data_ploted.Value - value) < 1e-3, :);
end
% open the corresponding image/ROI
% set patient
set(handles.MIA_data.MIA_name_list, 'Value',...
    find(ismember(get(handles.MIA_data.MIA_name_list, 'String'), char(match_data.Patient)) ==1 ));
MIA('MIA_name_list_Callback',hObject, eventdata,  findobj('Tag', 'MIA_GUI'));
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

% set timepoint
set(handles.MIA_data.MIA_time_points_list, 'Value',...
    find(ismember(get(handles.MIA_data.MIA_time_points_list, 'String'), char(match_data.Timepoint)) ==1 ));
MIA('MIA_time_points_list_Callback',  findobj('Tag', 'MIA_GUI'), eventdata, handles.MIA_data);
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

% set and open scan
set(handles.MIA_data.MIA_scan_VOIs_button, 'Value', 0)
MIA('MIA_scan_VOIs_button_Callback',  findobj('Tag', 'MIA_GUI'), eventdata, handles.MIA_data);
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

set(handles.MIA_data.MIA_scans_list, 'Value',...
    find(ismember(get(handles.MIA_data.MIA_scans_list, 'String'), char(match_data.Parameter)) ==1 ));
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

MIA('MIA_load_axes_Callback',  findobj('Tag', 'MIA_GUI'), eventdata, handles.MIA_data)
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));


% set and open VOI
set(handles.MIA_data.MIA_scan_VOIs_button, 'Value', 1)
MIA('MIA_scan_VOIs_button_Callback',  findobj('Tag', 'MIA_GUI'), eventdata, handles.MIA_data);
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

set(handles.MIA_data.MIA_scans_list, 'Value',...
    find(ismember(get(handles.MIA_data.MIA_scans_list, 'String'), char(match_data.VOI)) ==1 ));
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

MIA('MIA_load_axes_Callback',  findobj('Tag', 'MIA_GUI'), eventdata, handles.MIA_data)
handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));

% if plot type selected ==> display the corresponding slice
if get(handles.speedy_mean_by_slice, 'Value') == 1
    set(handles.MIA_data.MIA_slider_slice, 'Value',...
        find(handles.MIA_data.data_selected.z_offset == match_data.Coord_Z ==1 ));
    handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));
    
    MIA('MIA_load_axes_Callback',  findobj('Tag', 'MIA_GUI'), eventdata, handles.MIA_data)
    handles.MIA_data = guidata(findobj('Tag', 'MIA_GUI'));
end



% --- Executes on button press in speedy_mean_by_group.
function speedy_mean_by_group_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_mean_by_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(findobj('Tag', 'speedy_histogram_based_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_histogram_based_by_ROI'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_id'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_group'), 'Value', 0)

% --- Executes on button press in speedy_mean_by_subject.
function speedy_mean_by_subject_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_mean_by_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(findobj('Tag', 'speedy_histogram_based_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_histogram_based_by_ROI'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_id'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_group'), 'Value', 0)


% --- Executes on button press in speedy_mean_by_slice.
function speedy_mean_by_slice_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_mean_by_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(findobj('Tag', 'speedy_histogram_based_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_histogram_based_by_ROI'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_id'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_group'), 'Value', 0)


% --- Executes on button press in speedy_histogram_based_by_group.
function speedy_histogram_based_by_group_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_histogram_based_by_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedy_histogram_based_by_group

set(findobj('Tag', 'speedy_mean_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_subject'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_slice'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_id'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_group'), 'Value', 0)

% --- Executes on button press in speedy_histogram_based_by_ROI.
function speedy_histogram_based_by_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_histogram_based_by_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedy_histogram_based_by_ROI
set(findobj('Tag', 'speedy_mean_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_subject'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_slice'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_id'), 'Value', 0)
set(findobj('Tag', 'speedy_4Ddata_by_echo_group'), 'Value', 0)


% --- Executes on button press in speedy_4Ddata_by_echo_id.
function speedy_4Ddata_by_echo_id_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_4Ddata_by_echo_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedy_4Ddata_by_echo_id
set(findobj('Tag', 'speedy_histogram_based_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_histogram_based_by_ROI'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_subject'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_slice'), 'Value', 0)


% --- Executes on button press in speedy_4Ddata_by_echo_group.
function speedy_4Ddata_by_echo_group_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_4Ddata_by_echo_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedy_4Ddata_by_echo_group
set(findobj('Tag', 'speedy_histogram_based_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_histogram_based_by_ROI'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_group'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_subject'), 'Value', 0)
set(findobj('Tag', 'speedy_mean_by_slice'), 'Value', 0)

% --- Executes on button press in speedy_data_ploted_button.
function speedy_data_ploted_button_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_data_ploted_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


curve_obj = get(handles.speedy_graph1, 'Children');
for i=1:numel(curve_obj)
    data_x(i,:) = get(curve_obj(i), 'Xdata');
    data_y(i,:) = get(curve_obj(i), 'Ydata');
    data_name(i) = {get(curve_obj(i), 'DisplayName');};
end
data
% if ~isunix
% [filename,pathname] = uiputfile('.xls','Please enter le .xls filename');
% export(handles.data_ploted,'XLSfile',fullfile(pathname,filename))
% else
%     [filename,pathname] = uiputfile('.txt','Please enter .txt filename');
%     export(handles.data_ploted,'file',fullfile(pathname,filename),'delimiter','\t')
% end

% % --- Executes on button press in speedy_export_data_ploted_button.
% function speedy_export_data_ploted_button_Callback(hObject, eventdata, handles)
% % hObject    handle to speedy_export_data_ploted_button (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%



% --- Executes on button press in speedy_set_clip_button.
function speedy_set_clip_button_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_set_clip_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedy_set_clip_button

if ~isfield(handles, 'clips')
    return
end
if isempty(findobj('Tag', 'speedy_clip_table'))
    set(handles.speedy_set_clip_button, 'String', 'Update Clip');
    tmp(1:size(handles.clips(:,1),1),1:2) = 1;
    table_data(:,1) = cellstr(handles.clips.Parameter);
    table_data(:,2:3) = [num2cell(double(handles.clips.Clip_Min)) num2cell(double(handles.clips.Clip_Max))];
    uitable('Tag', 'speedy_clip_table', 'Units', 'normalized','Position', [0.60 0.015 0.44 0.60],...
        'ColumnWidth',{100},...
        'ColumnName',{'Parameters','Clip - Min','Clip - Max'},...
        'ColumnEditable', logical([0 1 1]),...
        'Data',table_data);
else
    set(handles.speedy_set_clip_button, 'String', 'Set Clip');
    new_data = get(findobj('Tag', 'speedy_clip_table'), 'Data');
    handles.clips.Clip_Min = cell2num(new_data(:,2));
    handles.clips.Clip_Max = cell2num(new_data(:,3));
    %clear table1
    delete(findobj('Tag', 'speedy_clip_table'));
    
    guidata(hObject, handles);
    
end


% --- Executes on button press in speedy_del_data.
function speedy_del_data_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_del_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'data_ploted')
    handles = rmfield(handles, 'data_ploted');
end
if ~isempty(get(handles.speedy_graph1, 'Children'))
    delete(get(handles.speedy_graph1, 'Children'))
    legend(handles.speedy_graph1,'off');
    hold(handles.speedy_graph1, 'off');
end
set(handles.speedy_table1, 'Data', {'', '', '', '', ''});

guidata(hObject, handles)


% if ~isempty( handles.MIA_data.database(1).databaseinfo.voxels_database_filename)
%     user_response = questdlg('Do you want to delete the voxel database file?', 'Warning', 'Yes', 'No', 'Cancel', 'Cancel');
%     if strcmp(user_response, 'Cancel') || strcmp(user_response, 'No')
%         return
%     end
%     if strcmp(user_response, 'Yes')
%         user_response = questdlg('Do you REALLY wantto delete the voxel database file?', 'Warning', 'No', 'Yes', 'Cancel', 'Cancel');
%     end
%     if strcmp(user_response, 'Cancel') || strcmp(user_response, 'No')
%         return
%     end
%     delete(fullfile(handles.MIA_data.database(1).databaseinfo.pathname,  handles.MIA_data.database(1).databaseinfo.voxels_database_filename))
%     handles.MIA_data.database(1).databaseinfo.voxels_database_filename = '';
%     guidata(hObject, handles);
%     speedy_getData_Callback(hObject, eventdata, handles)
% else
%     speedy_getData_Callback(hObject, eventdata, handles)
% end


% --- Executes on button press in speedy_loadData.
function speedy_loadData_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global voxels_database
logbook = {};
[voxels_database_file, voxels_database_directory] = uigetfile({'*.mat'}, 'Please select a voxels database file','MultiSelect','off');
set(handles.speedy_info_text, 'String', 'Current Prossess : loading...');
load(fullfile(voxels_database_directory, voxels_database_file));
if ~exist('voxels_database', 'var')
    warndlg(strcat('Somthing wrong with the voxels database file:', {' '},fullfile(voxels_database_directory, voxels_database_file)),'Warning');
    return
end
if  handles.MIA_data.database(1).databaseinfo.voxels_database_need_to_update == 0
    logbook{numel(logbook)+1} =[handles.MIA_data.database(1).databaseinfo.voxels_database_filename ' loaded'];
else % need to update the voxels_database
    logbook{numel(logbook)+1} =[handles.MIA_data.database(1).databaseinfo.voxels_database_filename ' is loaded but need to be updated.... not coded yet!!'];
end
handles.voxels_database_shortlist = unique(voxels_database(:,1:5)); %#ok<NODEF>

% update the handles
parameters_list = cellstr(unique(voxels_database.Parameter));
for i =1:size(parameters_list,1)
    min(i) = handles.MIA_data.clips(strcmp(handles.MIA_data.clips(:,1), parameters_list(i)),2);
    max(i) =  handles.MIA_data.clips(strcmp(handles.MIA_data.clips(:,1), parameters_list(i)),3);
end
handles.clips =  dataset({nominal(parameters_list), 'Parameter'},...
    {cell2num(min)', 'Clip_Min'},...
    {cell2num(max)', 'Clip_Max'});
guidata(hObject, handles);

set(handles.speedy_info_text, 'String', 'Current Prossess : Get/Update data done!');
%update display
IA_update_database_display(hObject, eventdata, handles)

% update clipped by another parameter popupmenu menu
set(handles.speedy_clipped_by_another_parameter_popupmenu, 'String', ['No', cellstr(handles.clips.Parameter)'], 'Value',  1);

% display info
if ~isempty(logbook)
    listdlg('ListString', logbook', 'ListSize',[250 350], 'Name', 'logbook');
end


% --- Executes on mouse motion over figure - except title and menu.
function speedy_GUI_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to speedy_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'data_ploted')
    return
end
set(handles.speedy_graph1_tooltip, 'Visible', 'off');
position_fig1 = get(handles.speedy_graph1, 'Position');

% currPt in percentage of the speedy_GUI size
currPt=get(handles.speedy_GUI,'CurrentPoint');
GUI_size = get(handles.speedy_GUI, 'Position');
currPt(1) = currPt(1)/GUI_size(3);
currPt(2) = currPt(2)/GUI_size(4);

if currPt(1) > position_fig1(1) && currPt(1) < position_fig1(1)+position_fig1(3) && ...
        currPt(2) > position_fig1(2) && currPt(2) < position_fig1(2)+position_fig1(4)
    speedy_graph1_currentpostion(hObject);
end

function speedy_graph1_currentpostion(hObject, eventdata, handles)
handles = guidata(hObject);

if isempty(get(handles.speedy_graph1, 'Children'))
    return
else
    curve_nb = get(handles.speedy_graph1, 'Children');
    x_label = [];
    y_label = [];
    for i = 1:numel(curve_nb)
        x_label = [x_label get(curve_nb(i), 'XData')];
        y_label = [y_label get(curve_nb(i), 'YData')];
    end
    position = get(handles.speedy_graph1, 'CurrentPoint');
    if get(handles.speedy_histogram_based_by_ROI, 'Value') == 1
        [~, y_index] = min(abs(y_label - position(1,2)));
        x_list = [];
        for i = 1:numel(x_label)
            if y_label(i)==y_label(y_index)
                x_list = [x_list x_label(i)];
            end
        end
        [~, x_index] = min(abs(x_list - position(1)));
        Curr_x_value = x_list(x_index);
        %find curve
        for i = 1:numel(curve_nb)
            x_labeltps = get(curve_nb(i), 'XData');
            y_labeltps = get(curve_nb(i), 'YData');
            for j = 1:numel(x_labeltps)
                if y_labeltps(j)==y_label(y_index) & x_labeltps(j) == x_list(x_index) %#ok<AND2>
                    caption_label = get(curve_nb(i), 'DisplayName');
                end
            end
        end
        
    else
        [~, x_index] = min(abs(x_label - position(1)));
        y_list = [];
        for i = 1:numel(y_label)
            if x_label(i)==x_label(x_index)
                y_list = [y_list y_label(i)];
            end
        end
        [~, y_index] = min(abs(y_list - position(1,2)));
        Curr_y_value = y_list(y_index);
        %find curve
        for i = 1:numel(curve_nb)
            x_labeltps = get(curve_nb(i), 'XData');
            y_labeltps = get(curve_nb(i), 'YData');
            for j = 1:numel(y_labeltps)
                if x_labeltps(j)==x_label(x_index) & y_labeltps(j) == y_list(y_index) %#ok<AND2>
                    caption_label = get(curve_nb(i), 'DisplayName');
                end
            end
        end
    end
    
    
    currPt=get(handles.speedy_GUI,'CurrentPoint');
    if ~exist('caption_label', 'var')
        txt = strcat(num2str(Curr_y_value));
        set(handles.speedy_graph1_tooltip, 'String', txt);
    else
        if get(handles.speedy_histogram_based_by_ROI, 'Value') == 1 || ...
                get(handles.speedy_histogram_based_by_group, 'Value') == 1
            txt = {caption_label};
            set(handles.speedy_graph1_tooltip, 'String', txt);
        else
            txt = strcat(num2str(Curr_y_value), {' '}, caption_label);
            set(handles.speedy_graph1_tooltip, 'String', txt);
        end
    end
    if strcmp(txt, 'NaN')
        return
    end
    length_size = length(txt{:});
    set(handles.speedy_graph1_tooltip, 'Position', [currPt(1)+0.5 currPt(2)+0.5 length_size+12 1]);
    set(handles.speedy_graph1_tooltip, 'Visible', 'on');
end
guidata(hObject, handles);


% --- Executes on button press in speedy_export_for_clustering_button.
function speedy_export_for_clustering_button_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_export_for_clustering_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Rat Group Tp ROI Offset X Y ADC CBV ...
global voxels_database

if isempty(voxels_database.Properties.UserData)
    warning('You choose the original resolution for each image. If those resolutions are differents, it will impact your clustering results.')
end

Normalisation_mode = questdlg('How would you normalize your data?', 'Warning', 'Patient-by-Patient', 'All database', 'None', 'None');
Clipping_mode = questdlg('Do you want to clip your data?', 'Warning', 'Yes', 'No', 'No');

% create the new dataset format
VOI_list = unique(handles.voxels_database_shortlist_subgroup.VOI);
Id_list = unique(handles.voxels_database_shortlist_subgroup.Patient);
Tp_list = unique(handles.voxels_database_shortlist_subgroup.Timepoint);
All_Pameter_list = unique(handles.voxels_database_shortlist_subgroup.Parameter);
data_in_table = dataset;
template_dataset= voxels_database(1,:);
template_dataset.Value = [];
template_dataset.Echo = [];
template_dataset.Expt = [];
template_dataset.Parameter = [];
for ii = 1:numel(All_Pameter_list)
    para_name = char(All_Pameter_list(ii));
    para_name = clean_variable_name(para_name);
    eval(['template_dataset.' para_name '=nan;']);
end
template_dataset(1,:) = [];

% complet the new dataset format
warning off
for i=1:numel(Id_list)
    i
    for j=1:numel(Tp_list)
        for x=1:numel(VOI_list)
            tmp_data = template_dataset;
            tmp_database = voxels_database(voxels_database.Patient == Id_list(i) &...
                voxels_database.Timepoint == Tp_list(j) &...
                voxels_database.VOI == VOI_list(x),:);
            if ~isempty(tmp_database)
                for ii = 1:numel(All_Pameter_list)
                    para_name = char(All_Pameter_list(ii));
                    para_name = clean_variable_name(para_name);
                    if sum(tmp_database.Parameter == All_Pameter_list(ii)) ~= 0
                        tmp = tmp_database(tmp_database.Parameter == All_Pameter_list(ii),:);
                        if isempty(tmp_data)
                            tmp_data.Group = tmp.Group;
                            tmp_data.Patient = tmp.Patient;
                            tmp_data.Timepoint = tmp.Timepoint;
                            tmp_data.VOI = tmp.VOI;
                            tmp_data.Coord_X = tmp.Coord_X;
                            tmp_data.Coord_Y = tmp.Coord_Y;
                            tmp_data.Coord_Z =tmp.Coord_Z;
                            eval(['tmp_data.' para_name '=tmp.Value;']);
                        else
                            if size(tmp_data,1) == size(tmp,1)
                                eval(['tmp_data.' para_name '=tmp.Value;']);
                            end
                        end
                    else
                        if isempty(tmp_data)
                            tmp = tmp_database(tmp_database.Parameter == tmp_database.Parameter(1),:);
                            tmp_data.Group = tmp.Group;
                            tmp_data.Patient = tmp.Patient;
                            tmp_data.Timepoint = tmp.Timepoint;
                            tmp_data.VOI = tmp.VOI;
                            tmp_data.Coord_X = tmp.Coord_X;
                            tmp_data.Coord_Y = tmp.Coord_Y;
                            tmp_data.Coord_Z =tmp.Coord_Z;
                            eval(['tmp_data.' para_name '=nan([size(tmp,1),1]);']);
                        else
                            eval(['tmp_data.' para_name '=nan([size(tmp_data,1),1]);']);  %join
                        end
                    end
                    clear tmp
                end
                %% Clip data
                if strcmp(Clipping_mode, 'Yes')
                    for iii=1:numel(All_Pameter_list)
                        para_name = char(All_Pameter_list(iii));
                        para_name = clean_variable_name(para_name);
                        % tmp_data.ADC(tmp_data.ADC< handles.clips.Clip_Min(handles.clips.Parameter == char(All_Pameter_list(iii))) = NaN;
                        eval(['tmp_data.' para_name '(tmp_data.' para_name '< handles.clips.Clip_Min(handles.clips.Parameter == char(All_Pameter_list(iii)))) = NaN;']);
                        eval(['tmp_data.' para_name '(tmp_data.' para_name '> handles.clips.Clip_Max(handles.clips.Parameter == char(All_Pameter_list(iii)))) = NaN;']);
                    end
                end
                if strcmp(Normalisation_mode, 'Patient-by-Patient')
                    for iii=1:numel(All_Pameter_list)
                        para_name = char(All_Pameter_list(iii));
                        para_name = clean_variable_name(para_name);
                        % Normalize data
                        % data_normalize = data - mean /SD
                        eval(['tmp_mean = nanmean(tmp_data.' para_name ');']);
                        eval(['tmp_std = nanstd(tmp_data.' para_name ');']);
                        eval(['tmp_data.' para_name '=(tmp_data.' para_name ' - tmp_mean )/ tmp_std;']);
                    end
                end
                if isempty(data_in_table)
                    data_in_table = tmp_data;
                else
                    data_in_table = union(data_in_table, tmp_data);
                end
                clear tmp_data
            end
        end
    end
end

% normalize data (data pooled)
% data_normalize = data - mean /SD
if strcmp(Normalisation_mode, 'All database')
    for i=1:numel(All_Pameter_list)
        para_name = char(All_Pameter_list(i));
        para_name = clean_variable_name(para_name);
        eval(['tmp_mean = nanmean(data_in_table.' para_name ');']);
        eval(['tmp_std = nanstd(data_in_table.' para_name ');']);
        eval(['data_in_table.' para_name '=(data_in_table.' para_name ' - tmp_mean )/ tmp_std;']);
    end
end
warning on %#ok<WNON>

% release memory
clear tmp_database

% export data in csv files
[filename,pathname] = uiputfile('.csv','Enter .csv filename');
export(data_in_table,'File',fullfile(pathname,filename), 'Delimiter', ',');

msgbox('Done', 'logbook') ;



% % prepare data for the clustering
% data_for_R = [];
% for i = 1:numel(All_Pameter_list)
%     para_name = char(All_Pameter_list(i));
%     para_name = clean_variable_name(para_name);
%     eval(['data_for_R = [data_for_R data_in_table.' para_name '];']);
% end
% % save NaN voxels and remove them from data_in_table dataset if NaNs
% nan_voxel_nbr = findn(isnan(data_for_R));
% if ~isempty(nan_voxel_nbr)
%     nan_voxel_nbr =unique(nan_voxel_nbr(:,1));
%     save_nan_voxels=data_in_table(nan_voxel_nbr,:);
%     % no clustering for the NaN data
%     save_nan_voxels.mc = nan(size(save_nan_voxels,1),1);
%
%     data_in_table(nan_voxel_nbr,:)=[];
%     % remove NaN voxels in data_for_R
%     data_for_R(nan_voxel_nbr,:)=[];
% end
%
%
% saveR('tmp.R','data_for_R')
%
% % create and save a new R_script for each clustering
% % create_mclust_script_for_R(handles)
%
%
%
% % run mclust_script
% old_pwd =pwd;
% cd 'C:/Program Files/R/R-3.0.2/bin';
% system('Rscript.exe script_mclust.r')
% cd(old_pwd)
%
% load('Results_tmp.mat')
% delete('tmp.R', 'Results_tmp.mat')
%
% data_in_table.mc = results;
%
% % merge both dataset (clutered one and NaNs one) if NaNs
% if ~isempty(nan_voxel_nbr)
%     data_in_table=cat(1,data_in_table, save_nan_voxels);
% end
% % patient_list = unique(data_in_table.Patient);
% patient_list =unique(cellstr(data_in_table.Patient));
% for p = 1:size(patient_list,1)
%     z=5;
%     tmp_patient = data_in_table(data_in_table.Patient == patient_list{p},:);
%     z_offset = unique(tmp_patient.Coord_Z);
%     if size(z_offset,1) < 3
%         z=size(z_offset,1);
%     end
%     tmp = tmp_patient(tmp_patient.Coord_Z == z_offset(z), :) ;
%
%     test = zeros([128 128]);
%     for i = 1:size(tmp.Coord_Y,1)
%         test(tmp.Coord_X(i),tmp.Coord_Y(i)) = tmp.mc(i);
%     end
%     figure('Name', patient_list{p});
%     imagesc(test)
% end




function speedy_clipped_by_another_parameter_min_value_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_clipped_by_another_parameter_min_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speedy_clipped_by_another_parameter_min_value as text
%        str2double(get(hObject,'String')) returns contents of speedy_clipped_by_another_parameter_min_value as a double


% --- Executes during object creation, after setting all properties.
function speedy_clipped_by_another_parameter_min_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_clipped_by_another_parameter_min_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_clipped_by_another_parameter_popupmenu.
function speedy_clipped_by_another_parameter_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_clipped_by_another_parameter_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_clipped_by_another_parameter_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_clipped_by_another_parameter_popupmenu

popupmenu_list = get(handles.speedy_clipped_by_another_parameter_popupmenu, 'String');
set(handles.speedy_clipped_by_another_parameter_min_value, 'String',...
    handles.clips.Clip_Min(handles.clips.Parameter == popupmenu_list(get(handles.speedy_clipped_by_another_parameter_popupmenu, 'Value'))))

set(handles.speedy_clipped_by_another_parameter_max_value, 'String',...
    handles.clips.Clip_Max(handles.clips.Parameter == popupmenu_list(get(handles.speedy_clipped_by_another_parameter_popupmenu, 'Value'))))

% --- Executes during object creation, after setting all properties.
function speedy_clipped_by_another_parameter_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_clipped_by_another_parameter_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function speedy_clipped_by_another_parameter_max_value_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_clipped_by_another_parameter_max_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speedy_clipped_by_another_parameter_max_value as text
%        str2double(get(hObject,'String')) returns contents of speedy_clipped_by_another_parameter_max_value as a double


% --- Executes during object creation, after setting all properties.
function speedy_clipped_by_another_parameter_max_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_clipped_by_another_parameter_max_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in speedy_PRM_ref_popupmenu.
function speedy_PRM_ref_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to speedy_PRM_ref_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns speedy_PRM_ref_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from speedy_PRM_ref_popupmenu


% --- Executes during object creation, after setting all properties.
function speedy_PRM_ref_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedy_PRM_ref_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in speedy_job_status.
function speedy_job_status_Callback(~, ~, ~)
system('/data/home/sandra/Documents/Stage_Sofiane/Code_R_cluster/job_status.sh');


% --- Executes on button press in speedy_delete_job.
function speedy_delete_job_Callback(~, ~, ~)
system('/data/home/sandra/Documents/Stage_Sofiane/Code_R_cluster/delete_job.sh');

% --- Executes on button press in speedy_copy_to_clipboard_button.
function copy_to_clipboard_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to speedy_copy_to_clipboard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedy_copy_to_clipboard_button
mat2clipboard(get(handles.speedy_table1, 'Data'));
msgbox('Data copied!', 'logbook') ;

