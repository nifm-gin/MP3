function varargout = compute_parametric_maps_menu(varargin)
% COMPUTE_PARAMETRIC_MAPS_MENU MATLAB code for compute_parametric_maps_menu.fig
%      COMPUTE_PARAMETRIC_MAPS_MENU, by itself, creates a new COMPUTE_PARAMETRIC_MAPS_MENU or raises the existing
%      singleton*.
%
%      H = COMPUTE_PARAMETRIC_MAPS_MENU returns the handle to a new COMPUTE_PARAMETRIC_MAPS_MENU or the handle to
%      the existing singleton*.
%
%      COMPUTE_PARAMETRIC_MAPS_MENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPUTE_PARAMETRIC_MAPS_MENU.M with the given input arguments.
%
%      COMPUTE_PARAMETRIC_MAPS_MENU('Property','Value',...) creates a new COMPUTE_PARAMETRIC_MAPS_MENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before compute_parametric_maps_menu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to compute_parametric_maps_menu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help compute_parametric_maps_menu

% Last Modified by GUIDE v2.5 06-Aug-2014 14:10:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @compute_parametric_maps_menu_OpeningFcn, ...
    'gui_OutputFcn',  @compute_parametric_maps_menu_OutputFcn, ...
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


% --- Executes just before compute_parametric_maps_menu is made visible.
function compute_parametric_maps_menu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to compute_parametric_maps_menu (see VARARGIN)

% Choose default command line output for compute_parametric_maps_menu
handles.output = hObject;

% Variables initiation
handles.MIA_data = varargin{1};
handles.CPM_parameters_list = varargin{2};
handles.CPM_list_of_maps = {'DSC', 'T2*pre','T2*post','deltaR2', 'deltaR2*', 'BVf', 'VSI', 'Density', 'T2map', 'T1map',...
    'T2starcorr3D', 'R2prim', 'SO2map', 'ASL_InvEff', 'CBF-ASL', 'ASL-dyn', 'CMRO2', 'DCE-phenomeno', 'T1map_from_MultiAngles','DCE-permeability', 'DCE-concentration', 'BVf/stO2/R'};
set(handles.CPM_maps_wanted_listbox, 'String', handles.CPM_list_of_maps');
set(handles.CPM_maps_name_listbox, 'String', [{'NOT COMPUTED YET'} handles.CPM_parameters_list]');
set(handles.CPM_additional_info_table, 'Data', '', 'ColumnName', '', 'columnFormat', '',  'columnEditable', false)
% patient selected name
handles.patient_selected.name_nbr = get(handles.MIA_data.MIA_name_list, 'Value');

data_selected = finddata_selected(handles.MIA_data);
set(handles.CPM_patient_selected, 'String', ['Patient selected:  ' char(handles.MIA_data.database.Patient(data_selected)) ' ' char(handles.MIA_data.database.Tp(data_selected))]);
set(handles.CPM_additional_info_table, 'Data', {});
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes compute_parametric_maps_menu wait for user response (see UIRESUME)
%  uiwait(handles.CPM_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = compute_parametric_maps_menu_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in CPM_maps_wanted_listbox.
function CPM_maps_wanted_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_maps_wanted_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CPM_maps_wanted_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CPM_maps_wanted_listbox

handles = guidata(hObject);
% reset all
set(handles.CPM_maps_needed_listbox, 'Value', 1);
set(handles.CPM_maps_name_listbox, 'Value', 1);
set(handles.CPM_additional_info_table, 'Data', '', 'ColumnName', '', 'columnFormat', '',  'columnEditable', false)
set(handles.CPM_control_box_listbox, 'String', '');
set(handles.CPM_go_button, 'BackgroundColor', [1 0 0]);
handles.ready_to_go.couple = [];
if isfield(handles, 'ready_to_go')
    handles = rmfield(handles, 'ready_to_go');
end
handles.ready_to_go.logical = 0;

handles.maps_wanted_selected=handles.CPM_list_of_maps(get(handles.CPM_maps_wanted_listbox,'Value'));

handles.maps_needed = [];
new_couple_nbr = 0;
set(handles.CPM_maps_name_listbox, 'Max', 1);
set(handles.CPM_maps_name_listbox, 'Min', 1);
set(handles.CPM_maps_needed_listbox,'Value', 1);
for i = 1:numel(handles.maps_wanted_selected)
    switch handles.maps_wanted_selected{i}
        case 'T2*pre'
            handles.maps_needed = [handles.maps_needed {'MGEpre'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T2*pre'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MGEpre'};
        case 'T2*post'
            handles.maps_needed = [handles.maps_needed {'MGEpost'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T2*post'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MGEpost'};
        case 'T2'
            handles.maps_needed = [handles.maps_needed {'RAW_T2'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T2'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'RAW_T2'};
        case 'deltaR2'
            handles.maps_needed = [handles.maps_needed {'MGESEpre'} {'MGESEpost'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'deltaR2'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MGESEpre', 'MGESEpost'};
        case 'deltaR2*'
            handles.maps_needed = [handles.maps_needed {'T2*pre'} {'T2*post'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'deltaR2*'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'T2*pre', 'T2*post'};
        case 'BVf'
            handles.maps_needed = [handles.maps_needed {'deltaR2*'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'BVf'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'deltaR2*'};
        case 'VSI'
            handles.maps_needed = [handles.maps_needed {'deltaR2*'} {'deltaR2'} {'ADC'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'VSI'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'deltaR2*', 'deltaR2', 'ADC'};
        case 'Density'
            handles.maps_needed = [handles.maps_needed {'deltaR2*'} {'deltaR2'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'Density'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'deltaR2*', 'deltaR2'};
        case 'T2map'
            handles.maps_needed = [handles.maps_needed {'MSE'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T2map'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MSE'};
        case 'T1map'
            handles.maps_needed = [handles.maps_needed {'MTI-T1'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T1map'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MTI-T1'};
        case 'T2starcorr3D'
            handles.maps_needed = [handles.maps_needed {'MGE3D'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T2starcorr3D'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MGE3D'};
        case 'R2prim'
            handles.maps_needed = [handles.maps_needed {'T2starcorr3D'} {'T2map'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'R2prim'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'T2starcorr3D', 'T2map'};
        case 'SO2map'
            handles.maps_needed = [handles.maps_needed {'R2prim'} {'BVf'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'SO2map'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'R2prim', 'BVf'};
        case 'ASL_InvEff'
            handles.maps_needed = [handles.maps_needed {'ASL_Carotids'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'ASL_InvEff'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'ASL_Carotids'};
        case 'CBF-ASL'
            handles.maps_needed = [handles.maps_needed {'ASL'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'CBF-ASL'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'ASL'};
        case 'ASL-dyn'
            handles.maps_needed = [handles.maps_needed {'ASL-dyn'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'ASL-dyn'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'ASL-dyn'};
        case 'CMRO2'
            handles.maps_needed = [handles.maps_needed {'SO2map'} {'CBF'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'CMRO2'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'SO2map', 'CBF'};
        case 'DCE-phenomeno'
            handles.maps_needed = [handles.maps_needed {'DCE-phenomeno'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'DCE-phenomeno'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'DCE-phenomeno'};
        case 'T1map_from_MultiAngles'
            set(handles.CPM_maps_name_listbox, 'Max', 30);
            set(handles.CPM_maps_name_listbox, 'Min', 1);
            handles.maps_needed = [handles.maps_needed {'Multi_Angle_Scans_(min_3)'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'T1map_from_MultiAngles'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'Multi_Angle_Scans_(min_3)'};
        case 'DCE-permeability'
            handles.maps_needed = [handles.maps_needed {'DCE', 'T1map'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'DCE-permeability'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'DCE', 'T1map'};
        case 'DCE-concentration'
            handles.maps_needed = [handles.maps_needed {'DCE-SE', 'T1map'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'DCE-concentration'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'DCE-SE', 'T1map'};
        case 'DSC'
            handles.maps_needed = [handles.maps_needed {'Perf'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'DSC'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'Perf'};
        case 'BVf/stO2/R'
            handles.maps_needed = [handles.maps_needed {'MGEFIDSEpre'} {'MGEFIDSEpost'} {'ADCmap'}];
            new_couple_nbr = new_couple_nbr+1;
            handles.ready_to_go.couple(new_couple_nbr).map_wanted = {'BVf/stO2/R'};
            handles.ready_to_go.couple(new_couple_nbr).map_needed = {'MGEFIDSEpre', 'MGEFIDSEpost', 'ADCmap'};
    end
end

handles.maps_needed= unique(handles.maps_needed);
set(handles.CPM_maps_needed_listbox, 'String', handles.maps_needed');
guidata(hObject, handles);
CPM_maps_needed_listbox_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function CPM_maps_wanted_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPM_maps_wanted_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CPM_maps_needed_listbox.
function CPM_maps_needed_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_maps_needed_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CPM_maps_needed_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CPM_maps_needed_listbox

handles = guidata(hObject);
handles.map_needed_selected=handles.maps_needed(get(handles.CPM_maps_needed_listbox,'Value'));

additional_info = {};


switch handles.map_needed_selected{1}
    % info for T2*map_pre
    case 'MGEpre'
        additional_info_map_wanted = 'T2*pre';
        additional_info_name = {'Threshold', '      Method      ',  'Trash below (ms)', 'Trash after (ms)', 'Mask?'};
        additional_info_data = {'5','Non lineraire', '0', 'Inf', 'None'};
        additional_info_format = {'numeric', {'Non lineraire' 'Lineaire'}, 'numeric', 'numeric', ['None' handles.MIA_data.VOIs(1:end-2)]};
        additional_info_editable = [1 1 1 1 1];
        % info for T2*map_post
    case 'MGEpost'
        additional_info_map_wanted = 'T2*post';
        additional_info_name = {'Threshold', '      Method      ', 'Trash below (ms)', 'Trash after (ms)', 'Mask?'};
        additional_info_data = {'5','Non lineraire', '0', 'Inf', 'None'};
        additional_info_format = {'numeric', {'Non lineraire' 'Lineaire'}, 'numeric', 'numeric', ['None' handles.MIA_data.VOIs(1:end-2)]};
        additional_info_editable = [1 1 1  1 1];
        % info for deltaR2*
    case 'T2*post'
        additional_info_map_wanted = 'deltaR2*';
        additional_info_name = {'Final Res'};
        additional_info_data = {'Original'};
        additional_info_format = {{'Original' '64' '128' '256' '512'}};
        additional_info_editable = 1;
        % info for T2map
    case 'MSE'
        additional_info_map_wanted = 'T2map';
        additional_info_name = {'Threshold', 'Trash below (ms)', 'Trash after (ms)','      Method      ','Remove last echo?', 'Mask?'};
        additional_info_data = {'5', '0', 'Inf', 'Mono expo', 'No', 'None'};
        additional_info_format = {'numeric', 'numeric', 'numeric', {'Mono expo' 'Multi expo' 'EPG'}, {'No' 'Yes'}, ['None' handles.MIA_data.VOIs(1:end-2)]};
        additional_info_editable = [1 1 1 1 1 1];
        % info for T1map
    case 'MTI-T1'
        additional_info_map_wanted = 'T1map';
        additional_info_name = {'First scan','M0 map','IE map', 'Mask?' };
        additional_info_data = {'1','No','No', 'None'};
        additional_info_format = {'numeric', {'No' 'Yes'}, {'No' 'Yes'}, ['None' handles.MIA_data.VOIs(1:end-2)]};
        additional_info_editable = [1 1 1 1];
        % info for T2*2Dfrom3D
    case 'MGE3D'
        additional_info_map_wanted = 'T2starcorr3D';
        additional_info_name = {'Threshold', '      Method      ', 'First slice', 'Nbr of slice to sum', 'Final Res', 'Trash below (ms)', 'Trash after (ms)', 'Trash echo(s) after echo number', 'Mask?'};
        additional_info_data = {'5','Non lineraire', '4', '4', '128', '10', 'Inf', 'end', 'None'};
        additional_info_format = {'numeric', {'Non lineraire' 'Lineaire'}, 'numeric', 'numeric', {'64' '112' '128' '256' '512'}, 'numeric', 'numeric', 'char', ['None' handles.MIA_data.VOIs(1:end-2)]};
        additional_info_editable = [1 1 1 1 1 1 1 1 1];
        % info for R2prim
    case 'T2map'
        additional_info_map_wanted = 'R2prim';
        additional_info_name = {'Final Res'};
        additional_info_data = {'Original'};
        additional_info_format = { {'Original' '64' '128' '256' '512'}};
        additional_info_editable = 1;
        % info for deltaR2
    case 'MGESEpost'
        additional_info_map_wanted = 'deltaR2';
        additional_info_name = {'Nbr of echo used',  'Final Res'};
        additional_info_data = {'1', 'Original'};
        additional_info_format = {{'1' '3'}, {'Original' '64' '128' '256' '512'}};
        additional_info_editable = [1 1];
        % info for BVf
    case 'deltaR2*'
        if sum(strcmp(handles.maps_wanted_selected, 'BVf')) == 1
            additional_info_map_wanted = 'BVf';
            additional_info_name = {'   B0   ', 'gamma (10^8)', 'deltaxi CA (10^-6)'};
            additional_info_data = {'4.7', '2.67502', '0.28'};
            additional_info_format = {{'4.7' '7' '9.4'}, 'numeric', {'0.84' '0.56' '0.28' '0.1867' '0.14' '0.0933'}};
            additional_info_editable = [1 1 1];
        else
            additional_info_map_wanted = '';
            additional_info_name = '';
            additional_info_data = '';
            additional_info_format = '';
            additional_info_editable = 0;
        end
        % info for VSI map
    case 'ADC'
        additional_info_map_wanted = 'VSI';
        additional_info_name = {'   B0   ', 'gamma (10^8)', 'deltaxi CA (10^-6)'};
        additional_info_data = {'4.7', '2.67502', '0.28'};
        additional_info_format = {{'4.7' '7' '9.4'}, 'numeric', {'0.84' '0.56' '0.28' '0.1867' '0.14' '0.0933'}};
        additional_info_editable = [1 1 1];
        % info for Density map
    case 'deltaR2'
        if sum(strcmp(handles.maps_wanted_selected, 'Density')) == 1
            additional_info_map_wanted = 'Density';
            additional_info_name = {'Final Res', 'ADCmap?'};
            additional_info_data = {'Original', 'None'};
            additional_info_format = {{'Original' '64' '128' '256' '512'}, ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list]};
            additional_info_editable = 1;
        else
            additional_info_map_wanted = '';
            additional_info_name = '';
            additional_info_data = '';
            additional_info_format = '';
            additional_info_editable = 0;
        end
        % info for SO2map
    case 'BVf'
        additional_info_map_wanted = 'SO2map';
        additional_info_name = {'   B0   ', 'gamma (10^8)', 'deltaxi HBC (10^-6)', 'HCT', 'Using T1map?', 'HCT_map'};
        additional_info_data = {'4.7', '2.67502', '0.264', '0.375', 'No',  'None'};
        additional_info_format = {{'3' '4.7' '7' '9.4'}, 'numeric', 'numeric', 'numeric', {'No' 'Yes'}, ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list]};
        additional_info_editable = [1 1 1 1 1 1];
        % info for ASL inversion efficiency
    case 'ASL_Carotids'
        additional_info_map_wanted = 'ASL_InvEff';
        additional_info_name = '';
        additional_info_data = '';
        additional_info_format = '';
        additional_info_editable = 0;
        % info for CBF-ASL
    case 'ASL'
        if sum(strcmp(handles.maps_wanted_selected, 'CBF-ASL')) == 1
            additional_info_map_wanted = 'CBF-ASL';
            additional_info_name = {'ASL_InvEff', 'InvEff ROI', 'InvEff value', ...
                'Method quantif',  'Mask?', 'Output: Mean or Dyn?', 'Outlier: Exclude?',  'Outlier: Plot?',...
                'M0 map', 'M0 value', 'T1 map', 'T1 value', 'Transit map', 'Transit value'};
            % initialize value according to magnetic field and method of
            % first selected scan
            if max(strcmp('T1map',handles.CPM_parameters_list))
                T1mapState = 'T1map';
            else
                T1mapState = 'None';
            end
            if max(strcmp('ASL_InvEff',handles.CPM_parameters_list))
                ASL_InvEffState = 'ASL_InvEff';
            else
                ASL_InvEffState = 'None';
            end
            
            additional_info_data = {ASL_InvEffState, 'None', 'Default',...
                'Buxton1998', 'No', 'Mean', 'Yes', 'Yes',...
                'None', 'Default', T1mapState, 'Default','None', 'Default',};
            additional_info_format = {['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], ['None' handles.MIA_data.VOIs(1:end-2)], 'numeric',...
                {'Buxton1998' 'ASL-diff' 'ASL-Hadamard-Bolus' 'HadamardFitTransitTime' 'Alsop1996'  'Parkes2002'}, ['None' handles.MIA_data.VOIs(1:end-2)],  {'Mean' 'Dyn'},...
                {'Yes' 'No'}, {'Yes' 'No'},...
                ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], 'numeric', ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], 'numeric',...
                ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], 'numeric'};
            additional_info_editable = [1 1 1 ...
                1 1 1 ...
                1 1 1 1 ...
                1 1 1 1];
        else
            additional_info_map_wanted = '';
            additional_info_name = '';
            additional_info_data = '';
            additional_info_format = '';
            additional_info_editable = 0;
        end
        % info for ASL-dynamique
    case 'ASL-dyn'
        ListVoi = handles.MIA_data.VOIs;
        if strcmp(ListVoi{1},'Other')
            ListVoi = 'char';
            EditVoi = 0;
        else
            ListVoi = ['All pixels' ListVoi(1:end-1) ];
            EditVoi = 1;
        end
        additional_info_map_wanted = 'ASL-dyn';
        additional_info_name = {'Mask','Inv Eff','T1b'};
        additional_info_data = {'None','0.89','2.4'};
        additional_info_format = {ListVoi,'numeric','numeric'};
        additional_info_editable = [EditVoi 1 1];
        % info for CMRO2
    case 'CBF'
        additional_info_map_wanted = 'CMRO2' ;
        additional_info_name = {'Constant (Ca)', 'Final Res'};
        additional_info_data = {'20.43', 'Original'};
        additional_info_format = {'numeric', {'Original' '64' '128' '256' '512'}};
        additional_info_editable = [1 1];
        % info for DCE-permeability
    case 'DCE'
        additional_info_map_wanted = 'DCE-permeability';
        additional_info_name = {'CA Relaxivity  L/(mmol.s)', 'AIF_x', 'AIF_y', 'AIF_z','Realign First DCE scan first?','DCE tardif?', 'Delay between scan (sec)'};
        additional_info_data = {'3.5', 'NaN', 'NaN', 'NaN', 'No', 'No', 'NaN'};
        additional_info_format = {'numeric', 'numeric', 'numeric', 'numeric', {'Yes' 'No'}, ['No' handles.CPM_parameters_list], 'numeric'};
        additional_info_editable = [0 1 1 1 1 1 1];
        
        % info for T1map_from_MutliAngle
    case 'Multi_Angle_Scans_(min_3)'
        additional_info_map_wanted = 'T1map_from_MultiAngles';
        additional_info_name = {'Threshold', 'Fit', 'Fit type', 'Realign Raw scan first?'};
        additional_info_data = {'5' 'non lineaire', 'A*[1-exp(-t/T1)]', 'No'};
        additional_info_format = {'numeric', 'char', 'char', {'Yes' 'No'}};
        additional_info_editable = [1 0 0 1 1];
        
        % info for DCE-phenomeno
    case 'DCE-phenomeno'
        additional_info_map_wanted = 'DCE-phenomeno';
        additional_info_name = {'Last dyn before bolus', 'end of the analysis (in sec)'};
        additional_info_data = {'Auto', '180'};
        additional_info_format = {'numeric', 'numeric'};
        additional_info_editable = [1 1];
        
        % info for DCE-concentration
    case 'DCE-SE'
        additional_info_map_wanted = 'DCE-concentration';
        additional_info_name = {'Beginning of the bolus', 'CA Relaxivity  L/(mmol.s)'};
        additional_info_data = {'4', '3.5'};
        additional_info_format = {'numeric', 'numeric'};
        additional_info_editable = [1 1];
        
        % info for DSC analysis
    case 'Perf'
        additional_info_map_wanted = 'DSC';
        additional_info_name = {'Mask?'};
        additional_info_data = {'None'};
        additional_info_format = {['None' handles.MIA_data.VOIs(1:end-2)]};
        additional_info_editable = 1;
        % info for BVf/stO2/R
    case 'ADCmap'
        additional_info_map_wanted = 'BVf/stO2/R';
        additional_info_name = {' B0 (T  ', 'Mask?', 'Final Res', 'HCT', 'Mask?', 'Dico', 'MaxTE (nbr)'};
        additional_info_data = {'4.7', 'None', 'Original', '0.375', 'No', 'TE60ms', '23'};
        additional_info_format = {{'4.7'}, ['None' handles.MIA_data.VOIs(1:end-2)],{'Original' '64' '128' '256' '512'},'numeric', {'Yes' 'No'},...
            {'TE60ms' 'dico_full_se_60ms', 'dico_full_se_60ms+ADC', 'dico_full_se_60ms_extendeddiff'},...
            'numeric'};
        additional_info_editable = [1 1 1 1 1 1 1];
    otherwise
        additional_info_map_wanted = '';
        additional_info_name = '';
        additional_info_data = '';
        additional_info_format = '';
        additional_info_editable = 0;
end

% set title
set(handles.CPM_additional_info_title, 'String', ['Additional info for: ' additional_info_map_wanted]);
% set ColumnWidth to auto
set(handles.CPM_additional_info_table, 'ColumnWidth' , 'auto');
% set names of the columns
set(handles.CPM_additional_info_table, 'ColumnName', additional_info_name);
% set data (default's parameters)
set(handles.CPM_additional_info_table, 'Data', additional_info_data);
% set columnFormat (option for each parameters)
set(handles.CPM_additional_info_table, 'columnFormat', additional_info_format);
% set each colomn editable
set(handles.CPM_additional_info_table, 'columnEditable', logical(additional_info_editable));


guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function CPM_maps_needed_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPM_maps_needed_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CPM_maps_name_listbox.
function CPM_maps_name_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_maps_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CPM_maps_name_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CPM_maps_name_listbox


% --- Executes during object creation, after setting all properties.
function CPM_maps_name_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPM_maps_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CPM_control_box_listbox.
function CPM_control_box_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_control_box_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CPM_control_box_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CPM_control_box_listbox

% --- Executes during object creation, after setting all properties.
function CPM_control_box_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPM_control_box_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CPM_add_new_couple_button.
function CPM_add_new_couple_button_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_add_new_couple_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

list = [{'NOT COMPUTED YET'} handles.CPM_parameters_list]';
handles.map_name_selected = list(get(handles.CPM_maps_name_listbox, 'Value'));
control_data = get(handles.CPM_control_box_listbox, 'String');

% Return if the used try to select the same couple twice
for i = 1:size(control_data,1)
    match = strfind(control_data{i}, '//');
    tmp = control_data{i};
    map_needed_already_selected = tmp(1:match(1)-1);
    if strcmp(handles.map_needed_selected, map_needed_already_selected)
        warndlg('Already added!','Warning');
        return
    end
end
if isempty(control_data)
    control_data = [];
end

handles.additional_info_selected = get(handles.CPM_additional_info_table, 'Data');
if ~isempty(handles.additional_info_selected)
    new_couple = strcat(handles.map_needed_selected, '//', handles.map_name_selected, '-->',  [handles.additional_info_selected{:}]);
else
    new_couple = strcat(handles.map_needed_selected, '//', handles.map_name_selected);
end

% save the "couple"
couple_nbr = [];
map_needed_nbr = [];
for i=1:numel(handles.ready_to_go.couple)
    if sum(strcmp(handles.map_needed_selected, handles.ready_to_go.couple(i).map_needed)) > 0
        couple_nbr = [couple_nbr i];
        map_needed_nbr = [map_needed_nbr find(strcmp(handles.map_needed_selected, handles.ready_to_go.couple(i).map_needed) == 1)];
    end
end
% if numel(couple_nbr) > 1
%     warndlg('Contact BL','Warning');
% end

for i = 1:numel(couple_nbr)
    handles.ready_to_go.couple(couple_nbr(i)).map_name{map_needed_nbr(i)} = handles.map_name_selected;
    if ~isempty(handles.additional_info_selected)
        if strcmp(['Additional info for: ' handles.ready_to_go.couple(couple_nbr(i)).map_wanted{:}],  get(handles.CPM_additional_info_title, 'String'))
            handles.ready_to_go.couple(couple_nbr(i)).additional_information{1} = handles.additional_info_selected;
        end
    end
end
% update control box data
control_data = [control_data' {new_couple}]';
set(handles.CPM_control_box_listbox, 'String', control_data);

% check if everything is ready to start
% numel(handles.ready_to_go.couple) == numel(handles.ready_to_go.couple.map_name)
for i=1:numel(handles.ready_to_go.couple)
    if isempty(handles.ready_to_go.couple(i).map_name)
        set(handles.CPM_go_button, 'BackgroundColor', [1 0 0]);
        handles.ready_to_go.logical = 0;
        guidata(hObject, handles)
        return
    else
        if strcmp(handles.ready_to_go.couple(i).map_needed, 'Multi_Angle_Scans_(min_3)')
            if numel(handles.ready_to_go.couple(i).map_name{:}) < 2
                set(handles.CPM_go_button, 'BackgroundColor', [1 0 0]);
                handles.ready_to_go.logical = 0;
                guidata(hObject, handles)
                warndlg('Please select more than 2 maps to calculate a mutli angle T1','Warning');
                guidata(hObject, handles)
                return
            else
                set(handles.CPM_go_button, 'BackgroundColor', [0 1 0]);
                handles.ready_to_go.logical = 1;
                guidata(hObject, handles)
            end
        else
            set(handles.CPM_go_button, 'BackgroundColor', [0 1 0]);
            handles.ready_to_go.logical = 1;
            guidata(hObject, handles)
        end
    end
end
% if numel(get(handles.CPM_maps_needed_listbox, 'String')) == numel(control_data)
%     set(handles.CPM_go_button, 'BackgroundColor', [0 1 0]);
%     handles.ready_to_go.logical = 1;
%
% else
%     set(handles.CPM_go_button, 'BackgroundColor', [1 0 0]);
%     handles.ready_to_go.logical = 0;
% end
% guidata(hObject, handles)


% --- Executes on button press in CPM_remove_new_couple_button.
function CPM_remove_new_couple_button_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_remove_new_couple_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

couple_data_selected = get(handles.CPM_control_box_listbox, 'Value');
control_data = get(handles.CPM_control_box_listbox, 'String');
if isempty(control_data)
    return
end
match = strfind(control_data(couple_data_selected), '//');
tmp = control_data{couple_data_selected};
map_needed_name = tmp(1:cell2num(match(1))-1);
for i=1:numel(handles.ready_to_go.couple)
    for j = 1:numel(handles.ready_to_go.couple(i).map_needed)
        if strcmp({map_needed_name}, handles.ready_to_go.couple(i).map_needed(j))
            couple_data_nbr = i;
            map_needed_nbr = j;
        end
    end
end

control_data(couple_data_selected) = [];
set(handles.CPM_control_box_listbox, 'Value', 1);
set(handles.CPM_control_box_listbox, 'String', control_data);

handles.ready_to_go.couple(couple_data_nbr).map_name{map_needed_nbr} = {};
if isfield(handles.ready_to_go.couple(couple_data_nbr), 'additional_information') && ~isempty(handles.ready_to_go.couple(couple_data_nbr).additional_information)
    handles.ready_to_go.couple(couple_data_nbr).additional_information = [];
end

% check if everything is ready to start
if numel(get(handles.CPM_maps_needed_listbox, 'String')) == numel(control_data)
    set(handles.CPM_go_button, 'BackgroundColor', [0 1 0]);
    handles.ready_to_go.logical = 1;
else
    set(handles.CPM_go_button, 'BackgroundColor', [1 0 0]);
    handles.ready_to_go.logical = 0;
end

set(handles.CPM_control_box_listbox, 'Value', 1);
guidata(hObject, handles)

% --- Executes on button press in CPM_go_button.
function CPM_go_button_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.ready_to_go.logical == 0
    warndlg('Please select every maps needed','Warning');
    return
end
options = get(get(handles.CPM_patient_option_uipanel, 'Children'), 'Value');
if  options{1} == 1
    patient = finddata_selected(handles.MIA_data);
    data_available = handles.MIA_data.database(handles.MIA_data.database.Patient == handles.MIA_data.database.Patient(patient) & ...
        handles.MIA_data.database.Tp == handles.MIA_data.database.Tp(patient),:);
    
elseif options{2} == 1
    patient = finddata_selected(handles.MIA_data);
    data_available = handles.MIA_data.database(handles.MIA_data.database.Patient == handles.MIA_data.database.Patient(patient),:);
else
    data_available =  handles.MIA_data.database;
end
% % Add additional map to the map needed listing
% for ii = 1:numel(handles.maps_wanted_selected)
%     switch handles.maps_wanted_selected{ii}
%         case  'SO2map'
%             match = strcmp('SO2map', [handles.ready_to_go.couple.map_wanted]') == 1;
%             add_parameters = handles.ready_to_go.couple(match).additional_information;
%             if strcmp(add_parameters{:}(5), 'Yes')
%                 [scan_name, ok] = listdlg('PromptString','SO2 map OPTION: Please select the T1 map''s name:',...
%                     'Name', 'Selecting...',...
%                     'SelectionMode','single',...
%                     'ListSize', [400 300],...
%                     'ListString',handles.CPM_parameters_list');
%                 if ok == 1
%                     handles.ready_to_go.couple(match).map_name{3} = handles.CPM_parameters_list(scan_name);
%                 end
%             end
%     end
% end


% % First loop to define specifique condition (patient by patient) if option
% % 2 or 3
file_name_extension = get(handles.CPM_file_name_extension, 'String');
% if options{2} == 1 || options{3} == 1
%     for i = 1:numel(patient)
%         if  options{1} == 1
%             time_point = handles.patient_selected.tp_nbr;
%         else
%             time_point = 1:numel(handles.MIA_data.database(patient(i)).day);
%         end
%         for j = 1:numel(time_point)
%             for ii = 1:numel(handles.maps_wanted_selected)
%                 switch handles.maps_wanted_selected{ii}
%                     case 'DCE-phenomeno'
%                         if  sum(strcmp(['DCE-Max'  file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%
%                             match = strcmp('DCE-phenomeno', [handles.ready_to_go.couple.map_wanted]') == 1;
%                             add_parameters = handles.ready_to_go.couple(match).additional_information;
%                             if ~strcmp(add_parameters{:}(1), 'Auto')
%                                 question = {'Beginning of the bolus :','End of the bolus :'};
%                                 answer  = inputdlg(question,...
%                                     ['DCE-phenomeno :' handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date],...
%                                     [1 50; 1 50],...
%                                     {'4', '20'});
%                                 if isempty(answer)
%                                     return
%                                 end
%                                 DCEphenomeno_debut(patient(i),j) = answer(1);
%                                 DCEphenomeno_fin(patient(i),j) = answer(2);
%
%                             end
%                         end
%                     case 'DCE-concentration'
%                         if  sum(strcmp(['DCE-concentration'  file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%                             question = strcat('Beginning of the bolus for :', handles.MIA_data.database(patient(i)).name, handles.MIA_data.database(patient(i)).day(time_point(j)).date);
%                             answer  = inputdlg(question,'DCE-concentration',1,{'4'});
%                             if isempty(answer)
%                                 return
%                             end
%                             DCEconcentration_debut(patient(i),j) = answer;
%
%                         end
%                     case 'SO2map'
%                         if  sum(strcmp(['SO2map'  file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%                             question = strcat('Hct:', handles.MIA_data.database(patient(i)).name, handles.MIA_data.database(patient(i)).day(time_point(j)).date);
%                             answer  = inputdlg(question,'SO2map',1,{'0.375'});
%                             if isempty(answer)
%                                 return
%                             end
%                             Hct_values(patient(i),j) = answer;
% %                             Hct_values(patient(i),j) = {'0.375'};
%                         end
%                     case 'DCE-permeability'
%                         match = strcmp('DCE-permeability', [handles.ready_to_go.couple.map_wanted]') == 1;
%                         map_name_selected = handles.ready_to_go.couple(match).map_name;
%                         DCE_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
%                          if  sum(strcmp(['Ktrans'  file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0 &&...
%                                sum(DCE_nbr) ==1
%                             question = {'AIF coordinates X:','AIF coordinates Y:' ,'AIF coordinates Z:'};
%                             answer  = inputdlg(question,['AIF coordinate for :' handles.MIA_data.database(patient(i)).name, handles.MIA_data.database(patient(i)).day(time_point(j)).date], [1 75], {'1','1','1'});
%                             if isempty(answer)
%                                 return
%                             end
%                             DCEpermeability_AIF_voxel(patient(i),j,:) = answer';
%                          end
%                 end
%             end
%         end
%     end
% end
%update handles.MIA_data
logbook = {};

% for i = 1:numel(unique(data_available.Patient))
%     for j = 1:numel(unique(data_available.Tp))
for ii = 1:numel(handles.maps_wanted_selected)
    
    %             handles.MIA_data = guidata(handles.MIA_data.MIA_GUI);
    %             set(handles.CPM_patient_selected, 'String', ['Current process:  ' handles.MIA_data.database(patient(i)).name ' ' handles.MIA_data.database(patient(i)).day(time_point(j)).date...
    %                 '...' handles.maps_wanted_selected{ii}]);
    %             drawnow;
    switch handles.maps_wanted_selected{ii}
        case 'T2*pre'
            match = strcmp('T2*pre', [handles.ready_to_go.couple.map_wanted]') == 1;
            map_name_selected = handles.ready_to_go.couple(match).map_name{:};
            add_parameters = handles.ready_to_go.couple(match).additional_information;
            raw_data_available = find(data_available.SequenceName == map_name_selected);
            for x=1:numel(raw_data_available)
                patient_analyzed = data_available.Patient(raw_data_available(x));
                tp_analyzed = data_available.Tp(raw_data_available(x));
                if sum(data_available.Patient == patient_analyzed & data_available.Tp == tp_analyzed & data_available.SequenceName == ['T2*pre' file_name_extension]) == 0
                    
                    % Load mask (if selected and/or exist)
                    %                 if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                    %                     Mask_filename = {''};
                    %                 else
                    %                     Mask_nbr = strcmp(add_parameters{:}(6), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    %                     Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                    %                     if isempty(Mask_filename)
                    Mask_filename = {''};
                    %                     end
                    %                 end
                    [new_data] = parametric_map_T2star_nii(data_available.nii(raw_data_available(x)), data_available.json(raw_data_available(x)), Mask_filename{:}, add_parameters);
                    % test matlab coder
                    % [T2map_a] = parametric_map_T2_testmex(MSE_map_filename{:}, Mask_filename{:}, add_parameters{1,1});
                    
                    if ~isempty(new_data)
                        for y = 1: size(new_data,1)
                            json_data = spm_jsonread(strrep(char(new_data(y)), '.nii', '.json'));
                            New_scan = {'Patient','Tp','nii','json','type', 'SequenceName', 'Group';
                                char(patient_analyzed),tp_analyzed,new_data(y),strrep(char(new_data(y)), '.nii', '.json'), 'scan', json_data.SequenceName, 'Undefined'};
                            handles.MIA_data.database = [handles.MIA_data.database;cell2dataset(New_scan)];
                        end
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        logbook{numel(logbook)+1} =[char(patient_analyzed) '-' char(tp_analyzed) '-T2map Done'];
                    else
                        logbook{numel(logbook)+1} =[char(patient_analyzed) '-' char(tp_analyzed) '-T2map Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[char(patient_analyzed) '-' char(tp_analyzed) '-T2map already computed'];
                end
            end
%             if  sum(strcmp(['T2*pre' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%                 match = strcmp('T2*pre', [handles.ready_to_go.couple.map_wanted]') == 1;
%                 map_name_selected = handles.ready_to_go.couple(match).map_name{:};
%                 add_parameters = handles.ready_to_go.couple(match).additional_information;
%                 % find the parameter position for this patient (if
%                 % exist)
%                 parameter_nbr = strcmp(map_name_selected, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
%                 MGE_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
%                 
%                 % Load mask (if selected and/or exist)
%                 if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
%                     Mask_filename = {''};
%                 else
%                     Mask_nbr = strcmp(add_parameters{:}(6), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
%                     Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
%                     if isempty(Mask_filename)
%                         Mask_filename = {''};
%                     end
%                 end
%                 
%                 if ~isempty(MGE_map_filename) && sum(parameter_nbr) == 1
%                     T2star_map = parametric_map_T2star(MGE_map_filename{:}, Mask_filename{:}, add_parameters);
%                     
%                     if ~isempty(T2star_map)
%                         % save the new map andupdate MIA database
%                         uvascim.image = T2star_map;
%                         % save file
%                         scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starpre' file_name_extension '.mat'];
%                         save(fullfile([path scan_name]),'uvascim');
%                         clear uvascim;
%                         
%                         % add to the database and update handles.parameters and handles.clips
%                         handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, T2star_map, patient(i), time_point(j), ['T2*pre' file_name_extension], scan_name);
%                         guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
%                         
%                         clear T2star_map;
%                         logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*pre Done'];
%                     else
%                         logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*pre Skiped'];
%                     end
%                 else
%                     logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*pre Skiped'];
%                 end
%             else
%                 logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*pre already computed'];
%             end
        case 'T2*post'
            if  sum(strcmp(['T2*post' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('T2*post', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name{:};
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                % find the parameter position for this patient (if
                % exist)
                parameter_nbr = strcmp(map_name_selected, handles.MIA_data.database(patient(i) ).day(time_point(j)).parameters) == 1;
                MGE_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                % Load mask (if selected and/or exist)
                if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                    Mask_filename = {''};
                else
                    Mask_nbr = strcmp(add_parameters{:}(6), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                    if isempty(Mask_filename)
                        Mask_filename = {''};
                    end
                end
                if ~isempty(MGE_map_filename) && sum(parameter_nbr) == 1
                    T2star_map = parametric_map_T2star(MGE_map_filename{:}, Mask_filename{:}, add_parameters);
                    if ~isempty(T2star_map)
                        % save the new map andupdate MIA database
                        uvascim.image = T2star_map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starpost' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, T2star_map, patient(i), time_point(j), ['T2*post' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear T2star_map;
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*post Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*post Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*post Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2*post already computed'];
            end
        case 'deltaR2*'
            % check if deltaR2* exist already
            if  sum(strcmp(['deltaR2*' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('deltaR2*', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                % find the pre and post T2* exist find position for this patient (if exist)
                % T2*pre
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    T2star_pre_nbr = strcmp(['T2*pre' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    T2star_pre_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_pre = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T2star_pre_nbr));
                % T2*post
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    T2star_post_nbr = strcmp(['T2*post' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    T2star_post_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_post = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T2star_post_nbr));
                if ~isempty(filename_pre) && ~isempty(filename_post) && sum(T2star_post_nbr) == 1 && sum(T2star_pre_nbr) == 1
                    deltaR2star_map = parametric_map_deltaR2star(filename_pre{:}, filename_post{:}, add_parameters);
                    if ~isempty(deltaR2star_map)
                        % save the new map andupdate MIA database
                        uvascim.image = deltaR2star_map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2star' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, deltaR2star_map, patient(i), time_point(j), ['deltaR2*' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear deltaR2*;
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2* Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2* Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2* Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2* already computed'];
            end
        case 'deltaR2'
            % check if deltaR2* exist already
            if  sum(strcmp(['deltaR2' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('deltaR2', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                % find the pre MGESE or MGESEFID scan
                RAW_data_pre_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                filename_pre = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(RAW_data_pre_nbr));
                % find the post MGESE or MGESEFID scan
                RAW_data_post_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                filename_post = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(RAW_data_post_nbr));
                if ~isempty(filename_pre) && ~isempty(filename_post) && sum(RAW_data_pre_nbr) == 1  && sum(RAW_data_post_nbr)== 1
                    deltaR2_map = parametric_map_deltaR2(filename_pre{:}, filename_post{:}, add_parameters);
                    if ~isempty(deltaR2_map)
                        % save the new map andupdate MIA database
                        uvascim.image = deltaR2_map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                            '-deltaR2' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, deltaR2_map, patient(i), time_point(j), ['deltaR2' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear deltaR2_map;
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2 Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2 Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2 Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-deltaR2 already computed'];
            end
        case 'BVf'
            if  sum(strcmp(['BVf' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('BVf', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                % find deltaR2* position for this patient (if exist)
                if  isempty(handles.ready_to_go.couple(match).map_name) || strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    deltaR2star_nbr = strcmp(['deltaR2*' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    deltaR2star_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_deltaR2star = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(deltaR2star_nbr));
                
                if ~isempty(filename_deltaR2star) && sum(deltaR2star_nbr) == 1
                    BVf_map = parametric_map_BVf(filename_deltaR2star{:}, add_parameters);
                    if ~isempty(BVf_map)
                        % save the new map andupdate MIA database
                        uvascim.image = BVf_map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, BVf_map, patient(i), time_point(j), ['BVf' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear BVf_map;
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf already computed'];
            end
        case 'VSI'
            if  sum(strcmp(['VSI' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('VSI', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                % load deltaR2star
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    deltaR2star_nbr = strcmp(['deltaR2*' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    deltaR2star_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                deltaR2star_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(deltaR2star_nbr));
                
                % load deltaR2_map
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    deltaR2_nbr = strcmp(['deltaR2' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    deltaR2_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                deltaR2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(deltaR2_nbr));
                
                % load ADC map
                ADC_nbr = strcmp(map_name_selected{3}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                ADC_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ADC_nbr));
                
                if ~isempty(deltaR2star_filename) && ~isempty(deltaR2_filename) && ~isempty(ADC_filename) &&...
                        sum(deltaR2star_nbr) == 1  && sum(deltaR2_nbr) == 1 && sum(ADC_nbr) == 1
                    
                    % compute the VSI map
                    VSImap = parametric_map_VSI(deltaR2star_filename{:}, deltaR2_filename{:},ADC_filename{:}, add_parameters);
                    if ~isempty(VSImap)
                        % save the new map andupdate MIA database
                        uvascim.image = VSImap;
                        % save file
                        VSI_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-VSI' file_name_extension '.mat'];
                        save(fullfile([path VSI_scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, VSImap, patient(i), time_point(j), ['VSI' file_name_extension], VSI_scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-VSI Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-VSI Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-VSI Skiped'];
                end
                clear VSI ADC deltaR2 deltaR2star
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-VSI already computed'];
            end
        case 'Density'
            if  sum(strcmp(['Density' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('Density', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                % load deltaR2star
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    deltaR2star_nbr = strcmp(['deltaR2*' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    deltaR2star_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                deltaR2star_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(deltaR2star_nbr));
                
                % load deltaR2_map
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    deltaR2_nbr = strcmp(['deltaR2' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    deltaR2_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                deltaR2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(deltaR2_nbr));
                
                % find ADC map if needed
                if ~strcmp(add_parameters{:}(2), 'None')
                    ADCmap_nbr = strcmp(add_parameters{:}(2), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    if isempty(strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ADCmap_nbr)))
                        add_parameters{:}(2) = {''};
                    else
                        add_parameters{:}(2) = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ADCmap_nbr));
                    end
                end
                
                if ~isempty(deltaR2star_filename) && ~isempty(deltaR2_filename) &&...
                        sum(deltaR2star_nbr) == 1  && sum(deltaR2_nbr) == 1
                    
                    % compute the Density map
                    Density_map = parametric_map_density(deltaR2star_filename{:}, deltaR2_filename{:},add_parameters);
                    if ~isempty(Density_map)
                        % save the new map andupdate MIA database
                        uvascim.image = Density_map;
                        % save file
                        Density_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Density' file_name_extension '.mat'];
                        save(fullfile([path Density_scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, Density_map, patient(i), time_point(j), ['Density' file_name_extension], Density_scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Density Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Density Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Density Skiped'];
                end
                clear Density deltaR2 deltaR2star
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Density already computed'];
            end
            
            
        case 'T2map'
            match = strcmp('T2map', [handles.ready_to_go.couple.map_wanted]') == 1;
            map_name_selected = handles.ready_to_go.couple(match).map_name{:};
            add_parameters = handles.ready_to_go.couple(match).additional_information;
            if strcmp(add_parameters{:}(4), 'Multi expo')
                file_name_mid_part = {'-max', '-mean'};
            else
                file_name_mid_part = {''};
            end
            raw_data_available = find(data_available.SequenceName == map_name_selected);
            for x=1:numel(raw_data_available)
                patient_analyzed = data_available.Patient(raw_data_available(x));
                tp_analyzed = data_available.Tp(raw_data_available(x));
                if sum(data_available.Patient == patient_analyzed & data_available.Tp == tp_analyzed & data_available.SequenceName == ['T2map' file_name_mid_part{1} file_name_extension]) == 0
                    %                 % find the parameter position for this patient (if
                    %                 % exist)
                    %                 parameter_nbr = strcmp(map_name_selected, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    %                 MSE_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                    %
                    % Load mask (if selected and/or exist)
                    %                 if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                    %                     Mask_filename = {''};
                    %                 else
                    %                     Mask_nbr = strcmp(add_parameters{:}(6), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    %                     Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                    %                     if isempty(Mask_filename)
                    Mask_filename = {''};
                    %                     end
                    %                 end
                    [new_data] = parametric_map_T2_nii(data_available.nii(raw_data_available(x)), data_available.json(raw_data_available(x)), Mask_filename{:}, add_parameters);
                    % test matlab coder
                    % [T2map_a] = parametric_map_T2_testmex(MSE_map_filename{:}, Mask_filename{:}, add_parameters{1,1});
                    
                    if ~isempty(new_data)
                        for y = 1: size(new_data,1)
                        json_data = spm_jsonread(strrep(char(new_data(y)), '.nii', '.json'));
                        New_scan = {'Patient','Tp','nii','json','type', 'SequenceName', 'Group';
                            char(patient_analyzed),tp_analyzed,new_data(y),strrep(char(new_data(y)), '.nii', '.json'), 'scan', json_data.SequenceName, 'Undefined'};
                       handles.MIA_data.database = [handles.MIA_data.database;cell2dataset(New_scan)];
                        end
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)

                        logbook{numel(logbook)+1} =[char(patient_analyzed) '-' char(tp_analyzed) '-T2map Done'];
                    else
                        logbook{numel(logbook)+1} =[char(patient_analyzed) '-' char(tp_analyzed) '-T2map Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[char(patient_analyzed) '-' char(tp_analyzed) '-T2map already computed'];
                end
            end
            
            
        case 'T1map'
            if   sum(strcmp(['T1map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0 %%|| sum(strcmp(['M0map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0 || sum(strcmp(['IEmap' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0)
                match = strcmp('T1map', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name{:};
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                parameter_nbr = strcmp(map_name_selected, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                MTI_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                
                % Load mask (if selected and/or exist)
                if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                    Mask_filename = {''};
                else
                    Mask_nbr = strcmp(add_parameters{:}(4), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                    if isempty(Mask_filename)
                        Mask_filename = {''};
                    end
                end
                
                if ~isempty(MTI_map_filename)  && sum(parameter_nbr) == 1
                    [T1map,M0map,IEmap] = parametric_map_T1(MTI_map_filename{:}, Mask_filename{:}, add_parameters);
                    if ~isempty(T1map)
                        % save the new map andupdate MIA database
                        uvascim.image = T1map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, T1map, patient(i), time_point(j), ['T1map' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear T1map
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Skiped'];
                    end
                    if ~isempty(M0map)
                        % save the new map andupdate MIA database
                        uvascim.image = M0map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, M0map, patient(i), time_point(j), ['M0map' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear M0map
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map Skiped'];
                    end
                    if ~isempty(IEmap)
                        % save the new map andupdate MIA database
                        uvascim.image = IEmap;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-IEmap' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, IEmap, patient(i), time_point(j), ['IEmap' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear IEmap
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-IEmap Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-IEmap Skiped'];
                    end
                else
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map already computed'];
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map already computed'];
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-IEmap already computed'];
            end
        case 'T2starcorr3D'
            if  sum(strcmp(['T2starcorr3D' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('T2starcorr3D', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name{:};
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                % find the parameter position for this patient (if
                % exist)
                
                % calculate MGE2Dfrom3D
                parameter_nbr = strcmp(map_name_selected, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                MGE3D_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                
                if ~isempty(MGE3D_map_filename) && sum(parameter_nbr) == 1
                    if sum(strcmp(['MGE2Dfrom3D' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        MGE2Dfrom3D = parametric_map_MGE2Dfrom3D(MGE3D_map_filename{:}, add_parameters);
                        if ~isempty(MGE2Dfrom3D)
                            % save the new map andupdate MIA database
                            uvascim.image = MGE2Dfrom3D;
                            % save file
                            MGE2Dfrom3D_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D' file_name_extension '.mat'];
                            save(fullfile([path MGE2Dfrom3D_scan_name]),'uvascim');
                            clear uvascim;
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, MGE2Dfrom3D, patient(i), time_point(j), ['MGE2Dfrom3D' file_name_extension], MGE2Dfrom3D_scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D Skiped'];
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starcorr3D Skiped'];
                            continue
                        end
                    else
                        logbook{numel(logbook)+1} = [handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D already computed'];
                        MGE2Dfrom3D_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D' file_name_extension '.mat'];
                    end
                    
                    % calculate T2starcorr3D
                    T2star_parameter = add_parameters;
                    T2star_parameter{:}(3) = {'MGE'};
                    T2star_parameter{:}(4) = T2star_parameter{:}(6);
                    T2star_parameter{:}(5) = T2star_parameter{:}(7);
                    T2star_parameter{:}(6:7) = [];
                    % Load mask (if selected and/or exist)
                    if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                        Mask_filename = {''};
                    else
                        Mask_nbr = strcmp(add_parameters{:}(9), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                        Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                        if isempty(Mask_filename)
                            Mask_filename = {''};
                        end
                    end
                    T2starcorr3D = parametric_map_T2star(strcat(path,MGE2Dfrom3D_scan_name), Mask_filename{:}, T2star_parameter);
                    
                    if ~isempty(T2starcorr3D)
                        % save the new map andupdate MIA database
                        uvascim.image = T2starcorr3D;
                        % save file
                        T2starcorr3D_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starcorr3D' file_name_extension '.mat'];
                        save(fullfile([path T2starcorr3D_scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, T2starcorr3D, patient(i), time_point(j), ['T2starcorr3D' file_name_extension], T2starcorr3D_scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        
                        clear MGE2Dfrom3D T2starcorr3D
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D Done'];
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starcorr3D Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D Done'];
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starcorr3D Skiped'];
                    end
                    
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-MGE2Dfrom3D Skiped'];
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starcorr3D Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T2starcorr3D already computed'];
            end
        case 'R2prim'
            if  sum(strcmp(['R2prim' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('R2prim', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                % find T2*map filename
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    T2star_nbr = strcmp(['T2starcorr3D' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    T2star_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                T2star_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T2star_nbr));
                
                % find T2map filename
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    T2map_nbr = strcmp(['T2map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    T2map_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                T2map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T2map_nbr));
                
                if ~isempty(T2star_filename) && ~isempty(T2map_filename) && sum(T2star_nbr) == 1  && sum(T2map_nbr) == 1
                    R2prim_map = parametric_map_R2prim(T2star_filename{:}, T2map_filename{:}, add_parameters);
                    if ~isempty(R2prim_map)
                        % save R2prim map andupdate MIA database
                        uvascim.image = R2prim_map;
                        R2prim_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R2prim_map' file_name_extension '.mat'];
                        save(fullfile([path R2prim_scan_name]),'uvascim');
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), ['R2prim' file_name_extension], R2prim_scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear  uvascim
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R2prim Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R2prim Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R2prim Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R2prim already computed'];
            end
            
        case 'SO2map'
            if  sum(strcmp(['SO2map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('SO2map', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                if options{2} == 1 || options{3} == 1
                    add_parameters{:}(4) =  Hct_values(patient(i),j);
                end
                % find R2prim
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    R2prim_nbr = strcmp(['R2prim' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    R2prim_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                R2prim_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(R2prim_nbr));
                
                % find BVf_map
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    BVf_map_nbr = strcmp(['BVf' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    BVf_map_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                BVf_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(BVf_map_nbr));
                
                % find T1 map if needed
                if strcmp(add_parameters{:}(5), 'Yes')
                    if  strcmp(map_name_selected{3}, 'NOT COMPUTED YET')
                        T1map_nbr = strcmp(['T1map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    else
                        T1map_nbr = strcmp(map_name_selected{3}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    end
                    T1map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T1map_nbr));
                else
                    T1map_filename ={''};
                end
                % find HCT map if needed
                if ~strcmp(add_parameters{:}(6), 'None')
                    HCTmap_nbr = strcmp(add_parameters{:}(6), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    if isempty(strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(HCTmap_nbr)))
                        add_parameters{:}(6) = {''};
                    else
                        add_parameters{:}(6) = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(HCTmap_nbr));
                    end
                end
                
                
                if ~isempty(R2prim_filename) && ~isempty(BVf_map_filename) && sum(R2prim_nbr) == 1  && sum(BVf_map_nbr) == 1
                    SO2map = parametric_map_SO2(R2prim_filename{:}, BVf_map_filename{:},T1map_filename{:},  add_parameters);
                    if ~isempty(SO2map)
                        % save the SO2 map andupdate MIA database
                        uvascim.image = SO2map;
                        % save file
                        SO2map_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-SO2map' file_name_extension '.mat'];
                        save(fullfile([path SO2map_scan_name]),'uvascim');
                        %                             clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, SO2map, patient(i), time_point(j), ['SO2map' file_name_extension], SO2map_scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-SO2map Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-SO2map Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-SO2map Skiped'];
                end
                clear SO2map BVf_map R2prim uvascim
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-SO2map already computed'];
            end
        case 'ASL_InvEff'
            % check if ASL_InvEff exist already
            if  sum(strcmp(['ASL_InvEff' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('ASL_InvEff', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                % find the ASL_Carotids map (i.e: the position for this patient (if exist)
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    ASL_nbr = strcmp(['ASL_Carotids' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    ASL_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_ASL = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ASL_nbr));
                
                if ~isempty(filename_ASL) && sum(ASL_nbr) == 1
                    ASL_InvEff = parametric_map_ASL_InvEff(filename_ASL{:});
                    
                    if ~isempty(ASL_InvEff)
                        % save the new map andupdate MIA database
                        uvascim.image = ASL_InvEff;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ASL_InvEff.reco.texte file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, ASL_InvEff, patient(i), time_point(j), [ASL_InvEff.reco.texte file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ASL_InvEff.reco.texte ' Done'];
                        clear ASL_InvEff ;
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ASL_InvEff.reco.texte ' Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-ASL_InvEff Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-ASL_InvEff already computed'];
            end
        case 'CBF-ASL'
            % check if CBF-ASL exist already
            if  sum(strcmp(['ASL_CBF' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('CBF-ASL', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                % find the T1map, M0map, TTmap and the ASL map (i.e: the position for this patient (if exist)
                % find T1 map if needed
                if ~strcmp(add_parameters{:}(11), 'None')
                    if  strcmp(add_parameters{:}(11), 'NOT COMPUTED YET')
                        T1map_nbr = strcmp(['T1map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    else
                        T1map_nbr = strcmp(add_parameters{:}(11), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    end
                    T1map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T1map_nbr));
                    if isempty(T1map_filename)
                        T1map_filename ={''};
                    end
                else
                    T1map_filename ={''};
                end
                % find M0 map if needed
                if ~strcmp(add_parameters{:}(9), 'None')
                    if  strcmp(add_parameters{:}(9), 'NOT COMPUTED YET')
                        M0map_nbr = strcmp(['M0map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    else
                        M0map_nbr = strcmp(add_parameters{:}(9), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    end
                    M0map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(M0map_nbr));
                    if isempty(M0map_filename)
                        M0map_filename ={''};
                    end
                else
                    M0map_filename ={''};
                end
                % find TT map if needed
                if ~strcmp(add_parameters{:}(13), 'None')
                    if  strcmp(add_parameters{:}(13), 'NOT COMPUTED YET')
                        TTmap_nbr = strcmp(['TTmap' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    else
                        TTmap_nbr = strcmp(add_parameters{:}(13), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    end
                    TTmap_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(TTmap_nbr));
                    if isempty(TTmap_filename)
                        TTmap_filename ={''};
                    end
                else
                    TTmap_filename ={''};
                end
                % ASL
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    ASL_nbr = strcmp(['ASL' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    ASL_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_ASL = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ASL_nbr));
                
                % VOI inversion efficiency
                if ~strcmp(add_parameters{:}(2), 'None')
                    InvEff_ROI_nbr = strcmp(add_parameters{:}(2), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    InvEff_ROI_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(InvEff_ROI_nbr));
                    if isempty(InvEff_ROI_filename) || sum(InvEff_ROI_nbr) == 0
                        InvEff_ROI_filename ={''};
                    end
                else
                    InvEff_ROI_filename ={''};
                end
                
                % VOI Mask
                if  ~strcmp(add_parameters{:}(5), 'No')
                    VoiMask_nbr = strcmp(add_parameters{:}(5), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    VoiMask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(VoiMask_nbr));
                else
                    VoiMask_filename ={''};
                end
                % find inversion efficiency scan if needed
                if ~strcmp(add_parameters{:}(1), 'None')
                    if  strcmp(add_parameters{:}(1), 'NOT COMPUTED YET')
                        ASL_InvEff_nbr = strcmp(['InvEff' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    else
                        ASL_InvEff_nbr = strcmp(add_parameters{:}(1), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    end
                    ASL_InvEff_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ASL_InvEff_nbr));
                    if isempty(ASL_InvEff_filename)
                        ASL_InvEff_filename ={''};
                    end
                else
                    ASL_InvEff_filename ={''};
                end
                
                if ~isempty(filename_ASL) && sum(ASL_nbr) == 1
                    ASL = parametric_map_CBF_ASL(InvEff_ROI_filename{:}, VoiMask_filename{:}, ASL_InvEff_filename{:}, T1map_filename{:}, M0map_filename{:}, TTmap_filename{:}, filename_ASL{:},  add_parameters);
                    if ~isempty(ASL)
                        maps = {'ASL_CBF','ASL_Transit','ASL_T1app'};
                        for jj = 1:numel(ASL)
                            % save the new map andupdate MIA database
                            uvascim.image = ASL{jj};
                            % save file
                            scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' maps{jj} file_name_extension '.mat'];
                            save(fullfile([path scan_name]),'uvascim');
                            
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image , patient(i), time_point(j), [maps{jj} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            clear uvascim ;
                            
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' maps{jj} '- Done'];
                        end
                        clear ASL
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBF-ASL Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBF-ASL Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBF-ASL already computed'];
            end
        case 'ASL-dyn'
            if  ( sum(strcmp(['Transitmap' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0 || sum(strcmp(['M0map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0 || sum(strcmp(['CBFmap' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0)
                match = strcmp('ASL-dyn', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name{:};
                %handles.ready_to_go.couple(match).additional_information =  {{'1'    'Yes'    'Yes' 'cerveau'}};
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                parameter_nbr = strcmp(map_name_selected, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                ASL_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                % VOI mask
                if  max(strcmp(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs, add_parameters{:}(1)))==1
                    Voi_nbr = strcmp(add_parameters{:}(1), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    Voi_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Voi_nbr));
                else
                    Voi_filename ={''};
                end
                
                
                if ~isempty(ASL_map_filename)  && sum(parameter_nbr) == 1
                    [Transitmap,M0map,CBFmap,T1map] = parametric_map_ASL_dyn(ASL_map_filename{:},Voi_filename{:}, add_parameters);
                    if ~isempty(Transitmap)
                        % save the new map andupdate MIA database
                        uvascim.image = Transitmap;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Transitmap' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, Transitmap, patient(i), time_point(j), ['Transitmap' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear Transitmap
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Transitmap Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Transitmap Skiped'];
                    end
                    if ~isempty(M0map)
                        % save the new map andupdate MIA database
                        uvascim.image = M0map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, M0map, patient(i), time_point(j), ['M0map' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear M0map
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map Skiped'];
                    end
                    if ~isempty(CBFmap)
                        % save the new map andupdate MIA database
                        uvascim.image = CBFmap;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBFmap' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, CBFmap, patient(i), time_point(j), ['CBFmap' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear CBFmap
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBFmap Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBFmap Skiped'];
                    end
                    if ~isempty(T1map)
                        % save the new map andupdate MIA database
                        uvascim.image = T1map;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, T1map, patient(i), time_point(j), ['T1map' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear T1map
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Skiped'];
                    end
                    
                else
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Transitmap already computed'];
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-M0map already computed'];
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CBFmap already computed'];
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map already computed'];
            end
        case 'CMRO2'
            % check if CMRO2 exist already
            if  sum(strcmp(['CMRO2' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('CMRO2', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                
                % find the SO2 and the CBF map (i.e: the position for this patient (if exist)
                % SO2map
                if  strcmp(map_name_selected{1}, 'NOT COMPUTED YET')
                    SO2map_nbr = strcmp(['SO2map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    SO2map_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_SO2map= strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(SO2map_nbr));
                % CBF
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    CBF_nbr = strcmp(['ASL_CBF' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    CBF_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                filename_CBF = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(CBF_nbr));
                
                if ~isempty(filename_SO2map) && ~isempty(filename_CBF) && sum(SO2map_nbr) == 1  && sum(CBF_nbr) == 1
                    CMRO2 = parametric_map_CMRO2(filename_SO2map{:}, filename_CBF{:},  add_parameters);
                    if ~isempty(CMRO2)
                        % save the new map andupdate MIA database
                        uvascim.image = CMRO2;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CMRO2' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, CMRO2, patient(i), time_point(j), ['CMRO2' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear CMRO2 ;
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CMRO2 Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CMRO2 Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CMRO2 Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-CMRO2 already computed'];
            end
            
        case 'DCE-phenomeno'
            % check if DCE-phenomeno exist already
            if  sum(strcmp(['DCE-Max' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('DCE-phenomeno', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                DCE_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                filename_DCE= strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DCE_nbr));
                
                if (options{2} == 1 || options{3} == 1) && ~strcmp(add_parameters{:}(1), 'Auto')
                    add_parameters{:}(1) =  DCEphenomeno_debut(patient(i),j);
                    add_parameters{:}(2) =  DCEphenomeno_fin(patient(i),j);
                end
                if ~isempty(filename_DCE) && sum(DCE_nbr) == 1
                    [maxi_struc, ttp_struc, rehaus_struc, AUC_struc, ~] = parametric_map_DCEphenomeno(filename_DCE{:}, add_parameters);
                    if ~isempty(maxi_struc)
                        % save the new maps andupdate MIA database
                        uvascim.image = maxi_struc;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-Max' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, maxi_struc, patient(i), time_point(j), ['DCE-Max' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear maxi_struc ;
                        
                        uvascim.image = ttp_struc;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-ttp' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, ttp_struc, patient(i), time_point(j), ['DCE-ttp' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear ttp_struc ;
                        
                        uvascim.image = rehaus_struc;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-rehaus' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, rehaus_struc, patient(i), time_point(j), ['DCE-rehaus' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear rehaus_struc ;
                        
                        uvascim.image = AUC_struc;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-AUC' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, AUC_struc, patient(i), time_point(j), ['DCE-AUC' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear AUC_struc ;
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-phenomeno Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-phenomeno Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-phenomeno Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-phenomeno already computed'];
            end
        case {'T1map_from_MultiAngles'}
            if  sum(strcmp(['T1map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('T1map_from_MultiAngles' , [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                NbAngles = numel(map_name_selected{:});
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                FA_nbr = false([NbAngles, numel(handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)]);
                FA_filenames = cell([1 NbAngles]);
                for angles=1:NbAngles
                    FA_nbr(angles,:) = strcmp(map_name_selected{:}{angles}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    FA_filenames{angles} = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(FA_nbr(angles,:)));
                end
                if ~isempty(FA_filenames) && sum(FA_nbr(:)) == NbAngles
                    
                    % compute the T1map_from_MutliAngle map
                    T1map = parametric_map_T1fromVariableMultiAngles(FA_filenames, add_parameters);
                    if ~isempty(T1map)
                        
                        % save the report (as .ps)
                        date_str = date;
                        if exist([path filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
                            movefile([path filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
                                [path, handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-T1mapRawData' file_name_extension] '_SPM_realign.ps']);
                        end
                        % save the new map andupdate MIA database
                        uvascim.image = T1map;
                        % save file
                        T10_scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map' file_name_extension '.mat'];
                        save(fullfile([path T10_scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, T1map, patient(i), time_point(j), ['T1map' file_name_extension], T10_scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map Skiped'];
                end
                clear T1map
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-T1map already computed'];
            end
            
        case 'DCE-permeability'
            if  sum(strcmp(['Ktrans' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('DCE-permeability', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                if options{2} == 1 || options{3} == 1
                    add_parameters{:}(2) =  DCEpermeability_AIF_voxel(patient(i),j,1);
                    add_parameters{:}(3) =  DCEpermeability_AIF_voxel(patient(i),j,2);
                    add_parameters{:}(4) = DCEpermeability_AIF_voxel(patient(i),j,3);
                end
                % load T1map
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    T1map_nbr = strcmp(['T1map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    T1map_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                T1map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T1map_nbr));
                
                % load DCE map
                DCE_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                DCE_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DCE_nbr));
                
                if  ~strcmp(add_parameters{:}(6), 'No')
                    DCEtardif_nbr = strcmp(add_parameters{:}(6), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    DCEtardif_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DCEtardif_nbr));
                else
                    DCEtardif_filename = {'None'};
                end
                
                if ~isempty(T1map_filename) && ~isempty(DCE_filename) &&...
                        sum(T1map_nbr) == 1  && sum(DCE_nbr) == 1
                    
                    % check if the AIF is defined if not as the
                    % user
                    AIF_x = str2double(add_parameters{:}(2));
                    AIF_y = str2double(add_parameters{:}(3));
                    AIF_z = str2double(add_parameters{:}(4));
                    if isnan(AIF_x) || isnan(AIF_y) || isnan(AIF_z)
                        question = {'AIF coordinates X:','AIF coordinates Y:' ,'AIF coordinates Z:'};
                        answer  = inputdlg(question,['AIF coordinate for :' handles.MIA_data.database(patient(i)).name, handles.MIA_data.database(patient(i)).day(time_point(j)).date], [1 75], {'1','1','1'});
                        if isempty(answer)
                            return
                        end
                        add_parameters{:}(2) = answer(1);
                        add_parameters{:}(3) = answer(2);
                        add_parameters{:}(4) = answer(3);
                    end
                    % compute the Permeabilities maps
                    [Vp, Kep, Ktrans, Ctfit, Ctmap, T1t, Err] =...
                        parametric_map_DCEpermeability(DCE_filename{:}, T1map_filename{:},DCEtardif_filename{:}, add_parameters); %#ok<ASGLU>
                    if ~isempty(Vp)
                        %                                 maps = {'Vp',  'Kep','Ktrans', 'Ctfit', 'Ctmap', 'T1t', 'Err'};
                        maps = {'Vp',  'Kep','Ktrans'};% 'Ctfit', 'Ctmap', 'T1t', 'Err'};
                        for xx = 1:numel(maps)
                            uvascim.image = eval(maps{xx});
                            uvascim.image.reco.globalmax = max(uvascim.image.reco.data(:));
                            uvascim.image.reco.globalmin = min(uvascim.image.reco.data(:));
                            uvascim.image.clip=[0 0.5 1];
                            uvascim.image.reco=orderfields(uvascim.image.reco);
                            
                            % save the new maps andupdate MIA database
                            % save file
                            scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' maps{xx}  file_name_extension '.mat'];
                            save(fullfile([path scan_name]),'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [maps{xx} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            
                            clear uvascim;
                        end
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-permeability Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-permeability Skiped'];
                    end
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-permeability already computed'];
            end
        case 'DCE-concentration'
            if  sum(strcmp(['DCE-concentration' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('DCE-concentration', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                if options{2} == 1 || options{3} == 1
                    add_parameters{:}(1) =  DCEconcentration_debut(patient(i),j);
                end
                
                % load T1map
                if  strcmp(map_name_selected{2}, 'NOT COMPUTED YET')
                    T1map_nbr = strcmp(['T1map' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                else
                    T1map_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                end
                T1map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(T1map_nbr));
                
                % load DCE map
                DCE_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                DCE_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DCE_nbr));
                
                if ~isempty(T1map_filename) && ~isempty(DCE_filename) &&...
                        sum(T1map_nbr) == 1  && sum(DCE_nbr) == 1
                    
                    % compute DCE-concentration from DCE-SE and
                    % T1map
                    [DCEconcentration] = parametric_map_DCEconcentration(DCE_filename{:}, T1map_filename{:}, add_parameters);
                    if ~isempty(DCEconcentration)
                        uvascim.image = DCEconcentration;
                        uvascim.image.reco.globalmax = max(uvascim.image.reco.data(:));
                        uvascim.image.reco.globalmin = min(uvascim.image.reco.data(:));
                        uvascim.image.clip=[0 uvascim.image.reco.globalmax 1];
                        uvascim.image.reco=orderfields(uvascim.image.reco);
                        
                        % save the new maps andupdate MIA database
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' 'DCE-concentration'  file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), ['DCE-concentration' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        
                        clear uvascim;
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-concentration Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-concentration Skiped'];
                    end
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DCE-concentration already computed'];
            end
            
        case 'DSC'
            % check if DSC analysis scan exist already
            if  sum(strcmp(['CBV' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                match = strcmp('DSC', [handles.ready_to_go.couple.map_wanted]') == 1;
                map_name_selected = handles.ready_to_go.couple(match).map_name;
                add_parameters = handles.ready_to_go.couple(match).additional_information;
                DSC_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                filename_DSC= strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DSC_nbr));
                
                if ~isempty(filename_DSC) && sum(DSC_nbr) == 1
                    
                    fid=fopen(filename_DSC{:} ,'r');
                    if fid>0
                        fclose(fid);
                        load(filename_DSC{:});
                        data = uvascim.image;
                        data.reco.data= double(data.reco.data);
                    end
                    if ~strcmp(add_parameters{1}, 'None')
                        if ~isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                            Mask_nbr = strcmp(add_parameters{:}(1), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                            Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                            if ~isempty(Mask_filename)
                                fid=fopen(Mask_filename{:} ,'r');
                                if fid>0
                                    fclose(fid);
                                    Mask = load(Mask_filename{:});
                                    Mask = Mask.uvascroi;
                                    if ~isempty(Mask)
                                        [~, ~, E, Z,repet] = size(data.reco.data);
                                        
                                        tmp_data = data.reco.data;
                                        for slice = 1:Z
                                            for repetition = 1:repet
                                                for roi_number = 1:numel(Mask)
                                                    if abs(data.reco.fov_offsets(3,1,slice,repetition) - Mask(roi_number).fov_offsets(3)) < 1e-5
                                                        data.reco.data(:,:,:,slice,repetition) = tmp_data(:,:,:,slice,repetition).*repmat(Mask(roi_number).value,[1 1 E 1]);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    else
                        
                    end
                    
                    % algo to detect the AIF
                    meansignal = mean(squeeze(data.reco.data),4);
                    volume_mask = meansignal>max(data.reco.data(:))*0.01;
                    [aif,scores] = extraction_aif_volume(squeeze(data.reco.data),volume_mask);
                    figure; plot(aif); legend(handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(DSC_nbr), 'Location','NorthEast');
                    if sum(cell2num(scores(:,5))) > 10
                        %skip the process because the AIF is not
                        %good enough
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC skip (No correct AIF)'];
                    else
                        
                        drawnow
                        % compute the parameteric maps using the AIF
                        % detected automatically
                        [~,~,~,TMAX,TTP,T0, CBV,CBF,MTT,~,~,~,~] = deconvolution_perfusion_gui(aif,squeeze(data.reco.data(:,:,1,:,:)),data.acq.tr*10^(-3),data.acq.echotime*10^(-3),handles); %#ok<ASGLU>
                        maps = {'CBV','CBF','MTT','TMAX','TTP','T0'};
                        if ~isempty(CBV)
                            for xx = 1:numel(maps)
                                uvascim.image = data;
                                uvascim.image.reco.texte=maps{xx};
                                uvascim.image.reco.no_expts = 1;
                                uvascim.image.reco = rmfield(uvascim.image.reco, 'data');
                                uvascim.image.reco.data(:,:,1,:) = eval(maps{xx});
                                uvascim.image.reco.globalmax = max(uvascim.image.reco.data(:));
                                uvascim.image.reco.globalmin = min(uvascim.image.reco.data(:));
                                uvascim.image.clip=[0 250 1];
                                uvascim.image.reco.unit = '%';
                                ParamConfig=sprintf('##$QuantifMethod=Automatique deconvolution\n##$DSC=%s\n##$AIF=\n%s\n##$AIFScores=\n%s\n%s\n%s\n%s\n##$Warning=\n%d %d %d %d% d\n##$Warning Message=\n%s\n%s\n%s\n%s',...
                                    filename_DSC{:},num2str(aif), num2str([scores{:,1}]), num2str([scores{:,2}]), num2str([scores{:,3}]), num2str([scores{:,4}]),...
                                    [scores{:,5}], [scores{:,6}], [scores{:,7}], [scores{:,8}], [scores{:,9}]);
                                uvascim.image.reco.paramQuantif = ParamConfig;
                                uvascim.image.reco=orderfields(uvascim.image.reco);
                                
                                %Adapt the fov offsets and orientation infos
                                uvascim.image.reco.fov_offsets=uvascim.image.reco.fov_offsets(:,:,:,1);
                                uvascim.image.reco.fov_orientation=uvascim.image.reco.fov_orientation(:,:,:,1);
                                uvascim.image.reco.label=uvascim.image.reco.label(:,:,1);
                                uvascim.image.reco.phaselabel=uvascim.image.reco.phaselabel(:,:,1);
                                uvascim.image.reco.fov_phase_orientation=uvascim.image.reco.fov_phase_orientation(:,:,1);
                                uvascim.image.reco.scaling_factor=uvascim.image.reco.scaling_factor(:,:,1);
                                uvascim.image.reco.scaling_offset=uvascim.image.reco.scaling_offset(:,:,1);
                                
                                % save the new maps andupdate MIA database
                                % save file
                                scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' maps{xx}  file_name_extension];
                                save(fullfile(path, [scan_name '.mat']),'uvascim');
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [maps{xx} file_name_extension], [scan_name '.mat']);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                
                                clear uvascim;
                            end
                            clear CBV CBF MTT TMAX TTP T0;
                            
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC Done'];
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC Skiped'];
                        end
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-DSC already computed'];
            end
        case 'BVf/stO2/R'
            
            match = strcmp('BVf/stO2/R', [handles.ready_to_go.couple.map_wanted]') == 1;
            map_name_selected = handles.ready_to_go.couple(match).map_name;
            add_parameters = handles.ready_to_go.couple(match).additional_information;
            if  sum(strcmp(['BVf-Num' file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                % load pre MGEFIDSE map
                MGEFIDSEpre_nbr = strcmp(map_name_selected{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                MGEFIDSEpre_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(MGEFIDSEpre_nbr));
                
                % load post MGEFIDSE map
                MGEFIDSEpost_nbr = strcmp(map_name_selected{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                MGEFIDSEpost_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(MGEFIDSEpost_nbr));
                
                % load ADCmap
                ADC_nbr = strcmp(map_name_selected{3}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                ADC_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(ADC_nbr));
                
                % Load mask (if selected and/or exist)
                if isempty(handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs)
                    Mask_filename = {''};
                else
                    Mask_nbr = strcmp(add_parameters{:}(2), handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                    Mask_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(Mask_nbr));
                    if isempty(Mask_filename)
                        Mask_filename = {''};
                    end
                end
                if ~isempty(MGEFIDSEpre_filename) && ~isempty(MGEFIDSEpost_filename) && ~isempty(ADC_filename) &&...
                        sum(MGEFIDSEpre_nbr) == 1  && sum(MGEFIDSEpost_nbr) == 1 && sum(ADC_nbr) == 1
                    
                    % compute the VSI map
                    [BVf_Num_struct, stO2_Num_struct,Vessel_orient_struct, ADC_Num_struct, R_Num_struct, R2_Num_struct]  = parametric_map_BVf_stO2_R(MGEFIDSEpre_filename{:},...
                        MGEFIDSEpost_filename{:},ADC_filename{:},Mask_filename{:}, add_parameters);
                    
                    if ~isempty(BVf_Num_struct)
                        % save the new maps andupdate MIA database
                        uvascim.image = BVf_Num_struct;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf-Num' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, BVf_Num_struct, patient(i), time_point(j), ['BVf-Num' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        MIA('MIA_update_database_display',  hObject, eventdata, findobj('Tag', 'MIA_GUI'));
                        clear BVf_Num_struct ;
                        
                        uvascim.image = stO2_Num_struct;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-stO2-Num' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, stO2_Num_struct, patient(i), time_point(j), ['stO2-Num' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear stO2_Num_struct ;
                        
                        uvascim.image = R_Num_struct;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R-Num' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, R_Num_struct, patient(i), time_point(j), ['R-Num' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear R_Num_struct ;
                        
                        uvascim.image = R2_Num_struct;
                        % save file
                        scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-R2-Num' file_name_extension '.mat'];
                        save(fullfile([path scan_name]),'uvascim');
                        clear uvascim;
                        
                        % add to the database and update handles.parameters and handles.clips
                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, R2_Num_struct, patient(i), time_point(j), ['R2-Num' file_name_extension], scan_name);
                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                        clear R2_Num_struct ;
                        
                        if ~isempty(Vessel_orient_struct)
                            uvascim.image = Vessel_orient_struct;
                            % save file
                            scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Vessel_orient-Num' file_name_extension '.mat'];
                            save(fullfile([path scan_name]),'uvascim');
                            clear uvascim;
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, Vessel_orient_struct, patient(i), time_point(j), ['Vessel_orient-Num' file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            clear Vessel_orient_struct ;
                        end
                        if ~isempty(ADC_Num_struct)
                            uvascim.image = ADC_Num_struct;
                            % save file
                            scan_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-ADC-Num' file_name_extension '.mat'];
                            save(fullfile([path scan_name]),'uvascim');
                            clear uvascim;
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, ADC_Num_struct, patient(i), time_point(j), ['ADC-Num' file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            clear ADC_Num_struct ;
                            
                        end
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf/stO2/R Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf/stO2/R Skiped'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf/stO2/R Skiped'];
                end
            else
                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-BVf/stO2/R already computed'];
            end
    end
    % save database after each patient analysis
    path_root = pwd;
    handles.MIA_data = guidata(handles.MIA_data.MIA_GUI);
    cd(handles.MIA_data.database(1).databaseinfo.pathname);
    handles.database(1).databaseinfo.VOIs = handles.MIA_data.VOIs;
    handles.database(1).databaseinfo.clips = handles.MIA_data.clips;
    database = handles.MIA_data.database; %#ok<NASGU>
    save(handles.MIA_data.database(1).databaseinfo.filename, 'database');
    cd(path_root);
end
%     end
%
% end
set(handles.CPM_patient_selected, 'String', 'Current process:  Done!');
%clear all nii files
nii_listing  = dir([handles.MIA_data.database(1).path '*.nii']);
if ~isempty(nii_listing)
    for x=1:size(nii_listing,1)
        delete([handles.MIA_data.database(1).path nii_listing(x).name])
    end
end

if ~isempty(logbook)
    listdlg('ListString', logbook', 'ListSize',[250 350], 'Name', 'logbook');
else
    msgbox('Done', 'logbook') ;
end
MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data)



% --- Executes on button press in CPM_close_button.
function CPM_close_button_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.CPM_GUI);


% --- Executes on button press in CPM_refresh_maps_name_listbox.
function CPM_refresh_maps_name_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_refresh_maps_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.MIA_data = guidata(handles.MIA_data.MIA_GUI);
parameters_list = [];
for i=1:numel(get(handles.MIA_data.MIA_name_list, 'String'))
    for j = 1:numel(handles.MIA_data.database(i).day)
        parameters_list = [parameters_list handles.MIA_data.database(i).day(j).parameters];
    end
end
handles.CPM_parameters_list = unique(parameters_list);
set(handles.CPM_maps_name_listbox, 'String', [{'NOT COMPUTED YET'} handles.CPM_parameters_list]');

% Update handles structure
guidata(hObject, handles);



function CPM_file_name_extension_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_file_name_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CPM_file_name_extension as text
%        str2double(get(hObject,'String')) returns contents of CPM_file_name_extension as a double


% --- Executes during object creation, after setting all properties.
function CPM_file_name_extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPM_file_name_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function CPM_help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CPM_Multiparametric_maps_figures_menu_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_Multiparametric_maps_figures_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = which('Compute_parametric_map_planche.pdf');
uiopen(filename,1)




% --- Executes on button press in CPM_save_protocol.
function CPM_save_protocol_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_save_protocol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'ready_to_go') ||handles.ready_to_go.logical == 0;
    warndlg('Please select every maps needed','Warning');
    return
end
cd(handles.MIA_data.database(1).databaseinfo.pathname)
[filename,pathname] = uiputfile('.mat','Please enter the name of the protocol');
if ischar(filename) == 0
    return
end
proto_to_save.maps_wanted_selected = handles.maps_wanted_selected;
proto_to_save.maps_needed = handles.maps_needed;
proto_to_save.map_needed_selected = handles.map_needed_selected;
proto_to_save.ready_to_go = handles.ready_to_go;
proto_to_save.map_name_selected = handles.map_name_selected;
proto_to_save.additional_info_selected = handles.additional_info_selected; %#ok<STRNU>

proto_to_save.CPM_control_box_listbox_Value = get(handles.CPM_control_box_listbox, 'Value');
proto_to_save.CPM_control_box_listbox_String = get(handles.CPM_control_box_listbox, 'String');

save(fullfile(pathname, filename), 'proto_to_save');



% --- Executes on button press in CPM_load_protocol.
function CPM_load_protocol_Callback(hObject, eventdata, handles)
% hObject    handle to CPM_load_protocol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd(handles.MIA_data.database(1).databaseinfo.pathname)
[filename,pathname] = uigetfile('.mat','Please select the protocol to load');
if ischar(filename) == 0
    return
end
load(fullfile(pathname, filename));

handles.maps_wanted_selected = proto_to_save.maps_wanted_selected;
handles.maps_needed = proto_to_save.maps_needed;
handles.map_needed_selected = proto_to_save.map_needed_selected;
handles.ready_to_go = proto_to_save.ready_to_go;
set(handles.CPM_go_button, 'BackgroundColor', [0 1 0]);

handles.map_name_selected = proto_to_save.map_name_selected;
handles.additional_info_selected = proto_to_save.additional_info_selected;

%set every display boxes
for i=1:numel(handles.maps_wanted_selected)
    tmp(i) = find(strcmp(get(handles.CPM_maps_wanted_listbox, 'String'), handles.maps_wanted_selected{i}), 1); %#ok<AGROW>
end
set(handles.CPM_maps_wanted_listbox, 'Value', tmp');

set(handles.CPM_maps_needed_listbox, 'String', handles.maps_needed');
set(handles.CPM_maps_needed_listbox, 'Value',  find(strcmp(get(handles.CPM_maps_needed_listbox, 'String'), handles.map_needed_selected), 1));

set(handles.CPM_maps_name_listbox, 'Value',  find(strcmp(get(handles.CPM_maps_name_listbox, 'String'),handles.map_name_selected ), 1));


set(handles.CPM_control_box_listbox, 'Value', proto_to_save.CPM_control_box_listbox_Value) ;
set(handles.CPM_control_box_listbox, 'String', proto_to_save.CPM_control_box_listbox_String) ;

% Update handles structure
guidata(hObject, handles);

