function varargout = Math_module(varargin)
% MATH_MODULE MATLAB code for Math_module.fig
%      MATH_MODULE, by itself, creates a new MATH_MODULE or raises the existing
%      singleton*.
%
%      H = MATH_MODULE returns the handle to a new MATH_MODULE or the handle to
%      the existing singleton*.
%
%      MATH_MODULE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATH_MODULE.M with the given input arguments.
%
%      MATH_MODULE('Property','Value',...) creates a new MATH_MODULE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Math_module_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Math_module_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Math_module

% Last Modified by GUIDE v2.5 22-Oct-2014 10:30:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Math_module_OpeningFcn, ...
    'gui_OutputFcn',  @Math_module_OutputFcn, ...
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


% --- Executes just before Math_module is made visible.
function Math_module_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Math_module (see VARARGIN)

% Choose default command line output for Math_module
handles.output = hObject;

% Variables initiation
handles.MIA_data = varargin{1};
handles.Math_parameters_list = varargin{2};
handles.Math_VOIs_list = varargin{3};
%% probl?me pour calculer les nouveau offset lor de transformation d'image (ratation/permutation)
handles.Math_operations_listing = {'Arithmetic', 'Mean slices', 'Smooth', 'Add slices', ...
    'SPM: Coreg',...
    'SPM: Realign&Coreg',...
    'SPM: Reslice','SPM: Realign', 'SPM: Realign (Over time)', 'Same registration as', 'Normalization', 'Repair outlier',...
    'Remove images', 'Shift images', 'Import Atlas (and ROI)','Export to Nifti'};
handles.Math_operations_listing = sort(handles.Math_operations_listing);
set(handles.Math_maps_listbox, 'String', handles.Math_parameters_list');
set(handles.Math_operation_listbox, 'String', [handles.Math_operations_listing]');
set(handles.Math_additional_info_table, 'Data', '', 'ColumnName', '', 'columnFormat', '',  'columnEditable', false)

% find out he time point list
% set popup menu
id_list = {}; tp_list = {}; group_list = {};
if isfield(handles.MIA_data, 'database')
    for i=1:numel(handles.MIA_data.database)
        for j = 1:numel(handles.MIA_data.database(i).day)
            tp_list = [tp_list {handles.MIA_data.database(i).day(j).date}];
        end
        id_list = [id_list {handles.MIA_data.database(i).name}];
        group_list = [group_list {handles.MIA_data.database(i).group}];
    end
end
handles.Math_id_list = sort(unique(id_list));
handles.Math_tp_list = sort(unique(tp_list));
handles.Math_group_list = sort(unique(group_list)); 

% patient selected name
handles.patient_selected.name_nbr = get(handles.MIA_data.MIA_name_list, 'Value');
handles.patient_selected.name = handles.MIA_data.database(handles.patient_selected.name_nbr).name;
handles.patient_selected.tp_nbr = get(handles.MIA_data.MIA_time_points_list, 'Value');
handles.patient_selected.tp = handles.MIA_data.database(handles.patient_selected.name_nbr).day(handles.patient_selected.tp_nbr).date;
set(handles.Math_patient_selected, 'String', ['Patient selected:  ' handles.patient_selected.name ' ' handles.patient_selected.tp]);
set(handles.Math_additional_info_table, 'Data', {});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Math_module wait for user response (see UIRESUME)
% uiwait(handles.Math_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = Math_module_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Math_file_name_extension_Callback(hObject, eventdata, handles)
% hObject    handle to Math_file_name_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Math_file_name_extension as text
%        str2double(get(hObject,'String')) returns contents of Math_file_name_extension as a double


% --- Executes during object creation, after setting all properties.
function Math_file_name_extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Math_file_name_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Math_maps_listbox.
function Math_maps_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Math_maps_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Math_maps_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Math_maps_listbox


% --- Executes during object creation, after setting all properties.
function Math_maps_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Math_maps_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Math_operation_listbox.
function Math_operation_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Math_operation_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Math_operation_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Math_operation_listbox
handles = guidata(hObject);
handles.operation_selected=handles.Math_operations_listing(get(handles.Math_operation_listbox,'Value'));

additional_info = {};


switch handles.operation_selected{1}
    case 'Arithmetic'
        additional_info_name = {'Operation', 'Day of the 2nd scan', 'Name of the 2nd scan'};
        additional_info_data = {'Addition', 'first', handles.Math_parameters_list{1}};
        additional_info_format = {{'Addition', 'Subtraction', 'Multiplication', 'Division', 'Percentage', 'Concentration'},...
            ['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list};
        additional_info_editable = [1 1 1];
    case 'Geometric'
        additional_info_name = {'Operation', 'If permute --> new order: x y echo z rep', 'If Rotate (degree of rotation)'};
        additional_info_data = {'Permute', '4 2 3 1 5', '90'};
        additional_info_format = {{'Permute', 'Rotate'}, 'numeric', 'numeric'};
        additional_info_editable = [1 1 1];
    case 'Mean slices'
        additional_info_name = {'First slice', 'Nbr of slice to sum'};
        additional_info_data = {'1','2'};
        additional_info_format = {'numeric', 'numeric'};
        additional_info_editable = [1 1];
    case 'Add slices'
        additional_info_name = {'First slice', 'Nbr of slice to sum'};
        additional_info_data = {'1','2'};
        additional_info_format = {'numeric', 'numeric'};
        additional_info_editable = [1 1];
    case 'Smooth'
        additional_info_name = {'Type', 'Hsize', 'Sigma'};
        additional_info_data = {'gaussian', '3', '1'};
        additional_info_format = {{'gaussian'}, 'numeric', 'numeric'};
        additional_info_editable = [1 1 1];
    case 'Registration'
        additional_info_name = {'Mode', 'iteration', 'Tolerance', 'Transformation', 'PyramidLevels', 'Interpolant', 'padmethod', 'Ref day', 'Ref Scan'};
        additional_info_data = {'multimodal', '250', '0.01', 'affine', '2', 'nearest', 'fill', 'first', handles.Math_parameters_list{1}};
        additional_info_format = {{'multimodal', 'monomodal'}, 'numeric', 'numeric', {'translation', 'rigid', 'similarity', 'affine'},...
            'numeric', {'cubic', 'linear', 'nearest'}, {'bound', 'circular', 'fill', 'replicate', 'symmetric'}, ['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list};
        additional_info_editable = [1 1 1 1 1 1 1 1 1];
    case 'Registration-Kroon'
        additional_info_name = {'Ref day', 'Ref Scan'};
        additional_info_data = {'first', handles.Math_parameters_list{1}};
        additional_info_format = {['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list};
        additional_info_editable = [1 1];
        
    case 'SPM: Coreg'
        additional_info_name = { 'Ref day', 'Ref Scan', '',  'Apply to other images?', 'Final resolution', 'Function',  'Separation', 'Tolerence','Hist Smooth', 'Interpolation', 'Warpping',...
            'Masking', 'Filename Prefix', 'Set the Origine?'};
        additional_info_data(1,1:14) = {'same', handles.Math_parameters_list{1}, 'None', false, 'Unchanged', 'mi', 'Auto= [slice thickness voxel_size voxel_size/2]', '0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001', '7 7' , '4th Degree B-Spline', 'No wrap', 'Dont mask images', 'Coreg', 'No'};
        additional_info_data(1:size(handles.Math_parameters_list,2),3) = handles.Math_parameters_list';
        additional_info_data(1:size(handles.Math_parameters_list,2),4) = {false}; 
        additional_info_format = {['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list, 'char', 'logical', {'Same as Ref', 'Unchanged', '64', '112', '128', '192', '256', '384', '512'}, {'ncc', 'mi', 'nmi','ecc'},  'numeric',...
            'numeric', 'numeric', {'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline',  '3rd Degree B-Spline',  '4th Degree B-Spline',  '5th Degree B-Spline',  '6th Degree B-Spline',  '7th Degree B-Spline'},...
            {'No wrap', 'Warp X', 'Warp Y', 'Warp X&Y', 'Warp Z', 'Warp X&Z', 'Warp Y&Z', 'Warp X,Y&Z'}, {'Dont mask images', 'Mask image'}, 'char', {'No', 'Yes'}};
        additional_info_editable = [1 1 0 1 1 1 1 1 1 1 1 1 1 1];
    case 'SPM: Realign&Coreg'
        additional_info_name = { 'Ref day', 'Ref Scan', 'Prefix'};
        additional_info_data(1,1:size(additional_info_name,2)) = {'same', handles.Math_parameters_list{1}, 'Rea-Coreg'};
        additional_info_format = {['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list, 'char' };
        additional_info_editable = [1 1 1];
        
    case 'SPM: Reslice'
        additional_info_name = { 'Ref day', 'Ref Scan', 'Final resolution', 'Interpolation', 'Warpping', 'Masking', 'Filename Prefix'};
        additional_info_data= {'same', handles.Math_parameters_list{1}, 'Unchanged', '4th Degree B-Spline', 'No wrap', 'Dont mask images', 'reslice'};
        additional_info_format = {['same', 'first', '-1', '+1', handles.Math_tp_list], handles.Math_parameters_list,...
            {'Same as Ref', 'Unchanged'}, {'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline',  '3rd Degree B-Spline',  '4th Degree B-Spline',  '5th Degree B-Spline',  '6th Degree B-Spline',  '7th Degree B-Spline'},...
            {'No wrap', 'Warp X', 'Warp Y', 'Warp X&Y', 'Warp Z', 'Warp X&Z', 'Warp Y&Z', 'Warp X,Y&Z'}, {'Dont mask images', 'Mask image'}, 'char'};
        additional_info_editable = [1 1 1 1 1 1 1];
        
    case 'Export to Nifti'
        additional_info_name = { ''};
        additional_info_data= {''};
        additional_info_format = {''};
        additional_info_editable = [];
        
    case 'SPM: Realign'
        additional_info_name = {'EstiOption: Quatlity', 'EstiOption: Sep', 'EstiOption: Smooth', 'EstiOption: Num Pass', 'EstiOption: Interpol', 'EstiOption: Warp','EstiOption: Weight',...
            'ResliceOption: Res Image', 'ResliceOption: Interpol', 'ResliceOption: Warp', 'ResliceOption: Mask', 'ResliceOption: Prefix', 'Is Bolus?'};
        additional_info_data= {'0.9', '4', '5', 'Register to mean', '2nd Degree B-Spline', 'No Wrap', '0 files',...
            'All Images + Mean Image', '4th Degree B-Spline', 'No Wrap', 'Mask images', 'realign', 'No'};
        additional_info_format = {'char', 'char', 'char', {'Register to mean', 'Register to first'},...
            {'Trilinear', '2nd Degree B-Spline',  '3rd Degree B-Spline',  '4th Degree B-Spline',  '5th Degree B-Spline',  '6th Degree B-Spline',  '7th Degree B-Spline'},...
            {'No Wrap', 'Warp X', 'Warp Y', 'Warp X&Y', 'Warp Z', 'Warp X&Z', 'Warp Y&Z', 'Warp X,Y&Z'},...
            ['0 files', handles.Math_parameters_list],...
            {'All files (1..n)', 'Images 2..n,', 'All Images + Mean Image', 'Mean Images Only'},...
            {'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline',  '3rd Degree B-Spline',  '4th Degree B-Spline',  '5th Degree B-Spline',  '6th Degree B-Spline',  '7th Degree B-Spline'},...
            {'No Wrap', 'Warp X', 'Warp Y', 'Warp X&Y', 'Warp Z', 'Warp X&Z', 'Warp Y&Z', 'Warp X,Y&Z'},...
            {'Mask images', 'Dont mask images'}, 'char', {'No' 'Yes'}};
        additional_info_editable = [1 1 1 1 1 1 1 1 1 1 1 1 1];
        
    case 'SPM: Realign (Over time)'
          additional_info_name = {'EstiOption: Quatlity', 'EstiOption: Sep', 'EstiOption: Smooth', 'EstiOption: Num Pass', 'EstiOption: Interpol', 'EstiOption: Warp','EstiOption: Weight',...
            'ResliceOption: Res Image', 'ResliceOption: Interpol', 'ResliceOption: Warp', 'ResliceOption: Mask', 'ResliceOption: Prefix', 'Is Bolus?'};
        additional_info_data= {'0.9', '4', '5', 'Register to mean', '2nd Degree B-Spline', 'No Wrap', '0 files',...
            'All Images + Mean Image', '4th Degree B-Spline', 'No Wrap', 'Mask images', 'realign', 'No'};
        additional_info_format = {'char', 'char', 'char', {'Register to mean', 'Register to first'},...
            {'Trilinear', '2nd Degree B-Spline',  '3rd Degree B-Spline',  '4th Degree B-Spline',  '5th Degree B-Spline',  '6th Degree B-Spline',  '7th Degree B-Spline'},...
            {'No Wrap', 'Warp X', 'Warp Y', 'Warp X&Y', 'Warp Z', 'Warp X&Z', 'Warp Y&Z', 'Warp X,Y&Z'},...
            ['0 files', handles.Math_parameters_list],...
            {'All files (1..n)', 'Images 2..n,', 'All Images + Mean Image', 'Mean Images Only'},...
            {'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline',  '3rd Degree B-Spline',  '4th Degree B-Spline',  '5th Degree B-Spline',  '6th Degree B-Spline',  '7th Degree B-Spline'},...
            {'No Wrap', 'Warp X', 'Warp Y', 'Warp X&Y', 'Warp Z', 'Warp X&Z', 'Warp Y&Z', 'Warp X,Y&Z'},...
            {'Mask images', 'Dont mask images'}, 'char', {'No' 'Yes'}};
        additional_info_editable = [1 1 1 1 1 1 1 1 1 1 1 1 1];
        
    case 'Same registration as'
        additional_info_name = {'Day', 'tform scan', 'Interpolant', 'padmethod'};
        additional_info_data = {'first', handles.Math_parameters_list{1}, 'linear', 'fill'};
        additional_info_format = {['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list, {'cubic', 'linear', 'nearest'}, {'bound', 'circular', 'fill', 'replicate', 'symmetric'}};
        additional_info_editable = [1 1 1 1];
    case 'Normalization'
        additional_info_name = {'VOI day', 'VOI used to normalized', 'Value of normalization'};
        additional_info_data = {'same', handles.Math_VOIs_list{1}, '1'};
        additional_info_format = {['same' handles.Math_tp_list], handles.Math_VOIs_list, 'char' };
        additional_info_editable = [1 1 1];
    case 'Repair outlier'
        additional_info_name = {'VOI used to identify outlier', 'Alpha Thomson', 'Interpolation', 'Plot diagnostique', 'Detrend', 'Mean map'};
        additional_info_data = {'All pixels', '1e-6', 'No', 'No', 'No', 'No'};
        additional_info_format = {['All pixels' handles.Math_VOIs_list], 'numeric', {'No' 'Yes'}, {'No' 'Yes'}, {'No' 'Yes'}, {'No' 'Yes'}};
        additional_info_editable = [1 1 1 1 1 1];
    case 'Remove images'
        additional_info_name = {'Echoes to remove', 'Slices to remove', 'Expts to remove'};
        additional_info_data = {'None', 'None', 'None'};
        additional_info_format = {'numeric', 'numeric', 'numeric'};
        additional_info_editable = [1 1 1];
    case 'Shift images'
        additional_info_name = {'X shift: Left (-) / Right (+)', 'Y shift: Down (-) / Up (+)', 'Modify offset?'};
        additional_info_data = {'0', '0', 'No'};
        additional_info_format = {'numeric', 'numeric', {'No' 'Yes'}};
        additional_info_editable = [1 1 1];
    case 'Import Atlas (and ROI)'
        additional_info_name = {'Scan list', 'Scan saved', 'ROI list', 'ROI saved', 'Threshold min', 'Threshold max' };
        soft_path =  fileparts(which('MIA.m'));
        ROI_listing = importdata(fullfile(soft_path, 'tools', 'Automatic_ROI', 'Atlas', 'label_key_atlas.txt'));
        additional_info_data(1:4,1) = {'Atlas', 'TPM_CSF', 'TPM_White', 'TPM_Gray' }';
        additional_info_data(1:4,2) = {true};
        additional_info_data(1:size(ROI_listing,1),3) = ROI_listing;
        additional_info_data(1:9,4) = {true};
        additional_info_data(10:size(ROI_listing,1),4) = {false};
        additional_info_data(1:9,5) = {'0','0', '0', '127.5', '127.5', '0', '48.5', '120', '250'};
        additional_info_data(1:9,6) = {'260', '260', '260', '260', '260', '48.5', '100', '130', '260'};
        additional_info_format = {'char', 'logical','char', 'logical', 'numeric', 'numeric'};
        additional_info_editable = [0 1 0 1 1 1];
        
    case 'Register to template'
        additional_info_name = {'',  'Apply to other images?'};
        additional_info_data(1:2,1) = {'None', false};
        additional_info_data(1:size(handles.Math_parameters_list,2),1) = handles.Math_parameters_list';
        additional_info_data(1:size(handles.Math_parameters_list,2),2) = {false};
        additional_info_format = {'char', 'logical'};
        additional_info_editable = [0 1];
        set(handles.Math_file_name_extension, 'String', '_rTemplate');
    otherwise
        additional_info_name = {''};
        additional_info_data = {''};
        additional_info_format = {''};
        additional_info_editable = 0;
end

% set title
set(handles.Math_additional_info_title, 'String', ['Additional info for: ' handles.operation_selected{1}]);
% set names of the columns
set(handles.Math_additional_info_table, 'ColumnName', additional_info_name);
% set data (default's parameters)
set(handles.Math_additional_info_table, 'Data', additional_info_data);
% set columnFormat (option for each parameters)
set(handles.Math_additional_info_table, 'columnFormat', additional_info_format);

% set ColumnWidth to auto

width= num2cell(max(cellfun(@length,[additional_info_data' additional_info_name']')).*7);
set(handles.Math_additional_info_table,'ColumnWidth',width)
% set each colomn editable
set(handles.Math_additional_info_table, 'columnEditable', logical(additional_info_editable));

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function Math_operation_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Math_operation_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Math_go_button.
function Math_go_button_Callback(hObject, eventdata, handles)
% hObject    handle to Math_go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


options = get(get(handles.Math_patient_option_uipanel, 'Children'), 'Value');
% list of Maps selected
map_selected_listing = get(handles.Math_maps_listbox, 'Value');
% Additional info
add_parameters = get(handles.Math_additional_info_table, 'Data');

logbook = {};

if  options{1} == 1 || options{2} == 1
    patient = handles.patient_selected.name_nbr;
else
    patient = 1:numel(handles.MIA_data.database);
end

% these Operation required user input
switch handles.operation_selected{1}
    case 'Import Atlas (and ROI)'
        for i = 1:numel(patient)
            if  options{1} == 1
                time_point = handles.patient_selected.tp_nbr;
            else
                time_point = 1:numel(handles.MIA_data.database(patient(i)).day);
            end
            if strcmp(handles.MIA_data.database(patient(i)).path(end), filesep) == 0
                path = strcat(handles.MIA_data.database(patient(i)).path, filesep);
            else
                path = handles.MIA_data.database(patient(i)).path;
            end
            for j = 1:numel(time_point)
                if numel(map_selected_listing) > 1
                    warndlg('Please select only 1 scan (T2w)!','Warning');
                    return
                end
                if  sum(strcmp('Atlas', handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                    % find the parameter position for this patient (if exist)
                    parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                    map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                    % Load maps
                    if ~isempty(map_filename) && sum(parameter_nbr) == 1
                        % load the parameter
                        fid=fopen(map_filename{:} ,'r');
                        if fid>0
                            fclose(fid);
                            data_loaded = load(map_filename{:});
                            
                        else
                            logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$',...
                                map_filename{:});
                            continue
                        end
                    else
                        continue
                    end
                else
                    continue
                end
                %Create the binary mask
%                 g = double(squeeze(data_loaded.uvascim.image.reco.data(:,:,1,:,1)));
%                 %     max_ref = max(g(list_ind));
%                 mask = g <= 3*nanmedian(g(:));
%                 % Connected component labelling and brain selection
%                 mask = ~mask;
%                 label = bwlabeln(mask,6);
%                 STATS = regionprops(label, 'Area');
%                 [~,plus_grand_label] = max([STATS.Area]);
%                 if ~isempty(plus_grand_label)
%                     mask (label ~= plus_grand_label) = 0;
%                 end
%                 mask = imerode(mask,strel('disk',3));
%                 label = bwlabeln(mask,6);
%                 STATS = regionprops(label, 'Area'); 
%                 [~,plus_grand_label] = max([STATS.Area]);
%                 if ~isempty(plus_grand_label)
%                     mask (label ~= plus_grand_label) = 0;
%                 end
%                 mask = imdilate(mask,strel('disk',3));
%                 data_loaded.uvascim.image.reco.data=data_loaded.uvascim.image.reco.data .* repmat(permute(mask, [1 2 4 3]), [1 1 size(data_loaded.uvascim.image.reco.data, 3), 1, size(data_loaded.uvascim.image.reco.data, 5)]);
%                 
                % convert Anat to .nii file
                file_name = handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr};
                convert_uvascim2nii(data_loaded,...
                    file_name(1:end-4),...
                    path);
%                 matlabbatch{1}.spm.util.reorient.srcfiles = {
%                     [path file_name(1:end-4) '.nii,1']
%                     };
                
%                 matlabbatch{1}.spm.util.reorient.transform.transprm = [312 -89 -24 0 0 0 10 10 10 0 0 0];
%                 matlabbatch{1}.spm.util.reorient.prefix = 's_';
%                 spm('defaults', 'FMRI');
%                 spm_jobman('initcfg');
%                 spm_jobman('run', matlabbatch(1));
%                 spm_image('Display', [path 's_' file_name(1:end-4) '.nii']);
%                 disp('press enter when the origin is updated')
%                 pause
%                 clear matlabbatch
               
                
                % Coregister: estimate to template (T2w data) + shift and
                soft_path =  fileparts(which('MIA.m'));
                matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 's_Template_Tohoku_TPM4D.nii') ',1']};

                matlabbatch{1}.spm.spatial.coreg.estimate.source = {[path file_name(1:end-4) '.nii,1']};
                matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};

                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [5 2 1 0.5 0.25 0.1];
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [3 3];
                [SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);

                spm('defaults', 'FMRI');
                spm_jobman('initcfg');
                spm_jobman('run', matlabbatch(1));
                clear matlabbatch
               
                close(SPMinter)
                close(SPMgraph)

            end
            
        end
%     case 'Register to template'
%         for i = 1:numel(patient)
%             if  options{1} == 1
%                 time_point = handles.patient_selected.tp_nbr;
%             else
%                 time_point = 1:numel(handles.MIA_data.database(patient(i)).day);
%             end
%             if strcmp(handles.MIA_data.database(patient(i)).path(end), filesep) == 0
%                 path = strcat(handles.MIA_data.database(patient(i)).path, filesep);
%             else
%                 path = handles.MIA_data.database(patient(i)).path;
%             end
%             for j = 1:numel(time_point)
%                 if numel(map_selected_listing) > 1
%                     warndlg('Please select only 1 scan (T2w)!','Warning');
%                     return
%                 end
%                 
%                 if  sum(strcmp([handles.Math_parameters_list{map_selected_listing} get(handles.Math_file_name_extension, 'String')], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
%                     % find the parameter position for this patient (if exist)
%                     parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
%                     map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
%                     
%                     % Load maps
%                     if ~isempty(map_filename) && sum(parameter_nbr) == 1
%                         % load the parameter
%                         fid=fopen(map_filename{:} ,'r');
%                         if fid>0
%                             fclose(fid);
%                             data_loaded = load(map_filename{:});
%                             
%                         else
%                             logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$',...
%                                 map_filename{:});
%                             continue
%                         end
%                     else
%                         continue
%                     end
%                 else
%                     continue
%                 end
%                 
%                 %Create the binary mask
%                 g = double(squeeze(data_loaded.uvascim.image.reco.data));
%                 %     max_ref = max(g(list_ind));
%                 mask = g <= 3*median(g(:));
%                 % Connected component labelling and brain selection
%                 mask = ~mask;
%                 label = bwlabeln(mask,6);
%                 STATS = regionprops(label, 'Area');
%                 [~,plus_grand_label] = max([STATS.Area]);
%                 if ~isempty(plus_grand_label)
%                     mask (label ~= plus_grand_label) = 0;
%                 end
%                 mask = imerode(mask,strel('disk',3));
%                 label = bwlabeln(mask,6);
%                 STATS = regionprops(label, 'Area');
%                 [~,plus_grand_label] = max([STATS.Area]);
%                 if ~isempty(plus_grand_label)
%                     mask (label ~= plus_grand_label) = 0;
%                 end
%                 mask = imdilate(mask,strel('disk',3));
%                 data_loaded.uvascim.image.reco.data=data_loaded.uvascim.image.reco.data .* permute(mask, [1 2 4 3]);
%                 
%                 % convert uvascim to .nii
%                 % check if Others images are selected
%                 increment = 1;
%                 other_map_filename = '';
%                 other_map_name = '';
%                 
%                 if sum(cell2mat(add_parameters(:,2)) == 1) ~= 0
%                     scan_checked = findn(cell2mat(add_parameters(:,2)) == 1 );
%                     for x = 1: size(scan_checked,1)
%                         others_map_nbr = strcmp( handles.Math_parameters_list(scan_checked(x,1)), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
%                         if sum(others_map_nbr) == 1
%                             other_map_filename{increment}= handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{others_map_nbr};
%                             other_map_name{increment} = handles.MIA_data.database(patient(i)).day(time_point(j)).parameters{others_map_nbr};
%                             increment = increment+ 1;
%                         end
%                     end
%                 end
%                 
%                 
%                 list_of_uvascim_file = [strcat(handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}) other_map_filename]';
%                 list_of_parameter = [handles.Math_parameters_list{map_selected_listing} other_map_name];
%                 
%                 nii_filename = '';
%                 for x = 1:size(list_of_parameter,2)
%                     map_filename = list_of_uvascim_file{x};
%                     if x ~= 1
%                         data_loaded = load([path  map_filename]);
%                     end
%                     convert_uvascim2nii(data_loaded,...
%                         map_filename(1:end-4),...
%                         path);
%                     nii_filename{x,1} = [map_filename(1:end-4) '.nii,1'];
%                 end
%                 
%                 matlabbatch{1}.spm.util.reorient.srcfiles = strcat(repmat({path}, size(nii_filename)), nii_filename);
%                 matlabbatch{1}.spm.util.reorient.transform.transprm = [312 -89 -24 0 0 0 10 10 10 0 0 0];
%                 matlabbatch{1}.spm.util.reorient.prefix = 's_';
%                 spm('defaults', 'FMRI');
%                 spm_jobman('initcfg');
%                 spm_jobman('run', matlabbatch(1));
%                 
%                 
%                 
%                 spm_image('Display', [path 's_' list_of_uvascim_file{1}(1:end-4) '.nii']);
%                 disp('press enter when the origin is updated')
%                 pause
%                 clear matlabbatch
%                 
%                 % Coregister: estimate to template (T2w data) + shift and
%                 soft_path =  fileparts(which('MIA.m'));
%                 matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 'Template_Tohoku_TPM4D.nii') ',1']};
%                 matlabbatch{1}.spm.spatial.coreg.estimate.source = {[path 's_' nii_filename{1}]};
%                 matlabbatch{1}.spm.spatial.coreg.estimate.other = strcat(repmat({strcat(path, 's_')}, size(nii_filename(2:end))), nii_filename(2:end));
%                 matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'ncc';
%                 matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
%                 matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%                 matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
%                 
%                 spm('defaults', 'FMRI');
%                 spm_jobman('initcfg');
%                 spm_jobman('run', matlabbatch(1));
%                 clear matlabbatch
%             end
%             
%         end
end



%% output's option
output_option =  get(get(handles.Math_output_option_uipanel, 'Children'), 'Value');
if  output_option{1} == 1
    if isempty(get(handles.Math_file_name_extension, 'String'))
        % The file name extension is not mendatory for this functions
        if  strcmp(handles.operation_selected{1}, 'Normalization') || strcmp(handles.operation_selected{1}, 'SPM: Coreg') ...
                || strcmp(handles.operation_selected{1}, 'SPM: Reslice') || strcmp(handles.operation_selected{1}, 'Import Atlas (and ROI)')...
                || strcmp(handles.operation_selected{1}, 'SPM: Realign') || strcmp(handles.operation_selected{1}, 'SPM: Realign&Coreg') ...
                 || strcmp(handles.operation_selected{1}, 'SPM: Realign (Over time)' )|| strcmp(handles.operation_selected{1}, 'Clustering: GMM')   
            file_name_extension = '';
        else
            warndlg('Please enter a file name extension!','Warning');
            return
        end
    else
        file_name_extension = get(handles.Math_file_name_extension, 'String');
    end
    
else
    user_response = questdlg('Do you want to overwrite your file?','Warning', 'Yes', 'No', 'No');
    if strcmp(user_response, 'No')
        return
    end
    user_response = questdlg('Do you REALLY want to overwrite your file?','Warning', 'No', 'Yes', 'No');
    % overwrite the file if requested
    if strcmp(user_response, 'Yes')
        file_name_extension = '';
    else
        return
    end
    
end

for i = 1:numel(patient)
    if  options{1} == 1
        time_point = handles.patient_selected.tp_nbr;
    else
        time_point = 1:numel(handles.MIA_data.database(patient(i)).day);
    end
    for j = 1:numel(time_point)
        if strcmp(handles.MIA_data.database(patient(i)).path(end), filesep) == 0
            path = strcat(handles.MIA_data.database(patient(i)).path, filesep);
        else
            path = handles.MIA_data.database(patient(i)).path;
        end
        handles.MIA_data = guidata(handles.MIA_data.MIA_GUI);
        set(handles.Math_patient_selected, 'String', ['Current process:  ' handles.MIA_data.database(patient(i)).name ' ' handles.MIA_data.database(patient(i)).day(time_point(j)).date...
            '...' handles.operation_selected{1}]);
        drawnow;
        switch handles.operation_selected{1}
            case 'Arithmetic'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        image1_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        image1_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(image1_nbr));
                        % Get info entered by the user
                        Operation = add_parameters{1};
                        image2_day = add_parameters{2};
                        image2_scan = add_parameters{3};
                        
                        switch image2_day;
                            case 'same'
                                image2_nbr = strcmp(image2_scan, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                                image2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(image2_nbr));
                            case 'first'
                                image2_nbr = strcmp(image2_scan, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
                                image2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(1).scans_file(image2_nbr));
                            case '-1'
                                if time_point(j) > 1
                                    image2_nbr = strcmp(image2_scan, handles.MIA_data.database(patient(i)).day(time_point(j)-1).parameters) == 1;
                                    image2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)-1).scans_file(image2_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-First_time_point'];
                                    continue
                                end
                            case '+1'
                                if time_point(j) < numel(handles.MIA_data.database(patient(i)).day)
                                    image2_nbr = strcmp(image2_scan, handles.MIA_data.database(patient(i)).day(time_point(j)+1).parameters) == 1;
                                    image2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)+1).scans_file(image2_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-Last_time_point selected'];
                                    continue
                                end
                            otherwise
                                match_image2_map = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',image2_day)==1);
                                if ~isempty(match_image2_map)
                                    image2_nbr = strcmp(image2_scan, handles.MIA_data.database(patient(i)).day(match_image2_map).parameters) == 1;
                                    image2_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_image2_map).scans_file(image2_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                        handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                        'No reference scance named ' image2_day '-' image2_scan];
                                    continue
                                end
                                
                        end
                        if  ~isempty(image1_filename) && ~isempty(image2_filename) &&...
                                sum(image1_nbr) == 1  && sum(image2_nbr) == 1    &&...
                                strcmp(image2_filename,image1_filename) == 0
                            
                            % Load data
                            image1 = load(image1_filename{:});
                            image2 = load(image2_filename{:});
                            % check if both scan have the same geometry
                            if ~(isequal(size(image1.uvascim.image.reco.data), size(image2.uvascim.image.reco.data)) &&...
                                    isequal( round(image2.uvascim.image.reco.fov_offsets(:,1,:)*100)/100, round(image2.uvascim.image.reco.fov_offsets(:,1,:)*100)/100))
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped- differences between geometries'];
                                continue
                                
                            end
                            result_map = image2.uvascim.image;
                            
                            % Perform the arithmetic operation!
                            switch add_parameters{1}
                                case 'Addition'
                                    result_map.reco.data = image1.uvascim.image.reco.data + image2.uvascim.image.reco.data;
                                    Methode_info = '';
                                case 'Subtraction'
                                    result_map.reco.data = image1.uvascim.image.reco.data - image2.uvascim.image.reco.data;
                                    Methode_info = '';
                                case 'Multiplication'
                                    result_map.reco.data = image1.uvascim.image.reco.data .* image2.uvascim.image.reco.data;
                                    Methode_info = '';
                                case'Division'
                                    result_map.reco.data = image1.uvascim.image.reco.data ./ image2.uvascim.image.reco.data;
                                    Methode_info = '';
                                case 'Percentage'
                                    result_map.reco.data = ((image1.uvascim.image.reco.data - image2.uvascim.image.reco.data) ./ image1.uvascim.image.reco.data) .* 100;
                                    Methode_info = '';
                                case 'Concentration'
                                    question = strcat('Please enter the relaxivity (mM-1.sec-1)');
                                    relaxivity  = inputdlg(question,'Relaxivity',1,{'6'});
                                    if isempty(relaxivity)
                                        return
                                    end
                                    result_map.reco.data = (1./(image1.uvascim.image.reco.data/1000) - 1/(image2.uvascim.image.reco.data/1000)) / str2double(relaxivity{:}) ;
                                    Methode_info = ['relaxivity = ' relaxivity{:} ' mM-1.sec-1']; 
                            end
                            % Save information into the structure of the
                            % new file
                            ParamConfig=sprintf('##$QuantifMethod=Arithmetic maps (math module) --> %s\n##$First scan=%s\n##$First scan file=%s\n##$Second scan=%s\n##$Second scan file=%s\n##$Other Info %s',...
                                add_parameters{1},...
                                handles.Math_parameters_list{map_selected_listing(ii)},...
                                image1_filename{:},...
                                image2_scan,...
                                image2_filename{:},...
                                Methode_info); 
                            result_map.reco.paramQuantif = ParamConfig;
                            
                            if ~isempty(result_map)
                                % save the new map and update IA database
                                uvascim.image = result_map;
                                % save file
                                result_map_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ...
                                    handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '.mat'];
                                result_map_name = strrep(result_map_name, '*', 'star');
                                save(fullfile([path result_map_name]),'uvascim');
                                clear uvascim;
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, result_map, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], result_map_name);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                            end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear result_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end
                
            case 'Geometric'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1
                            % load the file
                            fid=fopen(map_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(map_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the data\n##$%s\n##$',...
                                    map_filename);
                                continue
                            end
                            if strcmp(handles.operation_selected{1}, 'Permute')
                                permutation_order = str2num(add_parameters{2}); %#ok<ST2NM>
                                data_permuted = permute(uvascim.image.reco.data,permutation_order);
                                % figure; imshow3D(squeeze(data_permuted(:,:,1,:,1)))
                                % save info
                                if isfield(uvascim.image.reco, 'paramQuantif')
                                    ParamConfig=sprintf('##$Original file name=%s\n##$Original file info=%s\n##$Math-Operation= Permute dim\n##$New order= %s\n',...
                                        map_filename{:},...
                                        uvascim.image.reco.paramQuantif,...
                                        add_parameters{2});
                                else
                                    ParamConfig=sprintf('##$Original file name=%s\n##$Math-Operation= Permute dim\n##$New order= %s\n',...
                                        map_filename{:},...
                                        add_parameters{2});
                                end
                                uvascim.image.reco.no_samples = size(data_permuted, 1);
                                uvascim.image.reco.no_views = size(data_permuted, 2);
                                uvascim.image.reco.no_echoes = size(data_permuted, 3);
                                uvascim.image.reco.no_slices = size(data_permuted, 4);
                                uvascim.image.reco.no_expts = size(data_permuted, 5);
                                uvascim.image.reco.fov = [size(data_permuted, 1) size(data_permuted, 2) size(data_permuted, 4)];
                                uvascim.image.reco.phaselabel(1: uvascim.image.reco.no_slices) =  uvascim.image.reco.phaselabel(1);
                                uvascim.image.reco.label(1:uvascim.image.reco.no_slices) =  uvascim.image.reco.label(1);
                                uvascim.image.reco.scaling_factor(1:uvascim.image.reco.no_slices) =  uvascim.image.reco.scaling_factor(1);
                                uvascim.image.reco.scaling_offset(1:uvascim.image.reco.no_slices) =  uvascim.image.reco.scaling_offset(1);
                                
                                %% probleme pour calculer les nouveaux offset!!!!!!!
                                %                                 uvascim.image.reco.fov_offsets(:,uvascim.image.reco.no_slices,uvascim.image.reco.no_slices) =  uvascim.image.reco.fov_offsets(:,1,1)
                                %                                 uvascim.image.reco.fov_orientation =
                                %                                 uvascim.image.reco.fov_phase_orientation =
                                
                            elseif strcmp(handles.operation_selected{1}, 'Rotate')
                                
                            end
                            
                            
                            
                            uvascim.image.reco.paramQuantif = ParamConfig;
                            uvascim.image.reco.data = data_permuted;
                            
                            
                            % save file
                            [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                            scan_name = [name file_name_extension ext];
                            save(fullfile([path, scan_name]),'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            clear uvascim;
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                            
                            
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
            case {'Mean slices', 'Add slices'}
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1
                            % load the file
                            fid=fopen(map_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(map_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the data\n##$%s\n##$',...
                                    map_filename);
                                continue
                            end
                            [X, Y, echos3d, depths3d, expt3d]=size(uvascim.image.reco.data);
                            %                             nbr_of_final_slice = ceil((depths3d-str2double(add_parameters{1})+1)/str2double(add_parameters{2}))-1;
                            SliceAv = (str2double(add_parameters{1}):floor((depths3d-str2double(add_parameters{1})+1)/str2double(add_parameters{2}))*str2double(add_parameters{2})+str2double(add_parameters{1})-1);
                            SliceAv = reshape(SliceAv,[str2double(add_parameters{2}) numel(SliceAv)/str2double(add_parameters{2})])';
                            if isempty(SliceAv)
                                return
                            end
                            data_permuted = permute(uvascim.image.reco.data,[1 2 4 3 5]);
                            output_data = zeros([X Y size(SliceAv,1) echos3d expt3d]);
                            for xx = 1:size(SliceAv,1)
                                for k=1:echos3d
                                    if strcmp(handles.operation_selected{1}, 'Mean slices')
                                        output_data(:,:,xx,:) = nanmean(data_permuted(:,:,SliceAv(xx,:),:,:),3);
                                    else
                                        output_data(:,:,xx,:) = sum(data_permuted(:,:,SliceAv(xx,:),:,:),3);
                                    end
                                end
                            end
                            output_data =  permute(output_data,[1 2 4 3 5]);
                            % save info
                            if isfield(uvascim.image.reco, 'paramQuantif')
                                ParamConfig=sprintf('##$Original file name=%s\n##$Original file info=%s\n##$MathOperation=Mean slices\n##$First slice=%s\n##$Nbr of slice added=%s\n',...
                                    map_filename{:},...
                                    uvascim.image.reco.paramQuantif,...
                                    add_parameters{1},...
                                    add_parameters{2});
                            else
                                ParamConfig=sprintf('##$Original file name=%s\n##$MathOperation=Mean slices\n##$First slice=%s\n##$Nbr of slice added=%s\n',...
                                    map_filename{:},...
                                    add_parameters{1},...
                                    add_parameters{2});
                            end
                            uvascim.image.reco.paramQuantif = ParamConfig;
                            uvascim.image.reco.data = output_data;
                            clear temp temp2
                            
                            % update structure
                            uvascim.image.reco.no_slices = size(SliceAv,1);
                            if isfield(uvascim.image.reco, 'fov') && size(uvascim.image.reco.fov, 1) == 3
                                uvascim.image.reco.fov(3)=[];
                            end
                            
                            for xx = 1:size(SliceAv,1)
                                temp=squeeze(uvascim.image.reco.fov_offsets(3,1,str2double(add_parameters{1})+(xx*str2double(add_parameters{2}))-...
                                    str2double(add_parameters{2}):str2double(add_parameters{1})+(xx*str2double(add_parameters{2}))-1));
                                new_offset=median(temp);
                                uvascim.image.reco.fov_offsets(3,:,xx) = new_offset;
                            end
                            uvascim.image.reco.fov_offsets(:,:,size(SliceAv,1)+1:end)=[];
                            uvascim.image.reco.fov_orientation(:,:,size(SliceAv,1)+1:end)=[];
                            uvascim.image.reco.fov_phase_orientation(:,size(SliceAv,1)+1:end)=[];
                            
                            uvascim.image.reco.globalmax = max( uvascim.image.reco.data(:));
                            uvascim.image.reco.globalmin = min( uvascim.image.reco.data(:));
                            uvascim.image.reco.phaselabel(:,size(SliceAv,1)+1:end)=[];
                            if isfield(uvascim.image.reco, 'scaling_factor')
                                uvascim.image.reco.scaling_factor(:,size(SliceAv,1)+1:end)=[];
                                uvascim.image.reco.scaling_offset(:,size(SliceAv,1)+1:end)=[];
                            end
                            
                            uvascim.image.reco.thickness =  uvascim.image.reco.thickness*str2double(add_parameters{2});
                            uvascim.image.reco.unit(:,size(SliceAv,1)+1:end)=[];
                            uvascim.image.reco.no_slices = size(SliceAv,1);
                            uvascim.image.reco=orderfields( uvascim.image.reco);
                            
                            uvascim.image.clip = [uvascim.image.reco.globalmin  uvascim.image.reco.globalmax 1];
                            
                            % save file
                            [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                            scan_name = [name file_name_extension ext];
                            save(fullfile([path, scan_name]),'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
%                             MIA('MIA_update_database_display', hObject, eventdata, handles.MIA_data)
                            clear uvascim;
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                            
                            
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
            case 'Smooth'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1
                            % load the file
                            fid=fopen(map_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(map_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the data\n##$%s\n##$',...
                                    map_filename);
                                continue
                            end
                            h = fspecial(add_parameters{1},str2double(add_parameters{2}),str2double(add_parameters{3}));
                            [~, ~, E, Z] = size(uvascim.image.reco.data);
                            for s=1:Z
                                for e=1:E
                                    uvascim.image.reco.data(:,:,e,s) = imfilter(uvascim.image.reco.data(:,:,e,s),h,'replicate');
                                end
                            end
                            % save info
                            if isfield(uvascim.image.reco, 'paramQuantif')
                                ParamConfig=sprintf('##$Original file name=%s\n##$Original file info=%s\n##$MathOperation=Smoothing type %s\n##$Hsize=%s\n##$Sigma=%s\n',...
                                    map_filename{:},...
                                    uvascim.image.reco.paramQuantif,...
                                    add_parameters{1},...
                                    add_parameters{2},...
                                    add_parameters{3});
                            else
                                ParamConfig=sprintf('##$Original file name=%s\n##$MathOperation=Smoothing type %s\n##$Hsize=%s\n##$Sigma=%s\n',...
                                    map_filename{:},...
                                    add_parameters{1},...
                                    add_parameters{2},...
                                    add_parameters{3});
                            end
                            uvascim.image.reco.paramQuantif = ParamConfig;
                            
                            % save file
                            [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                            scan_name = [name file_name_extension ext];
                            save(fullfile([path, scan_name]),'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            clear uvascim;
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                            
                            
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
                
            case {'Registration', 'Registration-Kroon'}
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        moving_map_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        moving_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(moving_map_nbr));
                        if strcmp(handles.operation_selected{1}, 'Registration')
                            ref_day = add_parameters{8};
                            ref_scan = add_parameters{9};
                        else
                            ref_day = add_parameters{1};
                            ref_scan = add_parameters{2};
                        end
                        switch ref_day;
                            case 'same'
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(fixed_map_nbr));
                            case 'first'
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(1).scans_file(fixed_map_nbr));
                            case '-1'
                                if time_point(j) > 1
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)-1).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)-1).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-First_time_point'];
                                    continue
                                end
                            case '+1'
                                if time_point(j) < numel(handles.MIA_data.database(patient(i)).day)
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)+1).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)+1).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-Last_time_point selected'];
                                    continue
                                end
                            otherwise
                                match_fixed_map = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',ref_day)==1);
                                if ~isempty(match_fixed_map)
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(match_fixed_map).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_fixed_map).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                        handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                        'No reference scance named ' ref_day '-' ref_scan];
                                    continue
                                end
                                
                        end
                        if  ~isempty(moving_map_filename) && ~isempty(fixed_map_filename) &&...
                                sum(moving_map_nbr) == 1  && sum(fixed_map_nbr) == 1    &&...
                                strcmp(fixed_map_filename,moving_map_filename) == 0
                            
                            % Perform the registration 2 methods!
                            if strcmp(handles.operation_selected{1}, 'Registration')
                                moved_map = compute_registration_map(fixed_map_filename{:}, moving_map_filename{:}, add_parameters);
                            else
                                %                                 moved_map = compute_registration_map_SPM(fixed_map_filename{:}, moving_map_filename{:});
                                moved_map = compute_registration_map_Kroon(fixed_map_filename{:}, moving_map_filename{:});
                                
                                
                                %% code to add the signal of 2 scans
                                %                             un = load(fixed_map_filename{:});
                                %                             deux = load(moving_map_filename{:});
                                %                             moved_map = deux.uvascim.image;
                                %                             moved_map.reco.data = un.uvascim.image.reco.data + deux.uvascim.image.reco.data;
                                
                            end
                            if ~isempty(moved_map)
                                % save the new map and update IA database
                                uvascim.image = moved_map;
                                % save file
                                moved_map_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ...
                                    handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '.mat'];
                                moved_map_name = strrep(moved_map_name, '*', 'star');
                                save(fullfile([path moved_map_name]),'uvascim');
                                clear uvascim;
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, moved_map, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], moved_map_name);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                            end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear moved_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end
            case 'SPM: Coreg'
                for ii = 1:numel(map_selected_listing)
                    cd(handles.MIA_data.database(1).path);
                    if  sum(strcmp([add_parameters{1,size(add_parameters,2)-1} handles.Math_parameters_list{map_selected_listing(ii)}], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        moving_map_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        moving_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(moving_map_nbr));
                        ref_day = add_parameters{1,1};
                        ref_scan = add_parameters{1,2};
                        switch ref_day;
                            case 'same'
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(fixed_map_nbr));
                            case 'first'
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(1).scans_file(fixed_map_nbr));
                            case '-1'
                                if time_point(j) > 1
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)-1).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)-1).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-First_time_point'];
                                    continue
                                end
                            case '+1'
                                if time_point(j) < numel(handles.MIA_data.database(patient(i)).day)
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)+1).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)+1).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-Last_time_point selected'];
                                    continue
                                end
                            otherwise
                                match_fixed_map = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',ref_day)==1);
                                if ~isempty(match_fixed_map)
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(match_fixed_map).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_fixed_map).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                        handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                        'No reference scance named ' ref_day '-' ref_scan];
                                    
                                    continue
                                end
                                
                        end
                        if  ~isempty(moving_map_filename) && ~isempty(fixed_map_filename) &&...
                                sum(moving_map_nbr) == 1  && sum(fixed_map_nbr) == 1   % &&...
                            %strcmp(fixed_map_filename,moving_map_filename) == 0
                            
                            % check if Others images are selected
                            increment = 1;
                            other_map_filename = '';
                            other_map_name = '';
                            if sum(cell2mat(add_parameters(:,4)) == 1) ~= 0
                                scan_checked = findn(cell2mat(add_parameters(:,4)) == 1 );
                                for x = 1: size(scan_checked,1)
                                    others_map_nbr = strcmp( handles.Math_parameters_list(scan_checked(x,1)), handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                                    if sum(others_map_nbr) == 1
                                        other_map_filename{increment}= strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{others_map_nbr});
                                        other_map_name{increment} = handles.MIA_data.database(patient(i)).day(time_point(j)).parameters{others_map_nbr};
                                        increment = increment+ 1;
                                    end
                                end
                            end
                            
                           [Coreg_maps, data_to_save] = Math_module_SPM_coreg(fixed_map_filename{:}, moving_map_filename{:},other_map_filename,add_parameters);
                            
                           if Coreg_maps == 1 && ~isempty(data_to_save)
                               for z = 1:size(data_to_save,2)
                                   % save file
                                   if z==1
                                       coreg_map_name = [handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' add_parameters{1,size(add_parameters,2)-1} handles.Math_parameters_list{map_selected_listing(ii)} '.mat'];
                                   else
                                       coreg_map_name = [handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' add_parameters{1,size(add_parameters,2)-1} other_map_name{z-1} '.mat'];
                                   end
                                   coreg_map_name = strrep(coreg_map_name, '*', 'star'); coreg_map_name = {coreg_map_name};
                                   
                                   uvascim = data_to_save(z).uvascim;
                                   save(fullfile([path coreg_map_name{:}]),'uvascim');
                                   
                                   % add to the database and update handles.parameters and handles.clips
                                   if z==1
                                       handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j),...
                                           [add_parameters{1,size(add_parameters,2)-1} handles.MIA_data.database(patient(i)).day(time_point(j)).parameters{moving_map_nbr}], coreg_map_name);
                                   else
                                       handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j),...
                                           [add_parameters{1,size(add_parameters,2)-1} other_map_name{z-1}], coreg_map_name);
                                   end
                                   guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                    logbook{numel(logbook)+1} =[coreg_map_name{:} '- Done'];
                               end
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                               clear uvascim;
                           end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear moved_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end
            case 'SPM: Realign&Coreg'
                cd(handles.MIA_data.database(1).path);
                %                 for ii = 1:numel(map_selected_listing)
                % ii = 1 because all scan selected will be process togather
                ii=1;
                if  sum(strcmp([add_parameters{1,size(add_parameters,2)} handles.Math_parameters_list{map_selected_listing(ii)}], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                    % find the parameter position for this patient (if exist)
                    
                    ref_day = add_parameters{1,1};
                    ref_scan = add_parameters{1,2};
                    switch ref_day;
                        case 'same'
                            fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                            fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(fixed_map_nbr));
                        case 'first'
                            fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
                            fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(1).scans_file(fixed_map_nbr));
                        case '-1'
                            if time_point(j) > 1
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)-1).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)-1).scans_file(fixed_map_nbr));
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                    [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-First_time_point'];
                                continue
                            end
                        case '+1'
                            if time_point(j) < numel(handles.MIA_data.database(patient(i)).day)
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)+1).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)+1).scans_file(fixed_map_nbr));
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                    [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-Last_time_point selected'];
                                continue
                            end
                        otherwise
                            match_fixed_map = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',ref_day)==1);
                            if ~isempty(match_fixed_map)
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(match_fixed_map).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_fixed_map).scans_file(fixed_map_nbr));
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                    handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                    'No reference scance named ' ref_day '-' ref_scan];
                                continue
                            end    
                    end
                    
                    % scans to realign and co-reg
                    for z = 1:numel(map_selected_listing)
                        if ~isempty(find(strcmp(handles.Math_parameters_list{map_selected_listing(z)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)==1, 1))
                            Scan_to_realign_and_coreg_nbr(z) = find(strcmp(handles.Math_parameters_list{map_selected_listing(z)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)==1);
                            Scan_to_realign_and_coreg_filename(z,1) = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(Scan_to_realign_and_coreg_nbr(z)));
                        else
                            Scan_to_realign_and_coreg_nbr(z) = NaN;
                            Scan_to_realign_and_coreg_filename(z,1) = {'Do not exist'};        
                        end
                    end
                     Scan_to_realign_and_coreg_filename(isnan(Scan_to_realign_and_coreg_nbr)) =  [];
                    Scan_to_realign_and_coreg_nbr(isnan(Scan_to_realign_and_coreg_nbr)) =  [];
                    if  ~isempty(Scan_to_realign_and_coreg_filename) && ~isempty(fixed_map_filename) && sum(fixed_map_nbr) == 1
                        
                        [realign_an_coreg_maps, data_to_save] = Math_module_SPM_realign_and_coreg(fixed_map_filename{:}, Scan_to_realign_and_coreg_filename,add_parameters);
                        if realign_an_coreg_maps == 1 && ~isempty(data_to_save)
                            for z = 1:numel(map_selected_listing)
                                % save file
                                coregd_map_name = [handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' add_parameters{3} handles.Math_parameters_list{map_selected_listing(z)} '.mat'];
                                coregd_map_name = strrep(coregd_map_name, '*', 'star'); coregd_map_name = {coregd_map_name};
                                uvascim = data_to_save(z).uvascim;
                                save(fullfile([path coregd_map_name{:}]),'uvascim');
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j),...
                                    [add_parameters{3} handles.MIA_data.database(patient(i)).day(time_point(j)).parameters{Scan_to_realign_and_coreg_nbr(z)}], coregd_map_name);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            end
                            clear uvascim;
                            
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                    end
                    clear moved_map
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                end
                %                 end
                
                
            case 'SPM: Reslice'
                cd(handles.MIA_data.database(1).path);
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([add_parameters{1,size(add_parameters,2)} handles.Math_parameters_list{map_selected_listing(ii)}], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        
                        ref_day = add_parameters{1,1};
                        ref_scan = add_parameters{1,2};
                        switch ref_day;
                            case 'same'
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(fixed_map_nbr));
                            case 'first'
                                fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
                                fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(1).scans_file(fixed_map_nbr));
                            case '-1'
                                if time_point(j) > 1
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)-1).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)-1).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-First_time_point'];
                                    continue
                                end
                            case '+1'
                                if time_point(j) < numel(handles.MIA_data.database(patient(i)).day)
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(time_point(j)+1).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)+1).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-Last_time_point selected'];
                                    continue
                                end
                            otherwise
                                match_fixed_map = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',ref_day)==1);
                                if ~isempty(match_fixed_map)
                                    fixed_map_nbr = strcmp(ref_scan, handles.MIA_data.database(patient(i)).day(match_fixed_map).parameters) == 1;
                                    fixed_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_fixed_map).scans_file(fixed_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                        handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                        'No reference scance named ' ref_day '-' ref_scan];
                                    
                                    continue
                                end
                                
                        end
                        % scan to relice
                        Scan_to_reslice_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        Scan_to_reslice_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(Scan_to_reslice_nbr));
                        
                        if  ~isempty(Scan_to_reslice_filename) && ~isempty(fixed_map_filename) && sum(fixed_map_nbr) == 1  
                            [resliced_maps, uvascim] = Math_module_SPM_reslice(fixed_map_filename{:}, Scan_to_reslice_filename,add_parameters);
                            if resliced_maps == 1 && ~isempty(uvascim)
                                % save file
                                resliced_map_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ...
                                    add_parameters{1,size(add_parameters,2)} handles.Math_parameters_list{map_selected_listing(ii)} '.mat'];
                                resliced_map_name = strrep(resliced_map_name, '*', 'star');
                                save(fullfile([path resliced_map_name]),'uvascim');
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [add_parameters{1,size(add_parameters,2)} handles.MIA_data.database(patient(i)).day(time_point(j)).parameters{Scan_to_reslice_nbr}], resliced_map_name);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                clear uvascim;
                                
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                            end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear moved_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end
                
            case 'SPM: Realign'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([add_parameters{1,size(add_parameters,2)-1} handles.Math_parameters_list{map_selected_listing(ii)}], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % scan to relice
                        Scan_to_realign_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        Scan_to_realign_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(Scan_to_realign_nbr));
                        
                        if  ~isempty(Scan_to_realign_filename)
                            
                            [realign_maps, uvascim] = Math_module_SPM_realign( Scan_to_realign_filename,add_parameters);
                            if realign_maps == 1 && ~isempty(uvascim)
                                if size(uvascim,2) > 1
                                    for z= 1:size(uvascim,2)
                                        tmp(:,:, z,:,size(uvascim(1).uvascim.image.reco.data,5))=uvascim(z).uvascim.image.reco.data(:,:,1,:,1);
                                    end
                                    uvascim = uvascim.uvascim;
                                    uvascim.image.reco.data = tmp;
                                else
                                    uvascim = uvascim.uvascim;
                                end
                                
                                % save file
                                realign_map_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ...
                                    add_parameters{1,size(add_parameters,2)-1} handles.Math_parameters_list{map_selected_listing(ii)} '.mat'];
                                realign_map_name = strrep(realign_map_name, '*', 'star');
                                save(fullfile([path realign_map_name]),'uvascim');
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [add_parameters{1,size(add_parameters,2)-1} handles.MIA_data.database(patient(i)).day(time_point(j)).parameters{Scan_to_realign_nbr}], realign_map_name);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                clear uvascim;
                                
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                            end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear moved_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end
                
             case 'SPM: Realign (Over time)'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([add_parameters{1,size(add_parameters,2)-1} handles.Math_parameters_list{map_selected_listing(ii)}], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % scan to relice
%                         Scan_to_realign_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
%                         Scan_to_realign_filename(1,:) = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(Scan_to_realign_nbr));
                        Scan_to_realign_filename = {};
                        for tp_nbr = 1:numel(handles.MIA_data.database(patient(i)).day)
                            Scan_to_realign_nbr  = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(tp_nbr).parameters) == 1;
                            if sum(Scan_to_realign_nbr) == 1
                            Scan_to_realign_filename(size(Scan_to_realign_filename,1)+1,:) = strcat(path, handles.MIA_data.database(patient(i)).day(tp_nbr).scans_file(Scan_to_realign_nbr));
                            end
                        end
                        if  ~isempty(Scan_to_realign_filename)
                            
                            [Status, realign_maps] = Math_module_SPM_realign( Scan_to_realign_filename,add_parameters);
                            if Status == 1 && ~isempty(realign_maps)
                                map_nbr = 1;
                                for tp_nbr = 1:numel(handles.MIA_data.database(patient(i)).day)
                                    Scan_to_realign_nbr  = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(tp_nbr).parameters) == 1;
                                    if sum(Scan_to_realign_nbr) == 1
                                        % save file
                                        realign_map_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(tp_nbr).date '-' ...
                                            add_parameters{1,size(add_parameters,2)-1} handles.Math_parameters_list{map_selected_listing(ii)} '.mat'];
                                        realign_map_name = strrep(realign_map_name, '*', 'star');
                                        uvascim = realign_maps(map_nbr).uvascim;
                                        save(fullfile([path realign_map_name]),'uvascim');
                                        
                                        % add to the database and update handles.parameters and handles.clips
                                        handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), tp_nbr, [add_parameters{1,size(add_parameters,2)-1} handles.MIA_data.database(patient(i)).day(tp_nbr).parameters{Scan_to_realign_nbr}], realign_map_name);
                                        guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                        MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                        map_nbr = map_nbr+1;
                                        clear uvascim;
                                    end
                                end
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                            end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear moved_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end   
            case 'Same registration as'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        moving_map_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        moving_map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(moving_map_nbr));
                        
                        switch add_parameters{1};
                            case 'same'
                                tform_map_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                                tform_map_nbr_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(tform_map_nbr));
                            case 'first'
                                tform_map_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(1).parameters) == 1;
                                tform_map_nbr_filename = strcat(path, handles.MIA_data.database(patient(i)).day(1).scans_file(tform_map_nbr));
                            case '-1'
                                if time_point(j) > 1
                                    tform_map_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(time_point(j)-1).parameters) == 1;
                                    tform_map_nbr_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)-1).scans_file(tform_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-First_time_point'];
                                    continue
                                end
                            case '+1'
                                if time_point(j) < numel(handles.MIA_data.database(patient(i)).day)
                                    tform_map_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(time_point(j)+1).parameters) == 1;
                                    tform_map_nbr_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)+1).scans_file(tform_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date...
                                        [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped-Last_time_point selected'];
                                    continue
                                end
                            otherwise
                                match_fixed_map = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',add_parameters{1})==1);
                                if ~isempty(match_fixed_map)
                                    tform_map_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(match_fixed_map).parameters) == 1;
                                    tform_map_nbr_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_fixed_map).scans_file(tform_map_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                        handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                        'No reference scance named ' add_parameters{1} '-' add_parameters{2}];
                                    continue
                                end
                        end
                        
                        if  ~isempty(moving_map_filename) && ~isempty(tform_map_nbr_filename) &&...
                                sum(moving_map_nbr) == 1  && sum(tform_map_nbr) == 1
                            % Perform the registration
                            moved_map = same_registration_as_map(moving_map_filename{:}, tform_map_nbr_filename{:}, add_parameters);
                            if ~isempty(moved_map)
                                % save the new map and update IA database
                                uvascim.image = moved_map;
                                % save file
                                moved_map_name = [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' ...
                                    handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '.mat'];
                                moved_map_name = strrep(moved_map_name, '*', 'star');
                                save(fullfile([path moved_map_name]),'uvascim');
                                clear uvascim;
                                
                                % add to the database and update handles.parameters and handles.clips
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, moved_map, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], moved_map_name);
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Done'];
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                            end
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Skyped'];
                        end
                        clear moved_map
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension '- Already computed'];
                    end
                end
            case 'Normalization'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp(['r' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        
                        % find the VOI position for this patient (if exist)
                        switch add_parameters{1}
                            case 'same'
                                VOI_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                                VOI_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(VOI_nbr));
                            otherwise
                                match_tp = find(strcmp({handles.MIA_data.database(patient(i)).day.date}',add_parameters{1})==1);
                                if ~isempty(match_tp)
                                    VOI_nbr = strcmp(add_parameters{2}, handles.MIA_data.database(patient(i)).day(match_tp).VOIs) == 1;
                                    VOI_filename = strcat(path, handles.MIA_data.database(patient(i)).day(match_tp).VOIs_file(VOI_nbr));
                                else
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name '-'...
                                        handles.MIA_data.database(patient(i)).day(time_point(j)).date '-'...
                                        'No VOI ' add_parameters{1} '-' add_parameters{2} ' to normalize this scan' ];
                                    continue
                                end
                        end
                        
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1 && ...
                                ~isempty(VOI_filename) && sum(VOI_nbr) == 1
                            % load the parameter
                            fid=fopen(map_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(map_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$',...
                                    map_filename);
                                continue
                            end
                            % load the ROI
                            fid=fopen(VOI_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(VOI_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the VOI file\n##$%s\n##$',...
                                    VOI_filename);
                                continue
                            end
                            VOI_matrix = zeros(size(uvascim.image.reco.data));
                            for iii = 1:uvascim.image.reco.no_slices
                                for jjj = 1:numel(uvascroi)
                                    if abs(uvascim.image.reco.fov_offsets(3,1,iii) - uvascroi(jjj).fov_offsets(3)) < 0.0001
                                        VOI_matrix(:,:,:,iii,:)=  imresize(uvascroi(jjj).value(:,:),[size(VOI_matrix, 1) size(VOI_matrix,2)],'bilinear');
                                    end
                                end
                            end
                            VOI_matrix(VOI_matrix == 0) = NaN;
                            tmp = double(uvascim.image.reco.data).*VOI_matrix;
                            uvascim.image.reco.data = double(uvascim.image.reco.data)/nanmean(tmp(:));
                            uvascim.image.reco.data =  uvascim.image.reco.data .* str2double(add_parameters{3});
                            
                            % save info
                            if isfield(uvascim.image.reco, 'paramQuantif')
                                ParamConfig=sprintf('##$Original file info=%s\n%s\n##$MathOperation=Normalization\n##$VOIs used=%s\n',...
                                    map_filename{:},...
                                    uvascim.image.reco.paramQuantif,...
                                    add_parameters{1});
                            else
                                ParamConfig=sprintf('##$Original file info=%s\n##$MathOperation=Normalization\n##$VOIs used=%s\n',...
                                    map_filename{:},...
                                    add_parameters{1});
                            end
                            uvascim.image.reco.paramQuantif = ParamConfig;
                            
                            % save file
                            [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                            scan_name = [name '-norm' file_name_extension ext];
                            save(fullfile([path, scan_name]),'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), ['r' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            clear uvascim;
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-r' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                            
                            
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-r' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped -No scan or no VOI'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-r' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
                
            case 'Repair outlier'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp(['outlier_expt_removed-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        
                        % find the VOI position for this patient (if exist)
                        VOI_nbr = strcmp(add_parameters{1}, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs) == 1;
                        if VOI_nbr ~= 0
                            VOI_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).VOIs_file(VOI_nbr));
                        else
                            VOI_filename = cell(1);
                        end
                        
                        
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1
                            
                            alpha = str2double(add_parameters{2});
                            interpolation = char2logical(add_parameters{3});
                            plotDiagnostics = char2logical(add_parameters{4});
                            detrend = char2logical(add_parameters{5});
                            MeanMap = char2logical(add_parameters{6});
                            
                            output = repair_outlier_slices(map_filename,VOI_filename,alpha,interpolation,plotDiagnostics,detrend,MeanMap);
                            
                            if ~isempty(output)
                                % save info
                                infoOutlier = {'outlier_expt_removed','outlier_exclude_mean'};
                                for itOutput = 1:numel(output)
                                    if isfield(output{itOutput}.reco, 'paramQuantif')
                                        ParamConfig=sprintf('##$Original file info=%s\n%s\n##$MathOperation=Repair outlier\n##$Map output=%s\n##$VOIs used=%s\n##$alpha=%d\n##$interpolation=%s\n##$plotDiagnostics=%s\n##$detrend=%s\n##$MeanMap=%s\n',...
                                            map_filename{:},...
                                            output{1}.reco.paramQuantif,...
                                            infoOutlier{itOutput},...
                                            add_parameters{1},alpha,interpolation,plotDiagnostics,detrend,MeanMap);
                                    else
                                        ParamConfig=sprintf('##$Original file info=%s\n##$MathOperation=Repair outlier\n##$Map output=%s\n##$VOIs used=%s\n##$alpha=%d\n##$interpolation=%s\n##$plotDiagnostics=%s\n##$detrend=%s\n##$MeanMap=%s\n',...
                                            map_filename{:},...
                                            infoOutlier{itOutput},...
                                            add_parameters{1},alpha,interpolation,plotDiagnostics,detrend,MeanMap);
                                    end
                                    output{1}.reco.paramQuantif = ParamConfig;
                                    uvascim.image = output{itOutput};
                                    
                                    % save file
                                    [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                                    scan_name = [name '-' infoOutlier{itOutput} file_name_extension ext];
                                    save(fullfile([path, scan_name]),'uvascim');
                                    
                                    % add to the database and update handles.parameters and handles.clips
                                    handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [infoOutlier{itOutput} '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                                    guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                    MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                    clear uvascim;
                                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-' infoOutlier{itOutput} '-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                                end
                            else
                                logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-outlier_expt_removed-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped -No scan or no VOI'];
                            end
                            
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-outlier_expt_removed-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped -No scan or no VOI'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date ['-outlier_expt_removed-' handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
            case 'Remove images'
                for ii = 1:numel(map_selected_listing)
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        
                        % Load maps
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1
                            % load the parameter
                            fid=fopen(map_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(map_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$',...
                                    map_filename{:});
                                continue
                            end
                            
                            % Retrieve number of elements (echoes, slices,
                            % experiments) to remove
                            ImData = uvascim.image.reco.data;
                            KeepIm = cell(1,3);
                            ListIm = {'Echoes', 'Slices', 'Expts'};
                            test = 0;
                            for itParam = 1:size(add_parameters,2)
                                    switch add_parameters{itParam}
                                        case 'None'
                                            KeepIm{itParam} = 1:size(ImData,itParam+2);
                                        otherwise
%                                             KeepIm{itParam} = setdiff(1:size(ImData,itParam+2),eval(add_parameters{itParam}));
                                            KeepIm{itParam} = setdiff(1:size(ImData,itParam+2),eval(strrep(add_parameters{itParam}, 'end', num2str(size(ImData,itParam+2)))));
                                            if isempty(KeepIm{itParam})
                                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$All %s ask to removed so no %s removed \n##$',...
                                                    map_filename{:},ListIm{itParam},ListIm{itParam});
                                                KeepIm{itParam} = 1:size(ImData,itParam+2);
                                                test = 1;
                                            end
                                    end
                            end
                            if test == 1
                                continue
                            end
                            uvascim.image.reco.data = ImData(:,:,KeepIm{1},KeepIm{2},KeepIm{3});
%                             fov_offsets = [3 echoes slice repet];
%                             fov_orientation = [9 echoes slice repet];
                            % fov_phase_orientation = [echoes slice repet]
                            % scaling_factor = [echoes slice repet]
                            % scaling_offset = [echoes slice repet]
                            % phaselabel = [echoes slice repet]
                            % label = [echoes slice repet]
                            
                            uvascim.image.reco.fov_offsets = uvascim.image.reco.fov_offsets(:,KeepIm{1},KeepIm{2},KeepIm{3});
                            uvascim.image.reco.fov_orientation = uvascim.image.reco.fov_orientation(:,KeepIm{1},KeepIm{2},KeepIm{3});
                            uvascim.image.reco.fov_phase_orientation = uvascim.image.reco.fov_phase_orientation(:,KeepIm{1},KeepIm{2},KeepIm{3});
                            uvascim.image.reco.no_echoes = numel(KeepIm{1});
                            uvascim.image.reco.no_slices = numel(KeepIm{2});
                            uvascim.image.reco.no_expts = numel(KeepIm{3});
                            if numel(KeepIm{1}) >= numel(uvascim.image.reco.echotime)
                            else
                                uvascim.image.reco.echotime = uvascim.image.reco.echotime(KeepIm{1});
                            end
                            %% if the nb of expt is modified for a dynamic scan we need to update the ##$PVM_ScanTimeStr= in the "image.texte" 
                            
                            % need to be done!
                            
                            %%
                            uvascim.image.reco.label = uvascim.image.reco.label(KeepIm{1},KeepIm{2},KeepIm{3});
                            if numel(uvascim.image.reco.phaselabel) >= max([numel(KeepIm{1}) numel(KeepIm{2}) numel(KeepIm{3})])
                                uvascim.image.reco.phaselabel = uvascim.image.reco.phaselabel(KeepIm{1},KeepIm{2},KeepIm{3});
                            end
                            if numel(uvascim.image.reco.scaling_factor) >= max([numel(KeepIm{1}) numel(KeepIm{2}) numel(KeepIm{3})])
                                uvascim.image.reco.scaling_factor = uvascim.image.reco.scaling_factor(KeepIm{1},KeepIm{2},KeepIm{3});
                            end
                            if numel(uvascim.image.reco.scaling_offset) >= max([numel(KeepIm{1}) numel(KeepIm{2}) numel(KeepIm{3})])
                                uvascim.image.reco.scaling_offset = uvascim.image.reco.scaling_offset(KeepIm{1},KeepIm{2},KeepIm{3});
                            end
                            
                            % save file
                            [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                            scan_name = [name '-del' file_name_extension ext];
                            save(fullfile([path, scan_name]),'uvascim', '-v7.3');
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            clear uvascim;
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                                                 
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped -No scan'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
            case 'Shift images'
                for ii = 1:numel(map_selected_listing)
                    if ~exist('user_response','var')
                        user_response = 'No';
                    end
                    if  sum(strcmp([handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0 || strcmp('Yes',user_response)
                        % find the parameter position for this patient (if exist)
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing(ii)}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        
                        % Load maps
                        if ~isempty(map_filename) && sum(parameter_nbr) == 1
                            % load the parameter
                            fid=fopen(map_filename{:} ,'r');
                            if fid>0
                                fclose(fid);
                                load(map_filename{:});
                            else
                                logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$',...
                                    map_filename{:});
                                continue
                            end
                            
                            x_shift_in_pixel=str2double(add_parameters{1});
                            y_shift_in_pixel=str2double(add_parameters{2});
                            ChangeOffsetYN = char2logical(add_parameters{3});
                            x_shift_in_mm=x_shift_in_pixel*uvascim.image.acq.fov(2)/double(uvascim.image.reco.no_views);
                            y_shift_in_mm=y_shift_in_pixel*uvascim.image.acq.fov(1)/double(uvascim.image.reco.no_samples);
                            
                            for mexpt=1:uvascim.image.reco.no_expts
                                for mslice=1:uvascim.image.reco.no_slices
                                    for mecho=1:uvascim.image.reco.no_echoes
                                        %On stocke l'image a decaler dans un tampon
                                        tempim=squeeze(uvascim.image.reco.data(:,:,mecho,mslice,mexpt));
                                        %On decale
                                        tempim = imshift_hany( tempim, [-y_shift_in_pixel x_shift_in_pixel] );
                                        %On remet l'image decalee dans la pile
                                        uvascim.image.reco.data(:,:,mecho,mslice,mexpt)=tempim;
                                        %On modifie le cas echeant les offsets associes a l'image
                                        if ChangeOffsetYN
                                            uvascim.image.reco.fov_offsets(:,mecho,mslice,mexpt)=uvascim.image.reco.fov_offsets(:,mecho,mslice,mexpt) + double([x_shift_in_mm; y_shift_in_mm; 0]);
                                        end
                                    end
                                end
                            end
                            
                            % save file
                            [~, name, ext] = fileparts([path handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr}]);
                            scan_name = [name '-shift' file_name_extension ext];
                            save(fullfile([path, scan_name]),'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j), [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension], scan_name);
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            clear uvascim;
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Done'];
                        else
                            logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-Skiped -No scan'];
                        end
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date [handles.Math_parameters_list{map_selected_listing(ii)} file_name_extension] '-already computed'];
                    end
                end
            case 'Import Atlas (and ROI)'
                if numel(map_selected_listing) > 1
                    warndlg('Please select only 1 scan (T2w)!','Warning');
                    return
                end
                
                if  sum(strcmp('Atlas', handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                    if sum(strcmp(handles.Math_parameters_list{map_selected_listing}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 1
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        file_name = handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr};
                        
                        % check if the data exist
                         if exist([path file_name(1:end-4) '.nii'], 'file') ~= 2
%                         if exist([path 's_' file_name(1:end-4) '.nii'], 'file') ~= 2
                            continue
                        end
                        clear matlabbatch
                        % Job1 = Normalize the TPM maps and the atlas to the anat
                        
                        matlabbatch{1}.spm.tools.oldseg.data = {
                            [path file_name(1:end-4) '.nii,1']
%                             [path 's_' file_name(1:end-4) '.nii,1']
                            };
                        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 1];
                        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 1];
                        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 1];
                        matlabbatch{1}.spm.tools.oldseg.output.biascor = 1;
                        matlabbatch{1}.spm.tools.oldseg.output.cleanup = 0;
                        soft_path =  fileparts(which('MIA.m'));
                        matlabbatch{1}.spm.tools.oldseg.opts.tpm = {
                            [fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 's_Template_Tohoku_gray.nii') ',1']
                            [fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 's_Template_Tohoku_white.nii') ',1']
                            [fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 's_Template_Tohoku_csf.nii') ',1']};
                        matlabbatch{1}.spm.tools.oldseg.opts.ngaus = [2
                            2
                            2
                            4];
                        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'mni';
                        matlabbatch{1}.spm.tools.oldseg.opts.warpreg = 1;
                        matlabbatch{1}.spm.tools.oldseg.opts.warpco = 2.5;
                        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.0001;
                        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 6;
                        matlabbatch{1}.spm.tools.oldseg.opts.samp = 0.3;
                        matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};
                        
                        % Job3 = co-register the atlas to the anat
                        % (referencial of the rat)
%                         matlabbatch{2}.spm.tools.oldnorm.write.subj.matname = {[path 's_' file_name(1:end-4) '_seg_inv_sn.mat']};
                        matlabbatch{2}.spm.tools.oldnorm.write.subj.matname = {[path file_name(1:end-4) '_seg_inv_sn.mat']};
                        copyfile(fullfile(soft_path, 'tools', 'Automatic_ROI', 'Atlas', 'atlas-resized.nii'),...
                            [path 'Atlas-resized_' file_name(1:end-4) '.nii'], 'f');
                        matlabbatch{2}.spm.tools.oldnorm.write.subj.resample = {[path 'Atlas-resized_' file_name(1:end-4) '.nii,1']};
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.preserve = 0;
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.bb = [NaN NaN NaN
                            NaN NaN NaN];
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.vox = [NaN NaN NaN];
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.interp = 1;
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.prefix = 'w';
                        
                        
                        % run jobs
                        spm('defaults', 'FMRI');
                        spm_jobman('initcfg');
                        spm_jobman('run', matlabbatch);
                        
                        %save atlas
                        scan_prefix = {'wAtlas-resized_', 'c1', 'c2', 'c3' };
                        map_name = {'Atlas' 'TPM_Gray', 'TPM_White', 'TPM_CSF'};
                        
                        
                        % Load T2w map
                        data_loaded = load(map_filename{:});
                        
                        
                        scan_checked = find(cell2mat(add_parameters(1:4,2)) == 1 );
                        if ~isempty(scan_checked)
                            for ii = 1:numel(scan_checked)
                                
                                tmp = convert_nii2uvascim([path scan_prefix{scan_checked(ii)} file_name(1:end-4) '.nii']);
                                uvascim = data_loaded.uvascim;
                                uvascim.image.reco.data = tmp.image.reco.data;
                                
                                ParamConfig=sprintf('##Altas gave by Tohoku University and pipeline developped by Margaux Faucher');
                                uvascim.image.reco.paramQuantif = ParamConfig;
                                
                                % save file
                                moved_map_name = [path handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' map_name{scan_checked(ii)} '.mat'];
                                moved_map_name = strrep(moved_map_name, '*', 'star');
                                save(moved_map_name,'uvascim');
                                
                                % add to the database and update handles.parameters and handles.clips
                                %                        MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j),...
                                    map_name{scan_checked(ii)},...
                                    [handles.MIA_data.database(patient(i)).name '-' handles.MIA_data.database(patient(i)).day(time_point(j)).date '-' map_name{scan_checked(ii)} '.mat']);
                                
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                                
                                clear uvascim;
                            end
                        end
                        
                        % delete unused file
                          delete(fullfile(path, [file_name(1:end-4) '_seg_inv_sn.mat']))
                        delete(fullfile(path, [ file_name(1:end-4) '_seg_sn.mat']))
                        
                        ROI_checked = find(cell2mat(add_parameters(:,4)) == 1 );
                        
                        % threshold used to create the ROI
                        if ~isempty(ROI_checked)
                            for ii = 1:numel(ROI_checked)
                                if strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_Brain') || strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_Cortex-R') || ...
                                        strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_Fit') ||...
                                        strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_Cortex-L') || strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_Striat-R') || ...
                                        strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_Striat-L')
                                    map_to_load = 'Atlas';
                                    
                                elseif strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_CSF')
                                    map_to_load = 'TPM_CSF';
                                elseif strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_WM')
                                    map_to_load = 'TPM_White';
                                elseif strcmp(add_parameters(ROI_checked(ii),3), 'Atlas_GM')
                                    map_to_load = 'TPM_Gray';
                                end
                                
                                handles.MIA_data =  MIA('MIA_save_automatic_ROI_form_Map', hObject, eventdata, handles.MIA_data, map_to_load, patient(i),...
                                    time_point(j),...
                                    [str2double(add_parameters{ROI_checked(ii),5}) str2double(add_parameters{ROI_checked(ii),6})],...
                                    add_parameters{ROI_checked(ii),3});
                                
                                guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                                MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            end
                        end
                        guidata(hObject, handles)
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Done'];
                    else
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Skiped -No scan'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '_Atlas-already computed'];
                end
                
            case 'Register to template'
                if numel(map_selected_listing) > 1
                    warndlg('Please select only 1 scan (T2w)!','Warning');
                    return
                end
                if  sum(strcmp([handles.Math_parameters_list{map_selected_listing} get(handles.Math_file_name_extension, 'String')], handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 0
                    if sum(strcmp(handles.Math_parameters_list{map_selected_listing}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters)) == 1
                        parameter_nbr = strcmp(handles.Math_parameters_list{map_selected_listing}, handles.MIA_data.database(patient(i)).day(time_point(j)).parameters) == 1;
                        map_filename = strcat(path, handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file(parameter_nbr));
                        file_name = handles.MIA_data.database(patient(i)).day(time_point(j)).scans_file{parameter_nbr};
                        
                        % check if the data exist
                        if exist([path 's_' file_name(1:end-4) '.nii'], 'file') ~= 2
                            continue
                        end
                        clear matlabbatch
                        % Job1 = Normalize the TPM maps and the atlas to the anat
                        matlabbatch{1}.spm.tools.oldseg.data = {
                            [path 's_' file_name(1:end-4) '.nii,1']
                            };
                        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 1];
                        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 1];
                        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 1];
                        matlabbatch{1}.spm.tools.oldseg.output.biascor = 1;
                        matlabbatch{1}.spm.tools.oldseg.output.cleanup = 0;
                        soft_path =  fileparts(which('MIA.m'));
                        matlabbatch{1}.spm.tools.oldseg.opts.tpm = {
                            [fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 'Template_Tohoku_gray.nii') ',1']
                            [fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 'Template_Tohoku_white.nii') ',1']
                            [fullfile(soft_path, 'tools', 'Automatic_ROI', 'TPM', 'Template_Tohoku_csf.nii') ',1']};
                        matlabbatch{1}.spm.tools.oldseg.opts.ngaus = [2
                            2
                            2
                            4];
                        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'subj';
                        matlabbatch{1}.spm.tools.oldseg.opts.warpreg = 1;
                        matlabbatch{1}.spm.tools.oldseg.opts.warpco = 25;
                        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.001;
                        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 60;
                        matlabbatch{1}.spm.tools.oldseg.opts.samp = 3;
                        matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};
                        
                        % Job2 = co-register the maps to the template
                        
                        matlabbatch{2}.spm.tools.oldnorm.write.subj.matname = {[path 's_' file_name(1:end-4) '_seg_sn.mat']};
                        matlabbatch{2}.spm.tools.oldnorm.write.subj.resample = strcat(repmat({strcat(path, 's_')}, size(nii_filename)), nii_filename);
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.preserve = 0;
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.bb = [NaN NaN NaN
                            NaN NaN NaN];
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.vox = [NaN NaN NaN];
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.interp = 1; % 1 = trilinear et 0 = nearst neighour
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
                        matlabbatch{2}.spm.tools.oldnorm.write.roptions.prefix = 'w';
                        
                        % run jobs
                        spm('defaults', 'FMRI');
                        spm_jobman('initcfg');
                        spm_jobman('run', matlabbatch);
                        
                        % save files
                        ParamConfig=sprintf('##Altas gave by Tohoku University and pipeline developped by Margaux Faucher');
                        for x = 1:size(list_of_parameter,2)
                            map_filename = list_of_uvascim_file{x};
                            uvascim = convert_nii2uvascim([path 'ws_' map_filename(1:end-4) '.nii']);
                            
                            uvascim.image.reco.paramQuantif = ParamConfig;
                            save([path map_filename(1:end-4) file_name_extension '.mat'],'uvascim');
                            
                            % add to the database and update handles.parameters and handles.clips
                            handles.MIA_data =  MIA('MIA_add_new_scan_to_database', handles.MIA_data, uvascim.image, patient(i), time_point(j),...
                                [list_of_parameter{x} file_name_extension],...
                                [map_filename(1:end-4) file_name_extension '.mat']);
                            
                            guidata(findobj('Tag', 'MIA_GUI'), handles.MIA_data)
                            MIA('MIA_update_database_display', hObject, eventdata,handles.MIA_data) ;
                            
                            clear uvascim;
                        end
                        
                        % delete unused file
                        delete(fullfile(path, ['s_' file_name(1:end-4) '_seg_inv_sn.mat']))
                        delete(fullfile(path, ['s_' file_name(1:end-4) '_seg_sn.mat']))
                        
                        guidata(hObject, handles)
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Done'];
                    else
                        
                        logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '-Skiped -No scan'];
                    end
                else
                    logbook{numel(logbook)+1} =[handles.MIA_data.database(patient(i)).name handles.MIA_data.database(patient(i)).day(time_point(j)).date '_Atlas-already computed'];
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
        
        for_spm_file_listing  = dir(fullfile(handles.MIA_data.database(1).path, '*for_spm*'));
        if ~isempty(for_spm_file_listing)
            for x=1:size(for_spm_file_listing,1)
                delete([handles.MIA_data.database(1).path for_spm_file_listing(x).name])
            end
        end  
    end   
end

%clear all nii files
nii_listing  = dir([handles.MIA_data.database(1).path '*.nii']);
if ~isempty(nii_listing)
    for x=1:size(nii_listing,1)
        delete([handles.MIA_data.database(1).path nii_listing(x).name])
    end
end
% %clear all for_spm files
% for_spm_file_listing  = dir(fullfile(handles.MIA_data.database(1).path, '*for_spm*'));
% if ~isempty(for_spm_file_listing)
%     for x=1:size(for_spm_file_listing,1)
%         delete([handles.MIA_data.database(1).path for_spm_file_listing(x).name])
%     end
% end
if strcmp(char(eventdata.Source.Callback), '@(hObject,eventdata)Math_module(''Math_go_button_Callback'',hObject,eventdata,guidata(hObject))')
    if ~isempty(logbook)
        listdlg('ListString', logbook', 'ListSize',[350 350], 'Name', 'logbook');
    else
        msgbox('Done', 'logbook') ;
    end
end
set(handles.Math_patient_selected, 'String', 'Current process:  ');

% --- Executes on button press in Math_close_button.
function Math_close_button_Callback(hObject, eventdata, handles)
% hObject    handle to Math_close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.Math_GUI);

% --- Executes on button press in Math_refresh_maps_name_listbox.
function Math_refresh_maps_name_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Math_refresh_maps_name_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MIA_data = guidata(handles.MIA_data.MIA_GUI);
parameters_list = [];
for i=1:numel(get(handles.MIA_data.MIA_name_list, 'String'))
    for j = 1:numel(handles.MIA_data.database(i).day)
        parameters_list = [parameters_list handles.MIA_data.database(i).day(j).parameters];
    end
end
handles.Math_parameters_list = unique(parameters_list);
set(handles.Math_maps_listbox, 'String', handles.Math_parameters_list');

% patient selected name
handles.patient_selected.name_nbr = get(handles.MIA_data.MIA_name_list, 'Value');
handles.patient_selected.name = handles.MIA_data.database(handles.patient_selected.name_nbr).name;
handles.patient_selected.tp_nbr = get(handles.MIA_data.MIA_time_points_list, 'Value');
handles.patient_selected.tp = handles.MIA_data.database(handles.patient_selected.name_nbr).day(handles.patient_selected.tp_nbr).date;
set(handles.Math_patient_selected, 'String', ['Patient selected:  ' handles.patient_selected.name ' ' handles.patient_selected.tp]);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Math_operation_help.
function Math_operation_help_Callback(hObject, eventdata, handles)
% hObject    handle to Math_operation_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

operation_selected = handles.operation_selected;
switch  operation_selected{:}
    case 'SPM: Coreg'
        [PATHSTR,~,~] = fileparts(which ('spm'));
        open(fullfile(PATHSTR, 'man', filesep, 'manual.pdf'));
end

