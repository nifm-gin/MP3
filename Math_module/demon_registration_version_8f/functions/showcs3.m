function varargout = showcs3(varargin)
% Function showcs3(V) will display x,y,z cross sections of the volume V
% The planes can be rotated with the mouse, and slice number with the
% sliders.
%
% Example,
%   load mri; D=squeeze(D);
%   scales=[1 1 0.5];
%   showcs3(D,scales);
%
% Function is written by D.Kroon University of Twente (July 2008)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showcs3_OpeningFcn, ...
                   'gui_OutputFcn',  @showcs3_OutputFcn, ...
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


% --- Executes just before showcs3 is made visible.
function showcs3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showcs3 (see VARARGIN)

% Choose default command line output for showcs3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

data.HandleWindow=gcf;
set(data.HandleWindow, 'Renderer', 'opengl')
hold on;  axis equal; axis xy;
if((~isempty(varargin))&&(numel(varargin{1})>8))
    data.voxelvolume=varargin{1};
    if((length(varargin)>=2)&&(length(varargin{2})==3))
        data.scales=varargin{2};
    else
        data.scales=[1 1 1];
    end
    data.sizes=size(data.voxelvolume);
    data.posx=0.5; data.posy=0.5;  data.posz=0.5; 
    set(handles.slider_x,'value',data.posx);
    set(handles.slider_y,'value',data.posy);
    set(handles.slider_z,'value',data.posz);
    daspect(data.scales)
    data=showslices(data);
    setMyData(data);
    rotate3d on;
else
    disp('Voxel volume not found');
end


% UIWAIT makes showcs3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = showcs3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_x_Callback(hObject, eventdata, handles)
% hObject    handle to slider_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.posx=get(handles.slider_x,'value');
data=showslices(data);
setMyData(data);
    
    

% --- Executes during object creation, after setting all properties.
function slider_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_y_Callback(hObject, eventdata, handles)
% hObject    handle to slider_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.posy=get(handles.slider_y,'value');
data=showslices(data);
setMyData(data);
    

% --- Executes during object creation, after setting all properties.
function slider_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_z_Callback(hObject, eventdata, handles)
% hObject    handle to slider_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.posz=get(handles.slider_z,'value');
data=showslices(data);
setMyData(data);
    
% --- Executes during object creation, after setting all properties.
function slider_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function data=showslices(data)
try delete(data.h1); delete(data.h2); delete(data.h3); catch end
posx=round(data.sizes(1)*data.posx);
posy=round(data.sizes(2)*data.posy);
posz=round(data.sizes(3)*data.posz);

if(ndims(data.voxelvolume)==3)
    slicexg = im2uint8(squeeze(data.voxelvolume(posx,:,:)))';
    sliceyg = im2uint8(squeeze(data.voxelvolume(:,posy,:)))';
    slicezg = im2uint8(squeeze(data.voxelvolume(:,:,posz)))';
    slicex=ind2rgb(slicexg,gray(256));
    slicey=ind2rgb(sliceyg,gray(256));
    slicez=ind2rgb(slicezg,gray(256));
else
    slicexs = im2uint8(squeeze(data.voxelvolume(posx,:,:,:))); 
    sliceys = im2uint8(squeeze(data.voxelvolume(:,posy,:,:))); 
    slicezs = im2uint8(squeeze(data.voxelvolume(:,:,posz,:))); 
    slicex=uint8(zeros([size(slicexs,2) size(slicexs,1) 3]));
    slicey=uint8(zeros([size(sliceys,2) size(sliceys,1) 3]));
    slicez=uint8(zeros([size(slicezs,2) size(slicezs,1) 3]));
    for i=1:3,
        slicex(:,:,i)=slicexs(:,:,i)'; slicey(:,:,i)=sliceys(:,:,i)';  slicez(:,:,i)=slicezs(:,:,i)';
    end
end

slicex_x=[posx posx;posx posx]; slicex_y=[0 (data.sizes(2)-1);0 (data.sizes(2)-1)]; slicex_z=[0 0;(data.sizes(3)-1) (data.sizes(3)-1)];
slicey_x=[0 (data.sizes(1)-1);0 (data.sizes(1)-1)]; slicey_y=[posy posy;posy posy]; slicey_z=[0 0;(data.sizes(3)-1) (data.sizes(3)-1)];
slicez_x=[0 (data.sizes(1)-1);0 (data.sizes(1)-1)]; slicez_y=[0 0;(data.sizes(2)-1) (data.sizes(2)-1)]; slicez_z=[posz posz;posz posz];

data.h1=surface(slicex_x,slicex_y,slicex_z, slicex,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct','FaceAlpha',1);
data.h2=surface(slicey_x,slicey_y,slicey_z, slicey,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct','FaceAlpha',1);
data.h3=surface(slicez_x,slicez_y,slicez_z, slicez,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct','FaceAlpha',1);

function data=getMyData()
data=getappdata(gcf,'showcsdata');

function setMyData(data)
setappdata(data.HandleWindow,'showcsdata',data);

