function varargout = deconv_perf(varargin)
% DECONV_PERF MATLAB code for deconv_perf.fig
%      DECONV_PERF, by itself, creates a new DECONV_PERF or raises the existing
%      singleton*.
%
%      H = DECONV_PERF returns the handle to a new DECONV_PERF or the handle to
%      the existing singleton*.
%
%      DECONV_PERF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECONV_PERF.M with the given input arguments.
%
%      DECONV_PERF('Property','Value',...) creates a new DECONV_PERF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deconv_perf_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deconv_perf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deconv_perf

% Last Modified by GUIDE v2.5 07-Jun-2013 14:21:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deconv_perf_OpeningFcn, ...
                   'gui_OutputFcn',  @deconv_perf_OutputFcn, ...
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
end

% --- Executes just before deconv_perf is made visible.
function deconv_perf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to deconv_perf (see VARARGIN)
 
handles.inputdata.filename_DSC = varargin{3};
add_parameters = varargin{4};
handles.inputdata.wmaxr = str2double(add_parameters{:}(1));
handles.inputdata.nb_vox_cand = str2double(add_parameters{:}(2));
handles.inputdata.nb_vox = str2double(add_parameters{:}(3));
handles.clips_list = varargin{5};
handles.map_name = varargin{6};

% load the DSC file
fid=fopen(handles.inputdata.filename_DSC ,'r');
if fid>0
    fclose(fid);
    DSCstruct = load(handles.inputdata.filename_DSC);
    handles.data = DSCstruct.uvascim.image;
    handles.data.reco.data= double(handles.data.reco.data);
else
    warning_text = sprintf('##$ Can not calculate the DSC map because there is\n##$ Somthing wrong with the data\n##$Perf scan=%s\n##$',...
        handles.inputdata.filename_DSC);
    msgbox(warning_text, 'CBV map warning') ;
    CBV_map = [];
    return
end

% Empty memory


if handles.data.reco.no_expts == 0
    warning_text = sprintf('##$ Can not calculate the DSC map because there is\n##$ Somthing wrong with the data\n##$Perf scan=%s\n##$',...
        handles.inputdata.filename_DSC);
    return 
end

if strcmp(add_parameters{:}(4), 'Yes')
    [PATHSTR, NAME_scan_to_realign, ~] = fileparts(handles.inputdata.filename_DSC);
    NAME_scan_to_realign = [NAME_scan_to_realign '-for_spm'];
    
    
    NHdr=CreateNiiHdr_from_uvascim(handles.inputdata.filename_DSC, DSCstruct, NAME_scan_to_realign);
    data_4d = squeeze(DSCstruct.uvascim.image.reco.data(:,:,1,:,:));
    WriteNii(fullfile(PATHSTR, [NAME_scan_to_realign '.nii']),NHdr, data_4d)
    
    % do not realign images (dynamics) during the bolus
    mean_signal=squeeze(mean(mean(mean(data_4d))));
    mean_beg = mean(mean_signal(3:8)); sd_beg = std(mean_signal(3:8));
    mean_end = mean(mean_signal(end-5:end)); sd_end = std(mean_signal(end-5:end));
    
    allscans = 1:size(data_4d,4);
%     peak =  find(mean_signal==min(mean_signal));
%     % detect the bolus 
%     % 1) scan before the peak with a signal < the mean_beg-2*sd_beg
%     % 2) scan after the peak with a signal > the mean_end-12*sd_end
%     ignorescans = find(mean_signal(1:peak) < mean_beg-2*sd_beg, 1 ) : (size(data_4d,4) -(size(mean_signal(peak:end),1) - find(mean_signal(peak:end) > mean_end-12*sd_end, 1 )));
%     keepscans = setdiff(allscans,ignorescans);


    
    for xx = 1:numel(allscans)
        if ~exist('Scan_to_realign_nii_filename', 'var')
            Scan_to_realign_nii_filename{1,:} =  fullfile(PATHSTR,[NAME_scan_to_realign '.nii,1']);
        else
            Scan_to_realign_nii_filename{size(Scan_to_realign_nii_filename,1)+1,:} =  fullfile(PATHSTR, [NAME_scan_to_realign '.nii,' num2str(allscans(xx))]); %#ok<AGROW>
        end
    end
     matlabbatch{1}.spm.spatial.realign.estwrite.data  = {Scan_to_realign_nii_filename};
    %%
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    jobs = repmat(matlabbatch, 1, 1);
    inputs = cell(0, 1);
    for crun = 1:1
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
    clear matlabbatch 
     
    date_str = date;
    if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
        movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
            [PATHSTR, filesep, NAME_scan_to_realign(1:end-8) '_SPM_realign.ps']);
    end
 
    % Load resliced data
    V_RealignandReslice =spm_vol(fullfile(PATHSTR, ['r', NAME_scan_to_realign, '.nii']));
    Data_RealignandReslice = spm_read_vols(V_RealignandReslice);
    Data_RealignandReslice = permute(Data_RealignandReslice, [1 2 5 3 4]);
    Data_RealignandReslice = permute(Data_RealignandReslice, [2 1 3 4 5]);
    Data_RealignandReslice = flip(Data_RealignandReslice, 1);
    Data_RealignandReslice = flip(Data_RealignandReslice, 2);
     handles.data.reco.data = Data_RealignandReslice;
%     for x = 1:numel(keepscans)
%         DSCstruct.uvascim.image.reco.data(:,:,1,:,keepscans(x)) = Data_RealignandReslice(:,:,1,:,keepscans(x));
%     end
%     handles.data = DSCstruct.uvascim.image;
end
% set figure
expt = 1;
nb_expts =handles.data.reco.no_expts;
set(handles.expts_disp,'String',sprintf('Slice : %02d / %02d',expt,nb_expts));
set(handles.expts_slider,'SliderStep',[1/nb_expts 2/nb_expts]);
set(handles.expts_slider,'Min', 1);
set(handles.expts_slider,'Max', nb_expts);

sli = 1;
nb_sli = handles.data.reco.no_slices;
if nb_sli > 1
    set(handles.slices_disp,'String',sprintf('Slice : %02d / %02d',sli,nb_sli));
    set(handles.slice_slider,'SliderStep',[1/nb_sli 2/nb_sli]);
    set(handles.slice_slider,'Min', 1);
    set(handles.slice_slider,'Max', nb_sli);
else
    set(handles.slice_slider,'Visible','Off');
    set(handles.slices_disp,'Visible','Off');
end

set(handles.Image_used_min, 'String', '1');
set(handles.Image_used_max, 'String', num2str(nb_expts));


% Choose default command line output for deconv_perf
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes deconv_perf wait for user response (see UIRESUME)
% uiwait(handles.deconv_perf);

if nb_sli > 1
    addlistener(handles.slice_slider,'Value','PostSet',@(hObject,eventdata)cont_slice_slider(hObject,eventdata,handles));
end
addlistener(handles.expts_slider,'Value','PostSet',@(hObject,eventdata)cont_expts_slider(hObject,eventdata,handles));

display_image(hObject,eventdata,handles)
end

% --- Outputs from this function are returned to the command line.
function varargout = deconv_perf_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout = {handles.aif handles.score};
guidata(handles.deconv_perf, handles)
end

% --- Executes on slider movement.
function slice_slider_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function slice_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function expts_slider_Callback(hObject, eventdata, handles)
% hObject    handle to expts_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function expts_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expts_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global output
output= {[], []};
delete(handles.deconv_perf)
end

% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
% hObject    handle to compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global output
score =handles.score;
data_table = get(handles.list_vox,'Data');
num=0;
for i=1:size(data_table,1)
    stri = num2str(i);
    if data_table{i,5} == 0
        score(i-num,:) =[];
        num = num+1;
    end
end
output= {handles.aif score get(handles.Image_used_min, 'String') get(handles.Image_used_max, 'String')};
delete(handles.deconv_perf)
end

% --- Executes during object creation, after setting all properties.
function deconv_perf_CreateFcn(hObject, eventdata, handles, varargin)
% hObject    handle to deconv_perf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% load('uvascim_Pat4.mat')
% 
colormap(gray(1024))
end

function cont_slice_slider(hObject,eventdata,handles)


n = 1;%get(findobj('Tag','uvasc_list_box'),'Value');
sli = round(get(handles.slice_slider,'Value'));
nb_sli = handles.data.reco.no_slices;
set(handles.slices_disp,'String',sprintf('Slice : %02d / %02d',sli,nb_sli));
display_image(hObject,eventdata,handles);
end

function cont_expts_slider(hObject,eventdata,handles)


n = 1;%get(findobj('Tag','uvasc_list_box'),'Value');
expt = round(get(handles.expts_slider,'Value'));
nb_expt = handles.data.reco.no_expts;
set(handles.expts_disp,'String',sprintf('Expt : %02d / %02d',expt,nb_expt));
display_image_expts(hObject,eventdata,handles);
end

function display_image(hObject, eventdata, handles)

% handles = guidata(handles.deconv_perf);

nb_sli = handles.data.reco.no_slices;
if nb_sli > 1
    sli = round(get(handles.slice_slider,'Value'));
else
    sli = 1;
end
expt = round(get(handles.expts_slider,'Value'));

%%% Find voxels for AIF
[aif,scores] =extraction_aif(squeeze(handles.data.reco.data(:,:,1,sli,:)), true(size(squeeze(handles.data.reco.data(:,:,1,sli,expt)))), handles.inputdata);

handles.score = scores;
handles.aif = aif;
guidata(handles.deconv_perf,handles);

colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
dat = cell(size(scores,1),5);
bin = false(size(scores,1),size(handles.data.reco.data(:,:,1,sli,expt),1),size(handles.data.reco.data(:,:,1,sli,expt),2));
tablecolor = {'#FF0000';'#00FF00';'#0000FF';'#FFFF00';'#00FFFF'};
old_data_table = get(handles.list_vox,'Data');
for i=1:size(scores,1)
    dat(i,:) = {colorgen(tablecolor{mod(i,5)+1},num2str(scores{i,3})) colorgen(tablecolor{mod(i,5)+1},num2str(scores{i,2})) scores{i,1} scores{i,4} true};
    bin(i,scores{i,2},scores{i,3}) = true;
end

set(handles.list_vox,'Data',dat);

col_plot = [1 0 0;0 1 0;0 0 1;1 1 0;0 1 1];
plot(handles.vox_curves,aif,'LineWidth',3,'Color','black')
hold(handles.vox_curves,'on')
for v=1:size(scores,1)
    plot(handles.vox_curves,squeeze(handles.data.reco.data(scores{v,2},scores{v,3},1,sli,:)),'Color',col_plot(mod(v,5)+1,1:3),'LineWidth',1.5, 'Tag', ['AIF' num2str(v)])
end
hold(handles.vox_curves,'off')
set(handles.vox_curves,'Xtick',[]);
set(handles.vox_curves,'Ytick',[]);
line_dynamics(hObject,eventdata,handles);
set(handles.vox_curves,'ButtonDownFcn',@vox_curves_ButtonDownFcn);
handles.aff_vox = true(size(scores,1),1);
guidata(handles.deconv_perf,handles);

axes(handles.slice_aif)
voximg(1,:,:,1) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(1,:,:,2) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(1,:,:,3) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(2,:,:,1) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(2,:,:,2) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(2,:,:,3) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(3,:,:,1) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(3,:,:,2) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(3,:,:,3) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(4,:,:,1) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(4,:,:,2) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(4,:,:,3) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(5,:,:,1) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(5,:,:,2) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(5,:,:,3) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));


img = squeeze(handles.data.reco.data(:,:,1,sli,expt));
% if ~isfield(handles.data,'clip')
%     NumScanClip = find(strcmp(handles.clips_list(:,1), handles.map_name{1})==1);
%     clip_min = handles.clips_list{NumScanClip,2};
%     clip_max = handles.clips_list{NumScanClip,3};
%     handles.data.clip = [clip_min clip_max 1];
% end
if sum(isnan([min(min(img)) max(max(img)) ])) > 1
    imagesc(img,'Parent',handles.slice_aif);
else
imagesc(img,'Parent',handles.slice_aif, [min(min(img)) max(max(img)) ])%,imgadjustlim(img));
end
hold(handles.slice_aif,'on')
for i=1:size(scores,1)
    image(squeeze(voximg(mod(i,5)+1,:,:,:)),'AlphaData',squeeze(bin(i,:,:)),'Parent',handles.slice_aif);
end
handles.aif_only = false;
guidata(handles.deconv_perf,handles);

if isfield(handles,'ROI')
    roi_contour=contourc(double(handles.ROI),1);
else
    roi_contour = [];
end
hold(handles.slice_aif,'off')

set(handles.slice_aif,'Xtick',[]);
set(handles.slice_aif,'Ytick',[]);

end

function display_image_expts(hObject, eventdata, handles)

handles = guidata(handles.deconv_perf);
score = handles.score;

n = 1;%get(findobj('Tag','uvasc_list_box'),'Value');
sli = round(get(handles.slice_slider,'Value'));
expt = round(get(handles.expts_slider,'Value'));

bin = false(5,size(handles.data.reco.data(:,:,1,sli,expt),1),size(handles.data.reco.data(:,:,1,sli,expt),2));
for i=1:5
    bin(i,score{i,2},score{i,3}) = true;
end

axes(handles.slice_aif)
voximg(1,:,:,1) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(1,:,:,2) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(1,:,:,3) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(2,:,:,1) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(2,:,:,2) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(2,:,:,3) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(3,:,:,1) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(3,:,:,2) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(3,:,:,3) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(4,:,:,1) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(4,:,:,2) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(4,:,:,3) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(5,:,:,1) = zeros(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(5,:,:,2) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));
voximg(5,:,:,3) = ones(size(handles.data.reco.data(:,:,1,sli,expt)));


img = squeeze(handles.data.reco.data(:,:,1,sli,expt));
imagesc(img,'Parent',handles.slice_aif,[min(min(img)) max(max(img)) ]);%imgadjustlim(img);
hold(handles.slice_aif,'on')
for i=1:5
    image(squeeze(voximg(i,:,:,:)),'AlphaData',squeeze(bin(i,:,:)),'Parent',handles.slice_aif);
end

if isfield(handles,'ROI')
    roi_contour=contourc(double(handles.ROI),1);
else
    roi_contour = [];
end
            
if ~isempty(roi_contour)
    % JW plot lines in contour one by one
    curr_idx = 1;
    while curr_idx < size(roi_contour,2),
        npts = roi_contour(2,curr_idx);
        plot_idx = curr_idx + (1:npts);
        curr_idx = curr_idx + 1+npts;
        line(roi_contour(1,plot_idx),roi_contour(2,plot_idx),...
            'parent',handles.slice_aif,...
            'color',[0 1 0],...
            'LineStyle', '-',...
            'LineWidth', 2,...
            'tag','ROI_contour',...
            'visible','on');
    end;
end

hold(handles.slice_aif,'off')

set(handles.slice_aif,'Xtick',[]);
set(handles.slice_aif,'Ytick',[]);

line_dynamics(hObject,eventdata,handles);
end

% --- Executes during object creation, after setting all properties.
function slices_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slices_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end

% --- Executes during object creation, after setting all properties.
function expts_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expts_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


end




% --- Executes on mouse press over axes background.
function vox_curves_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to vox_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
aff_vox = handles.aff_vox;
col_curves = [1 0 0;0 1 0;0 0 1;1 1 0;0 1 1];
if handles.aif_only
    for v=1:5
        if aff_vox(v)
            set(findobj(handles.vox_curves,'Type','line','Color',col_curves(v,1:3)),'Visible','on');
        end
    end
    handles.aif_only = false;
else
    for v=1:5
        set(findobj(handles.vox_curves,'Type','line','-not','Color','black'),'Visible','off')
    end
    handles.aif_only = true;
end
guidata(hObject,handles);
line_dynamics(hObject,eventdata,handles);
refreshdata(handles.vox_curves);
end


% --- Executes when entered data in editable cell(s) in list_vox.
function list_vox_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to list_vox (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(handles.deconv_perf);

sli = round(get(handles.slice_slider,'Value'));

scores = handles.score;
tmpdata = squeeze(handles.data.reco.data(:,:,1,sli,:));

data_table = get(handles.list_vox,'Data');
handles.aff_vox  = [data_table{:,5}]';
handles.aif  = zeros(1,length(handles.aif));
for i = 1:size(data_table,1)
    if data_table{i,5}
        if ~handles.aif_only
            set(findobj('Tag', ['AIF' num2str(i)]),'Visible','on');
        end
        handles.aif  = handles.aif  + squeeze(tmpdata(scores{i,2},scores{i,3},:)).';
    else
        set(findobj('Tag', ['AIF' num2str(i)]),'Visible','off');
    end
    
end
handles.aif  = handles.aif ./numel(find(handles.aff_vox ));
set(findobj(handles.vox_curves,'Type','line','Color','black'),'YData',handles.aif);

guidata(handles.deconv_perf,handles);
line_dynamics(hObject,eventdata,handles);
refreshdata(handles.vox_curves);

end



function Image_used_min_Callback(hObject, eventdata, handles)
% hObject    handle to Image_used_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Image_used_min as text
%        str2double(get(hObject,'String')) returns contents of Image_used_min as a double
end

% --- Executes during object creation, after setting all properties.
function Image_used_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Image_used_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Image_used_max_Callback(hObject, eventdata, handles)
% hObject    handle to Image_used_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Image_used_max as text
%        str2double(get(hObject,'String')) returns contents of Image_used_max as a double
end

% --- Executes during object creation, after setting all properties.
function Image_used_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Image_used_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function line_dynamics(hObject,eventdata,handles)

expt = round(get(handles.expts_slider,'Value'));
y_lims = get(handles.vox_curves,'YLim');
child_vox_curves = get(handles.vox_curves,'Children');
vert_dyn_line = findall(child_vox_curves,'Tag','vert_line');
if ~isempty(vert_dyn_line)
    set(vert_dyn_line,'XData',[expt expt]);
    set(vert_dyn_line,'YData',y_lims);
else
    axes(handles.vox_curves);
    vert_dyn_line = line([expt expt],y_lims,'Color','black','LineWidth',1);
    set(vert_dyn_line,'Tag','vert_line');
end
end
