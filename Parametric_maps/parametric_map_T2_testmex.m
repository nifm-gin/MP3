%function [varargout{1}, varargout{2}] = parametric_map_T2_testmex(MSE_map_filename, Mask_filename, add_parameters)
function [varargout] = parametric_map_T2_testmex(varargin)
%% (xdata,ydata,n,seuil_du_fit)
% generate a T2_mapfrom a RAW MSE scan
% this code come from the SO2_fig function
% additional_parameters correspond to the threshold used and the echo used

%#codegen
tic
if nargin == 3
    MSE_map_filename= varargin{1};
    Mask_filename= varargin{2};
    add_parameters= varargin{3};
    seuil_du_fit=str2double(add_parameters(1));
    trash_below=str2double(add_parameters(2));
    trash_after = str2double(add_parameters(3));
    remove_last_echo = add_parameters(5) ;
end
% load the MSME file
fid=fopen(MSE_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(MSE_map_filename);
    reco = data.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the T2 map because there is\n##$ Somthing wrong with the data\n##$MSE=%s\n##$',...
        MSE_map_filename);
    msgbox(warning_text, 'T2 map warning') ;
    varargout{1} = [];
    varargout{2} = [];
    return
end
% Load mask (if selected and/or exist)
if ~isempty(Mask_filename)
    fid=fopen(Mask_filename ,'r');
    if fid>0
        fclose(fid);
        Mask = load(Mask_filename);
        Mask = Mask.uvascroi;
    else
        Mask = []; 
    end
else
    Mask = [];
end


% save imformation
varargout{1}.acq=data.uvascim.image.acq;
varargout{1}.filename=data.uvascim.image.filename;

% Empty memory
clear data

%  human data
if strcmp(remove_last_echo, 'Yes')
    reco.data= reco.data(:,:,1:end-1,:);
    reco.echo_label = reco.echo_label(1:end-1,:);
    reco.echotime = reco.echotime(1:end-1);
    reco.fov_offsets = reco.fov_offsets(:,1:end-1,1:end-1);
    reco.fov_orientation = reco.fov_orientation(:,1:end-1,1:end-1);
    reco.fov_phase_orientation = reco.fov_phase_orientation(1:end-1,:);
    reco.label = reco.label(1:end-1,:);
    reco.no_echoes = reco.no_echoes-1;
    reco.phaselabel = reco.phaselabel(1:end-1,:);
    reco.scaling_factor = reco.scaling_factor(1:end-1,:);
    reco.scaling_offset = reco.scaling_offset(1:end-1,:);
    reco.unit = reco.unit(1:end-1,1:end-1);
end

% initializing variables

echotime=reco.echotime(:);
% Selectet echoes used
if trash_after == Inf
    lastecho = length(echotime);
else
    lastecho = sum(echotime < trash_after);
end

firstecho = sum(echotime < trash_below) + 1;
varargout{1}.reco='';
varargout{1}.reco.data=NaN(reco.no_samples,reco.no_views,1,reco.no_slices);     %
data_MSME=squeeze(reco.data(:,:,:,:));
maxim=max(data_MSME(:)) * seuil_du_fit / 100;
echotimes_used = echotime(firstecho:lastecho)';

tmp_data = reco.data(:,:,firstecho:lastecho,:);

if ~isempty(Mask)
    [~, ~, E, Z] = size(tmp_data);
    for i = 1:Z
        for j = 1:numel(Mask)
            if abs(reco.fov_offsets(3,1,i) - Mask(j).fov_offsets(3)) < 1e-5
                 tmp_data(:,:,:,i) = tmp_data(:,:,:,i).*repmat(Mask(j).value,[1 1 E 1]);          
            end
        end
    end
end

tmp_data = permute(tmp_data, [1 2 4 3]);
data_in_vector = reshape(tmp_data, [size(reco.data,1)*size(reco.data,2)*size(reco.data,4), size(tmp_data,4)]);
fit_result_a = NaN(size(data_in_vector,1),1);
voxel_to_fit=findn(sum(data_in_vector,2) >1);

for voxel_nbr = 1:size(voxel_to_fit,1)
    tmp_voxel_data=data_in_vector(voxel_to_fit(voxel_nbr),:);
    if max(tmp_voxel_data(:))>= maxim
        [t2,~,~]=fit_exp(echotimes_used,tmp_voxel_data,1);
        fit_result_a(voxel_to_fit(voxel_nbr))=t2;
    end
end


% Complete the T2star_map.reco structure
varargout{1}.reco.no_samples=reco.no_samples;
varargout{1}.reco.no_views=reco.no_views;

varargout{1}.reco.no_echoes=1;
varargout{1}.reco.no_slices=reco.no_slices;
varargout{1}.reco.no_expts=1;

%Adapt the fov offsets and orientation infos
tmpfov=reco.fov_offsets(:,1,:,1);
tmpori=reco.fov_orientation(:,1,:,1);
for m_slice=1:reco.no_slices,
    tmp=squeeze(tmpfov(:,m_slice));
    varargout{1}.reco.fov_offsets(:,:,m_slice)=repmat(tmp,1,varargout{1}.reco.no_echoes);
    tmp=squeeze(tmpori(:,m_slice));
    varargout{1}.reco.fov_orientation(:,:,m_slice)=repmat(tmp,1,varargout{1}.reco.no_echoes);
    varargout{1}.reco.label=repmat(reco.label(1,m_slice,1),1,varargout{1}.reco.no_echoes);
    varargout{1}.reco.phaselabel=repmat(reco.phaselabel(1,m_slice,1),1,varargout{1}.reco.no_echoes);
end
varargout{1}.reco.fov_phase_orientation=reco.fov_phase_orientation;
% T2map.reco.fitmethod='levenbergmarquardt';
% T2map.reco.seuil_du_fit=seuil_du_fit;
varargout{1}.reco.mask_to_use='';
varargout{1}.reco.texte='T2map';
varargout{1}.reco.unit='ms';
varargout{1}.reco.date=date;
varargout{1}.reco.thickness = reco.thickness;
varargout{1}.reco.displayedecho=1;
varargout{1}.reco.displayedslice=1;
varargout{1}.reco.displayedexpt=1;
varargout{1}.reco.echo_label ='';
% ParamConfig=sprintf('##$QuantifMethod=levenbergmarquardt\n##$Fit threshold=%s\n#$Echo trash below (ms)=%s\n##$Echo trash after (ms)=%s\n##$Raw scan used=%s\n##$Method used=%s\n',...
%     add_parameters{:}{1},...
%     add_parameters{:}{2},...
%     add_parameters{:}{3},...
%     add_parameters{:}{4},...
%     MSE_map_filename);
% varargout{1}.reco.paramQuantif = ParamConfig;

tmp_a=reshape(fit_result_a,[size(reco.data,1),size(reco.data,2),size(reco.data,4)]);
varargout{1}.reco.data=permute(tmp_a, [1 2 4 3]);
varargout{1}.reco.globalmax=max(varargout{1}.reco.data(:));
varargout{1}.reco.globalmin=min(varargout{1}.reco.data(:));

varargout{1}.scan_number = 104;
varargout{1}.clip = [0 200 1];
varargout{1}.reco_number = 1;
varargout{1}.reco=orderfields(varargout{1}.reco);
toc


