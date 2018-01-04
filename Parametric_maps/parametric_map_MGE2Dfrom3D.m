function MGE2Dfrom3D = parametric_map_MGE2Dfrom3D(MGE3D_map_filename, add_parameters)
%% (xdata,ydata,n,seuil_du_fit)
% generate a MGE2D map from a MGE3D scan by addition the signal of several slices


% load the MGE3D file
fid=fopen(MGE3D_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(MGE3D_map_filename);
    reco = data.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the MGE2Dfrom3D map because there is\n##$ Somthing wrong with the data\n##$MGE3D=%s\n##$',...
        MGE3D_map_filename);
    msgbox(warning_text, 'MGE2Dfrom3D map warning') ;
    MGE2Dfrom3D = [];
    return
end
MGE2Dfrom3D = data.uvascim.image;

% Empty memory
clear data

fist_slice =str2double(add_parameters{:}(3));
nbr_of_slice = str2double(add_parameters{:}(4));
resolution = str2double(add_parameters{:}(5));
echo_to_trash = add_parameters{:}(8);

% Some echoes need to be remove --> human data
if ~strcmp(echo_to_trash, 'end')
    reco.data= reco.data(:,:,1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.echo_label = reco.echo_label(1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.echotime = reco.echotime(1:str2double(echo_to_trash));
    reco.fov_offsets = reco.fov_offsets(:,1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.fov_orientation = reco.fov_orientation(:,1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.fov_phase_orientation = reco.fov_phase_orientation(1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.label = reco.label(1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.no_echoes = str2double(echo_to_trash);
    MGE2Dfrom3D.reco.phaselabel = reco.phaselabel(1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.scaling_factor = reco.scaling_factor(1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.scaling_offset = reco.scaling_offset(1:str2double(echo_to_trash),:);
    MGE2Dfrom3D.reco.unit = reco.unit(1:str2double(echo_to_trash),1:str2double(echo_to_trash));
end


[rows3d, cols3d, echos3d, depths3d]=size(reco.data);
% nbr_of_final_slice = round((depths3d-fist_slice+1)/nbr_of_slice)-1;
nbr_of_final_slice = fix((depths3d-fist_slice+1)/nbr_of_slice);
if nbr_of_final_slice ==0
    return
end
temp2  =zeros([size(reco.data, 1) size(reco.data, 2) size(reco.data, 3) nbr_of_final_slice]);
for i = 1:nbr_of_final_slice
    for k=1:echos3d
        temp=(reco.data(:,:,k,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
        temp=sum(abs(temp),4);
        temp2(:,:,k,i)=temp;
    end
end

xfactor = rows3d/resolution;
yfactor = cols3d/resolution;

x = 1:xfactor:reco.no_samples;
y = 1:yfactor:reco.no_views;
MGE2Dfrom3D.reco.data = NaN(size(reco.data((x+xfactor-1)/2,(y+yfactor-1)/2,:,1:nbr_of_final_slice)));

for i = 1:nbr_of_final_slice
    MGE2Dfrom3D.reco.data((x+xfactor-1)/2,(y+yfactor-1)/2,:,i)=(temp2(x,y,:,i)+...
        (temp2(x+xfactor-1,y,:,i)+...
        temp2(x,y+yfactor-1,:,i))+...
         temp2(x+xfactor-1,y+yfactor-1,:,i));
end
if length(MGE2Dfrom3D.reco.fov) == 3
    MGE2Dfrom3D.reco.fov(3)=[];
end
MGE2Dfrom3D.reco.fov_offsets(:,:,nbr_of_final_slice+1:end)=[];
for i = 1:nbr_of_final_slice
    temp=squeeze(reco.fov_offsets(:,1,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
    new_offset=median(temp,2);
    MGE2Dfrom3D.reco.fov_offsets(:,:,i) = repmat(new_offset, [1 MGE2Dfrom3D.reco.no_echoes]);
end
for i = 1:nbr_of_final_slice
    temp=squeeze(reco.fov_orientation(:,1,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
    new_offset=median(temp,2);
    MGE2Dfrom3D.reco.fov_orientation(:,:,i) = repmat(new_offset, [1 MGE2Dfrom3D.reco.no_echoes]);
end
for i = 1:nbr_of_final_slice
    temp=squeeze(reco.fov_phase_orientation(:,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
    MGE2Dfrom3D.reco.fov_phase_orientation(:,i)=median(temp,2);
end
 

MGE2Dfrom3D.reco.globalmax = max(MGE2Dfrom3D.reco.data(:));
MGE2Dfrom3D.reco.globalmin = min(MGE2Dfrom3D.reco.data(:));
MGE2Dfrom3D.reco.iminfos = 'MGE2Dfrom3D';
MGE2Dfrom3D.reco.no_samples = resolution;
MGE2Dfrom3D.reco.no_views = resolution;
MGE2Dfrom3D.reco.phaselabel(:,nbr_of_final_slice+1:end)=[];
MGE2Dfrom3D.reco.reco_meth = 'MGE2Dfrom3D';
for i = 1:nbr_of_final_slice
    temp=squeeze(reco.scaling_factor(:,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
    MGE2Dfrom3D.reco.scaling_factor(:,i)=median(temp,2);
end
for i = 1:nbr_of_final_slice
    temp=squeeze(reco.scaling_offset(:,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
    MGE2Dfrom3D.reco.scaling_offset(:,i)=median(temp,2);
end

MGE2Dfrom3D.reco.thickness = MGE2Dfrom3D.reco.thickness*nbr_of_slice;
MGE2Dfrom3D.reco.unit(:,nbr_of_final_slice+1:end)=[];
MGE2Dfrom3D.reco.no_slices = nbr_of_final_slice;

MGE2Dfrom3D.acq.thickness = MGE2Dfrom3D.reco.thickness;
MGE2Dfrom3D.acq.pix_spacing = [MGE2Dfrom3D.reco.fov(1)/MGE2Dfrom3D.reco.no_samples MGE2Dfrom3D.reco.fov(2)/MGE2Dfrom3D.reco.no_views]';
MGE2Dfrom3D.acq.fov_orientation =   median(squeeze(MGE2Dfrom3D.reco.fov_orientation(1:3,1,:)),2)';
MGE2Dfrom3D.acq.fov_offsets =  median(squeeze(MGE2Dfrom3D.reco.fov_offsets(1:3,1,:)),2)';


ParamConfig=sprintf('##$QuantifMethod=Merge slices\n##$First slice used=%s\n##$Number of slice merged=%s\n##$New resolution=%s\n##$Raw scan used=%s\n',...
    add_parameters{:}{3},...
    add_parameters{:}{4},...
    add_parameters{:}{5},...
    MGE3D_map_filename);
MGE2Dfrom3D.reco.paramQuantif = ParamConfig;
MGE2Dfrom3D.reco=orderfields(MGE2Dfrom3D.reco);

MGE2Dfrom3D.clip = [MGE2Dfrom3D.reco.globalmin MGE2Dfrom3D.reco.globalmax 1];
