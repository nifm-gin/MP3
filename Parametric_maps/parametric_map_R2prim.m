function R2prim_map = parametric_map_R2prim(T2star_filename, T2map_filename, add_parameters)
% generate a R2prim.reco from a T2* and a T2 maps

final_res = add_parameters{:}(1);
if strcmp(final_res, 'Original')
   rescale = 0;
else
    rescale = 1;
    final_res = str2double(add_parameters{:}(1));
end

fid=fopen(T2star_filename ,'r');
if fid>0
    fclose(fid);
    data_T2star = load(T2star_filename);
    RAW_T2star = data_T2star.uvascim.image.reco;
else
   warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ Somthing wrong with the data\n##$T2starcorr3D=%s\n##$',...
        T2star_filename);
    msgbox(warning_text, 'R2prim map warning') ;
    R2prim_map = [];
    return
end
fid=fopen(T2map_filename ,'r');
if fid>0
    fclose(fid);
    data_T2map= load(T2map_filename);
    RAW_T2map = data_T2map.uvascim.image.reco;
else
     warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ Somthing wrong with the data\n##$T2map=%s\n##$',...
        T2map_filename);
    msgbox(warning_text, 'R2prim map warning') ;
    R2prim_map = [];
    return
end

% check data compatibility (slice thickness and slice number)
if RAW_T2star.thickness ~= RAW_T2map.thickness
    warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ a slice thickness missmatch between\n##$T2*map=%s\n##$ and \n##$Tmap=%s',...
        T2star_filename,T2map_filename);
    msgbox(warning_text, 'R2prim map warning') ;
    R2prim_map = [];
    return
end
R2prim_map = data_T2star;
if RAW_T2star.no_slices ~= RAW_T2map.no_slices
    R2prim_slice_nbr = 0;
    R2prim_map.reco.fov_offsets = [];
    R2prim_map.reco.fov_orientation = [];
    R2prim_map.reco.label = {};
    R2prim_map.reco.phaselabel = {};
    R2prim_map.reco.fov_phase_orientation = [];
    for i = 1:size(RAW_T2star.data, 4)
        for j = 1:size(RAW_T2map.data, 4)
            if abs(RAW_T2star.fov_offsets(3,1,i) - RAW_T2map.fov_offsets(3,1,j)) < 1e-5
                R2prim_slice_nbr = R2prim_slice_nbr+1;
                % Compute the CMRO2 map each slice with the same offset
                R2prim_map.reco.data(:,:,1,R2prim_slice_nbr)=1./(RAW_T2star.data(:,:,1,i)*10^-3)-1./(RAW_T2map.data(:,:,1,j)*10^-3);
                
                % Update the SO2map structure 
                R2prim_map.reco.fov_offsets(:,1,R2prim_slice_nbr,1) = RAW_T2star.fov_offsets(:,1,i,1);
                R2prim_map.reco.fov_orientation(:,1,R2prim_slice_nbr,1) = RAW_T2star.fov_orientation(:,1,i,1);
                R2prim_map.reco.label(1,R2prim_slice_nbr,1) = RAW_T2star.label(1,i,1);
                R2prim_map.reco.phaselabel(1,R2prim_slice_nbr,1) = RAW_T2star.phaselabel(1,i,1);
                R2prim_map.reco.fov_phase_orientation(1,R2prim_slice_nbr,1) = RAW_T2star.fov_phase_orientation(1,i,1);
                
            end
        end
    end
    if R2prim_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ no slice offset match between\n##$T2*map=%s\n##$ and \n##$T2map=%s',...
            T2star_filename,T2map_filename);
        msgbox(warning_text, 'R2prim map warning') ;
        return
    end
    R2prim_map.reco.no_slices=R2prim_slice_nbr;
else
    R2prim_map.reco='';
    if rescale == 1
        if size(RAW_T2star.data,1) ~= final_res || size(RAW_T2star.data,2) ~= final_res
            for i = 1:RAW_T2star.no_slices
            tmp(:,:,1,i)= imresize(RAW_T2star.data(:,:,1,i),[final_res final_res],'bilinear');
            end
            RAW_T2star.data = tmp;
        end
        if size(RAW_T2map.data,1) ~= final_res || size(RAW_T2map.data,2) ~= final_res
             for i = 1:RAW_T2map.no_slices
            tmp(:,:,1,i)= imresize(RAW_T2map.data(:,:,1,i),[final_res final_res],'bilinear');
            end
            RAW_T2map.data = tmp;
        end
        clear tmp
        
    end
    % Compute the R2prim map for all slices
    R2prim_map.reco.data=1./(RAW_T2star.data*10^-3)-1./(RAW_T2map.data*10^-3);
    R2prim_map.reco.no_slices=RAW_T2star.no_slices;
    %Adapt the fov offsets and orientations infos
    R2prim_map.reco.fov_offsets=RAW_T2star.fov_offsets(:,1,:,1);
    R2prim_map.reco.fov_orientation=RAW_T2star.fov_orientation(:,1,:,1);
    R2prim_map.reco.label=RAW_T2star.label(1,:,1);
    R2prim_map.reco.phaselabel=RAW_T2star.phaselabel(1,:,1);
    R2prim_map.reco.fov_phase_orientation=RAW_T2star.fov_phase_orientation;
    
end


%Set the dimensions of the proc structure
if rescale == 0
    R2prim_map.reco.no_samples=RAW_T2star.no_samples;
    R2prim_map.reco.no_views=RAW_T2star.no_views;
else
    R2prim_map.reco.no_samples = final_res;
    R2prim_map.no_views = final_res;
end

%replace negative value of R2prim by the interpolation of their neibourgh 
R2prim_map_m = convn(R2prim_map.reco.data,ones(3,3,3)./26,'same');
R2prim_map.reco.data(R2prim_map.reco.data(:)<=0) =R2prim_map_m(R2prim_map.reco.data(:)<=0);


R2prim_map.reco.no_echoes=1; %Number of parameters stored here
R2prim_map.reco.no_expts=1;

%Complete the reco structure
R2prim_map.reco.texte='R2prim';
R2prim_map.reco.unit='s-1';
R2prim_map.reco.date=date;
R2prim_map.reco.displayedecho=1;
R2prim_map.reco.displayedslice=1;
R2prim_map.reco.displayedexpt=1;
R2prim_map.reco_number = 1;
R2prim_map.scan_number = 104;
R2prim_map.clip=[0 50 1];
R2prim_map.reco.globalmax = max(R2prim_map.reco.data(:));
R2prim_map.reco.globalmin = min(R2prim_map.reco.data(:));
R2prim_map.reco.thickness = RAW_T2star.thickness;
if isfield(RAW_T2map, 'paramQuantif') & isfield(RAW_T2map, 'paramQuantif')
ParamConfig=sprintf('##$QuantifMethod=''1./(T2star*10^-3)-1./(T2map*10^-3)''\n##$T2*map=%s\n##$T2map=%s\n\n##$Sequence T2*map info\n%s\n\n##$Sequence T2map info\n%s',...
   T2star_filename, T2map_filename, RAW_T2map.paramQuantif, RAW_T2map.paramQuantif);
else
 ParamConfig=sprintf('##$QuantifMethod=''1./(T2star*10^-3)-1./(T2map*10^-3)''\n##$T2*map=%s\n##$T2map=%s\n',...
   T2star_filename, T2map_filename);   
end
R2prim_map.reco.paramQuantif = ParamConfig;
R2prim_map.reco=orderfields(R2prim_map.reco);

%complete the structure
R2prim_map.acq = data_T2star.uvascim.image.acq;
R2prim_map.filename = data_T2star.uvascim.image.filename;


