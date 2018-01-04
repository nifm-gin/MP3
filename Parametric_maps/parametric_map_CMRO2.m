function CMRO2 = parametric_map_CMRO2(filename_SO2map, filename_CBF,  add_parameters)
% generate a CMRO2 from a SO2map and a CBF scan
% this code come from the the ASL module coded by C. Debacker


% load the SO2map file
fid=fopen(filename_SO2map ,'r');
if fid>0
    fclose(fid);
    SO2map = load(filename_SO2map);
    SO2map = SO2map.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ Somthing wrong with the data\n##$SO2map=%s\n##$',...
        filename_SO2map);
    msgbox(warning_text, 'CMRO2 map warning') ;
    CMRO2 = [];
    return
end

% load the CBF file
fid=fopen(filename_CBF ,'r');
if fid>0
    fclose(fid);
    CBF = load(filename_CBF);
    CBF = CBF.uvascim.image;
    % If CBF map has multiple echoes (for excluded data) keep the CBF data
    % only
    if size(CBF.reco.data, 3) == 2
       CBF.reco.data= CBF.reco.data(:,:,1,:);
    end
else
    warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ Somthing wrong with the data\n##$CBFmap=%s\n##$',...
        filename_CBF);
    msgbox(warning_text, 'CMRO2 map warning') ;
    CMRO2 = [];
    return
end
CMRO2 = SO2map;
% data = CBF.reco.data.*(1/SO2map.reco.data).*Ca;

% formule Julien
% data = Ca*(CBF./100).*(OEF./100);

Ca = str2double(add_parameters{:}(1)); %14.7*1.39;  % 14.7g/dl d'Hb      1g Hb transporte 1.39mL d'O2
final_res = add_parameters{:}(2);
if strcmp(final_res, 'Original')
   rescale = 0;
else
    rescale = 1;
    final_res = str2double(add_parameters{:}(2));
end


% check data compatibility (slice thickness and slice number)
if CBF.reco.thickness ~= SO2map.reco.thickness
    warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ a slice thickness missmatch between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
        filename_SO2map,filename_CBF);
    msgbox(warning_text, 'CMRO2 map warning') ;
     CMRO2 = [];
    return
end
if CBF.reco.no_samples ~= SO2map.reco.no_samples && rescale == 0
     warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ a resolution missmatch between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
        filename_SO2map,filename_CBF);
    msgbox(warning_text, 'CMRO2 map warning') ;
     CMRO2 = [];
    return
end
if CBF.reco.no_slices ~= SO2map.reco.no_slices
    CMRO2_slice_nbr = 0;
    CMRO2.reco.fov_offsets = [];
    CMRO2.reco.fov_orientation = [];
    CMRO2.reco.label = {};
    CMRO2.reco.phaselabel = {};
    CMRO2.reco.fov_phase_orientation = [];
    for i = 1:size(CBF.reco.data, 4)
        for j = 1:size(SO2map.reco.data, 4)
            if abs(CBF.reco.fov_offsets(3,1,i) - SO2map.reco.fov_offsets(3,1,j)) < 1e-5
                CMRO2_slice_nbr = CMRO2_slice_nbr+1;
                % Compute the CMRO2 map each slice with the same offset
                %                 data(:,:,:,CMRO2_slice_nbr)=Ca.*(imresize(squeeze(CBF.reco.data(:,:,:,i)), 0.5)./100).*((100-SO2map.reco.data(:,:,:,j))./100);
                if rescale == 0
                     data(:,:,:,CMRO2_slice_nbr)=Ca.*(CBF.reco.data(:,:,:,i)./100).*((100-SO2map.reco.data(:,:,:,j))./100);
                else
                     data(:,:,:,CMRO2_slice_nbr)=Ca.*(imresize(CBF.reco.data(:,:,:,i),[final_res final_res],'bilinear')./100).*...
                         ((100-imresize(SO2map.reco.data(:,:,:,j),[final_res final_res],'bilinear'))./100);
                end
                % Update the SO2map structure
                CMRO2.reco.fov_offsets(:,1,CMRO2_slice_nbr,1) = SO2map.reco.fov_offsets(:,1,j,1);
                CMRO2.reco.fov_orientation(:,1,CMRO2_slice_nbr,1) = SO2map.reco.fov_orientation(:,1,j,1);
                CMRO2.reco.label(1,CMRO2_slice_nbr,1) = SO2map.reco.label(1,j,1);
                CMRO2.reco.phaselabel(1,CMRO2_slice_nbr,1) = SO2map.reco.phaselabel(1,j,1);
                CMRO2.reco.fov_phase_orientation(1,CMRO2_slice_nbr,1) = SO2map.reco.fov_phase_orientation(1,j,1);
                
            end
        end
    end
    if CMRO2_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ no slice offset match between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
            filename_SO2map,filename_CBF);
        msgbox(warning_text, 'CMRO2 map warning') ;
        return
    end
    CMRO2.reco.no_slices=CMRO2_slice_nbr;
else
    % Compute the CMRO2 map for all slices
    if rescale == 0
        data=Ca.*(CBF.reco.data./100).*((100-SO2map.reco.data)./100);
    else
        for i = 1:CBF.reco.no_slices
            data(:,:,:,i)=Ca.*(imresize(CBF.reco.data(:,:,1,i),[final_res final_res],'bilinear')./100).*...
                ((100-imresize(SO2map.reco.data(:,:,1,i),[final_res final_res],'bilinear'))./100);
        end
    end
    
end


% copy and update the SO2map structure
CMRO2.reco = rmfield(CMRO2.reco, 'data');
CMRO2.reco.data = data;
if rescale == 1
CMRO2.reco.no_samples = final_res;
CMRO2.reco.no_views = final_res;
end
CMRO2.reco.texte = 'CMRO2';
CMRO2.reco.globalmax = max(data(:));
CMRO2.reco.globalmin = min(data(:));
if isfield(SO2map.reco, 'paramQuantif')
    ParamConfig=sprintf('##$QuantifMethod=''Ca*(CBF/100).*((100-SO2map)/100)''\n##$final_res=%s\n##$Ca=%s\n##$SO2map=%s\n##$CBFmap=%s\n\n##$SO2map info\n%s\n##$CBFmap info\n%s',...
       add_parameters{:}{2}, add_parameters{:}{1}, filename_SO2map,filename_CBF, SO2map.reco.paramQuantif, CBF.reco.paramQuantif);
else
    ParamConfig=sprintf('##$QuantifMethod=''Ca*(CBF/100).*((100-SO2map)/100)''\n##$final_res=%s##$Ca=%d\s##$SO2map=%s\n##$CBFmap=%s',...
        add_parameters{:}{2},add_parameters{:}{1}, filename_SO2map,filename_CBF);
end
CMRO2.reco.paramQuantif = ParamConfig;
if isfield(CMRO2.reco, 'wrongpix')
CMRO2.reco = rmfield(CMRO2.reco, 'wrongpix');
end
if isfield(CMRO2.reco, 'err')
CMRO2.reco = rmfield(CMRO2.reco, 'err');
end
CMRO2.reco=orderfields(CMRO2.reco);
