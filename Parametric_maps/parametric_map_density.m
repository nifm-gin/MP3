function Density_map = parametric_map_density(deltaR2star_filename, deltaR2_filename, add_parameters)
% generate a Density map from a deltaR2star map and a deltaR2 map

final_res = add_parameters{:}(1);
ADCmap_filename = add_parameters{:}(2);

% load the deltaR2star file
fid=fopen(deltaR2star_filename ,'r');
if fid>0
    fclose(fid);
    deltaR2star = load(deltaR2star_filename);
    deltaR2star = deltaR2star.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ Somthing wrong with the data\n##$deltaR2*=%s\n##$',...
        deltaR2star_filename);
    msgbox(warning_text, 'VSI map warning') ;
    Density_map = [];
    return
end
% load the deltaR2 file
fid=fopen(deltaR2_filename ,'r');
if fid>0
    fclose(fid);
    deltaR2 = load(deltaR2_filename);
    deltaR2 = deltaR2.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ Somthing wrong with the data\n##$deltaR2=%s\n##$',...
        deltaR2_filename);
    msgbox(warning_text, 'Density map warning') ;
    Density_map = [];
    return
end

% load the ADC_map if needed
if ~strcmp(ADCmap_filename, 'None')
    if ~isempty(ADCmap_filename{:})
        fid=fopen(ADCmap_filename{:} ,'r');
        if fid>0
            fclose(fid);
            ADCmap = load(ADCmap_filename{:});
            ADCmap = ADCmap.uvascim.image;
        else
            warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ Somthing wrong with the data\n##$ADCmap=%s\n##$',...
                ADCmap_filename);
            msgbox(warning_text, 'Density map warning') ;
            Density_map = [];
            return
        end
    else
        warning_text = sprintf('##$ Can not calculate the Density map because \n##$ No ADCmap has been found\n');
        msgbox(warning_text, 'Density map warning') ;
        ADCmap = [];
        return
    end
end


if strcmp(final_res, 'Original')
   rescale = 0;
else
    rescale = 1;
    final_res = str2double(add_parameters{:}(1));
end




Density_map = deltaR2star;

% check data compatibility (slice thickness and slice number)
if deltaR2star.reco.thickness ~= deltaR2.reco.thickness
    warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ a slice thickness missmatch between\n##$deltaR2star map=%s\n##$ and \n##$deltaR2 map=%s',...
        deltaR2star_filename,deltaR2_filename);
    msgbox(warning_text, 'Density map warning') ;
     Density_map = [];
    return
end
if deltaR2star.reco.no_samples ~= deltaR2.reco.no_samples && rescale == 0
     warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ a resolution missmatch between\n##$deltaR2star map=%s\n##$ and \n##$deltaR2 map=%s',...
        deltaR2star_filename,deltaR2_filename);
    msgbox(warning_text, 'Density map warning') ;
     Density_map = [];
    return
end
if deltaR2star.reco.no_slices ~= deltaR2.reco.no_slices
    Density_map_slice_nbr = 0;
    Density_map.reco.fov_offsets = [];
    Density_map.reco.fov_orientation = [];
    Density_map.reco.label = {};
    Density_map.reco.phaselabel = {};
    Density_map.reco.fov_phase_orientation = [];
    for i = 1:size(deltaR2.reco.data, 4)
        for j = 1:size(deltaR2star.reco.data, 4)
            if abs(deltaR2.reco.fov_offsets(3,1,i) - deltaR2star.reco.fov_offsets(3,1,j)) < 1e-5
                Density_map_slice_nbr = Density_map_slice_nbr+1;
                if strcmp(ADCmap_filename, 'None')
                    %ref lemasson et al 2013
                    if rescale == 0
                        data(:,:,:,Density_map_slice_nbr)=329000.*((deltaR2.reco.data(:,:,:,i).^3)./(deltaR2star.reco.data(:,:,:,j).^2));
                    else
                        data(:,:,:,Density_map_slice_nbr)=329000.*((imresize(deltaR2.reco.data(:,:,1,i),[final_res final_res],'bilinear').^3)./...
                            (imresize(deltaR2star.reco.data(:,:,1,j),[final_res final_res],'bilinear').^2));
                    end
                else
                    %ref Boehm-Sturm et al 2012 Contrast media molecular imaging
                     if rescale == 0
                        data(:,:,:,Density_map_slice_nbr)= (218456000./ADCmap.reco.data(:,:,:,i)).*((deltaR2.reco.data(:,:,:,i).^3)./(deltaR2star.reco.data(:,:,:,j).^2));
                    else
                        data(:,:,:,Density_map_slice_nbr)=(218456000./imresize(ADCmap.reco.data(:,:,1,i),[final_res final_res],'bilinear')).*((imresize(deltaR2.reco.data(:,:,1,i),[final_res final_res],'bilinear').^3)./...
                            (imresize(deltaR2star.reco.data(:,:,1,j),[final_res final_res],'bilinear').^2));
                    end
                end
                
                % Update the SO2map structure
                Density_map.reco.fov_offsets(:,1,Density_map_slice_nbr,1) = deltaR2star.reco.fov_offsets(:,1,j,1);
                Density_map.reco.fov_orientation(:,1,Density_map_slice_nbr,1) = deltaR2star.reco.fov_orientation(:,1,j,1);
                Density_map.reco.label(1,Density_map_slice_nbr,1) = deltaR2star.reco.label(1,j,1);
                Density_map.reco.phaselabel(1,Density_map_slice_nbr,1) = deltaR2star.reco.phaselabel(1,j,1);
                Density_map.reco.fov_phase_orientation(1,Density_map_slice_nbr,1) = deltaR2star.reco.fov_phase_orientation(1,j,1);
                
            end
        end
    end
    if Density_map_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ no slice offset match between\n##$DeltaR2=%s\n##$ and \n##$DeltaR2*=%s',...
            filename_SO2map,filename_CBF);
        msgbox(warning_text, 'Density map warning') ;
        return
    end
    Density_map.reco.no_slices=Density_map_slice_nbr;
else
    % Compute theDensity_map for all slices
    if strcmp(ADCmap_filename, 'None')
        %ref lemasson et al 2013
        if rescale == 0
            data=329000.*((deltaR2.reco.data.^3)./(deltaR2star.reco.data.^2));
        else
            for i = 1:deltaR2star.reco.no_slices
                data(:,:,:,i)=329000.*((imresize(deltaR2.reco.data(:,:,1,i),[final_res final_res],'bilinear').^3)./...
                    (imresize(deltaR2star.reco.data(:,:,1,i),[final_res final_res],'bilinear').^2));
            end
        end
    else
        %ref Boehm-Sturm et al 2012 Contrast media molecular imaging
        if rescale == 0
            data=(218456000./ADCmap.reco.data).*((deltaR2.reco.data.^3)./(deltaR2star.reco.data.^2));
        else
            for i = 1:deltaR2star.reco.no_slices
                data(:,:,:,i)=(218456000./imresize(ADCmap.reco.data(:,:,1,i),[final_res final_res],'bilinear')).*((imresize(deltaR2.reco.data(:,:,1,i),[final_res final_res],'bilinear').^3)./...
                    (imresize(deltaR2star.reco.data(:,:,1,i),[final_res final_res],'bilinear').^2));
            end
        end
    end
    
    
end


% copy and update the SO2map structure
Density_map.reco = rmfield(Density_map.reco, 'data');
Density_map.reco.data = data;
if rescale == 1
Density_map.reco.no_samples = final_res;
Density_map.reco.no_views = final_res;
end
Density_map.reco.texte = 'Density';
Density_map.reco.globalmax = max(data(:));
Density_map.reco.globalmin = min(data(:));
if strcmp(ADCmap_filename, 'None')
    if isfield(deltaR2star.reco, 'paramQuantif')
        ParamConfig=sprintf('##$QuantifMethod=''328999*((deltaR2^3)/(deltaR2*^2))''\n##$final_res=%s\n##$DeltaR2* map=%s\n##$DeltaR2 map=%s\n\n##$DeltaR2* scan info\n%s\n\n##$DeltaR2 scan info\n%s\n',...
            add_parameters{:}{1}, deltaR2star_filename,deltaR2_filename, deltaR2star.reco.paramQuantif, deltaR2.reco.paramQuantif);
    else
        ParamConfig=sprintf('##$QuantifMethod=''329000*((deltaR2^3)/(deltaR2*^2))''\n##$final_res=%s\n##$DeltaR2* map=%s\n##$DeltaR2 map=%s\n\',...
            add_parameters{:}{1}, deltaR2star_filename,deltaR2_filename);
    end
else
     if isfield(deltaR2star.reco, 'paramQuantif')
        ParamConfig=sprintf('##$QuantifMethod=''218456000/ADC *((deltaR2^3)/(deltaR2*^2))''\n##$final_res=%s\n##$ADCmap map=%s\n##$DeltaR2* map=%s\n##$DeltaR2 map=%s\n\n##$DeltaR2* scan info\n%s\n\n##$DeltaR2 scan info\n%s\n',...
            add_parameters{:}{1}, ADCmap_filename{:}, deltaR2star_filename,deltaR2_filename, deltaR2star.reco.paramQuantif, deltaR2.reco.paramQuantif);
    else
        ParamConfig=sprintf('##$QuantifMethod=''218456000/ADC*((deltaR2^3)/(deltaR2*^2))''\n##$final_res=%s\n##$ADCmap map=%s\n##$DeltaR2* map=%s\n##$DeltaR2 map=%s\n\',...
            add_parameters{:}{1}, ADCmap_filename{:}, deltaR2star_filename,deltaR2_filename);
    end
end
Density_map.reco.paramQuantif = ParamConfig;
if isfield(Density_map.reco, 'wrongpix')
Density_map.reco = rmfield(Density_map.reco, 'wrongpix');
end
if isfield(Density_map.reco, 'err')
Density_map.reco = rmfield(Density_map.reco, 'err');
end
Density_map.reco=orderfields(Density_map.reco);
Density_map.clip = [0 1000 1];


