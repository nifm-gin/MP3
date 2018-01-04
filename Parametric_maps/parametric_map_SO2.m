function SO2map = parametric_map_SO2(R2prim_filename, BVf_map_filename, T1map_filename,  add_parameters)
% generate a SO2 from a R2prim and a T2*corr3D scan
% this code come from the the S02-fig coded by T. Christen

B0=str2double(add_parameters{:}(1));
gamma=str2double(add_parameters{:}(2)); gamma = gamma*10^8;
deltaxi_HBC=str2double(add_parameters{:}(3)); deltaxi_HBC = deltaxi_HBC*10^-6;
HCT = str2double(add_parameters{:}(4));
HCTmap_filename = add_parameters{:}(6);


% load the R2prim file
fid=fopen(R2prim_filename ,'r');
if fid>0
    fclose(fid);
    R2prim = load(R2prim_filename);
    R2prim = R2prim.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ Somthing wrong with the data\n##$R2prim=%s\n##$',...
        R2prim_filename);
    msgbox(warning_text, 'SO2 map warning') ;
    SO2map = [];
    return
end

% load the BVf file
fid=fopen(BVf_map_filename ,'r');
if fid>0
    fclose(fid);
    BVf_map = load(BVf_map_filename);
    BVf_map = BVf_map.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ Somthing wrong with the data\n##$BVf=%s\n##$',...
        BVf_map_filename);
    msgbox(warning_text, 'SO2 map warning') ;
    SO2map = [];
    return
end

% load the T1map if needed
if ~isempty(T1map_filename)
    fid=fopen(T1map_filename ,'r');
    if fid>0
        fclose(fid);
        T1map = load(T1map_filename);
        T1map = T1map.uvascim.image;
    else
        warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ Somthing wrong with the data\n##$T1map=%s\n##$',...
            R2prim_filename);
        msgbox(warning_text, 'SO2 map warning') ;
        SO2map = [];
        return
    end
end
% load the HCT_map if needed
if ~strcmp(HCTmap_filename, 'None')
    if ~isempty(HCTmap_filename{:})
        fid=fopen(HCTmap_filename{:} ,'r');
        if fid>0
            fclose(fid);
            HCTmap = load(HCTmap_filename{:});
            HCTmap = HCTmap.uvascim.image;
        else
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ Somthing wrong with the data\n##$HCTmap=%s\n##$',...
                HCTmap_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
    else
        warning_text = sprintf('##$ Can not calculate the SO2 map because \n##$ No HCTmap has been found\n');
        msgbox(warning_text, 'SO2 map warning') ;
        SO2map = [];
        return
    end
end

% 1st case ==> apply the same DeltaKhi everywhere
% otherwise using the T1map to segment the brain and apply 3 differents
% DeltaKhi: gray, white and CF

if isempty(T1map_filename)
    % Apply the same DeltaKhi AND use HCT variable (using a HCTmap)
    if ~strcmp(HCTmap_filename, 'None')
        if ~isfield( HCTmap.reco, 'no_slices')
            HCTmap.reco.no_slices = 1;
        end
        temp_slice_nbr = [R2prim.reco.no_slices BVf_map.reco.no_slices  HCTmap.reco.no_slices];
        HCTmap_resized = HCTmap;
        
        for a=1:size(HCTmap.reco.data,4)
            tmp(:,:,1,a) = imresize(HCTmap.reco.data(:,:,1,a),[size(R2prim.reco.data, 1) size(R2prim.reco.data,2)],'bilinear');
        end
        HCTmap_resized.reco.data = tmp;
        SO2map = BVf_map;
        SO2map.reco.data = [];
        if length(find(temp_slice_nbr == R2prim.reco.no_slices)) ~= numel(temp_slice_nbr)
            SO2map_slice_nbr = 0;
            SO2map.reco.fov_offsets = [];
            SO2map.reco.fov_orientation = [];
            SO2map.reco.label = {};
            SO2map.reco.phaselabel = {};
            SO2map.reco.fov_phase_orientation = [];
            for i = 1:size(R2prim.reco.data, 4)
                for j = 1:size(BVf_map.reco.data, 4)
                    if abs(R2prim.reco.fov_offsets(3,1,i) - BVf_map.reco.fov_offsets(3,1,j)) < 1e-5
                        for x = 1:size(HCTmap_resized.reco.data, 4)
                            if abs(R2prim.reco.fov_offsets(3,1,i) - HCTmap_resized.reco.fov_offsets(3,1,x)) < 1e-5
                                SO2map_slice_nbr = SO2map_slice_nbr+1;
                                % Compute the SO2map each slice with the same offset
                                SO2map.reco.data(:,:,1,SO2map_slice_nbr)=(1-(R2prim.reco.data(:,:,1,i))./...
                                    ((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*deltaxi_HBC.*(HCTmap_resized.reco.data(:,:,1,x)./100)*B0))*100;%en %
                                % Update the SO2map structure
                                SO2map.reco.fov_offsets(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_offsets(:,1,j,1);
                                SO2map.reco.fov_orientation(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_orientation(:,1,j,1);
                                SO2map.reco.label(1,SO2map_slice_nbr,1) = BVf_map.reco.label(1,j,1);
                                SO2map.reco.phaselabel(1,SO2map_slice_nbr,1) = BVf_map.reco.phaselabel(1,j,1);
                                SO2map.reco.fov_phase_orientation(1,SO2map_slice_nbr,1) = BVf_map.reco.fov_phase_orientation(1,j,1);
                            end
                        end
                    end
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s ##$T1map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
            SO2map.reco.no_slices=SO2map_slice_nbr;
        else
            % calculate SO2 map for all slices
            %             SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
        end
        
        
        % Apply the same DeltaKhi AND use constant HCT value
    else
        % check data compatibility (slice thickness and slice number)
        if R2prim.reco.thickness ~= BVf_map.reco.thickness
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a slice thickness missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s',...
                R2prim_filename,BVf_map_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        SO2map = BVf_map;
        if R2prim.reco.no_slices ~= BVf_map.reco.no_slices
            SO2map_slice_nbr = 0;
            SO2map.reco.fov_offsets = [];
            SO2map.reco.fov_orientation = [];
            SO2map.reco.label = {};
            SO2map.reco.phaselabel = {};
            SO2map.reco.fov_phase_orientation = [];
            for i = 1:size(R2prim.reco.data, 4)
                for j = 1:size(BVf_map.reco.data, 4)
                    if abs(R2prim.reco.fov_offsets(3,1,i) - BVf_map.reco.fov_offsets(3,1,j)) < 1e-5
                        SO2map_slice_nbr = SO2map_slice_nbr+1;
                        % Compute the SO2map map each slice with the same offset
                        SO2map.reco.data(:,:,1,SO2map_slice_nbr)=(1-(R2prim.reco.data(:,:,1,i))./((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
                        % Update the SO2map structure
                        SO2map.reco.fov_offsets(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_offsets(:,1,j,1);
                        SO2map.reco.fov_orientation(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_orientation(:,1,j,1);
                        SO2map.reco.label(1,SO2map_slice_nbr,1) = BVf_map.reco.label(1,j,1);
%                         SO2map.reco.phaselabel(1,SO2map_slice_nbr,1) = BVf_map.reco.phaselabel(1,j,1);
                        SO2map.reco.fov_phase_orientation(1,SO2map_slice_nbr,1) = BVf_map.reco.fov_phase_orientation(1,j,1);
                    end
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
            SO2map.reco.no_slices=SO2map_slice_nbr;
        else
            % calculate SO2 map for all slices
            SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
        end
    end
else
    % case using the T1map to segment the brain and apply 3 differents
    % DeltaKhi: gray, white and CF   AND use HCT variable (using a HCTmap)
    if ~strcmp(HCTmap_filename, 'None')
        % check data compatibility (slice thickness and slice number)
        temp_thickness = [R2prim.reco.thickness BVf_map.reco.thickness  T1map.reco.thickness];
        temp_slice_nbr = [R2prim.reco.no_slices BVf_map.reco.no_slices  T1map.reco.no_slices];
        temp_resolution = [R2prim.reco.no_samples BVf_map.reco.no_samples  T1map.reco.no_samples];
        if length(find(temp_thickness == R2prim.reco.thickness)) ~= numel(temp_thickness)
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a slice thickness missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s and \n##$T1map=%s',...
                R2prim_filename,BVf_map_filename, T1map_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        if length(find(temp_resolution == R2prim.reco.no_samples)) ~= numel(temp_resolution)
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a resolution missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s\n##$ and\n##$T1map=%s\n',...
                deltaR2star_filename,deltaR2_filename,ADC_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        
        if ~isfield( HCTmap.reco, 'no_slices')
            HCTmap.reco.no_slices = 1;
        end
        temp_slice_nbr = [R2prim.reco.no_slices BVf_map.reco.no_slices  HCTmap.reco.no_slices];
        HCTmap_resized = HCTmap;
        
        for a=1:size(HCTmap.reco.data,4)
            tmp(:,:,1,a) = imresize(HCTmap.reco.data(:,:,1,a),[size(R2prim.reco.data, 1) size(R2prim.reco.data,2)],'bilinear');
        end
        HCTmap_resized.reco.data = tmp;
        SO2map = BVf_map;
        SO2map.reco.data = [];
        if length(find(temp_slice_nbr == R2prim.reco.no_slices)) ~= numel(temp_slice_nbr)
            SO2map_slice_nbr = 0;
            SO2map.reco.fov_offsets = [];
            SO2map.reco.fov_orientation = [];
            SO2map.reco.label = {};
            SO2map.reco.phaselabel = {};
            SO2map.reco.fov_phase_orientation = [];
            for i = 1:size(R2prim.reco.data, 4)
                for j = 1:size(BVf_map.reco.data, 4)
                    if abs(R2prim.reco.fov_offsets(3,1,i) - BVf_map.reco.fov_offsets(3,1,j)) < 1e-5
                        for x = 1:size(HCTmap_resized.reco.data, 4)
                            if abs(R2prim.reco.fov_offsets(3,1,i) - HCTmap_resized.reco.fov_offsets(3,1,x)) < 1e-5
                                SO2map_slice_nbr = SO2map_slice_nbr+1;
                                % Compute the SO2map each slice with the same offset
                                
                                mask_white = T1map.reco.data(:,:,1,j) <= 1300;
                                mask_gray = logical(T1map.reco.data(:,:,1,j) > 1300) & logical(T1map.reco.data(:,:,1,j) < 1900);
                                mask_CF = T1map.reco.data(:,:,1,j) >= 1900;
                                %
                                % calculate SO2 map for all slices
                                tmp=(1-(R2prim.reco.data(:,:,1,i))./((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*5.36e-07.*(HCTmap_resized.reco.data(:,:,1,x)./100)*B0))*100;%en %
                                tmp(isnan(tmp)) = 0;
                                mask_white2 = mask_white.*tmp;
                                %
                                tmp=(1-(R2prim.reco.data(:,:,1,i))./((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*2.6e-07.*(HCTmap_resized.reco.data(:,:,1,x)./100)*B0))*100;%en %
                                tmp(isnan(tmp)) = 0;
                                mask_gray2 = mask_gray.*tmp;
                                
                                tmp=(1-(R2prim.reco.data(:,:,1,i))./((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*0.33e-07.*(HCTmap_resized.reco.data(:,:,1,x)./100)*B0))*100;%en %
                                tmp(isnan(tmp)) = 0;
                                mask_CF2 = mask_CF.*tmp;
                                
                                SO2map.reco.data(:,:,1,SO2map_slice_nbr)= mask_white2 + mask_gray2 + mask_CF2;
                                
                                
                                %                                 SO2map.reco.data(:,:,1,SO2map_slice_nbr)=(1-(R2prim.reco.data(:,:,1,i))./...
                                %                                     ((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*deltaxi_HBC.*(HCTmap_resized.reco.data(:,:,1,x)./100)*B0))*100;%en %
                                % Update the SO2map structure
                                SO2map.reco.fov_offsets(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_offsets(:,1,j,1);
                                SO2map.reco.fov_orientation(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_orientation(:,1,j,1);
                                SO2map.reco.label(1,SO2map_slice_nbr,1) = BVf_map.reco.label(1,j,1);
                                SO2map.reco.phaselabel(1,SO2map_slice_nbr,1) = BVf_map.reco.phaselabel(1,j,1);
                                SO2map.reco.fov_phase_orientation(1,SO2map_slice_nbr,1) = BVf_map.reco.fov_phase_orientation(1,j,1);
                            end
                        end
                    end
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s ##$T1map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
            SO2map.reco.no_slices=SO2map_slice_nbr;
            %         else
            %             % calculate SO2 map for all slices
            % %             SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
            %
            %             % classification of the brain in 3 classes (Gray, White and
            %             % Cerebrospinal fluid (CF))based on the T1 value
            %             mask_white = T1map.reco.data <= 1300;
            %             mask_gray = logical(T1map.reco.data > 1300) & logical(T1map.reco.data < 1900);
            %             mask_CF = T1map.reco.data >= 1900;
            %             %
            %             % calculate SO2 map for all slices
            %             tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*5.36e-07*HCT*B0))*100;%en %
            %             tmp(isnan(tmp)) = 0;
            %             mask_white2 = mask_white.*tmp;
            %             %
            %             tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*2.6e-07*HCT*B0))*100;%en %
            %             tmp(isnan(tmp)) = 0;
            %             mask_gray2 = mask_gray.*tmp;
            %
            %             tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*0.33e-07*HCT*B0))*100;%en %
            %             tmp(isnan(tmp)) = 0;
            %             mask_CF2 = mask_CF.*tmp;
            %
            %             SO2map.reco.data = mask_white2 + mask_gray2 + mask_CF2;
        end
        
        % case using the T1map to segment the brain and apply 3 differents
        % DeltaKhi: gray, white and CF   AND use constant HCT value
    else
        % check data compatibility (slice thickness and slice number)
        temp_thickness = [R2prim.reco.thickness BVf_map.reco.thickness  T1map.reco.thickness];
        temp_slice_nbr = [R2prim.reco.no_slices BVf_map.reco.no_slices  T1map.reco.no_slices];
        temp_resolution = [R2prim.reco.no_samples BVf_map.reco.no_samples  T1map.reco.no_samples];
        if length(find(temp_thickness == R2prim.reco.thickness)) ~= numel(temp_thickness)
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a slice thickness missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s and \n##$T1map=%s',...
                R2prim_filename,BVf_map_filename, T1map_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        if length(find(temp_resolution == R2prim.reco.no_samples)) ~= numel(temp_resolution)
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a resolution missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s\n##$ and\n##$T1map=%s\n',...
                deltaR2star_filename,deltaR2_filename,ADC_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        
        SO2map = BVf_map;
        if length(find(temp_slice_nbr == R2prim.reco.no_slices)) ~= numel(temp_slice_nbr)
            SO2map_slice_nbr = 0;
            SO2map.reco.fov_offsets = [];
            SO2map.reco.fov_orientation = [];
            SO2map.reco.label = {};
            SO2map.reco.phaselabel = {};
            SO2map.reco.fov_phase_orientation = [];
            for i = 1:size(R2prim.reco.data, 4)
                for j = 1:size(BVf_map.reco.data, 4)
                    if abs(R2prim.reco.fov_offsets(3,1,i) - BVf_map.reco.fov_offsets(3,1,j)) < 1e-5
                        for x = 1:size(T1map.reco.data, 4)
                            if abs(R2prim.reco.fov_offsets(3,1,i) - T1map.reco.fov_offsets(3,1,x)) < 1e-5
                                SO2map_slice_nbr = SO2map_slice_nbr+1;
                                % Compute the SO2map each slice with the same offset
                                SO2map.reco.data(:,:,1,SO2map_slice_nbr)=(1-(R2prim.reco.data(:,:,1,i))./((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
                                % Update the SO2map structure
                                SO2map.reco.fov_offsets(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_offsets(:,1,j,1);
                                SO2map.reco.fov_orientation(:,1,SO2map_slice_nbr,1) = BVf_map.reco.fov_orientation(:,1,j,1);
                                SO2map.reco.label(1,SO2map_slice_nbr,1) = BVf_map.reco.label(1,j,1);
                                SO2map.reco.phaselabel(1,SO2map_slice_nbr,1) = BVf_map.reco.phaselabel(1,j,1);
                                SO2map.reco.fov_phase_orientation(1,SO2map_slice_nbr,1) = BVf_map.reco.fov_phase_orientation(1,j,1);
                            end
                        end
                    end
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s ##$T1map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
            SO2map.reco.no_slices=SO2map_slice_nbr;
        else
            % classification of the brain in 3 classes (Gray, White and
            % Cerebrospinal fluid (CF))based on the T1 value
            mask_white = T1map.reco.data <= 1300;
            mask_gray = logical(T1map.reco.data > 1300) & logical(T1map.reco.data < 1900);
            mask_CF = T1map.reco.data >= 1900;
            %
            % calculate SO2 map for all slices
            tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*5.36e-07*HCT*B0))*100;%en %
            tmp(isnan(tmp)) = 0;
            mask_white2 = mask_white.*tmp;
            %
            tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*2.6e-07*HCT*B0))*100;%en %
            tmp(isnan(tmp)) = 0;
            mask_gray2 = mask_gray.*tmp;
            
            tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*0.33e-07*HCT*B0))*100;%en %
            tmp(isnan(tmp)) = 0;
            mask_CF2 = mask_CF.*tmp;
            
            SO2map.reco.data = mask_white2 + mask_gray2 + mask_CF2;
            
            % test if a relationship between the T1 value and the deltaKhi
            % exist
            %         deltakhi_map = ((T1map.reco.data.*-0.0041)+9.23).*0.0000001;
            %         SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi.*deltakhi_map*HCT*B0))*100;%en %
            
        end
    end
end


% complete structure
SO2map.reco.globalmax = max(SO2map.reco.data(:));
SO2map.reco.globalmin = min(SO2map.reco.data(:));
SO2map.reco.texte = 'SO2';
SO2map.reco.unit = '%';

ParamConfig=sprintf('##$QuantifMethod=''(1 - (1/(R2prim) / ((BVf_map/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0)) * 100''\n##$B0=%s\n##$gamma=%s\n##$deltaxi_HBC=%s\n##$Hct=%s\n##$R2prim file=%s\n##$BVf_map file=%s\n##$R2prim scan info\n%s\n\n##$BVf_map scan info\n%s',...
    add_parameters{:}{1},...
    strcat(add_parameters{:}{2}, '.10^8'),strcat(add_parameters{:}{3}, '.10^-6'), add_parameters{:}{4},...
    R2prim_filename,BVf_map_filename, R2prim.reco.paramQuantif, BVf_map.reco.paramQuantif);
SO2map.reco.paramQuantif = ParamConfig;
SO2map.reco=orderfields(SO2map.reco);

SO2map.clip=[0 100 1];
