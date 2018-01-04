function VSImap = parametric_map_VSI(deltaR2star_filename, deltaR2_filename, ADC_filename, add_parameters)
% generate a SO2 from a R2prim and a T2*corr3D scan
% this code come from the the S02-fig coded by T. Christen

B0=str2double(add_parameters{:}(1));
gamma=str2double(add_parameters{:}(2)); gamma = gamma*10^8;
deltaxi=str2double(add_parameters{:}(3)); deltaxi = deltaxi*10^-6;

% load the deltaR2star file
fid=fopen(deltaR2star_filename ,'r');
if fid>0
    fclose(fid);
    deltaR2star = load(deltaR2star_filename);
    deltaR2star = deltaR2star.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ Somthing wrong with the data\n##$deltaR2*=%s\n##$',...
        deltaR2star_filename);
    msgbox(warning_text, 'VSI map warning') ;
    VSImap = [];
    return
end
% load the deltaR2 file
fid=fopen(deltaR2_filename ,'r');
if fid>0
    fclose(fid);
    deltaR2 = load(deltaR2_filename);
    deltaR2 = deltaR2.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ Somthing wrong with the data\n##$deltaR2=%s\n##$',...
        deltaR2_filename);
    msgbox(warning_text, 'VSI map warning') ;
    VSImap = [];
    return
end
% load the ADC file
fid=fopen(ADC_filename ,'r');
if fid>0
    fclose(fid);
    ADC = load(ADC_filename);
    ADC = ADC.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ Somthing wrong with the data\n##$ADC=%s\n##$',...
        ADC_filename);
    msgbox(warning_text, 'VSI map warning') ;
    VSImap = [];
    return
end
temp_thickness = [deltaR2star.reco.thickness deltaR2.reco.thickness  ADC.reco.thickness];
temp_slice_nbr = [deltaR2star.reco.no_slices deltaR2.reco.no_slices  ADC.reco.no_slices];
temp_resolution = [deltaR2star.reco.no_samples deltaR2.reco.no_samples  ADC.reco.no_samples];

% check data compatibility (slice thickness and slice number)
if  length(find(temp_thickness == deltaR2star.reco.thickness)) ~= numel(temp_thickness)
    warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ a slice thickness missmatch between\n##$deltaR2star=%s\n##$ and \n##$deltaR2=%s\n##$ and\n##$ADC=%s\n',...
        deltaR2star_filename,deltaR2_filename,ADC_filename);
    msgbox(warning_text, 'VSI map warning') ;
     VSImap = [];
    return
end
if length(find(temp_resolution == deltaR2star.reco.no_samples)) ~= numel(temp_resolution)
     warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ a resolution missmatch between\n##$deltaR2star=%s\n##$ and \n##$deltaR2=%s\n##$ and\n##$ADC=%s\n',...
       deltaR2star_filename,deltaR2_filename,ADC_filename);
    msgbox(warning_text, 'VSI map warning') ;
     VSImap = [];
    return
end

VSImap = deltaR2;

if length(find(temp_slice_nbr == deltaR2star.reco.no_slices)) ~= numel(temp_slice_nbr)
    VSImap_slice_nbr = 0;
    VSImap.reco.fov_offsets = [];
    VSImap.reco.fov_orientation = [];
    VSImap.reco.label = {};
    VSImap.reco.phaselabel = {};
    VSImap.reco.fov_phase_orientation = [];
    for i = 1:size(deltaR2star.reco.data, 4)
        for j = 1:size(deltaR2.reco.data, 4)
            if abs(deltaR2star.reco.fov_offsets(3,1,i) - deltaR2.reco.fov_offsets(3,1,j)) < 1e-5
                for x = 1:size(ADC.reco.data, 4)
                    if abs(deltaR2star.reco.fov_offsets(3,1,i) - ADC.reco.fov_offsets(3,1,x)) < 1e-5
                        VSImap_slice_nbr = VSImap_slice_nbr+1;
                        % Compute the VSI map each slice with the same offset
                        ratio = deltaR2star.reco.data(:,:,1,i) ./ deltaR2.reco.data(:,:,1,j); % ratio : no unit
                        index_ratiofinite= find(~isfinite(ratio));
                        index_rationan= find(isnan(ratio));
                        if ~isempty(index_ratiofinite)
                            ratio(index_ratiofinite) = 0;
                        end
                        if ~isempty(index_rationan)
                            ratio(index_rationan)    = 0;
                        end
                        index_neg= find(ratio<0);
                        if ~isempty(index_neg)
                            ratio(index_neg)    = 0;
                        end
                        ratio = ratio.^1.5;
                          slice_data = 1.77^(-1.5) * sqrt(squeeze(ADC.reco.data(:,:,1,x)) ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
                       %%%% Rui data
%                        test = imresize(squeeze(ADC.reco.data(:,:,1,x)), [128 128]);
%                        index_inferieurzero=find(test < 0);
%                        test(index_inferieurzero) = 0;
%                        slice_data = 1.77^(-1.5) * sqrt(test ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
                        %%%
                        
                        index_vsifinite=find(~isfinite(slice_data));
                        slice_data(index_vsifinite)= nan;
                        index_vsinan=find(isnan(slice_data));
                        slice_data(index_vsinan)= nan;
                        index_infzero=find(slice_data < 0);
                        slice_data(index_infzero)= nan;
                        VSImap.reco.data(:,:,1,VSImap_slice_nbr) = slice_data;

                        
                        %Compute the wrongpix map
                        % initiation  
                        wrongpix_deltaR2star=zeros(size(deltaR2star.reco.data, 1),size(deltaR2star.reco.data, 2));
                        wrongpix_deltaR2 = wrongpix_deltaR2star;
                        wrongpix_index_ratiofinite = wrongpix_deltaR2star;
                        wrongpix_index_rationan = wrongpix_deltaR2star;
                        wrongpix_index_infzero = wrongpix_deltaR2star;
                        wrongpix_index_vsifinite = wrongpix_deltaR2star;
                        wrongpix_index_vsinan = wrongpix_deltaR2star;
                        
                        wrongpix_deltaR2star(find(deltaR2star.reco.wrongpix(:,:,i)>0)) = 1;%#ok<FNDSB>
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_deltaR2star; %Propagation of existing wronpix
                        wrongpix_deltaR2(find(deltaR2.reco.wrongpix(:,:,j)>0)) = 1;%#ok<FNDSB>
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_deltaR2; %Propagation of existing wronpi
                        
                        wrongpix_index_ratiofinite(index_ratiofinite) = 1; %Pixels where proc.var is infinite
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_index_ratiofinite;
                        
                        wrongpix_index_rationan(index_rationan) = 1; %Pixels where proc.var is NaN
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_index_rationan;
                        
                        wrongpix_index_infzero(index_infzero) = 1; %Pixels where VSI < 0
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_index_infzero;
                        wrongpix_index_vsifinite(index_vsifinite) = 1; %Pixels where VSI is infinite
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_index_vsifinite;
                        wrongpix_index_vsinan(index_vsinan) = 1; %Pixels where VSI is nan
                        VSImap.reco.wrongpix(:,:,1,VSImap_slice_nbr)=wrongpix_index_vsinan;
                        
                        % Update the VSI structure
                        VSImap.reco.fov_offsets(:,1,VSImap_slice_nbr,1) = deltaR2.reco.fov_offsets(:,1,j,1);
                        VSImap.reco.fov_orientation(:,1,VSImap_slice_nbr,1) = deltaR2.reco.fov_orientation(:,1,j,1);
                        VSImap.reco.label(1,VSImap_slice_nbr,1) = deltaR2.reco.label(1,j,1);
                        VSImap.reco.phaselabel(1,VSImap_slice_nbr,1) = deltaR2.reco.phaselabel(1,j,1);
                        VSImap.reco.fov_phase_orientation(1,VSImap_slice_nbr,1) = deltaR2.reco.fov_phase_orientation(1,j,1);
                        
                    end
                end
            end
        end
    end
    if VSImap_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ no slice offset match between\n##$deltaR2star=%s\n##$ and \n##$deltaR2=%s\n##$ and\n##$ADC=%s\n',...
            deltaR2star_filename,deltaR2_filename,ADC_filename);
        msgbox(warning_text, 'VSI map warning') ;
        return
    end
    VSImap.reco.no_slices=VSImap_slice_nbr;
else
    % calculate VSI map for all slices
    ratio = deltaR2star.reco.data ./ deltaR2.reco.data; % ratio : no unit
    index_ratiofinite= find(~isfinite(ratio));
    index_rationan= find(isnan(ratio));
    if ~isempty(index_ratiofinite)
        ratio(index_ratiofinite) = 0;
    end
    if ~isempty(index_rationan)
        ratio(index_rationan)    = 0;
    end
     index_neg= find(ratio<0);
    if ~isempty(index_neg)
        ratio(index_neg)    = 0;
    end
    ratio = ratio.^1.5;
    VSImap.reco.data(:,:,1,:) = 1.77^(-1.5) * sqrt(squeeze(ADC.reco.data) ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
    
    index_vsifinite=find(~isfinite(VSImap.reco.data));
    VSImap.reco.data(index_vsifinite)= 0;
    index_vsinan=find(isnan(VSImap.reco.data));
    VSImap.reco.data(index_vsinan)= 0;
    index_infzero=find(VSImap.reco.data < 0);
    VSImap.reco.data(index_infzero)= 0;
    
    %Compute the wrongpix map
    VSImap.reco.wrongpix(find(deltaR2star.reco.wrongpix>0))=1;%#ok<FNDSB> %Propagation of existing wronpix
    VSImap.reco.wrongpix(find(deltaR2.reco.wrongpix>0))=1; %#ok<FNDSB>
    VSImap.reco.wrongpix(index_ratiofinite)=1;					%Pixels where proc.var is infinite
    VSImap.reco.wrongpix(index_rationan)=1;					%Pixels where proc.var is NaN
    VSImap.reco.wrongpix(index_infzero)=1;						%Pixels where VSI < 0
    VSImap.reco.wrongpix(index_vsifinite)=1;					%Pixels where VSI is infinite
    VSImap.reco.wrongpix(index_vsinan)=1;						%Pixels where VSI is nan
    
end


% complete structure
VSImap.reco.globalmax = max(VSImap.reco.data(:));
VSImap.reco.globalmin = min(VSImap.reco.data(:));
VSImap.reco.texte = 'VSI';
VSImap.reco.unit = 'um';
ParamConfig=sprintf('##$QuantifMethod=''1.77^(-1.5) * sqrt(ADC / (gamma * deltaxi * B0)) * (deltaR2* / DeltaR2))\n##$B0=%s\n##$gamma=%s\n##$deltaxi_CA=%s\n##$ADC file=%s\n##$DeltaR2* file=%s\n##$DeltaR2 file=%s\n\n##$ADC scan info\n%s\n\n##$DeltaR2* scan info\n%s\n\n##$DeltaR2 scan info\n%s\n\n',...
    strcat(add_parameters{:}{2}, '.10^8'),strcat(add_parameters{:}{3}, '.10^-6'), '2*pi*gamma*deltaxi_CA*B0',...
    ADC_filename,deltaR2star_filename, deltaR2_filename, [ADC.reco.iminfos{:}], deltaR2star.reco.paramQuantif, deltaR2.reco.paramQuantif);
VSImap.reco.paramQuantif = ParamConfig;
VSImap.reco=orderfields(VSImap.reco);

VSImap.clip=[0 50 1];

