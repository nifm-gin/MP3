function deltaR2_map = parametric_map_deltaR2(filename_pre, filename_post, add_parameters)
%generate a deltaR2_map.reco from a pre and post contrast agent T2* maps 
% this code come form the VSI_compute_deltaR2 function -EB 2005

% Load data (Pre/post contrast agent)
% load the MGESE or MGEFIDSE file




fid=fopen(filename_pre ,'r');
if fid>0
    fclose(fid);
    data_pre = load(filename_pre);
    RAW_data_pre = data_pre.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ Somthing wrong with the data\n##$MGESEpre=%s\n##$',...
        filename_pre);
    msgbox(warning_text, 'deltaR2 map warning') ;
    deltaR2_map = [];
    return
end
fid=fopen(filename_post ,'r');
if fid>0
    fclose(fid);
    data_post = load(filename_post);
    RAW_data_post = data_post.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ Somthing wrong with the data\n##$MGESEpost=%s\n##$',...
        filename_post);
    msgbox(warning_text, 'deltaR2 map warning') ;
    deltaR2_map = [];
    return
end
% %% smooth data
% for i=1:size(RAW_data_pre,1)
%     if sum(data_in_vector(i,:)) ~= 0
%         data_in_vector(i,:) = smooth(data_in_vector(i,:), 'lowess')';
%     end
% end
% %% smooth data
% for i=1:size(data_in_vector,1)
%     if sum(data_in_vector(i,:)) ~= 0
%         data_in_vector(i,:) = smooth(data_in_vector(i,:), 'lowess')';
%     end
% end


se_echotime = scan_acqp('##$SpinEchoTime=',data_pre.uvascim.image.texte,1);
if isnan(se_echotime)
    se_echotime = scan_acqp('##$PVM_EchoTime=',data_pre.uvascim.image.texte,1);
end
se_echo_pos = abs(RAW_data_pre.echotime - se_echotime) <= 2;
se_echo_pos = find(se_echo_pos == 1);

nb_echo_used=str2double(add_parameters{:}(1));
final_res = add_parameters{:}(2);
if strcmp(final_res, 'Original')
   rescale = 0;
else
    rescale = 1;
    final_res = str2double(add_parameters{:}(2));
end


% check data compatibility (slice thickness and slice number)
if RAW_data_pre.thickness ~= RAW_data_post.thickness
    warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ a slice thickness missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
        filename_pre,filename_post);
    msgbox(warning_text, 'deltaR2 map warning') ;
    deltaR2_map = [];
    return
end
if RAW_data_pre.no_samples ~= RAW_data_post.no_samples && rescale == 0
     warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ a resolution missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
        filename_pre,filename_post);
    msgbox(warning_text, 'deltaR2 map warning') ;
    deltaR2_map = [];
    return
end
if RAW_data_pre.no_slices ~= RAW_data_post.no_slices
    deltaR2_slice_nbr = 0;
    deltaR2_map.reco.fov_offsets = [];
    deltaR2_map.reco.fov_orientation = [];
    deltaR2_map.reco.label = {};
    deltaR2_map.reco.phaselabel = {};
    deltaR2_map.reco.fov_phase_orientation = [];
    for i = 1:size(RAW_data_pre.data, 4)
        for j = 1:size(RAW_data_post.data, 4)
            if abs(RAW_data_pre.fov_offsets(3,1,i) - RAW_data_post.fov_offsets(3,1,j)) < 1e-5
                deltaR2_slice_nbr = deltaR2_slice_nbr+1;
                % Compute the CMRO2 map each slice with the same offset
                if nb_echo_used == 1
                    if rescale == 0
                        temp_avant=squeeze(RAW_data_pre.data(:,:,se_echo_pos,i));
                        temp_apres=squeeze(RAW_data_post.data(:,:,se_echo_pos,j));
                    else
                        temp_avant=imresize(squeeze(RAW_data_pre.data(:,:,se_echo_pos,i)),[final_res final_res],'bilinear');
                        temp_apres=imresize(squeeze(RAW_data_post.data(:,:,se_echo_pos,j)),[final_res final_res],'bilinear');
                    end
                else
                    if rescale == 0
                        temp_avant=squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,i),3));
                        temp_apres=squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,j),3));
                    else
                        temp_avant=imresize(squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,i),3)),[final_res final_res],'bilinear');
                        temp_apres=imresize(squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,j),3)),[final_res final_res],'bilinear');
                    end
                end
                index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
                index_apres=find(temp_apres<(1e-3*mean(temp_apres(:))));
                
                %To avoid division by small numbers
                temp_apres(index_apres)=1;
                warning off %#ok<WNOFF>
                deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
                warning on %#ok<WNON>
                deltaR2(index_avant)=0;
                deltaR2(index_apres)=0;
                deltaR2_map.reco.data(:,:,1,deltaR2_slice_nbr)=deltaR2;
                
                %Compute the wrongpix map : wrongpix=1 if temp_apres or temp_avant < 0.001*mean(image)
                if rescale == 0
                    tmpim=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views);
                else
                    tmpim=zeros(final_res,final_res);
                end
                tmpim(index_avant)=1;
                tmpim(index_apres)=1;
                tmpim(deltaR2<0)=1;
                deltaR2_map.reco.wrongpix(:,:,1,deltaR2_slice_nbr)=tmpim;
                % Update the SO2map structure
                deltaR2_map.reco.fov_offsets(:,1,deltaR2_slice_nbr,1) = RAW_data_pre.fov_offsets(:,1,i,1);
                deltaR2_map.reco.fov_orientation(:,1,deltaR2_slice_nbr,1) = RAW_data_pre.fov_orientation(:,1,i,1);
                deltaR2_map.reco.label(1,deltaR2_slice_nbr,1) = RAW_data_pre.label(1,i,1);
                deltaR2_map.reco.phaselabel(1,deltaR2_slice_nbr,1) = RAW_data_pre.phaselabel(1,i,1);
                deltaR2_map.reco.fov_phase_orientation(1,deltaR2_slice_nbr,1) = RAW_data_pre.fov_phase_orientation(1,i,1);    
            end
        end
    end
    if deltaR2_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the deltaR2 map because there is\n##$ no slice offset match between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
            filename_SO2map,filename_CBF);
        msgbox(warning_text, 'deltaR2 map warning') ;
        return
    end
    deltaR2_map.reco.no_slices=deltaR2_slice_nbr;
else
    % Compute the deltaR2 map for all slices
    deltaR2_map.reco='';
    if rescale == 0
        deltaR2_map.reco.data=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views,1,RAW_data_pre.no_slices);     % DeltaR2 map
        deltaR2_map.reco.wrongpix=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views,1,RAW_data_pre.no_slices);  % excluded voxels
    else
        deltaR2_map.reco.data=zeros(final_res,final_res,1,RAW_data_pre.no_slices);     
        deltaR2_map.reco.wrongpix=zeros(final_res,final_res,1,RAW_data_pre.no_slices);  
    end
    
    for m_slice=1:RAW_data_pre.no_slices
        if nb_echo_used == 1
            if rescale == 0
                
                temp_avant=squeeze(RAW_data_pre.data(:,:,se_echo_pos,m_slice));
                temp_apres=squeeze(RAW_data_post.data(:,:,se_echo_pos,m_slice));
            else
                
                temp_avant=imresize(squeeze(RAW_data_pre.data(:,:,se_echo_pos,m_slice)),[final_res final_res],'bilinear');
                temp_apres=imresize(squeeze(RAW_data_post.data(:,:,se_echo_pos,m_slice)),[final_res final_res],'bilinear');
            end
        else
            if rescale == 0
                temp_avant=squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3));
                temp_apres=squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3));
            else
                temp_avant=imresize(squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3)),[final_res final_res],'bilinear');
                temp_apres=imresize(squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,m_slice),3)),[final_res final_res],'bilinear');
            end
        end
        index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
        index_apres=find(temp_apres<(1e-3*mean(temp_apres(:))));
        
        %To avoid division by small numbers
        temp_apres(index_apres)=1;
        warning off %#ok<WNOFF>
        deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
        warning on %#ok<WNON>
        deltaR2(index_avant)=0;
        deltaR2(index_apres)=0;
        deltaR2_map.reco.data(:,:,1,m_slice)=deltaR2;
        
        %Compute the err map
        %proc.err(:,:,1,m_slice)=tmpim;
        
        %Compute the wrongpix map : wrongpix=1 if temp_apres or temp_avant < 0.001*mean(image)
        if rescale == 0
            tmpim=zeros(RAW_data_pre.no_samples,RAW_data_pre.no_views);
        else
            tmpim=zeros(final_res,final_res);
        end
        tmpim(index_avant)=1;
        tmpim(index_apres)=1;
        tmpim(deltaR2<0)=1;
        deltaR2_map.reco.wrongpix(:,:,1,m_slice)=tmpim;
    end
    %Adapt the fov offsets and orientations infos
    deltaR2_map.reco.fov_offsets=RAW_data_pre.fov_offsets(:,1,:,1);
    deltaR2_map.reco.fov_orientation=RAW_data_pre.fov_orientation(:,1,:,1);
    deltaR2_map.reco.label=RAW_data_pre.label(1,:,1);
    deltaR2_map.reco.phaselabel=RAW_data_pre.phaselabel(1,:,1);
    deltaR2_map.reco.fov_phase_orientation=RAW_data_pre.fov_phase_orientation(1,:,1);
    deltaR2_map.reco.no_slices=RAW_data_pre.no_slices;
end

%Set the dimensions of the proc structure
if rescale == 0
    deltaR2_map.reco.no_samples=RAW_data_pre.no_samples;
    deltaR2_map.reco.no_views=RAW_data_pre.no_views;
else
    deltaR2_map.reco.no_samples = final_res;
    deltaR2_map.no_views = final_res;
end
deltaR2_map.reco.no_echoes=1; %Number of parameters stored here
deltaR2_map.reco.no_expts=1;

%Complete the reco structure
deltaR2_map.reco.texte='DeltaR2';
deltaR2_map.reco.unit='ms-1';
deltaR2_map.reco.date=date;
deltaR2_map.reco.displayedecho=1;
deltaR2_map.reco.displayedslice=1;
deltaR2_map.reco.displayedexpt=1;
deltaR2_map.reco_number = 1;
deltaR2_map.scan_number = 104;
deltaR2_map.clip=[0 0.4 1];
deltaR2_map.reco.globalmax = max(deltaR2_map.reco.data(:));
deltaR2_map.reco.globalmin = min(deltaR2_map.reco.data(:));
deltaR2_map.reco.thickness = data_pre.uvascim.image.reco.thickness;

ParamConfig=sprintf('##$QuantifMethod=''-(1/se_echotime) * log(temp_apres ./ temp_avant)''\n##$Spin Echo time=%d\n##$Spin Echo position=%d\n##$Number of echo used=%s\n##$T2pre=%s\n##$T2post=%s\n\n##$Sequence pre info\n%s\n\n##$Sequence post info\n%s',...
    se_echotime, se_echo_pos, add_parameters{:}{1},filename_pre,filename_post, [RAW_data_pre.iminfos{:}], [RAW_data_post.iminfos{:}]);
deltaR2_map.reco.paramQuantif = ParamConfig;
deltaR2_map.reco=orderfields(deltaR2_map.reco);

%complete the structure
deltaR2_map.acq = data_pre.uvascim.image.acq;
deltaR2_map.filename = data_pre.uvascim.image.filename;


