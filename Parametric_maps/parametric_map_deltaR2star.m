function deltaR2star_map = parametric_map_deltaR2star(filename_pre, filename_post,add_parameters)
%generate a deltaR2star_map.reco from a pre and post contrast agent T2* maps
% this code come form the VSI_compute_deltaR2s function -EB 2005

% Load data (Pre/post contrast agent)
% load the MGE file
fid=fopen(filename_pre ,'r');
if fid>0
    fclose(fid);
    data_pre = load(filename_pre);
    RAW_T2star_pre = data_pre.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ Somthing wrong with the data\n##$T2*pre=%s\n##$',...
        filename_pre);
    msgbox(warning_text, 'deltaR2star map warning') ;
    deltaR2star_map = [];
    return
end
fid=fopen(filename_post ,'r');
if fid>0
    fclose(fid);
    data_post = load(filename_post);
    RAW_T2star_post = data_post.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ Somthing wrong with the data\n##$T2*post=%s\n##$',...
        filename_pre);
    msgbox(warning_text, 'deltaR2star map warning') ;
    deltaR2star_map = [];
    return
end
final_res = add_parameters{:}(1);
if strcmp(final_res, 'Original')
   rescale = 0;
else
    rescale = 1;
    final_res = str2double(add_parameters{:}(1));
end

% check data compatibility (slice thickness and slice number)
if RAW_T2star_pre.thickness ~= RAW_T2star_post.thickness
    warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ a slice thickness missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
        filename_pre,filename_post);
    msgbox(warning_text, 'deltaR2star map warning') ;
     deltaR2star_map = [];
    return
end
if RAW_T2star_pre.no_samples ~= RAW_T2star_post.no_samples && rescale == 0
     warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ a resolution missmatch between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
        filename_pre,filename_post);
    msgbox(warning_text, 'deltaR2star map warning') ;
     deltaR2star_map = [];
    return
end
if RAW_T2star_pre.no_slices ~= RAW_T2star_post.no_slices 
    deltaR2star_slice_nbr = 0;
    deltaR2star_map.reco.fov_offsets = [];
    deltaR2star_map.reco.fov_orientation = [];
    deltaR2star_map.reco.label = {};
    deltaR2star_map.reco.phaselabel = {};
    deltaR2star_map.reco.fov_phase_orientation = [];
   
    for i = 1:size(RAW_T2star_pre.data, 4)
        for j = 1:size(RAW_T2star_post.data, 4)
            if abs(RAW_T2star_pre.fov_offsets(3,1,i) - RAW_T2star_post.fov_offsets(3,1,j)) < 1e-5
                deltaR2star_slice_nbr = deltaR2star_slice_nbr+1; 
                if rescale == 0
                     r2savant=RAW_T2star_pre.data(:,:,1,i);
                     r2sapres=RAW_T2star_post.data(:,:,1,j);
                     wrongpix_pre=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,1);  
                    wrongpix_post=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,1);  
                else
                    r2savant=imresize(RAW_T2star_pre.data(:,:,:,i),[final_res final_res],'bilinear');
                    r2sapres=imresize(RAW_T2star_post.data(:,:,:,i),[final_res final_res],'bilinear');
                    wrongpix_pre=zeros(final_res,final_res,1,1);  
                    wrongpix_post=zeros(final_res,final_res,1,1); 
                end

                index_avant=find(r2savant<10e-2);% looking for the T2*pre < 100µs
                warning off					%To avoid division by zero warning
                r2savant=1./r2savant;       % R2*pre (ms-1)
                warning on
                r2savant(index_avant)=0;

                index_apres=find(r2sapres<10e-2);% olooking for the T2*post  < 100µs
                warning off					%To avoid division by zero warning
                r2sapres=1./r2sapres;	% R2*post (ms-1)
                warning on %#ok<*WNON>
                r2sapres(index_apres)=0;
                
                deltaR2star_map.reco.data(:,:,1,deltaR2star_slice_nbr)=r2sapres-r2savant;
                index_nan=find(RAW_T2star_pre.data(:,:,1,j)==NaN); %#ok<FNAN>
                deltaR2star_map.reco.data(index_nan)=0;
                %Compute the error map
                %proc.err=;
                
                %Compute the wrongpix map
                if rescale == 0
                    wrongpix_pre(find(RAW_T2star_pre.wrongpix(:,:,j)>0)) = 1;%#ok<FNDSB>
                    wrongpix_post(find(RAW_T2star_post.wrongpix(:,:,j)>0)) = 1;%#ok<FNDSB>
                else
                    wrongpix_pre(find(imresize(RAW_T2star_pre.wrongpix(:,:,j),[final_res final_res],'bilinear')>0)) = 1;%#ok<FNDSB>
                    wrongpix_post(find(imresize(RAW_T2star_post.wrongpix(:,:,j),[final_res final_res],'bilinear')>0)) = 1;%#ok<FNDSB>    
                end
                deltaR2star_map.reco.wrongpix(:,:,1,deltaR2star_slice_nbr)=or(wrongpix_pre,wrongpix_post) ; %Propagation of existing wronpix %#ok<FNDSB>
                deltaR2star_map.reco.wrongpix(index_avant)=1;						%Pixels where r2savant has been set to zero
                deltaR2star_map.reco.wrongpix(index_apres)=1;						%Pixels where r2sapres has been set to zero
                deltaR2star_map.reco.wrongpix(index_nan)=1;						%Pixels where proc.data=NaN
                
                % Update the deltaR2star_map structure
                deltaR2star_map.reco.fov_offsets(:,1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.fov_offsets(:,1,i,1);
                deltaR2star_map.reco.fov_orientation(:,1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.fov_orientation(:,1,i,1);
                deltaR2star_map.reco.label(1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.label(1,i,1);
%                 deltaR2star_map.reco.phaselabel(1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.phaselabel(1,i,1);
                deltaR2star_map.reco.fov_phase_orientation(1,deltaR2star_slice_nbr,1) = RAW_T2star_pre.fov_phase_orientation(1,i,1);
                
            end
        end
    end
    if deltaR2star_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the deltaR2star map because there is\n##$ no slice offset match between\n##$T2*pre=%s\n##$ and \n##$T2*post=%s',...
           filename_pre,filename_post);
        msgbox(warning_text, 'deltaR2star map warning') ;
        return
    end
    deltaR2star_map.reco.no_slices=deltaR2star_slice_nbr;
else
    %Initialisation des variables
    deltaR2star_map.reco='';
    if rescale == 0
        r2savant=RAW_T2star_pre.data(:,:,1,:);
        r2sapres=RAW_T2star_post.data(:,:,1,:);
        deltaR2star_map.reco.data=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,RAW_T2star_pre.no_slices);     %paramètres (DeltaR2*)
        deltaR2star_map.reco.wrongpix=zeros(RAW_T2star_pre.no_samples,RAW_T2star_pre.no_views,1,RAW_T2star_pre.no_slices);  %carte de pixel à exclure
    else
        r2savant=imresize(RAW_T2star_pre.data(:,:,:,:),[final_res final_res],'bilinear');
        r2sapres=imresize(RAW_T2star_post.data(:,:,:,:),[final_res final_res],'bilinear');
        deltaR2star_map.reco.data=zeros(final_res,final_res,1,RAW_T2star_pre.no_slices);     %paramètres (DeltaR2*)
        deltaR2star_map.reco.wrongpix=zeros(final_res,final_res,1,RAW_T2star_pre.no_slices);  %carte de pixel à exclure
    end
    
    index_avant=find(r2savant<10e-2);% on cherche les T2* < 100µs
    warning off					%To avoid division by zero warning
    r2savant=1./r2savant;	% R2*avant en ms-1
    warning on
    r2savant(index_avant)=0;
    
    index_apres=find(r2sapres<10e-2);% on cherche les T2* < 100µs
    warning off					%To avoid division by zero warning
    r2sapres=1./r2sapres;	% R2*apres en ms-1
    warning on %#ok<*WNON>
    r2sapres(index_apres)=0;
    
    deltaR2star_map.reco.data(:,:,1,:)=r2sapres-r2savant;
    index_nan=find(RAW_T2star_pre.data==NaN); %#ok<FNAN>
    deltaR2star_map.reco.data(index_nan)=0;
    %Compute the error map
    %proc.err=;
    
    %Compute the wrongpix map
    if rescale == 0
        deltaR2star_map.reco.wrongpix(find(RAW_T2star_pre.wrongpix>0)) = 1;%#ok<FNDSB>
        deltaR2star_map.reco.wrongpix(find(RAW_T2star_post.wrongpix>0)) = 1;%#ok<FNDSB>
    else
        for i = 1:RAW_T2star_pre.no_slices
            wrongpix_pre = imresize(RAW_T2star_pre.wrongpix(:,:,i),[final_res final_res],'bilinear');
            wrongpix_pre(find(wrongpix_pre>0)) = 1;%#ok<FNDSB>
            index_nan_pre=isnan(wrongpix_pre);
            wrongpix_pre(index_nan_pre)=0;
             wrongpix_post = imresize(RAW_T2star_post.wrongpix(:,:,i),[final_res final_res],'bilinear');
            wrongpix_post(find(wrongpix_post>0)) = 1;%#ok<FNDSB>
            index_nan_post=isnan(wrongpix_post); 
            wrongpix_post(index_nan_post)=0;
            deltaR2star_map.reco.wrongpix(:,:,1,i)=or(wrongpix_pre,wrongpix_post) ;
        end
    end

    deltaR2star_map.reco.wrongpix(index_avant)=1;						%Pixels where r2savant has been set to zero
    deltaR2star_map.reco.wrongpix(index_apres)=1;						%Pixels where r2sapres has been set to zero
    deltaR2star_map.reco.wrongpix(index_nan)=1;						%Pixels where proc.data=NaN
    
    deltaR2star_map.reco.no_slices=RAW_T2star_pre.no_slices;
    % Adapt the fov offsets and orientations infos
    deltaR2star_map.reco.fov_offsets=RAW_T2star_pre.fov_offsets(:,1,:,1);
    deltaR2star_map.reco.fov_orientation=RAW_T2star_pre.fov_orientation(:,1,:,1);
    deltaR2star_map.reco.label=RAW_T2star_pre.label(1,:,1);
    deltaR2star_map.reco.phaselabel=RAW_T2star_pre.phaselabel(1,:,1);
    deltaR2star_map.reco.fov_phase_orientation=RAW_T2star_pre.fov_phase_orientation(1,:,1);
end

%Set the dimensions of the proc structure
if rescale == 0
    deltaR2star_map.reco.no_samples=RAW_T2star_pre.no_samples;
    deltaR2star_map.reco.no_views=RAW_T2star_pre.no_views;
else
    deltaR2star_map.reco.no_samples=final_res;
    deltaR2star_map.reco.no_views=final_res;
end
deltaR2star_map.reco.no_echoes=1; %Number of parameters stored here
deltaR2star_map.reco.no_expts=1;

%Complete the reco structure
deltaR2star_map.reco.texte='DeltaR2*';
deltaR2star_map.reco.unit='ms-1';
deltaR2star_map.reco.date=date;
deltaR2star_map.reco.displayedecho=1;
deltaR2star_map.reco.displayedslice=1;
deltaR2star_map.reco.displayedexpt=1;
deltaR2star_map.reco_number = 1;
deltaR2star_map.scan_number = 104;
deltaR2star_map.clip=[0 0.2 1];
deltaR2star_map.reco.globalmax = max(deltaR2star_map.reco.data(:));
deltaR2star_map.reco.globalmin = min(deltaR2star_map.reco.data(:));
deltaR2star_map.reco.thickness = data_pre.uvascim.image.reco.thickness;
deltaR2star_map.reco=orderfields(deltaR2star_map.reco);
if iscell(RAW_T2star_pre.paramQuantif)
    RAW_T2star_pre.paramQuantif = RAW_T2star_pre.paramQuantif{1};
end
if iscell(RAW_T2star_post.paramQuantif)
    RAW_T2star_post.paramQuantif = RAW_T2star_post.paramQuantif{1};
end
ParamConfig=sprintf('##$QuantifMethod=''(1/T2*post)-(1/T2*pre)''\n##$final_res=%s\n##$T2*pre=%s\n##$T2*post=%s\n\n##$T2*pre scan info\n%s\n\n##$T2*post scan info\n%s',...
    add_parameters{:}{1},filename_pre,filename_post, RAW_T2star_pre.paramQuantif, RAW_T2star_post.paramQuantif);
deltaR2star_map.reco.paramQuantif = ParamConfig;
deltaR2star_map.reco=orderfields(deltaR2star_map.reco);


%complete the structure
deltaR2star_map.acq = data_pre.uvascim.image.acq;
deltaR2star_map.filename = data_pre.uvascim.image.filename;


