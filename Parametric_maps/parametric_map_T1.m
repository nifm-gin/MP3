function [T1map,M0map,IEmap] = parametric_map_T1(MTI_map_filename, Mask_filename, add_parameters)
% generate a T1_map from a RAW MTI-T1 scan (FAIR-EPI kind)
% this code come from the the ASL module coded by C. Debacker

% retieve calculation parameter
first_scan=str2double(add_parameters{:}(1));
M0mapYN=add_parameters{:}(2);
M0map = [];
IEmapYN=add_parameters{:}(3);
IEmap = [];

%% load the MTI file
fid=fopen(MTI_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(MTI_map_filename);
    data = data.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the T1 map because there is\n##$ Somthing wrong with the data\n##$MTI=%s\n##$',...
        MTI_map_filename);
    msgbox(warning_text, 'T1 map warning') ;
    T1map = [];
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

% apply mask
if ~isempty(Mask)
    [~, ~, E, Z] = size(data.reco.data);
    for i = 1:Z
        for j = 1:numel(Mask)
            if abs(data.reco.fov_offsets(3,1,i) - Mask(j).fov_offsets(3)) < 1e-5
                 data.reco.data(:,:,:,i) = data.reco.data(:,:,:,i).*repmat(Mask(j).value,[1 1 E 1]);           
            end
        end
    end
end


%% retreive FAIR method
T1map_method = scan_acqp('##$Method=',data.texte,0);
if(~isempty(regexpi(T1map_method,'\w*fair\w*')))
    PvVersion = scan_acqp('##TITLE=Parameter List, ParaVision ',data.texte,0);
    if ~isnan(PvVersion(1)) && ~isempty(regexpi(PvVersion,'\w*6.\w*'))
        FairMode=scan_acqp('##$PVM_FairMode=',data.texte,0); % SELECTIVE, NONSELECTIVE, INTERLEAVED or INTERLEAVED2
    else
        FairMode=scan_acqp('##$FairMode=',data.texte,0); % SELECTIVE, NONSELECTIVE, INTERLEAVED or INTERLEAVED2
    end
    if(max(strcmp(FairMode,{'SELECTIVE','NONSELECTIVE'}))==1)
        NumFairMode = 1;
    elseif(max(strcmp(FairMode,{'INTERLEAVED','INTERLEAVED2'}))==1)
        NumFairMode = 2;
    end
else
    disp('scan not implemented')
end

%% slice acquisition order for correct inversion time
ObjOrderListName = sprintf('##$PVM_ObjOrderList=( %d )',size(data.reco.data,4));
PVM_ObjOrderList=scan_acqp(ObjOrderListName,data.texte,1)';
if(size(data.reco.data,4)==1)
    PVM_ObjOrderList = 0;
end

%% retrieve inversion time and inter slice time for different acquisition
% module of FAIR method (EPI, RARE, fisp, segm...)
if ( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ~isempty(regexpi(T1map_method,'\w*epi\w*')) )
    % inter slice time
    interSliceTime = InterSliceTimeFairEpi(data);
    
    % inversion time
    FairTIR_Arr=(scan_acqp('FairTIR_Arr=',data.texte,1))';
    InvTimeRaw = FairTIR_Arr(first_scan:end);   % liste  des temps inversion pour la methode fair
    InvTimeRaw = repmat(InvTimeRaw,NumFairMode,1);            % premiere ligne selective deuxieme nonselective si les deux mesures sont faites
elseif( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ( ~isempty(regexpi(T1map_method,'\w*segm\w*')) ||...
        ~isempty(regexpi(T1map_method,'\w*fisp\w*')) ) )
    % inversion time
    InvTimeSel = (scan_acqp('##$InvTimeSel=',data.texte,1))';
    InvTimeGlo=(scan_acqp('##$InvTimeGlob=',data.texte,1))';
    if(strcmp(FairMode,'SELECTIVE'))
        InvTimeRaw=InvTimeSel(first_scan:end);
    elseif(strcmp(FairMode,'NONSELECTIVE'))
        InvTimeRaw=InvTimeGlo(first_scan:end);
    elseif(strcmp(FairMode,'INTERLEAVED'))
        InvTimeRaw(1,:)=InvTimeSel(first_scan:end);
        InvTimeRaw(2,:)=InvTimeGlo(first_scan:end);
    end
    % inter slice time
    interSliceTime = 50; % ms
    fprintf('Correction temps interslice par defaut egale a %ims\n',interSliceTime)
end

%% Init variable
T1map.acq=data.acq;
T1map.filename=data.filename;
T1map.texte=data.texte;
if strcmp(M0mapYN,'Yes')
    M0map.acq=data.acq;
    M0map.filename=data.filename;
    M0map.texte=data.texte;
end
if strcmp(IEmapYN,'Yes')
    IEmap.acq=data.acq;
    IEmap.filename=data.filename;
    IEmap.texte=data.texte;
end
if(NumFairMode==1)
    data_in_vector = reshape(data.reco.data, [size(data.reco.data,1)*size(data.reco.data,2),...
        size(data.reco.data,4), size(data.reco.data,5)]);
elseif(NumFairMode==2)
    data_in_vector = reshape(data.reco.data, [size(data.reco.data,1)*size(data.reco.data,2),...
        size(data.reco.data,4), size(data.reco.data,5)/NumFairMode,NumFairMode]);
    data_in_vector(:,:,:,1) = reshape(data.reco.data(:,:,:,:,NumFairMode*first_scan-1:NumFairMode:end), [size(data.reco.data,1)*size(data.reco.data,2),...
        size(data.reco.data,4), size(data.reco.data,5)/NumFairMode]);
    data_in_vector(:,:,:,2) = reshape(data.reco.data(:,:,:,:,NumFairMode*first_scan:NumFairMode:end), [size(data.reco.data,1)*size(data.reco.data,2),...
        size(data.reco.data,4), size(data.reco.data,5)/NumFairMode]);
end
fit_T1_result = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
fit_T1_err= NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
if strcmp(M0mapYN,'Yes')
    fit_M0_result = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
    fit_M0_err = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
end
if strcmp(IEmapYN,'Yes')
    fit_IE_result = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
    fit_IE_err = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
end

%% Fitting
for it_mode = 1 : NumFairMode
    for it_slice = 1 : size(data.reco.data,4)
        data_in_vectorTemp = squeeze(data_in_vector(:,it_slice,:,it_mode));
        %     indT1 = find((data_in_vector(:,it_slice,1) ~= 0) & (isnan(data_in_vector(:,it_slice,1)) ~= 1));
        %     data_in_vectorTemp = squeeze(data_in_vector(indT1,it_slice,:));
        NumSlice = find(PVM_ObjOrderList==it_slice-1)-1;
        
        % retrieve parameters
        echotime_used = InvTimeRaw(it_mode,:) + NumSlice * interSliceTime;
        
        %         for voxel_nbr=1:size(data_in_vectorTemp,1)
        for voxel_nbr=1:size(data_in_vectorTemp,1)
            if sum(isnan(data_in_vectorTemp(voxel_nbr,:))) ~= length(data_in_vectorTemp(voxel_nbr,:))
                [~,voxel_min_nbr]=min(data_in_vectorTemp(voxel_nbr,:));
                [aaa, bbb,  ccc]=levenbergmarquardt('fit_T1_3param',echotime_used,data_in_vectorTemp(voxel_nbr,:),[echotime_used(voxel_min_nbr)/0.693*1.2 max(data_in_vectorTemp(voxel_nbr,:)) 1]);
                if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 & ccc==-1 %#ok<AND2>
                    fit_T1_result(voxel_nbr,it_slice,it_mode)=aaa(1);
                    fit_T1_err(voxel_nbr,it_slice,it_mode)=bbb(1);
                    if strcmp(M0mapYN,'Yes')
                        fit_M0_result(voxel_nbr,it_slice,it_mode)=aaa(2);
                        fit_M0_err(voxel_nbr,it_slice,it_mode)=bbb(2);
                    end
                    if strcmp(IEmapYN,'Yes')
                        fit_IE_result(voxel_nbr,it_slice,it_mode)=aaa(3);
                        fit_IE_err(voxel_nbr,it_slice,it_mode)=bbb(3);                        
                    end
                end
            end
        end
        
    end
end

%% reshape and formatting of data
tmp=reshape(fit_T1_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4),NumFairMode]);
T1map.reco.data=permute(tmp, [1 2 4 3]);
% correction Kober 2004 on apparent T1
if ( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ( ~isempty(regexpi(T1map_method,'\w*segm\w*')) ||...
        ~isempty(regexpi(T1map_method,'\w*fisp\w*')) ) )
    PVM_ExcPulseAngle=scan_acqp('##$PVM_ExcPulseAngle=',data.texte,1);
    Seg_time=scan_acqp('##$Seg_time=',data.texte,1);
%    T1map.reco.data = T1map.reco.data ./ ( 1 + T1map.reco.data * log( cos( PVM_ExcPulseAngle * pi / 180 ) ) / Seg_time );
end
tmp=reshape(fit_T1_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4),NumFairMode]);
T1map.reco.err=permute(tmp, [1 2 4 3]);
if strcmp(M0mapYN,'Yes')
    tmp=reshape(fit_M0_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4),NumFairMode]);
    M0map.reco.data=permute(tmp, [1 2 4 3]);
    tmp=reshape(fit_M0_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4),NumFairMode]);
    M0map.reco.err=permute(tmp, [1 2 4 3]);
end
if strcmp(IEmapYN,'Yes')
    tmp=reshape(fit_IE_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4),NumFairMode]);
    IEmap.reco.data=permute(tmp, [1 2 4 3]);
    tmp=reshape(fit_IE_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4),NumFairMode]);
    IEmap.reco.err=permute(tmp, [1 2 4 3]);
end

%% T1 map formating
T1map.reco.echo_label(1,1) = {'fit_T1'};
T1map.reco.unit(1,1) = {'ms'};
if(NumFairMode==2)
    T1map.reco.echo_label(1,1) = {'fit_T1 selective'};
    T1map.reco.echo_label(1,2) = {'fit_T1 non selective'};
    T1map.reco.unit(1,2) = {'ms'};
end
T1map.reco.texte = 'T1map';
T1map.reco.date = date;
T1map.reco.no_echoes = size(T1map.reco.echo_label,2);
T1map.reco.no_expts  = size(T1map.reco.data,5);
T1map.reco.no_slices = T1map.acq.no_slices;
T1map.reco.globalmin=min(T1map.reco.data(:));
T1map.reco.globalmax=max(T1map.reco.data(:));
for m_expt=1:T1map.reco.no_expts,
    for m_slice=1:T1map.reco.no_slices,
        for m_echo=1:T1map.reco.no_echoes
            T1map.reco.fov_offsets(:,m_echo,m_slice,m_expt) = data.reco.fov_offsets(:,1,m_slice,m_expt);
            T1map.reco.fov_orientation(:,m_echo,m_slice,m_expt) = data.reco.fov_orientation(:,1,m_slice,m_expt);
            T1map.reco.label(m_echo,m_slice,m_expt) = data.reco.label(1,m_slice,m_expt);
            T1map.reco.phaselabel(m_echo,m_slice,m_expt) = data.reco.phaselabel(1,m_slice,m_expt);
            T1map.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = data.reco.fov_phase_orientation(1,m_slice,m_expt);
            T1map.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
            T1map.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
        end
    end
end
T1map.reco = orderfields(T1map.reco);

T1map.reco.displayedecho=T1map.reco.no_echoes;
T1map.reco.displayedslice=data.reco.displayedslice;
T1map.reco.displayedexpt=1;
T1map.reco.thickness=data.reco.thickness;
T1map.reco.no_views=data.reco.no_views;
T1map.reco.no_samples=data.reco.no_samples;
T1map.reco.angAP=data.reco.angAP;
T1map.reco.angFH=data.reco.angFH;
T1map.reco.angRL=data.reco.angRL;
T1map.reco.angulation=data.reco.angulation;
T1map.reco.bitpix=data.reco.bitpix;
T1map.reco.fov=data.reco.fov;

ParamConfig=sprintf('##$QuantifMethod=levenbergmarquardt with function AB_expt1_v1.m\n##$First scan used=%s\n##$Raw scan used=%s\n##$interSliceTime=%0.2f\n##END=',...
    add_parameters{:}{1},...
    MTI_map_filename,...
    interSliceTime);
T1map.reco.paramQuantif = ParamConfig;
T1map.reco=orderfields(T1map.reco);

T1map.scan_number = 104;
T1map.reco_number = 1;
T1map.clip=[0 3500 1];

%% M0 map formating
if strcmp(M0mapYN,'Yes')
    M0map.reco.echo_label(1,1) = {'fit_M0'};
    M0map.reco.unit(1,1) = {'ms'};
    if(NumFairMode==2)
        M0map.reco.echo_label(1,1) = {'fit_M0 selective'};
        M0map.reco.echo_label(1,2) = {'fit_M0 non selective'};
        M0map.reco.unit(1,2) = {'a.u.'};
    end
    M0map.reco.texte = 'M0map';
    M0map.reco.date = date;
    M0map.reco.no_echoes = size(M0map.reco.echo_label,2);
    M0map.reco.no_expts  = size(M0map.reco.data,5);
    M0map.reco.no_slices = M0map.acq.no_slices;
    M0map.reco.globalmin=min(M0map.reco.data(:));
    M0map.reco.globalmax=max(M0map.reco.data(:));
    for m_expt=1:M0map.reco.no_expts,
        for m_slice=1:M0map.reco.no_slices,
            for m_echo=1:M0map.reco.no_echoes
                M0map.reco.fov_offsets(:,m_echo,m_slice,m_expt) = data.reco.fov_offsets(:,1,m_slice,m_expt);
                M0map.reco.fov_orientation(:,m_echo,m_slice,m_expt) = data.reco.fov_orientation(:,1,m_slice,m_expt);
                M0map.reco.label(m_echo,m_slice,m_expt) = data.reco.label(1,m_slice,m_expt);
                M0map.reco.phaselabel(m_echo,m_slice,m_expt) = data.reco.phaselabel(1,m_slice,m_expt);
                M0map.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = data.reco.fov_phase_orientation(1,m_slice,m_expt);
                M0map.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
                M0map.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
            end
        end
    end
    M0map.reco = orderfields(M0map.reco);
    
    M0map.reco.displayedecho=M0map.reco.no_echoes;
    M0map.reco.displayedslice=M0map.reco.no_slices;
    M0map.reco.displayedexpt=1;
    M0map.reco.thickness=data.reco.thickness;
    M0map.reco.no_views=data.reco.no_views;
    M0map.reco.no_samples=data.reco.no_samples;
    M0map.reco.angAP=data.reco.angAP;
    M0map.reco.angFH=data.reco.angFH;
    M0map.reco.angRL=data.reco.angRL;
    M0map.reco.angulation=data.reco.angulation;
    M0map.reco.bitpix=data.reco.bitpix;
    M0map.reco.fov=data.reco.fov;
    
    
    M0map.reco.paramQuantif = ParamConfig;
    M0map.reco=orderfields(M0map.reco);
    
    M0map.scan_number = 105;
    M0map.reco_number = 1;
end

%% IE map formating
if strcmp(IEmapYN,'Yes')
    IEmap.reco.echo_label(1,1) = {'fit_IE'};
    IEmap.reco.unit(1,1) = {'ms'};
    if(NumFairMode==2)
        IEmap.reco.echo_label(1,1) = {'fit_IE selective'};
        IEmap.reco.echo_label(1,2) = {'fit_IE non selective'};
        IEmap.reco.unit(1,2) = {'a.u.'};
    end
    IEmap.reco.texte = 'IEmap';
    IEmap.reco.date = date;
    IEmap.reco.no_echoes = size(IEmap.reco.echo_label,2);
    IEmap.reco.no_expts  = size(IEmap.reco.data,5);
    IEmap.reco.no_slices = IEmap.acq.no_slices;
    IEmap.reco.globalmin=min(IEmap.reco.data(:));
    IEmap.reco.globalmax=max(IEmap.reco.data(:));
    for m_expt=1:IEmap.reco.no_expts,
        for m_slice=1:IEmap.reco.no_slices,
            for m_echo=1:IEmap.reco.no_echoes
                IEmap.reco.fov_offsets(:,m_echo,m_slice,m_expt) = data.reco.fov_offsets(:,1,m_slice,m_expt);
                IEmap.reco.fov_orientation(:,m_echo,m_slice,m_expt) = data.reco.fov_orientation(:,1,m_slice,m_expt);
                IEmap.reco.label(m_echo,m_slice,m_expt) = data.reco.label(1,m_slice,m_expt);
                IEmap.reco.phaselabel(m_echo,m_slice,m_expt) = data.reco.phaselabel(1,m_slice,m_expt);
                IEmap.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = data.reco.fov_phase_orientation(1,m_slice,m_expt);
                IEmap.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
                IEmap.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
            end
        end
    end
    IEmap.reco = orderfields(IEmap.reco);
    
    IEmap.reco.displayedecho=IEmap.reco.no_echoes;
    IEmap.reco.displayedslice=IEmap.reco.no_slices;
    IEmap.reco.displayedexpt=1;
    IEmap.reco.thickness=data.reco.thickness;
    IEmap.reco.no_views=data.reco.no_views;
    IEmap.reco.no_samples=data.reco.no_samples;
    IEmap.reco.angAP=data.reco.angAP;
    IEmap.reco.angFH=data.reco.angFH;
    IEmap.reco.angRL=data.reco.angRL;
    IEmap.reco.angulation=data.reco.angulation;
    IEmap.reco.bitpix=data.reco.bitpix;
    IEmap.reco.fov=data.reco.fov;
    IEmap.reco.paramQuantif = ParamConfig;
    IEmap.reco=orderfields(IEmap.reco);
    
    IEmap.scan_number = 105;
    IEmap.reco_number = 1;
end