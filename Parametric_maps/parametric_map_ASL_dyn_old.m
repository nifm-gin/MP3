function [Transitmap,M0map,CBFmap,T1map] = parametric_map_ASL_dyn_old(ASL_map_filename, VOI_filename, add_parameters)
% generate a Transit_time_map, CBF map, M0 map, T1 map from a dynamic pCASL scan
% this code was coded by Lydiane Hirschler

% retrieve calculation parameter
%first_scan=str2double(add_parameters{:}(1));
M0map = [];
CBFmap = [];
T1map = [];

lambda = 0.9;

%T = str2double(add_parameters{:}(7)); %s
%step = str2double(add_parameters{:}(8));%s
%nb_cycles = str2double(add_parameters{:}(9));
alpha = str2double(add_parameters{:}(2));
T1blood= str2double(add_parameters{:}(3));
%CBFi = str2double(add_parameters{:}(4));
%T1ti = str2double(add_parameters{:}(5));%s


%% load the ASL file
fid=fopen(ASL_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(ASL_map_filename);
    data = data.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the maps because there is\n##$ Something wrong with the data\n##$ASL=%s\n##$',...
        ASL_map_filename);
    msgbox(warning_text, 'Dasl map warning') ;
    Transitmap = [];
    return
end

[Voi, VoiYn] = extract_map(VOI_filename,'ASL-dyn','VOI');
if isempty(Voi)
    VoiYn=0;
end
%% Lecture parametres fichier methode
step = scan_acqp('##$PCASL_TempResolution=',data.texte,1)/1000;
T = scan_acqp('##$PCASL_TotalLabelTime=',data.texte,1)/1000;

if isnan(step)
    subboliDuration = scan_acqp('PCASL_SubboliTime=',data.texte,1)/1000;
    step = subboliDuration(1);
    T = scan_acqp('PCASL_LabelTime=',data.texte,1)/1000;
end

%% Init variable
Transitmap.acq=data.acq;
Transitmap.filename=data.filename;
Transitmap.texte=data.texte;

M0map.acq=data.acq;
M0map.filename=data.filename;
M0map.texte=data.texte;

CBFmap.acq=data.acq;
CBFmap.filename=data.filename;
CBFmap.texte=data.texte;

T1map.acq=data.acq;
T1map.filename=data.filename;
T1map.texte=data.texte;


for it_sl = 1 : size(data.reco.data,4)
    if VoiYn
        factorResize = size(data.reco.data,1) / size(Voi(it_sl).value,1);
        VoiTemp = imresize(Voi(it_sl).value,factorResize,'nearest');
        % VoiTemp(VoiTemp==0)=NaN;
        data.reco.data(:,:,:,it_sl,:) = data.reco.data(:,:,:,it_sl,:).*repmat(VoiTemp,[1 1 size(data.reco.data,3) 1 size(data.reco.data,5)]);
        % data.reco.data(data.reco.data(:,:,:,it_sl,:) == 0) = NaN;
    end
end
dummyScans = size(data.reco.data,3)-numel(subboliDuration);
data_in_vector = reshape(data.reco.data, [size(data.reco.data,1)*size(data.reco.data,2),...
    size(data.reco.data,3), size(data.reco.data,4),size(data.reco.data,5)]);

fit_transit_result = NaN(size(data_in_vector,1),size(data_in_vector,3));
fit_transit_err= NaN(size(data_in_vector,1),size(data_in_vector,3));

fit_M0_result = NaN(size(data_in_vector,1),size(data_in_vector,3));
fit_M0_err = NaN(size(data_in_vector,1),size(data_in_vector,3));

fit_CBF_result = NaN(size(data_in_vector,1),size(data_in_vector,3));
fit_CBF_err = NaN(size(data_in_vector,1),size(data_in_vector,3));
fit_T1_result = NaN(size(data_in_vector,1),size(data_in_vector,3));
fit_T1_err = NaN(size(data_in_vector,1),size(data_in_vector,3));


%%
%% FITTING

for it_slice = 1 : size(data.reco.data,4)
%     %%TESTS
%     load('DASLdata.mat')
%     data_in_vector=y_data;
    
    %% A decommenter si donnees a l envers ou plusieurs cycles en une repetition
    %         data_in_vectorTemp_endroit = squeeze(data_in_vector(:,2:end,it_slice));
    %         data_in_vectorTemp_moy=(data_in_vectorTemp_endroit(:,1:T/step) + data_in_vectorTemp_endroit(:,T/step+1:nb_cycles*T/step))/nb_cycles; %moyenne sur les cycles
    %         for i = 1:T/step;
    %             data_in_vectorTemp(:,i) = data_in_vectorTemp_moy(:,T/step-i+1); %donnees dyn lues a l'envers
    %         end
    %%
    data_in_vectorTemp_cycles = data_in_vector(:,dummyScans+1:end,it_slice,:); % ignorer les premieres images (dummy scans "simules")
    data_in_vectorTemp = squeeze(mean(data_in_vectorTemp_cycles,4)); %moyenne sur les repetitions (ie. cycles)
    
    %     indT1 = find((data_in_vector(:,it_slice,1) ~= 0) & (isnan(data_in_vector(:,it_slice,1)) ~= 1));
    %     data_in_vectorTemp = squeeze(data_in_vector(indT1,it_slice,:));
    %NumSlice = find(PVM_ObjOrderList==it_slice-1)-1;
    
    % retrieve parameters
    %        echotime_used = InvTimeRaw(it_mode,:) + NumSlice * interSliceTime;
    N = size(data_in_vectorTemp,1);
    parfor_progress(N); % Initialize
    parfor_progress;
    %         for voxel_nbr=1:size(data_in_vectorTemp,1)
    for voxel_nbr=1:N
        %for voxel_nbr=1:N
        % Count
        % fprintf('voxel : %d / %d \n', voxel_nbr,N )
        if sum(isnan(data_in_vectorTemp(voxel_nbr,:))) ~= length(data_in_vectorTemp(voxel_nbr,:)) && sum(data_in_vectorTemp(voxel_nbr,:) == 0) ~= length(data_in_vectorTemp(voxel_nbr,:))
            
            %% FILTRAGE
            %                 dataFT = fft(data_in_vectorTemp(voxel_nbr,:));
            %
            %                 for i=3:2:T/step
            %                     dataFT(1,i) = 0;
            %                 end
            %
            %                 data_in_vectorTemp_filtre(voxel_nbr,:) = ifft(dataFT(1,:));
            data_in_vectorTemp_filtre(voxel_nbr,:) = data_in_vectorTemp(voxel_nbr,:);
            %% REDRESSEMENT COURBE
            debut = mean(data_in_vectorTemp_filtre(voxel_nbr,1:3));
            fin = mean(data_in_vectorTemp_filtre(voxel_nbr,T/step-3 : end));
            pente = (fin - debut)/(T/step);
            for point = 1:T/step
                data_in_vectorTemp_filtre_redres(voxel_nbr,point) = data_in_vectorTemp_filtre(voxel_nbr,point) - pente * point;
            end
            %      data_in_vectorTemp_filtre_redres(voxel_nbr,:) = data_in_vectorTemp_filtre(voxel_nbr,:);
            % Conditions initiales fit
            [minx, miny] = min(data_in_vectorTemp_filtre_redres(voxel_nbr,:)); %temps de transit
            Ttransiti = (miny*step-T/2);
            maxData = max(data_in_vectorTemp_filtre_redres(voxel_nbr,:));
            cstT1 = 0.63*(debut - minx) + minx;
            [~, abscisseT1] = min(abs(data_in_vectorTemp_filtre_redres(voxel_nbr,miny:end)- cstT1));
            T1ti = abscisseT1  * step + Ttransiti;
            %A =  -1/(abscisseT1  * step);
            %T1ti = -lambda / (lambda * A + CBFi);
            CBFi = (maxData - minx) *lambda / (2* maxData * T1ti * alpha * exp(-Ttransiti/T1blood) *(1- exp(-T/(2*T1ti))));
            A= -1/T1ti- CBFi/lambda ;
            eATs2=expm(A*T/2);
            coeff = -2 * CBFi/lambda *alpha*exp(-Ttransiti/T1blood);
            
            Finit=eATs2*(eATs2-1)*coeff/A; %A negatif donc signe inverse et A=-1/T1sat
            Xfinit=Finit/(1-expm(A*T)); %T=2*delta (publi)
            M0init = maxData / (1 + Xfinit);
            
            % FIT
            [aaa, bbb,  ccc]=levenbergmarquardtASL('adpki2',[lambda alpha T step T/step T1blood],data_in_vectorTemp_filtre_redres(voxel_nbr,:),[M0init CBFi Ttransiti T1ti]);
            % fprintf('CBF= %f, transit = %f, T1 = %f \n',aaa(2)*6000,aaa(3),aaa(4));
            if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 & ccc==-1 %#ok<AND2>
                fit_transit_result(voxel_nbr,it_slice)=aaa(3);
                fit_transit_err(voxel_nbr,it_slice)=bbb(3);
                
                fit_M0_result(voxel_nbr,it_slice)=aaa(1);
                fit_M0_err(voxel_nbr,it_slice)=bbb(1);
                
                fit_CBF_result(voxel_nbr,it_slice)=aaa(2);
                fit_CBF_err(voxel_nbr,it_slice)=bbb(2);
                fit_T1_result(voxel_nbr,it_slice)=aaa(4);
                fit_T1_err(voxel_nbr,it_slice)=bbb(4);
                
            end
            
        end
    end
    parfor_progress(0); % Clean up
end
%end

%% reshape and formatting of data
tmp=reshape(fit_transit_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%Transitmap.reco.data=permute(tmp, [1 2 4 3]);
Transitmap.reco.data=tmp*1000;
%Transitmap.reco.data=interpolation_Nans(Transitmap.reco.data, VoiTemp);

tmp=reshape(fit_transit_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%Transitmap.reco.err=permute(tmp, [1 2 4 3]);
Transitmap.reco.err=tmp*1000;

tmp=reshape(fit_M0_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%M0map.reco.data=permute(tmp, [1 2 4 3]);
M0map.reco.data=tmp;
%M0map.reco.data=interpolation_Nans(M0map.reco.data, VoiTemp);

tmp=reshape(fit_M0_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%M0map.reco.err=permute(tmp, [1 2 4 3]);
M0map.reco.err=tmp;

tmp=reshape(fit_CBF_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%CBFmap.reco.data=permute(tmp, [1 2 4 3]);
CBFmap.reco.data=tmp*6000;
%CBFmap.reco.data=interpolation_Nans(CBFmap.reco.data, VoiTemp);
tmp=reshape(fit_CBF_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%CBFmap.reco.err=permute(tmp, [1 2 4 3]);
CBFmap.reco.err=tmp*6000;

tmp=reshape(fit_T1_result,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%CBFmap.reco.data=permute(tmp, [1 2 4 3]);
T1map.reco.data=tmp;
%T1map.reco.data=interpolation_Nans(T1map.reco.data, VoiTemp);
tmp=reshape(fit_T1_err,[size(data.reco.data,1),size(data.reco.data,2),size(data.reco.data,4)]);
%CBFmap.reco.err=permute(tmp, [1 2 4 3]);
T1map.reco.err=tmp;


%% T1 map formating
Transitmap.reco.echo_label(1,1) = {'fit_transit'};
Transitmap.reco.unit(1,1) = {'ms'};
% if(NumFairMode==2)
%     Transitmap.reco.echo_label(1,1) = {'fit_T1 selective'};
%     Transitmap.reco.echo_label(1,2) = {'fit_T1 non selective'};
%     Transitmap.reco.unit(1,2) = {'ms'};
% end
Transitmap.reco.texte = 'TransitTimemap';
Transitmap.reco.date = date;
Transitmap.reco.no_echoes = 1;
Transitmap.reco.no_expts  = size(Transitmap.reco.data,5);
Transitmap.reco.no_slices = Transitmap.acq.no_slices;
Transitmap.reco.globalmin=min(Transitmap.reco.data(:));
Transitmap.reco.globalmax=max(Transitmap.reco.data(:));
for m_expt=1:Transitmap.reco.no_expts,
    for m_slice=1:Transitmap.reco.no_slices,
        for m_echo=1:Transitmap.reco.no_echoes
            Transitmap.reco.fov_offsets(:,m_echo,m_slice,m_expt) = data.reco.fov_offsets(:,1,m_slice,m_expt);
            Transitmap.reco.fov_orientation(:,m_echo,m_slice,m_expt) = data.reco.fov_orientation(:,1,m_slice,m_expt);
            Transitmap.reco.label(m_echo,m_slice,m_expt) = data.reco.label(1,m_slice,m_expt);
            Transitmap.reco.phaselabel(m_echo,m_slice,m_expt) = data.reco.phaselabel(1,m_slice,m_expt);
            Transitmap.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = data.reco.fov_phase_orientation(1,m_slice,m_expt);
            Transitmap.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
            Transitmap.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
        end
    end
end
Transitmap.reco = orderfields(Transitmap.reco);

Transitmap.reco.displayedecho=Transitmap.reco.no_echoes;
Transitmap.reco.displayedslice=Transitmap.reco.no_slices;
Transitmap.reco.displayedexpt=1;
Transitmap.reco.thickness=data.reco.thickness;
Transitmap.reco.no_views=data.reco.no_views;
Transitmap.reco.no_samples=data.reco.no_samples;

ParamConfig=sprintf('##$QuantifMethod=levenbergmarquardt with function AB_expt1_v1.m\n##$Raw scan used=%s\n##END=',...
    ASL_map_filename);
Transitmap.reco.paramQuantif = ParamConfig;
Transitmap.reco=orderfields(Transitmap.reco);

Transitmap.scan_number = 104;
Transitmap.reco_number = 1;
Transitmap.clip=[0 1000 1];

%% M0 map formating

M0map.reco.echo_label(1,1) = {'fit_M0'};
M0map.reco.unit(1,1) = {'ms'};
%     if(NumFairMode==2)
%         M0map.reco.echo_label(1,1) = {'fit_M0 selective'};
%         M0map.reco.echo_label(1,2) = {'fit_M0 non selective'};
%         M0map.reco.unit(1,2) = {'a.u.'};
%     end
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

M0map.reco.paramQuantif = ParamConfig;
M0map.reco=orderfields(M0map.reco);

M0map.scan_number = 105;
M0map.reco_number = 1;


%% CBF map formating

CBFmap.reco.echo_label(1,1) = {'fit_CBF'};
CBFmap.reco.unit = {'mL/100g/min'};
%     if(NumFairMode==2)
%         CBFmap.reco.echo_label(1,1) = {'fit_IE selective'};
%         CBFmap.reco.echo_label(1,2) = {'fit_IE non selective'};
%         CBFmap.reco.unit(1,2) = {'a.u.'};
%     end
CBFmap.reco.texte = 'CBFfitmap';
CBFmap.reco.date = date;
CBFmap.reco.no_echoes = size(CBFmap.reco.echo_label,2);
CBFmap.reco.no_expts  = size(CBFmap.reco.data,5);
CBFmap.reco.no_slices = CBFmap.acq.no_slices;
CBFmap.reco.globalmin=min(CBFmap.reco.data(:));
CBFmap.reco.globalmax=max(CBFmap.reco.data(:));
for m_expt=1:CBFmap.reco.no_expts,
    for m_slice=1:CBFmap.reco.no_slices,
        for m_echo=1:CBFmap.reco.no_echoes
            CBFmap.reco.fov_offsets(:,m_echo,m_slice,m_expt) = data.reco.fov_offsets(:,1,m_slice,m_expt);
            CBFmap.reco.fov_orientation(:,m_echo,m_slice,m_expt) = data.reco.fov_orientation(:,1,m_slice,m_expt);
            CBFmap.reco.label(m_echo,m_slice,m_expt) = data.reco.label(1,m_slice,m_expt);
            CBFmap.reco.phaselabel(m_echo,m_slice,m_expt) = data.reco.phaselabel(1,m_slice,m_expt);
            CBFmap.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = data.reco.fov_phase_orientation(1,m_slice,m_expt);
            CBFmap.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
            CBFmap.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
        end
    end
end
CBFmap.reco = orderfields(CBFmap.reco);

CBFmap.reco.displayedecho=CBFmap.reco.no_echoes;
CBFmap.reco.displayedslice=CBFmap.reco.no_slices;
CBFmap.reco.displayedexpt=1;
CBFmap.reco.thickness=data.reco.thickness;
CBFmap.reco.no_views=data.reco.no_views;
CBFmap.reco.no_samples=data.reco.no_samples;

CBFmap.reco.paramQuantif = ParamConfig;
CBFmap.reco=orderfields(CBFmap.reco);

CBFmap.scan_number = 105;
CBFmap.reco_number = 1;


%% T1 map formating

T1map.reco.echo_label(1,1) = {'fit_T1'};
T1map.reco.unit = {'s'};
%     if(NumFairMode==2)
%         CBFmap.reco.echo_label(1,1) = {'fit_IE selective'};
%         CBFmap.reco.echo_label(1,2) = {'fit_IE non selective'};
%         CBFmap.reco.unit(1,2) = {'a.u.'};
%     end
T1map.reco.texte = 'T1fitmap';
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
T1map.reco.displayedslice=T1map.reco.no_slices;
T1map.reco.displayedexpt=1;
T1map.reco.thickness=data.reco.thickness;
T1map.reco.no_views=data.reco.no_views;
T1map.reco.no_samples=data.reco.no_samples;

T1map.reco.paramQuantif = ParamConfig;
T1map.reco=orderfields(T1map.reco);

T1map.scan_number = 105;
T1map.reco_number = 1;
