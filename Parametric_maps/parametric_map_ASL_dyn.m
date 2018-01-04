function [output1,output2,output3] = parametric_map_ASL_dyn(it_sl,ASL,VoiOutlier_filename,PVM_ObjOrderList, interSliceTime,alpha,T1blood,lambda,fit_CBF_result,fit_transit_result,fit_T1app_result)
% generate a Transit_time_map, CBF map, M0 map, T1 map from a dynamic pCASL scan
% this code was coded by Lydiane Hirschler
T1blood = T1blood/1000;
%% Init variables
[Voi, VoiYn] = extract_map(VoiOutlier_filename,'ASL-dyn','VOI');
if isempty(VoiOutlier_filename)
    VoiYn = 0;
end

if VoiYn
    factorResize = size(ASL.reco.data,1) / size(Voi(it_sl).value,1);
    if size(ASL.reco.data,1)~= size(ASL.reco.data,2)
        factorResize = [size(ASL.reco.data,1) size(ASL.reco.data,2)];
    end
    VoiTemp = imresize(Voi(it_sl).value,factorResize,'nearest');
    ASL.reco.data(:,:,:,it_sl,:) = ASL.reco.data(:,:,:,it_sl,:).*repmat(VoiTemp,[1 1 size(ASL.reco.data,3) 1 size(ASL.reco.data,5)]);
end

%% Lecture parametres fichier methode
step = scan_acqp('##$PCASL_TempResolution=',ASL.texte,1)/1000;
T = scan_acqp('##$PCASL_TotalLabelTime=',ASL.texte,1)/1000;
PostLabelTime=scan_acqp('ASL_PostLabelTime=',ASL.texte,1)/1000;

if isnan(step)
    subboliDuration = scan_acqp('PCASL_SubboliTime=',ASL.texte,1)/1000;
    step = subboliDuration(1);
    T = scan_acqp('PCASL_LabelTime=',ASL.texte,1)/1000;
end

NumSlice = find(PVM_ObjOrderList==it_sl-1)-1;
if size(PVM_ObjOrderList) == 1 % cas 1 slice ou 3D
    NumSlice=0;
end
Effective_PLDslice = PostLabelTime + interSliceTime/1000 * NumSlice;

dummyScans = size(ASL.reco.data,3)-numel(subboliDuration);

data_in_vector = reshape(ASL.reco.data, [size(ASL.reco.data,1)*size(ASL.reco.data,2),...
    size(ASL.reco.data,3), size(ASL.reco.data,4),size(ASL.reco.data,5)]);

%%
%% FITTING

%     %%TESTS
%     load('DASLdata.mat')
%     data_in_vector=y_data;
    
    %% A decommenter si donnees a l envers ou plusieurs cycles en une repetition
    %         data_in_vectorTemp_endroit = squeeze(data_in_vector(:,2:end,it_sl));
    %         data_in_vectorTemp_moy=(data_in_vectorTemp_endroit(:,1:T/step) + data_in_vectorTemp_endroit(:,T/step+1:nb_cycles*T/step))/nb_cycles; %moyenne sur les cycles
    %         for i = 1:T/step;
    %             data_in_vectorTemp(:,i) = data_in_vectorTemp_moy(:,T/step-i+1); %donnees dyn lues a l'envers
    %         end
    %%
    data_in_vectorTemp_cycles = data_in_vector(:,dummyScans+1:end,it_sl,:); % ignorer les premieres images (dummy scans "simules")
    data_in_vectorTemp = squeeze(mean(data_in_vectorTemp_cycles,4)); %moyenne sur les repetitions (ie. cycles)
    
    %     indT1 = find((data_in_vector(:,it_sl,1) ~= 0) & (isnan(data_in_vector(:,it_sl,1)) ~= 1));
    %     data_in_vectorTemp = squeeze(data_in_vector(indT1,it_sl,:));
    %NumSlice = find(PVM_ObjOrderList==it_sl-1)-1;
    
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
           % if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 & ccc==-1 %#ok<AND2>
                fit_transit_result(voxel_nbr,it_sl)=aaa(3);
                fit_transit_err(voxel_nbr,it_sl)=bbb(3);
                
                fit_M0_result(voxel_nbr,it_sl)=aaa(1);
                fit_M0_err(voxel_nbr,it_sl)=bbb(1);
                
                fit_CBF_result(voxel_nbr,it_sl)=aaa(2);
                fit_CBF_err(voxel_nbr,it_sl)=bbb(2);
                fit_T1app_result(voxel_nbr,it_sl)=aaa(4);
                fit_T1_err(voxel_nbr,it_sl)=bbb(4);                
          %  end            
        end
    end
    parfor_progress(0); % Clean up


output1 = fit_CBF_result(:,it_sl);
output2 = fit_transit_result(:,it_sl);
output3 = fit_T1app_result(:,it_sl);


