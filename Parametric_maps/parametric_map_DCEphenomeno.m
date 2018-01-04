function [maxi_struc, ttp_struc, rehaus_struc, AUC_struc, exclus_struc] = parametric_map_DCEphenomeno(filename_DCE, add_parameters)

% generate a 4 phenomenological maps from a DCE scan: 
%  i) max enhance 
%  ii)time-to-peak
%  iii) precentage of enhancement 
%  iv) exclued pixels 
%
% this code come from the permeability module 
% additional_parameters correspond to the size of the windows used for the loating mean and the beggining of the bolus 

% load the DCE file
fid=fopen(filename_DCE ,'r');
if fid>0
    fclose(fid);
    data = load(filename_DCE);
    reco = data.uvascim.image.reco;
else
     warning_text = sprintf('##$ Can not calculate the DCE-phenomeno map because there is\n##$ Somthing wrong with the data\n##$DCE=%s\n##$',...
        filename_DCE);
    msgbox(warning_text, 'DCE-phenomeno map warning') ;
    maxi_struc = []; ttp_struc = []; rehaus_struc = []; AUC_struc = []; exclus_struc = [];
    return
end


% check if human data (par/rec data)
if strcmp(data.uvascim.image.filename(end-3:end), '.PAR')
    TR = data.uvascim.image.acq.tr;
else
    timing = sscanf(scan_acqp('##$PVM_ScanTimeStr=',data.uvascim.image.texte,2),'%dh%dm%ds%dms');
    duree_tot=timing(1)*3600+timing(2)*60+timing(3)+timing(4)/1000; %en seconde
    TR = duree_tot / data.uvascim.image.reco.no_expts;
end
repetition_nbr = size(reco.data,5);
if ~strcmp(add_parameters{:}(1), 'Auto')
    debut = str2double(add_parameters{:}(1));
    fin = str2double(add_parameters{:}(2));
else
    mean_signal =max(reshape(data.uvascim.image.reco.data, [size(data.uvascim.image.reco.data,1)*size(data.uvascim.image.reco.data,2)*size(data.uvascim.image.reco.data,4) size(data.uvascim.image.reco.data,5)]),[], 1);
    mean_baseline = nanmean(mean_signal(1:3));
    sd_baseline = nanstd(mean_signal(1:3));
    debut = find(mean_signal>(mean_baseline+2*sd_baseline), 1)-1;
    repetion_time = 1:repetition_nbr;
    repetion_time = repetion_time*TR;
    [~, fin] = min(abs(repetion_time -(debut*TR + str2double(add_parameters{:}(2)))));
end
% tmp_data = permute(reco.data, [1 2 4 3]);
data_in_vector = reshape(reco.data, [size(reco.data,1)*size(reco.data,2)*size(reco.data,4), size(reco.data,5)]);


maxi = NaN(size(data_in_vector,1),1);
ttp  = NaN(size(data_in_vector,1),1);
rehaus  = NaN(size(data_in_vector,1),1);
AUC  = NaN(size(data_in_vector,1),1);
exclus  = NaN(size(data_in_vector,1),1);


for voxel_nbr=1:size(data_in_vector,1)% 
%     if data_in_vector(voxel_nbr,1) ~= 0 && data_in_vector(voxel_nbr,2) ~= 0
        tmp = squeeze(data_in_vector(voxel_nbr,:)); %
        moy_ss_gd = mean(tmp(1:debut));
        std_ss_gd = std(tmp(1:debut));
        
        %if moy_ss_gd && sum(tmp>(moy_ss_gd+2*std_ss_gd))>0
        % moving average -------------------------------
        for m=2:repetition_nbr-1
            tmp(m) = (data_in_vector(voxel_nbr,m-1)+data_in_vector(voxel_nbr,m)+data_in_vector(voxel_nbr,m+1))/3;
        end
        % maxi et ttp -------------------------------------
        [maxi(voxel_nbr),ttp(voxel_nbr)] = max(tmp(debut:fin));
        ttp(voxel_nbr) = ttp(voxel_nbr) * TR;
        % rehaussement ------------------------------------
        rehaus(voxel_nbr) = ((max(tmp)-moy_ss_gd)/moy_ss_gd)*100;
        % maxi(voxel_nbr) = data(i,j,k,l,ttp(i,j,k,l)); % Max intensity with
        % the filter
        
        % AUC Area under the curve ----------------------------
        AUC(voxel_nbr) = 0;
        for m=debut:fin
            AUC(voxel_nbr) = AUC(voxel_nbr) + (data_in_vector(voxel_nbr,m)-moy_ss_gd);
        end
%     else
%         exclus(voxel_nbr) = 1;
%     end
end

% maxi = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% ttp = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% rehaus = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% exclus = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% AUC = zeros(size(reco.data,1),size(reco.data,2),size(reco.data,3),size(reco.data,4),'single');
% 
% h = waitbar(0, 'DCE-phenomeno maps in progress...');
% steps_tot= reco.no_slices*reco.no_samples*reco.no_views;
% current_step = 1;
% 
% 
% for i=1:size(reco.data,1) % lignes
%     for j=1:size(reco.data,2) % colonnes
%         for k=1:size(reco.data,3) % echo
%             for l=1:size(reco.data,4) % coupes
%                 waitbar(current_step/steps_tot, h)
%                 tmp = squeeze(reco.data(i,j,k,l,:)); % donn�es temporelles
%                 if tmp(1) ~= 0 && tmp(2) ~= 0
%                     moy_ss_gd = mean(tmp(1:debut));
%                     std_ss_gd = std(tmp(1:debut));
%                     %if moy_ss_gd && sum(tmp>(moy_ss_gd+2*std_ss_gd))>0
%                     % moyenne flottante -------------------------------
%                     for m=2:size(tmp)-1
%                         tmp(m) = (reco.data(i,j,k,l,m-1)+reco.data(i,j,k,l,m)+reco.data(i,j,k,l,m+1))/3;
%                     end
%                     % maxi et ttp -------------------------------------
%                     [maxi(i,j,k,l),ttp(i,j,k,l)] = max(tmp(debut:fin));
%                     ttp(i,j,k,l) = ttp(i,j,k,l) * TR;
%                     % rehaussement ------------------------------------
%                     rehaus(i,j,k,l) = ((max(tmp)-moy_ss_gd)/moy_ss_gd)*100;
%                     % maxi(i,j,k,l) = data(i,j,k,l,ttp(i,j,k,l)); % pour r�cup�rer l'intensit� max et pas celle filtr�e
%                     % AUC sous la courbe ----------------------------
%                     AUC(i,j,k,l) = 0;
%                     for m=debut:fin
%                         AUC(i,j,k,l) = AUC(i,j,k,l) + (reco.data(i,j,k,l,m)-moy_ss_gd);
%                     end
%                     %else
%                     %   exclus(i,j,k,l) = 1;
%                     %end
%                 else
%                     exclus(i,j,k,l) = 1;
%                 end
%                 current_step = current_step+1;
%             end
%         end
%     end
% end

% save imformation


tempstruct=data.uvascim.image;
tempstruct.reco.data = zeros(size(reco.data,1),...
    size(reco.data,2),1,size(reco.data,4),1,'single');
tempstruct.reco.no_echoes = 1;
tempstruct.reco.displayedecho = 1;
tempstruct.reco.displayedslice = 1;
tempstruct.reco.displayedexpt = 1;
tempstruct.reco.no_expts = 1;
tempstruct.reco.echotime = NaN('single');

%Adapt the fov offsets and orientation infos
tmpfov=reco.fov_offsets(:,1,:,1);
tmpori=reco.fov_orientation(:,1,:,1);
tmpphaselabel=reco.phaselabel(1,:,1);
tempstruct.reco.fov_offsets = [];
tempstruct.reco.fov_orientation = [];
tempstruct.reco.scaling_factor = [];
tempstruct.reco.scaling_offset = [];
tempstruct.reco.label = {};
tempstruct.reco.phaselabel = {};

for m_slice=1:reco.no_slices,
    tempstruct.reco.fov_offsets(:,:,m_slice)=squeeze(tmpfov(:,m_slice));
    tempstruct.reco.fov_orientation(:,:,m_slice)=squeeze(tmpori(:,m_slice));
    tempstruct.reco.label(1,m_slice)=reco.label(1,m_slice,1);
    tempstruct.reco.phaselabel(1,m_slice)=reco.phaselabel(1,m_slice,1);
    tempstruct.reco.scaling_factor(1,m_slice) = reco.scaling_factor(1,m_slice,1);
    tempstruct.reco.scaling_offset(1,m_slice) = reco.scaling_offset(1,m_slice,1);
end
tempstruct.reco.fov_phase_orientation=reco.fov_phase_orientation(:,:,1);


%delete useless field
tempstruct.reco = rmfield(tempstruct.reco, 'angAP');
tempstruct.reco = rmfield(tempstruct.reco, 'angFH');
tempstruct.reco = rmfield(tempstruct.reco, 'angRL');
tempstruct.reco = rmfield(tempstruct.reco, 'angulation');
tempstruct.reco = rmfield(tempstruct.reco, 'echotime');
tempstruct.reco = rmfield(tempstruct.reco, 'iminfos');
tempstruct.reco = rmfield(tempstruct.reco, 'reco_meth');

ParamConfig=sprintf('##$QuantifMethod=phenomenological maps from a DCE scan\n##$Beggining of the injection=%s\n##$End of the analysis=%s\n##$Raw scan used=%s',...
        num2str(debut),...
        num2str(fin),...
        filename_DCE);
tempstruct.reco.paramQuantif = ParamConfig;

tempstruct.reco = rmfield(tempstruct.reco, 'data');
maxi_struc = tempstruct;
ttp_struc = tempstruct;
rehaus_struc = tempstruct;
AUC_struc = tempstruct;
exclus_struc = tempstruct;

% maxi
maxi_struc.reco.data = reshape(maxi,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
maxi_struc.reco.unit = {'a.u.'};
maxi_struc.reco.echo_label = {'DCE-Max'};
% ttp (time-to-peak(
ttp_struc.reco.data =  reshape(ttp,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
ttp_struc.reco.unit = {'s'};
ttp_struc.reco.echo_label = {'DCE-ttp'};
% rehaus = (intensite-MoySansGd)/MoySansGd
rehaus_struc.reco.data = reshape(rehaus,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
rehaus_struc.reco.unit = {'%'};
rehaus_struc.reco.echo_label = {'DCE-rehaus'};
% Area under the curve (AUC)
AUC_struc.reco.data = reshape(AUC,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]); 
AUC_struc.reco.unit = {'a.u.'};
AUC_struc.reco.echo_label = {'DCE-AUC'};
% Exclued voxels
exclus_struc.reco.data = reshape(exclus,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
exclus_struc.reco.unit = {''};
exclus_struc.reco.echo_label = {'DCE-exclus'};

% close(h);
