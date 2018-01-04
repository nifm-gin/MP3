function [DCEconcentration] = parametric_map_DCEconcentration(filename_DCE, filename_T1map, add_parameters)

% generate a DCE scan in concentration of Gd


debut = str2double(add_parameters{:}(1));
r1 = str2double(add_parameters{:}(2));

fid=fopen(filename_DCE ,'r');
if fid>0
    fclose(fid);
    dataDCE = load(filename_DCE);
    DCE.reco = dataDCE.uvascim.image.reco;
else
     warning_text = sprintf('##$ Can not calculate the DCE-phenomeno map because there is\n##$ Somthing wrong with the data\n##$DCE=%s\n##$',...
        filename_DCE);
    msgbox(warning_text, 'DCE-phenomeno map warning') ;
    DCEconcentration = [];
    return
end

fid=fopen(filename_T1map ,'r');
if fid>0
    fclose(fid);
    dataT1 = load(filename_T1map);
    T1map.reco = dataT1.uvascim.image.reco;
else
     warning_text = sprintf('##$ Can not calculate the DCE-phenomeno map because there is\n##$ Somthing wrong with the data\n##$DCE=%s\n##$',...
        filename_T1map);
    msgbox(warning_text, 'filename_T1map map warning') ;
    DCEconcentration = [];
    return
end

TR = dataDCE.uvascim.image.acq.tr/1000;

% tmp_data = permute(reco.data, [1 2 4 3]);
data_in_vector = reshape(DCE.reco.data, [size(DCE.reco.data,1)*size(DCE.reco.data,2)*size(DCE.reco.data,4), size(DCE.reco.data,5)]);
T1map_in_vector = reshape(T1map.reco.data, [size(T1map.reco.data,1)*size(T1map.reco.data,2)*size(T1map.reco.data,4) 1]);
T1map_in_vector = T1map_in_vector/1000;
repetition_nbr = size(DCE.reco.data,5);

DCEconcentration_data = NaN(size(data_in_vector,1),size(DCE.reco.data,5));
ratio_T1map = 1./T1map_in_vector;
ratio_r1 = 1/r1;
if size(data_in_vector, 1) ~= size(T1map_in_vector,1)
    DCEconcentration = [];
    return
end
for voxel_nbr=1:size(data_in_vector,1)% 
    if ~isnan(T1map_in_vector(voxel_nbr))
        S = squeeze(data_in_vector(voxel_nbr,:));
        %% Compute the base line (S0)
        S0 = mean(S(1:debut));
        %% Compute T1(t) as equal to T1(t) = TR/-ln(S(t)/A - 1) and A = S0/1-e-(TR/T10)
        A = S0 / (1 - exp(-(TR/T1map_in_vector(voxel_nbr))));
        for i = 1 : size(S,2)
            tmpb=1-S(i)/A;
            if tmpb <0.01,
                T1_temps(i)=T1map_in_vector(voxel_nbr)*5;
            else
                T1_temps(i) = TR / -(log(abs(1-S(i)/A)));
            end
 
        end
        % compute the DCE in concentration
        DCEconcentration_data(voxel_nbr,:) = (1./T1_temps - ratio_T1map(voxel_nbr)) .* ratio_r1;
    else
        exclus(voxel_nbr) = 1;
    end
end
DCEconcentration_data(DCEconcentration_data<0) = 0;
DCEconcentration_data=reshape(DCEconcentration_data,[size(DCE.reco.data,1),size(DCE.reco.data,2),size(DCE.reco.data,4),size(DCE.reco.data,5)]);
DCEconcentration_data = permute(DCEconcentration_data, [1 2 5 3 4]);

% save imformation


DCEconcentration=dataDCE.uvascim.image;


ParamConfig=sprintf('##$QuantifMethod=DCEconcentration maps from a DCE-SE scan\n##$Beggining of the injection=%s\n##$Relaxivity of th CA=%s\n##$DCE-SE scan used=%s\n##$T1map used=%s',...
        add_parameters{:}{1},...
        add_parameters{:}{2},...
        filename_DCE,filename_T1map);
DCEconcentration.reco.paramQuantif = ParamConfig;

DCEconcentration.reco = rmfield(DCEconcentration.reco, 'data');
DCEconcentration.reco.data = DCEconcentration_data;
DCEconcentration.reco.unit = {'mmol'};
DCEconcentration.reco.echo_label = {'DCEconcentration'};

% % Exclued voxels
% exclus_struc.reco.data = reshape(exclus,[size(reco.data,1),size(reco.data,2),1,size(reco.data,4),1]);
% exclus_struc.reco.unit = {''};
% exclus_struc.reco.echo_label = {'DCE-exclus'};

% close(h);
