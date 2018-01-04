function BVf_map = parametric_map_BVf(deltaR2star_map_filename, add_parameters)
% generate a BVf map from a deltaR2* map
% additional_parameters correspond to :
%           - B0, gamma and the deltaxi of the contrast agent
% threshold used


B0=str2double(add_parameters{:}(1));
gamma=str2double(add_parameters{:}(2)); gamma = gamma*10^8;
deltaxi_CA=str2double(add_parameters{:}(3)); deltaxi_CA = deltaxi_CA*10^-6;

deltaomega=2*pi*gamma*deltaxi_CA*B0;

% load the MGE file
fid=fopen(deltaR2star_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(deltaR2star_map_filename);
    reco = data.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the BVf map because there is\n##$ Somthing wrong with the data\n##$deltaR2star=%s\n##$',...
        deltaR2star_map_filename);
    msgbox(warning_text, 'BVf map warning') ;
    BVf_map = [];
    return
end
% save imformation
BVf_map=data.uvascim.image;

% Empty memory
clear data

for i = 1:BVf_map.reco.no_slices
    BVf_map.reco.data(:,:,1,i)=(3/(2*deltaomega)*BVf_map.reco.data(:,:,1,i))*100000;
end
BVf_map.reco.texte='BVf';
BVf_map.reco.globalmax = max(BVf_map.reco.data(:));
BVf_map.reco.globalmin = min(BVf_map.reco.data(:));
BVf_map.clip=[0 50 1];
BVf_map.reco.unit = '%';
if isfield(BVf_map.reco, 'paramQuantif')
ParamConfig=sprintf('##$QuantifMethod=''(3/(2*deltaomega)*BVf)*100000''\n##$DeltaR2*=%s\n##$B0=%s\n##$gamma=%s\n##$deltaxi_CA=%s\n##$deltaomega=%s\n\n##$DeltaR2* info\n%s',...
    deltaR2star_map_filename,add_parameters{:}{1}, strcat(add_parameters{:}{2}, '10^8'),strcat(add_parameters{:}{3}, '10^-6'), '2*pi*gamma*deltaxi_CA*B0', reco.paramQuantif);
else
   ParamConfig=sprintf('##$QuantifMethod=''(3/(2*deltaomega)*BVf)*100000''\n##$DeltaR2*=%s\n##$B0=%s\n##$gamma=%s\n##$deltaxi_CA=%s\n##$deltaomega=%s',...
    deltaR2star_map_filename,add_parameters{:}{1}, strcat(add_parameters{:}{2}, '10^8'),strcat(add_parameters{:}{3}, '10^-6'), '2*pi*gamma*deltaxi_CA*B0');
end
BVf_map.reco.paramQuantif = ParamConfig;
BVf_map.reco=orderfields(BVf_map.reco);










