function moved_map = compute_registration_map_Kroon(fixed_map_filename, moving_map_filename)
% generate a moved_map from a moving map registered to a fixed map using
% the algorithum developped by Kroon
% Function is written by D.Kroon University of Twente (March 2009)


% load the fixed_map_filename file
fid=fopen(fixed_map_filename ,'r');
if fid>0
    fclose(fid);
    fixed = load(fixed_map_filename);
    fixed = fixed.uvascim.image;
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
        fixed_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    moved_map = [];
    return
end

% load the file to move
fid=fopen(moving_map_filename ,'r');
if fid>0
    fclose(fid);
    moving = load(moving_map_filename);
    moving = moving.uvascim.image;
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the moving scan\n##$Moving scan=%s\n##$',...
        moving_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    moved_map = [];
    return
end
moved_map = moving;


% perform the registration
% need to determine the echo or expt num (imregister_tform works only
% wi
Imoving =double(squeeze(moving.reco.data(:,:,1,:)));
Istatic= double(squeeze(fixed.reco.data(:,:,1,:)));

% [moved,Bx,By,Bz,~,~,~] = register_volumes(Imoving,Istatic);
[moved,~,~,~,~,~,~] = register_volumes(Imoving,Istatic);
tform = [];
% tform(:,:,:,1) = Bx;
% tform(:,:,:,2) = By;
% tform(:,:,:,3) = Bz;


if isfield(moving.reco, 'texte')
    moved_map.reco.texte =  [moving.reco.texte '-reg'];
end

moved_map.reco.data = permute(moved, [1 2 4 3]);


moved_map.reco.fov_offsets = fixed.reco.fov_offsets;
moved_map.reco.fov_orientation = fixed.reco.fov_orientation;
moved_map.reco.no_slices = fixed.reco.no_slices;
moved_map.reco.phaselabel = fixed.reco.phaselabel;
moved_map.reco.scaling_factor = fixed.reco.scaling_factor;
moved_map.reco.scaling_offset = fixed.reco.scaling_offset;

moved_map.reco.no_samples = fixed.reco.no_samples;
moved_map.reco.no_views = fixed.reco.no_views;
moved_map.reco.thickness = fixed.reco.thickness;

% save image registered
ParamConfig=sprintf('##$QuantifMethod=%s\n##$Reference_file=%s\n',...
    'Kroon algo with default parameters',...
    fixed_map_filename);

moved_map.reco.paramQuantif = ParamConfig;
moved_map.reco.tform = tform;

