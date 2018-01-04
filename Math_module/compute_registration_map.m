function moved_map = compute_registration_map(fixed_map_filename, moving_map_filename,  add_parameters)
% generate a moved_map from a moving map registered to a fixed map using
% several parameters defined by the user


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
    msgbox(moving_map_filename, 'Registration map warning') ;
    moved_map = [];
    return
end
moved_map = moving;
[optimizer,metric] = imregconfig('multimodal');
% [~,metric] = imregconfig(add_parameters{1});
%  metric = registration.metric.MattesMutualInformation;
% optimizer = registration.optimizer.RegularStepGradientDescent;
% optimizer.GradientMagnitudeTolerance = str2double(add_parameters{3});
% optimizer.MaximumIterations = str2double(add_parameters{2});
% optimizer.MinimumStepLength = 0.000001;
% optimizer.MaximumStepLength =  0.01;
% perform the registration
% need to determine the echo or expt num (imregister_tform works only
% wi
m =double(squeeze(moving.reco.data(:,:,1,:)));
f= double(squeeze(fixed.reco.data(:,:,1,:)));

% remove slice to the moving image to match with the fixed one
% if size(f,3) < size(m,3)
%     m = squeeze(m(:,:,1:size(f,3)));
% end
% if size(f,3) > size(m,3)
%     f = squeeze(f(:,:,1:size(m,3)));
% end
%check if NaN and inf in the matrix and replace by 0 because the imregister
%function do not like NaN

f(isnan(f)==1)=0;
m(isnan(m)==1)=0;
f(f==inf)=0;f(f==-inf)=0;
m(m==inf)=0;m(m==-inf)=0;
m(isnan(m)==1)=0;


% % resize the moving image to math with the fixed image
if size(m,1) > size(f,1)
    f_resize = imresize(f,[size(m,1) size(m,2)],'bilinear');
else
    m_resize = imresize(m,[size(f,1) size(f,2)],'bilinear');
end
% f_resize = imresize(f,[size(m,1) size(m,2)],'bilinear');
% 
% % [moved, tform] = tform_imregister(m, f, add_parameters{4} , optimizer, metric,'DisplayOptimization', true,...
% %     'PyramidLevels',str2double(add_parameters{5}), add_parameters{6},add_parameters{7});
% used with matlab 2012b
[moved, tform] = tform_imregister(m_resize, f, add_parameters{4} , optimizer, metric,'DisplayOptimization', true,...
    'PyramidLevels',str2double(add_parameters{5}), add_parameters{6},add_parameters{7});

%% matlab 2013b
% if size(m,1) > size(f,1)
%     [moved, tform] = imregister(m, f_resize, add_parameters{4} , optimizer, metric,'DisplayOptimization', true,...
%         'PyramidLevels',str2double(add_parameters{5}));
% else
%     [moved, tform] = imregister(m_resize, f, add_parameters{4} , optimizer, metric,'DisplayOptimization', true,...
%         'PyramidLevels',str2double(add_parameters{5}));
% end
if isfield(moving.reco, 'texte')
    moved_map.reco.texte =  [moving.reco.texte '-reg'];
end

moved_map.reco.data = permute(moved, [1 2 4 3]);


moved_map.reco.fov_offsets = fixed.reco.fov_offsets;
moved_map.reco.fov_orientation = fixed.reco.fov_orientation;
moved_map.reco.no_slices = fixed.reco.no_slices;
moved_map.reco.phaselabel = fixed.reco.phaselabel;
if isfield(fixed.reco, 'scaling_factor')
moved_map.reco.scaling_factor = fixed.reco.scaling_factor;
end
if isfield(fixed.reco, 'scaling_offset')
moved_map.reco.scaling_offset = fixed.reco.scaling_offset;
end
if size(m,1) <= size(f,1)
    moved_map.reco.no_samples = fixed.reco.no_samples;
    moved_map.reco.no_views = fixed.reco.no_views;
end
moved_map.reco.thickness = fixed.reco.thickness;

% save image registered
ParamConfig=sprintf('##$QuantifMethod=%s\n##$Iterations=%s\n##$Tolerance=%s\n##$Transformation=%s\n##$PyramidLevels=%s\n##$Interpolant=%s\n##$Pad Method=%s\n##$Reference_file=%s\n',...
    add_parameters{1},...
    add_parameters{2},...
    add_parameters{3},...
    add_parameters{4},...
    add_parameters{5},...
    add_parameters{6},...
    add_parameters{7},...
    fixed_map_filename);

moved_map.reco.paramQuantif = ParamConfig;
moved_map.reco.tform = tform;

