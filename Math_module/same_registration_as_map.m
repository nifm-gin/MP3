function moved_map = same_registration_as_map(moving_map_filename, tform_map_nbr_filename,  add_parameters)
% generate a moved_map from a matrix of transformation already computed
% and save in the .mat file


% load the fixed_map_filename file
fid=fopen(tform_map_nbr_filename ,'r');
if fid>0
    fclose(fid);
    fixed = load(tform_map_nbr_filename);
    fixed = fixed.uvascim.image;
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
        fixed_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    moved_map = [];
    return
end

% load the the matrix of tranformation
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

if ~isfield(fixed.reco, 'tform')
     warning_text = sprintf('##$ Can not perform the registration because \n##$ no transformation matrix has been found for the scan\n##$Moving scan=%s\n##$',...
        tform_map_nbr_filename);
    msgbox(warning_text, 'Registration map warning') ;
    moved_map = [];
    return
end
   
moved_map = moving;
if numel(size(fixed.reco.tform)) == 4
    Bx =single(squeeze(fixed.reco.tform(:,:,:,1)));
    By =single(squeeze(fixed.reco.tform(:,:,:,2)));
    Bz =single(squeeze(fixed.reco.tform(:,:,:,3)));
    registration_type = 'Kroon';
else
    registration_type = 'matlab';
end

for j=1:size(moving.reco.data, 3)
    % determine the number of expt
    for jj=1:size(moving.reco.data, 5)
        % 2 kinds of registration matrix: 1) from Kroon's method and 2)
        % from matlab method
        if strcmp(registration_type, 'Kroon')
            if sum(size(squeeze(moving.reco.data(:,:,1,:))) == size(Bx)) == 3
                tmp = squeeze(moving.reco.data(:,:,j,:,jj));
                moved(:,:,j,:,jj) =movepixels(tmp,Bx,By,Bz,3);
            else
                msgbox('probleme in the same_registration_as_map function', 'different matrix sizes for Kroon registration') ;
                return
            end
        else
            m =squeeze(moving.reco.data(:,:,j,:,jj));
            f= squeeze(fixed.reco.data(:,:,1,:,jj));
            nDims      = ndims(f);
            resampler  = makeresampler(add_parameters{3},add_parameters{4});
            fillValue  = 0;
            tdims_a    = 1:nDims;
            tdims_b    = 1:nDims;
            outSize    = size(f);
            moved(:,:,j,:,jj) = tformarray(m,fixed.reco.tform,resampler,...
                tdims_a, tdims_b,...
                outSize, [], fillValue);
        end
    end
end
%  if strcmp(registration_type, 'Kroon')
% if isfield(moving.reco, 'texte')
%     moved_map.reco.texte =  [moving.reco.texte '-reg'];
% end
%  end
% Set the class of output to input class
% if(strcmpi(Iclass,'uint8')), moved=uint8(moved*((2^8)-1)); end
% if(strcmpi(Iclass,'uint16')), moved=uint16(moved*((2^16)-1)); end
% if(strcmpi(Iclass,'uint32')), moved=uint32(moved*((2^32)-1)); end
% if(strcmpi(Iclass,'int8')), moved=int8(moved*((2^7)-1)); end
% if(strcmpi(Iclass,'int16')), moved=int16(moved*((2^15)-1)); end
% if(strcmpi(Iclass,'int32')), moved=int32(moved*((2^31)-1)); end

moved_map.reco.data = moved;

moved_map.reco.fov_offsets = fixed.reco.fov_offsets;
moved_map.reco.fov_orientation = fixed.reco.fov_orientation;
moved_map.reco.no_samples = fixed.reco.no_samples;
moved_map.reco.no_views = fixed.reco.no_views;
moved_map.reco.no_slices = fixed.reco.no_slices;
moved_map.reco.phaselabel = fixed.reco.phaselabel;
moved_map.reco.scaling_factor = fixed.reco.scaling_factor;
moved_map.reco.scaling_offset = fixed.reco.scaling_offset;
moved_map.reco.thickness = fixed.reco.thickness;

% save image registered
if isfield(fixed.reco, 'paramQuantif')
ParamConfig=sprintf('##$Interpolant=%s\n##$Pad Method=%s\n\n##$Matrix of tranformation origine=\n%s\n\n%s\n',...
    add_parameters{3},...
    add_parameters{4},...
    tform_map_nbr_filename,...
    fixed.reco.paramQuantif);
else
   ParamConfig=sprintf('##$Interpolant=%s\n##$Pad Method=%s\n\n##$Matrix of tranformation origine=\n%s\n\n',...
    add_parameters{3},...
    add_parameters{4},...
    tform_map_nbr_filename); 
end

moved_map.reco.paramQuantif = ParamConfig;
if numel(size(fixed.reco.tform)) == 4
    moved_map.reco.tform = [];
else
    moved_map.reco.tform = fixed.reco.tform;
end

