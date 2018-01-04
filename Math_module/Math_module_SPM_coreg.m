function [Coreg_maps, uvascim] = Math_module_SPM_coreg(fixed_map_filename, moving_map_filename, other_map_filename, add_parameters)



% load the file to move
fid=fopen(moving_map_filename ,'r');
if fid>0
    fclose(fid);
    moving = load(moving_map_filename);
    [PATHSTR, NAME_scan_to_coreg, ~] = fileparts(moving_map_filename);
    % convert moving to nii
    NAME_scan_to_coreg = [NAME_scan_to_coreg '-for_spm'];
    if numel(size(squeeze(moving.uvascim.image.reco.data))) < 5
        NHdr=CreateNiiHdr_from_uvascim(moving_map_filename, moving, NAME_scan_to_coreg);
        if  size(moving.uvascim.image.reco.data,3) >  size(moving.uvascim.image.reco.data,5)
            data_4d = permute(moving.uvascim.image.reco.data, [1 2 4 3]);
        else
            data_4d = squeeze(moving.uvascim.image.reco.data(:,:,1,:,:));
        end
        %data_4d(data_4d == 0) = NaN;
        WriteNii(fullfile(PATHSTR, [NAME_scan_to_coreg '.nii']),NHdr, data_4d)
    else
        % create 4D nii from 5D scan by adding the repetition after the echo
        data_4d = permute(moving.uvascim.image.reco.data, [1 2 4 3 5]);
        data_4d = reshape(data_4d, [size(data_4d, 1), size(data_4d, 2), size(data_4d, 3), size(data_4d, 4)*size(data_4d, 5)]);
        %
        NHdr=CreateNiiHdr_from_uvascim(moving_map_filename, moving, NAME_scan_to_coreg);
        %data_4d(data_4d == 0) = NaN;
        WriteNii(fullfile(PATHSTR, [NAME_scan_to_coreg '.nii']),NHdr, data_4d)
    end
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the moving scan\n##$Moving scan=%s\n##$',...
        moving_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    Coreg_maps = [];
    return
end
final_resolution = add_parameters(1,5);

% load the fixed_map_filename file
fid=fopen(fixed_map_filename ,'r');
if fid>0
    fclose(fid);
    fixed = load(fixed_map_filename);
    % convert fixed to nii
    [PATHSTR, NAME_fixed, ~] = fileparts(fixed_map_filename);
    NAME_fixed = [NAME_fixed '-for_spm'];
    if ~strcmp(final_resolution, 'Same as Ref')
        if strcmp(final_resolution, 'Unchanged')
            resolution = [size(moving.uvascim.image.reco.data,1) size(moving.uvascim.image.reco.data,2)];
        else
            resolution = [str2num(final_resolution{:}) str2num(final_resolution{:})];
        end
        for echo=1:size(fixed.uvascim.image.reco.data,3)
            for slice = 1:size(fixed.uvascim.image.reco.data,4)
                for repet = 1:size(fixed.uvascim.image.reco.data,5)
                    tmp (:,:,echo,slice,repet)= imresize(squeeze(fixed.uvascim.image.reco.data(:,:,echo,slice,repet)), resolution);
                end
            end
        end
        fixed.uvascim.image.reco = rmfield(fixed.uvascim.image.reco, 'data');
        fixed.uvascim.image.reco.data = tmp;
        % update uvasc structure
        fixed.uvascim.image.acq.no_samples = resolution(1);
        fixed.uvascim.image.acq.no_views = resolution(2);
        if ~isfield(fixed.uvascim.image.acq, 'pix_spacing')
            % maybe needs to performe this calcul for everyone!
            fixed.uvascim.image.acq.pix_spacing = [fixed.uvascim.image.acq.fov(1)/resolution(1) fixed.uvascim.image.acq.fov(2)/resolution(2)]';
        else
            fixed.uvascim.image.acq.pix_spacing = fixed.uvascim.image.acq.pix_spacing;
        end
        fixed.uvascim.image.reco.no_samples = resolution(1);
        fixed.uvascim.image.reco.no_views = resolution(2);
        
    end
    NHdr=CreateNiiHdr_from_uvascim(fixed_map_filename, fixed, NAME_fixed);
    if numel(size(squeeze(fixed.uvascim.image.reco.data))) < 5
        if  size(fixed.uvascim.image.reco.data,3) >  size(fixed.uvascim.image.reco.data,5)
            data_4d = permute(fixed.uvascim.image.reco.data, [1 2 4 3]);
        else
            data_4d = squeeze(fixed.uvascim.image.reco.data(:,:,1,:,:));
        end
        %data_4d(data_4d == 0) = NaN;
        WriteNii(fullfile(PATHSTR, [NAME_fixed '.nii']), NHdr, data_4d)
    else
        % create 4D nii from 5D scan
            % create 4D nii from 5D scan by adding the repetition after the echo
        data_4d = permute(fixed.uvascim.image.reco.data, [1 2 4 3 5]);
        data_4d = reshape(data_4d, [size(data_4d, 1), size(data_4d, 2), size(data_4d, 3), size(data_4d, 4)*size(data_4d, 5)]);
        %
%         NHdr=CreateNiiHdr_from_uvascim(fixed_map_filename, fixed, NAME_fixed);
    %data_4d(data_4d == 0) = NaN;
        WriteNii(fullfile(PATHSTR, [NAME_fixed '.nii']),NHdr, data_4d)
    end
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
        fixed_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    Coreg_maps = [];
    return
end



% load other images (if not empty)
if ~isempty(other_map_filename)
    for i=1:size(other_map_filename,2)
        fid=fopen(other_map_filename{:,i},'r');
        if fid>0
            fclose(fid);
            Other_maps(i)= load(other_map_filename{:,i});
            
            % convert Other_maps to nii
            [PATHSTR, tmp_name, ~] = fileparts(other_map_filename{:,i});
            NAME_scan_other(i,1) = {[tmp_name '-for_spm']}; %#ok<AGROW>
            if numel(size(squeeze(Other_maps(i).uvascim.image.reco.data))) < 5
                NHdr=CreateNiiHdr_from_uvascim(other_map_filename{:,i}, Other_maps(i), NAME_scan_other{i,1});
                if  size(Other_maps(i).uvascim.image.reco.data,3) >  size(Other_maps(i).uvascim.image.reco.data,5)
                    data_4d = permute(Other_maps(i).uvascim.image.reco.data, [1 2 4 3]);
                else
                    data_4d = squeeze(Other_maps(i).uvascim.image.reco.data(:,:,1,:,:));
                end
                WriteNii(fullfile(PATHSTR, [NAME_scan_other{i,1} '.nii']),NHdr, data_4d)
                
            else
                % create 4D nii from 5D scan by adding the repetition after the echo
                data_4d = permute(Other_maps(i).uvascim.image.reco.data, [1 2 4 3 5]);
                data_4d = reshape(data_4d, [size(data_4d, 1), size(data_4d, 2), size(data_4d, 3), size(data_4d, 4)*size(data_4d, 5)]);
                %
                NHdr=CreateNiiHdr_from_uvascim(other_map_filename{:,i}, Other_maps(i), NAME_scan_other{i,1});
                WriteNii(fullfile(PATHSTR, [NAME_scan_other{i,1} '.nii']),NHdr, data_4d)      
            end
            for xx = 1:size(data_4d,4)
                tmp_name2{xx,:} =  [ NAME_scan_other{i,1} '.nii,' num2str(xx)];
            end
            estwrite_data(i) = {tmp_name2};
            clear tmp_name2  
        else
            Other_maps(i).data = '';
        end
    end
end

if strcmp(add_parameters(1,12), 'Mask image')
    % % test with mask for image fixed
    % % create binary mask
    g = double(squeeze(fixed.uvascim.image.reco.data));
    %     max_ref = max(g(list_ind));
    mask = g <= 3*median(g(:));
    % Connected component labelling and brain selection
    mask = ~mask;
    label = bwlabeln(mask,6);
    STATS = regionprops(label, 'Area');
    [~,plus_grand_label] = max([STATS.Area]);
    if ~isempty(plus_grand_label)
        mask (label ~= plus_grand_label) = 0;
    end
    mask = imerode(mask,strel('disk',3));
    label = bwlabeln(mask,6);
    STATS = regionprops(label, 'Area');
    [~,plus_grand_label] = max([STATS.Area]);
    if ~isempty(plus_grand_label)
        mask (label ~= plus_grand_label) = 0;
    end
    mask = imdilate(mask,strel('disk',3));
    fixed.uvascim.image.reco.data=fixed.uvascim.image.reco.data .* permute(mask, [1 2 4 3]);
    
    % test with mask for image moving
    % create binary mask
    g = double(squeeze(moving.uvascim.image.reco.data));
    %     max_ref = max(g(list_ind));
    mask = g <= 3*median(g(:));
    % Connected component labelling and brain selection
    mask = ~mask;
    label = bwlabeln(mask,6);
    STATS = regionprops(label, 'Area');
    [~,plus_grand_label] = max([STATS.Area]);
    if ~isempty(plus_grand_label)
        mask (label ~= plus_grand_label) = 0;
    end
    mask = imerode(mask,strel('disk',3));
    label = bwlabeln(mask,6);
    STATS = regionprops(label, 'Area');
    [~,plus_grand_label] = max([STATS.Area]);
    if ~isempty(plus_grand_label)
        mask (label ~= plus_grand_label) = 0;
    end
    mask = imdilate(mask,strel('disk',3));
    moving.uvascim.image.reco.data=moving.uvascim.image.reco.data .* permute(mask, [1 2 4 3]);
    
end

% set the origine manually
if strcmp(add_parameters(1,14), 'Yes') == 1;
    spm_image('Display', strcat(PATHSTR, filesep,  NAME_scan_to_coreg, '.nii'));
    disp('press enter when the origin is updated')
    pause
    spm_image('Display', strcat(PATHSTR, filesep,  NAME_fixed, '.nii'))
    disp('press enter when the origin is updated')
    pause
end

%% perform co-reg alogrithum (using SPM 12)
% NAME_fixed = 'SO2_MSME-map';
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(PATHSTR, [NAME_fixed '.nii,1'])};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(PATHSTR,[NAME_scan_to_coreg '.nii,1'])};
other_files = '';
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};

aa = spm_select('ExtList', PATHSTR, [NAME_scan_to_coreg '.nii'], 1:999);
if size(aa,1) >1
    other_files = cellstr(aa)';
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(aa);
end
if ~isempty(other_map_filename)
    for i=1:numel(NAME_scan_other)
        aa = spm_select('ExtList', PATHSTR, [NAME_scan_other{i} '.nii'], 1:999);
        other_files = [other_files cellstr(aa)'];
    end
    matlabbatch{1}.spm.spatial.coreg.estwrite.other =  other_files';
end
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = add_parameters{1,6};
if strcmp(add_parameters{1,7}, 'Auto= [slice thickness voxel_size voxel_size/2]') 
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [fixed.uvascim.image.reco.thickness fixed.uvascim.image.reco.fov(1)/fixed.uvascim.image.reco.no_samples  fixed.uvascim.image.reco.fov(1)/fixed.uvascim.image.reco.no_samples/2];
else
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = str2num(add_parameters{1,7});
end
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = str2num(add_parameters{1,8}); %#ok<*ST2NM>
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = str2num(add_parameters{1,9});
%define options
% Type of interpolation
switch add_parameters{1,10}
    case 'Nearest neighbour'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
    case 'Trilinear'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7;
end
%  Type of Warpping
switch add_parameters{1,11}
    case 'No wrap'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    case 'Warp X'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [1 1 1];
end
%  Mask?
switch add_parameters{1,12}
    case 'Dont mask images'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask= 0;
    case 'Mask image'
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 1;
end
%% always set mask to 0
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask= 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = add_parameters{1,13};

[SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});

date_str = date;
if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
    movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
        [PATHSTR, filesep, NAME_scan_to_coreg(1:end-8) '_SPMcoreg.ps'], 'f'); 
end
close(SPMinter)
close(SPMgraph)


% create output structure
% update first output structure and in particular the reco field:
% based on the original structure (the moving structure) but with the info
% form the fixed structure
uvascim(1) = moving;

uvascim.uvascim.image.reco.angAP = fixed.uvascim.image.reco.angAP;
uvascim.uvascim.image.reco.angFH = fixed.uvascim.image.reco.angFH;
uvascim.uvascim.image.reco.angRL = fixed.uvascim.image.reco.angRL;
uvascim.uvascim.image.reco.angulation = fixed.uvascim.image.reco.angulation;
uvascim.uvascim.image.reco.fov = fixed.uvascim.image.reco.fov;
uvascim.uvascim.image.reco.no_slices = fixed.uvascim.image.reco.no_slices;
uvascim.uvascim.image.reco.no_samples = fixed.uvascim.image.reco.no_samples;
uvascim.uvascim.image.reco.no_views = fixed.uvascim.image.reco.no_views;
uvascim.uvascim.image.reco.thickness = fixed.uvascim.image.reco.thickness;

% fov_offsets = [3 echoes slice repet]
uvascim.uvascim.image.reco.fov_offsets=repmat(fixed.uvascim.image.reco.fov_offsets(:,1,:,1),[1,uvascim.uvascim.image.reco.no_echoes,1,uvascim.uvascim.image.reco.no_expts]);
% fov_orientation = [9 echoes slice repet]
uvascim.uvascim.image.reco.fov_orientation = repmat(fixed.uvascim.image.reco.fov_orientation(:,1,:,1),[1,uvascim.uvascim.image.reco.no_echoes,1,uvascim.uvascim.image.reco.no_expts]);
% fov_phase_orientation = [echoes slice repet]
uvascim.uvascim.image.reco.fov_phase_orientation = repmat(fixed.uvascim.image.reco.fov_phase_orientation(1,:,1),[uvascim.uvascim.image.reco.no_echoes,1,uvascim.uvascim.image.reco.no_expts]);
% scaling_factor = [echoes slice repet]
uvascim.uvascim.image.reco.scaling_factor = repmat(fixed.uvascim.image.reco.scaling_factor(:,1,1),[uvascim.uvascim.image.reco.no_echoes,uvascim.uvascim.image.reco.no_slices,uvascim.uvascim.image.reco.no_expts]);
% scaling_offset = [echoes slice repet]
uvascim.uvascim.image.reco.scaling_offset = repmat(fixed.uvascim.image.reco.scaling_offset(:,1,1),[uvascim.uvascim.image.reco.no_echoes,uvascim.uvascim.image.reco.no_slices,uvascim.uvascim.image.reco.no_expts]);
% phaselabel = [echoes slice repet]
uvascim.uvascim.image.reco.phaselabel = repmat(fixed.uvascim.image.reco.phaselabel(:,1,1),[uvascim.uvascim.image.reco.no_echoes,uvascim.uvascim.image.reco.no_slices,uvascim.uvascim.image.reco.no_expts]);
% label = [echoes slice repet]
uvascim.uvascim.image.reco.label = repmat(fixed.uvascim.image.reco.label(:,1,1),[uvascim.uvascim.image.reco.no_echoes,uvascim.uvascim.image.reco.no_slices,uvascim.uvascim.image.reco.no_expts]);


ParamConfig=sprintf('##$SPM_coreg_write algorithum=\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n##$Name of each parameters=\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n##$Parameters set up%s\n',...
    'Ref day', 'Ref Scan', 'Interpolation', 'Warpping', 'Masking', 'Filename Prefix',  add_parameters{1,:});
uvascim.uvascim.image.reco.paramQuantif = ParamConfig;

% Update now the acq field: based on the fixed structure but with info
% form the original file (the moving structure)
uvascim.uvascim.image.acq= fixed.uvascim.image.acq;
uvascim.uvascim.image.acq.ppl_name = moving.uvascim.image.acq.ppl_name;
if isfield(uvascim.uvascim.image.acq, 'flip_angle')
uvascim.uvascim.image.acq.flip_angle = moving.uvascim.image.acq.flip_angle;
end
uvascim.uvascim.image.acq.echotime = moving.uvascim.image.acq.echotime;
uvascim.uvascim.image.acq.champ_magnetique = moving.uvascim.image.acq.champ_magnetique;
uvascim.uvascim.image.acq.study = moving.uvascim.image.acq.study;
uvascim.uvascim.image.acq.subject = moving.uvascim.image.acq.subject;
uvascim.uvascim.image.acq.tr = moving.uvascim.image.acq.tr;

% load nii file(s), reorient them and store them in the uvascim structure
V =spm_vol(fullfile([PATHSTR,filesep, [add_parameters{1,size(add_parameters,2)-1} NAME_scan_to_coreg], '.nii']));
tmp = spm_read_vols(V);
if numel(size(squeeze(moving.uvascim.image.reco.data))) < 5
    if  size(moving.uvascim.image.reco.data,3) >  size(moving.uvascim.image.reco.data,5)
        uvascim.uvascim.image.reco.data = permute(tmp, [1 2 4 3]);
        uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4]);
    else
        uvascim.uvascim.image.reco.data = permute(tmp, [1 2 5 3 4]);
        uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4 5]);
    end
else
    % 4D nii file to 5D matlab matrix
    uvascim.uvascim.image.reco.data = reshape(tmp, [size(tmp,1), size(tmp,2),size(tmp,3),...
        size(moving.uvascim.image.reco.data,3), size(moving.uvascim.image.reco.data,5)]);
    uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [1 2 4 3 5]);
    uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4 5]);
    
end
% end

uvascim.uvascim.image.reco.data = flip(uvascim.uvascim.image.reco.data, 1);
uvascim.uvascim.image.reco.data = flip(uvascim.uvascim.image.reco.data, 2);
if exist('Other_maps', 'var')
    uvascim(2:size(Other_maps,2)+1) = Other_maps;
    for i=1:numel(NAME_scan_other)
        uvascim(i+1).uvascim.image.reco.angAP = fixed.uvascim.image.reco.angAP;
        uvascim(i+1).uvascim.image.reco.angFH = fixed.uvascim.image.reco.angFH;
        uvascim(i+1).uvascim.image.reco.angRL = fixed.uvascim.image.reco.angRL;
        uvascim(i+1).uvascim.image.reco.angulation = fixed.uvascim.image.reco.angulation;
        uvascim(i+1).uvascim.image.reco.fov = fixed.uvascim.image.reco.fov;
        uvascim(i+1).uvascim.image.reco.no_slices = fixed.uvascim.image.reco.no_slices;
        uvascim(i+1).uvascim.image.reco.no_samples = fixed.uvascim.image.reco.no_samples;
        uvascim(i+1).uvascim.image.reco.no_views = fixed.uvascim.image.reco.no_views;
        uvascim(i+1).uvascim.image.reco.thickness = fixed.uvascim.image.reco.thickness;
        
        % fov_offsets = [3 echoes slice repet]
        uvascim(i+1).uvascim.image.reco.fov_offsets=repmat(fixed.uvascim.image.reco.fov_offsets(:,1,:,1),[1,uvascim(i+1).uvascim.image.reco.no_echoes,1,uvascim(i+1).uvascim.image.reco.no_expts]);
        % fov_orientation = [9 echoes slice repet]
        uvascim(i+1).uvascim.image.reco.fov_orientation = repmat(fixed.uvascim.image.reco.fov_orientation(:,1,:,1),[1,uvascim(i+1).uvascim.image.reco.no_echoes,1,uvascim(i+1).uvascim.image.reco.no_expts]);
        % fov_phase_orientation = [echoes slice repet]
        uvascim(i+1).uvascim.image.reco.fov_phase_orientation = repmat(fixed.uvascim.image.reco.fov_phase_orientation(1,:,1),[uvascim(i+1).uvascim.image.reco.no_echoes,1,uvascim(i+1).uvascim.image.reco.no_expts]);
        % scaling_factor = [echoes slice repet]
        uvascim(i+1).uvascim.image.reco.scaling_factor = repmat(fixed.uvascim.image.reco.scaling_factor(:,1,1),[uvascim(i+1).uvascim.image.reco.no_echoes,uvascim(i+1).uvascim.image.reco.no_slices,uvascim(i+1).uvascim.image.reco.no_expts]);
        % scaling_offset = [echoes slice repet]
        uvascim(i+1).uvascim.image.reco.scaling_offset = repmat(fixed.uvascim.image.reco.scaling_offset(:,1,1),[uvascim(i+1).uvascim.image.reco.no_echoes,uvascim(i+1).uvascim.image.reco.no_slices,uvascim(i+1).uvascim.image.reco.no_expts]);
        % phaselabel = [echoes slice repet]
        uvascim(i+1).uvascim.image.reco.phaselabel = repmat(fixed.uvascim.image.reco.phaselabel(:,1,1),[uvascim(i+1).uvascim.image.reco.no_echoes,uvascim(i+1).uvascim.image.reco.no_slices,uvascim(i+1).uvascim.image.reco.no_expts]);
        % label = [echoes slice repet]
        uvascim(i+1).uvascim.image.reco.label = repmat(fixed.uvascim.image.reco.label(:,1,1),[uvascim(i+1).uvascim.image.reco.no_echoes,uvascim(i+1).uvascim.image.reco.no_slices,uvascim(i+1).uvascim.image.reco.no_expts]);
        
        
        ParamConfig=sprintf('##$SPM_coreg_write algorithum=\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n##$Name of each parameters=\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n##$Parameters set up%s\n',...
            'Ref day', 'Ref Scan', 'Interpolation', 'Warpping', 'Masking', 'Filename Prefix',  add_parameters{1,:});
        uvascim(i+1).uvascim.image.reco.paramQuantif = ParamConfig;
        
        % Update now the acq field: based on the fixed structure but with info
        % form the original file (the Other_maps structure)
        uvascim(i+1).uvascim.image.acq= fixed.uvascim.image.acq;
        uvascim(i+1).uvascim.image.acq.ppl_name = Other_maps(i).uvascim.image.acq.ppl_name;
        if isfield(uvascim(i+1).uvascim.image.acq, 'flip_angle')
        uvascim(i+1).uvascim.image.acq.flip_angle = Other_maps(i).uvascim.image.acq.flip_angle;
        end
        uvascim(i+1).uvascim.image.acq.echotime = Other_maps(i).uvascim.image.acq.echotime;
        uvascim(i+1).uvascim.image.acq.champ_magnetique = Other_maps(i).uvascim.image.acq.champ_magnetique;
        uvascim(i+1).uvascim.image.acq.study = Other_maps(i).uvascim.image.acq.study;
        uvascim(i+1).uvascim.image.acq.subject = Other_maps(i).uvascim.image.acq.subject;
        uvascim(i+1).uvascim.image.acq.tr = Other_maps(i).uvascim.image.acq.tr;
        
        % load nii file(s), reorient them and store them in the uvascim structure
        V =spm_vol(fullfile([ PATHSTR,filesep, [add_parameters{1,size(add_parameters,2)-1}  NAME_scan_other{i}], '.nii']));
        tmp = spm_read_vols(V);
        if numel(size(squeeze(Other_maps(i).uvascim.image.reco.data))) < 5
            if  size(Other_maps(i).uvascim.image.reco.data,3) >  size(Other_maps(i).uvascim.image.reco.data,5)
                uvascim(i+1).uvascim.image.reco.data = permute(tmp, [1 2 4 3]);
                uvascim(i+1).uvascim.image.reco.data = permute(uvascim(i+1).uvascim.image.reco.data, [2 1 3 4]);
            else
                uvascim(i+1).uvascim.image.reco.data = permute(tmp, [1 2 5 3 4]);
                uvascim(i+1).uvascim.image.reco.data = permute(uvascim(i+1).uvascim.image.reco.data, [2 1 3 4 5]);
            end
        else
            % 4D nii file to 5D matlab matrix
            uvascim(i+1).uvascim.image.reco.data = reshape(tmp, [size(tmp,1), size(tmp,2),size(tmp,3),...
                size(Other_maps(i).uvascim.image.reco.data,3), size(Other_maps(i).uvascim.image.reco.data,5)]);
            uvascim(i+1).uvascim.image.reco.data = permute(uvascim(i+1).uvascim.image.reco.data, [1 2 4 3 5]);
            uvascim(i+1).uvascim.image.reco.data = permute(uvascim(i+1).uvascim.image.reco.data, [2 1 3 4 5]);
        end
        % end
        
        uvascim(i+1).uvascim.image.reco.data = flip(uvascim(i+1).uvascim.image.reco.data, 1);
        uvascim(i+1).uvascim.image.reco.data = flip(uvascim(i+1).uvascim.image.reco.data, 2);
        
    end
end
Coreg_maps = 1;




