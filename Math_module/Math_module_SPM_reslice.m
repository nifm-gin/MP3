function [moved_maps, uvascim]= Math_module_SPM_reslice(fixed_map_filename, Scan_to_reslice_filename, add_parameters)

%% use the Coregister function of SPM12 to reslice scans
% load the  scan to reslice (if not empty)
Scan_to_reslice_filename = Scan_to_reslice_filename{:};
fid=fopen(Scan_to_reslice_filename ,'r');
if fid>0
    fclose(fid);
    moved = load(Scan_to_reslice_filename);
    [PATHSTR, NAME_scan_to_reslice, ~] = fileparts(Scan_to_reslice_filename);
    NAME_scan_to_reslice = [NAME_scan_to_reslice '-for_spm'];
    if numel(size(squeeze(moved.uvascim.image.reco.data))) < 5
         NHdr=CreateNiiHdr_from_uvascim(Scan_to_reslice_filename, moved, NAME_scan_to_reslice);
        if  size(moved.uvascim.image.reco.data,3) >  size(moved.uvascim.image.reco.data,5)
            data_4d_to_reslice = permute(moved.uvascim.image.reco.data, [1 2 4 3]);
        else
            data_4d_to_reslice = squeeze(moved.uvascim.image.reco.data(:,:,1,:,:));
        end
        WriteNii(fullfile(PATHSTR, [NAME_scan_to_reslice '.nii']),NHdr, data_4d_to_reslice)
    else
        % create 4D nii from 5D scan by adding the repetition after the echo
         data_4d_to_reslice = permute(moved.uvascim.image.reco.data, [1 2 4 3 5]);
         data_4d_to_reslice = reshape(data_4d_to_reslice, [size(data_4d_to_reslice, 1), size(data_4d_to_reslice, 2), size(data_4d_to_reslice, 3), size(data_4d_to_reslice, 4)*size(data_4d_to_reslice, 5)]);
        % 
         NHdr=CreateNiiHdr_from_uvascim(Scan_to_reslice_filename, moved, NAME_scan_to_reslice);
         WriteNii(fullfile(PATHSTR, [NAME_scan_to_reslice '.nii']),NHdr, data_4d_to_reslice)
    end
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
        fixed_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    moved_maps = [];
    return
end

final_resolution = add_parameters(1,3);
% load the fixed_map_filename file
fid=fopen(fixed_map_filename ,'r');
if fid>0
    fclose(fid);
    fixed = load(fixed_map_filename);
    [PATHSTR, NAME_fixed, ~] = fileparts(fixed_map_filename);
    NAME_fixed = [NAME_fixed '-for_spm'];
    if strcmp(final_resolution, 'Unchanged')
        for echo=1:size(fixed.uvascim.image.reco.data,3)
            for slice = 1:size(fixed.uvascim.image.reco.data,4)
                for repet = 1:size(fixed.uvascim.image.reco.data,5)
                    tmp(:,:,echo,slice,repet)= imresize(squeeze(fixed.uvascim.image.reco.data(:,:,echo,slice,repet)),...
                        [size(moved.uvascim.image.reco.data,1) size(moved.uvascim.image.reco.data,2)]);
                end
            end
        end
        fixed.uvascim.image.reco = rmfield(fixed.uvascim.image.reco, 'data');
        fixed.uvascim.image.reco.data = tmp;
        % update uvasc structure
        fixed.uvascim.image.acq.no_samples = moved.uvascim.image.acq.no_samples;
        fixed.uvascim.image.acq.no_views = moved.uvascim.image.acq.no_views;
        if ~isfield(fixed.uvascim.image.acq, 'pix_spacing')
            % maybe needs to perform this calcul for everyone!
            fixed.uvascim.image.acq.pix_spacing = [fixed.uvascim.image.acq.fov(1)/size(moved.uvascim.image.reco.data,2) fixed.uvascim.image.acq.fov(2)/size(moved.uvascim.image.reco.data,1)]';
        else
             if ~isfield(moved.uvascim.image.acq, 'pix_spacing')
            fixed.uvascim.image.acq.pix_spacing = [moved.uvascim.image.acq.fov(1)/size(moved.uvascim.image.reco.data,2) moved.uvascim.image.acq.fov(2)/size(moved.uvascim.image.reco.data,1)]';
             else
            fixed.uvascim.image.acq.pix_spacing = moved.uvascim.image.acq.pix_spacing;
             end
        end
        fixed.uvascim.image.reco.no_samples = moved.uvascim.image.reco.no_samples;
        fixed.uvascim.image.reco.no_views = moved.uvascim.image.reco.no_views;
    end
    NHdr=CreateNiiHdr_from_uvascim(fixed_map_filename, fixed, NAME_fixed);
    if numel(size(squeeze(fixed.uvascim.image.reco.data))) < 5
        if  size(fixed.uvascim.image.reco.data,3) >  size(fixed.uvascim.image.reco.data,5)
            data_4d_of_fixed_map = permute(fixed.uvascim.image.reco.data, [1 2 4 3]);
        else
            data_4d_of_fixed_map = squeeze(fixed.uvascim.image.reco.data(:,:,1,:,:));
        end
        WriteNii(fullfile(PATHSTR, [NAME_fixed '.nii']), NHdr, data_4d_of_fixed_map)
    else
        % create 4D nii from 5D scan
    end
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
        fixed_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    moved_maps = [];
    return
end

%% perform co-reg alogrithum (using SPM 12)
matlabbatch{1}.spm.spatial.coreg.write.ref = {[NAME_fixed '.nii,1']};
if ~isempty(Scan_to_reslice_filename)
    for xx = 1:size(data_4d_to_reslice,4)
        if ~exist('Scan_to_reslice_nii_filename', 'var')
            Scan_to_reslice_nii_filename{1,:} =  [NAME_scan_to_reslice '.nii,1'];
        else
            Scan_to_reslice_nii_filename{size(Scan_to_reslice_nii_filename,1)+1,:} =  [NAME_scan_to_reslice '.nii,' num2str(xx)]; %#ok<AGROW>
        end
    end
    matlabbatch{1}.spm.spatial.coreg.write.source = Scan_to_reslice_nii_filename;
end
%define options
% Type of interpolation
switch add_parameters{1,4}
    case 'Nearest neighbour'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    case 'Trilinear'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 7;
end
%  Type of Warpping
switch add_parameters{1,5}
    case 'No wrap'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    case 'Warp X'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [1 1 1];
end
%  Mask?
switch add_parameters{1,6}
    case 'Dont mask images'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask= 0;
    case 'Mask image'
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 1;
end

matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = add_parameters{1,size(add_parameters,2)};

[SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch

date_str = date;
if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
    movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
        [PATHSTR, filesep, NAME_scan_to_reslice(1:end-8) '_SPM_reslice.ps']); 
end
close(SPMinter)
close(SPMgraph)

% create output structure
% update first output structure and in particular the reco field: 
% based on the original structure (the moved structure) but with the info
% form the fixed structure
uvascim = moved.uvascim;
uvascim.image.reco.angAP = fixed.uvascim.image.reco.angAP;
uvascim.image.reco.angFH = fixed.uvascim.image.reco.angFH;
uvascim.image.reco.angRL = fixed.uvascim.image.reco.angRL;
uvascim.image.reco.angulation = fixed.uvascim.image.reco.angulation;
uvascim.image.reco.fov = fixed.uvascim.image.reco.fov;
uvascim.image.reco.no_slices = fixed.uvascim.image.reco.no_slices;
uvascim.image.reco.no_samples = fixed.uvascim.image.reco.no_samples;
uvascim.image.reco.no_views = fixed.uvascim.image.reco.no_views;
uvascim.image.reco.thickness = fixed.uvascim.image.reco.thickness;

% fov_offsets = [3 echoes slice repet]
uvascim.image.reco.fov_offsets=repmat(fixed.uvascim.image.reco.fov_offsets(:,1,:,1),[1,uvascim.image.reco.no_echoes,1,uvascim.image.reco.no_expts]);
% fov_orientation = [9 echoes slice repet]
uvascim.image.reco.fov_orientation = repmat(fixed.uvascim.image.reco.fov_orientation(:,1,:,1),[1,uvascim.image.reco.no_echoes,1,uvascim.image.reco.no_expts]);
% fov_phase_orientation = [echoes slice repet]
uvascim.image.reco.fov_phase_orientation = repmat(fixed.uvascim.image.reco.fov_phase_orientation(1,:,1),[uvascim.image.reco.no_echoes,1,uvascim.image.reco.no_expts]);
% scaling_factor = [echoes slice repet]
uvascim.image.reco.scaling_factor = repmat(fixed.uvascim.image.reco.scaling_factor(:,1,1),[uvascim.image.reco.no_echoes,uvascim.image.reco.no_slices,uvascim.image.reco.no_expts]);
% scaling_offset = [echoes slice repet]
uvascim.image.reco.scaling_offset = repmat(fixed.uvascim.image.reco.scaling_offset(:,1,1),[uvascim.image.reco.no_echoes,uvascim.image.reco.no_slices,uvascim.image.reco.no_expts]);
% phaselabel = [echoes slice repet]
uvascim.image.reco.phaselabel = repmat(fixed.uvascim.image.reco.phaselabel(:,1,1),[uvascim.image.reco.no_echoes,uvascim.image.reco.no_slices,uvascim.image.reco.no_expts]);
% label = [echoes slice repet]
uvascim.image.reco.label = repmat(fixed.uvascim.image.reco.label(:,1,1),[uvascim.image.reco.no_echoes,uvascim.image.reco.no_slices,uvascim.image.reco.no_expts]);


ParamConfig=sprintf('##$SPM_coreg_write algorithum=\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n##$Name of each parameters=\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n##$Parameters set up%s\n',...
    'Ref day', 'Ref Scan', 'Interpolation', 'Warpping', 'Masking', 'Filename Prefix',  add_parameters{1,:});
uvascim.image.reco.paramQuantif = ParamConfig;

% Update now the acq field: based on the fixed structure but with info
% form the original file (the moved structure)
uvascim.image.acq= fixed.uvascim.image.acq;
uvascim.image.acq.ppl_name = moved.uvascim.image.acq.ppl_name;
if isfield(uvascim.image.acq, 'flip_angle')
    uvascim.image.acq.flip_angle = moved.uvascim.image.acq.flip_angle;
end
uvascim.image.acq.echotime = moved.uvascim.image.acq.echotime;
uvascim.image.acq.champ_magnetique = moved.uvascim.image.acq.champ_magnetique;
uvascim.image.acq.study = moved.uvascim.image.acq.study;
uvascim.image.acq.subject = moved.uvascim.image.acq.subject;
uvascim.image.acq.tr = moved.uvascim.image.acq.tr;

% load nii file(s), reorient them and store them in the uvascim structure
V =spm_vol(fullfile([ PATHSTR,filesep, [add_parameters{1,size(add_parameters,2)} NAME_scan_to_reslice], '.nii']));
tmp = spm_read_vols(V);
if numel(size(squeeze(moved.uvascim.image.reco.data))) < 5
    if  size(moved.uvascim.image.reco.data,3) >  size(moved.uvascim.image.reco.data,5)
        uvascim.image.reco.data = permute(tmp, [1 2 4 3]);
        uvascim.image.reco.data = permute(uvascim.image.reco.data, [2 1 3 4]);
    else
        uvascim.image.reco.data = permute(tmp, [1 2 5 3 4]);
        uvascim.image.reco.data = permute(uvascim.image.reco.data, [2 1 3 4 5]);
    end
else
    % 4D nii file to 5D matlab matrix
    uvascim.image.reco.data = reshape(tmp, [size(tmp,1), size(tmp,2),size(tmp,3),...
        size(moved.uvascim.image.reco.data,3), size(moved.uvascim.image.reco.data,5)]);
    uvascim.image.reco.data = permute(uvascim.image.reco.data, [1 2 4 3 5]);
    uvascim.image.reco.data = permute(uvascim.image.reco.data, [2 1 3 4 5]);
end
% end

uvascim.image.reco.data = flip(uvascim.image.reco.data, 1);
uvascim.image.reco.data = flip(uvascim.image.reco.data, 2);


moved_maps = 1;
 
   


