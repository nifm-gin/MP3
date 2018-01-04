function [realign_and_coreg_maps, uvascim]= Math_module_SPM_realign_and_coreg(fixed_map_filename, Scan_to_realign_and_coreg_filename, add_parameters)

%% use the Coregister function of SPM12 to coreg scans

% load the fixed_map_filename file
fid=fopen(fixed_map_filename ,'r');
if fid>0   
    fclose(fid);
    fixed = load(fixed_map_filename);  
    [PATHSTR, NAME_fixed, ~] = fileparts(fixed_map_filename);
    NAME_fixed = [NAME_fixed '-for_spm'];
    NHdr=CreateNiiHdr_from_uvascim(fixed_map_filename, fixed, NAME_fixed);
    if numel(size(squeeze(fixed.uvascim.image.reco.data))) < 5
        if  size(fixed.uvascim.image.reco.data,3) >  size(fixed.uvascim.image.reco.data,5)
            data_4d = permute(fixed.uvascim.image.reco.data, [1 2 4 3]);
        else
            data_4d = squeeze(fixed.uvascim.image.reco.data(:,:,1,:,:));
        end
        WriteNii(fullfile(PATHSTR, [NAME_fixed '.nii']), NHdr, data_4d)
    else
        % create 4D nii from 5D scan
    end
else
    warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
        fixed_map_filename);
    msgbox(warning_text, 'Registration map warning') ;
    realign_and_coreg_maps = [];
    return
end

% load the  scan to coreg (if not empty)
for i=1:numel(Scan_to_realign_and_coreg_filename)
    tpm_filename = Scan_to_realign_and_coreg_filename{i};
    fid=fopen(tpm_filename ,'r');
    if fid>0
        fclose(fid);
        data = load(tpm_filename);
        moved(i) = data;
        [PATHSTR, tmp_name, ~] = fileparts(tpm_filename);
        NAME_scan_to_coreg(i,1) = {[tmp_name '-for_spm']}; %#ok<AGROW>
        if numel(size(squeeze(moved(i).uvascim.image.reco.data))) < 5
            NHdr=CreateNiiHdr_from_uvascim(tpm_filename, moved(i), NAME_scan_to_coreg{i,1});
            if  size(moved(i).uvascim.image.reco.data,3) >  size(moved(i).uvascim.image.reco.data,5)
                data_4d = permute(moved(i).uvascim.image.reco.data, [1 2 4 3]);
            else
                data_4d = squeeze(moved(i).uvascim.image.reco.data(:,:,1,:,:));
            end
            WriteNii(fullfile(PATHSTR, [NAME_scan_to_coreg{i,1} '.nii']),NHdr, data_4d)
            
        else
            % create 4D nii from 5D scan by adding the repetition after the echo
            data_4d = permute(moved(i).uvascim.image.reco.data, [1 2 4 3 5]);
            data_4d = reshape(data_4d, [size(data_4d, 1), size(data_4d, 2), size(data_4d, 3), size(data_4d, 4)*size(data_4d, 5)]);
            %
            NHdr=CreateNiiHdr_from_uvascim(tpm_filename, moved(i), NAME_scan_to_coreg{i,1});
            WriteNii(fullfile(PATHSTR, [NAME_scan_to_coreg{i,1} '.nii']),NHdr, data_4d)
            
        end
        for xx = 1:size(data_4d,4)
                tmp_name2{xx,:} =  [ NAME_scan_to_coreg{i,1} '.nii,' num2str(xx)];
        end
        estwrite_data(i) = {tmp_name2};
        clear tmp_name2
    else
        warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
            fixed_map_filename);
        msgbox(warning_text, 'Registration map warning') ;
        realign_and_coreg_maps = [];
        return
    end
end


%define options
%%
matlabbatch{1}.spm.spatial.realign.estwrite.data = estwrite_data;

%%
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

[SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm defaults fmri
spm_jobman initcfg
% spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch 

date_str = date;
if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
    movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
        [PATHSTR, filesep, NAME_fixed(1:end-8) '_SPM_realign_and_coreg(1).ps']); 
end

% Co-reg (estimage/write)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(PATHSTR, [NAME_fixed '.nii,1'])};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(PATHSTR,['mean' estwrite_data{1}{1}])};
other_files = {};
for i=1:numel(NAME_scan_to_coreg)
    aa = spm_select('ExtList', PATHSTR, ['r' NAME_scan_to_coreg{i} '.nii'], 1:999);
    other_files = [other_files cellstr(aa)'];
end
matlabbatch{1}.spm.spatial.coreg.estwrite.other = other_files';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch

if exist([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'], 'file') == 2
    movefile([PATHSTR filesep 'spm_' date_str(end-3:end) date_str(end-7:end-5) date_str(1:end-9) '.ps'],...
         [PATHSTR, filesep, NAME_fixed(1:end-8) '_SPM_realign_and_coreg(2).ps']); 
end
close(SPMinter)
close(SPMgraph)

% create output structure
% update first output structure and in particular the reco field: 
% based on the original structure (the moved structure) but with the info
% form the fixed structure
uvascim = moved;
for i=1:numel(Scan_to_realign_and_coreg_filename)
    uvascim(i).uvascim.image.reco.angAP = fixed.uvascim.image.reco.angAP;
    uvascim(i).uvascim.image.reco.angFH = fixed.uvascim.image.reco.angFH;
    uvascim(i).uvascim.image.reco.angRL = fixed.uvascim.image.reco.angRL;
    uvascim(i).uvascim.image.reco.angulation = fixed.uvascim.image.reco.angulation;
    uvascim(i).uvascim.image.reco.fov = fixed.uvascim.image.reco.fov;
    uvascim(i).uvascim.image.reco.no_slices = fixed.uvascim.image.reco.no_slices;
    uvascim(i).uvascim.image.reco.no_samples = fixed.uvascim.image.reco.no_samples;
    uvascim(i).uvascim.image.reco.no_views = fixed.uvascim.image.reco.no_views;
    uvascim(i).uvascim.image.reco.thickness = fixed.uvascim.image.reco.thickness;
    
    % fov_offsets = [3 echoes slice repet]
    uvascim(i).uvascim.image.reco.fov_offsets=repmat(fixed.uvascim.image.reco.fov_offsets(:,1,:,1),[1,uvascim(i).uvascim.image.reco.no_echoes,1,uvascim(i).uvascim.image.reco.no_expts]);
    % fov_orientation = [9 echoes slice repet]
    uvascim(i).uvascim.image.reco.fov_orientation = repmat(fixed.uvascim.image.reco.fov_orientation(:,1,:,1),[1,uvascim(i).uvascim.image.reco.no_echoes,1,uvascim(i).uvascim.image.reco.no_expts]);
    % fov_phase_orientation = [echoes slice repet]
    uvascim(i).uvascim.image.reco.fov_phase_orientation = repmat(fixed.uvascim.image.reco.fov_phase_orientation(1,:,1),[uvascim(i).uvascim.image.reco.no_echoes,1,uvascim(i).uvascim.image.reco.no_expts]);
    % scaling_factor = [echoes slice repet]
    uvascim(i).uvascim.image.reco.scaling_factor = repmat(fixed.uvascim.image.reco.scaling_factor(:,1,1),[uvascim(i).uvascim.image.reco.no_echoes,uvascim(i).uvascim.image.reco.no_slices,uvascim(i).uvascim.image.reco.no_expts]);
    % scaling_offset = [echoes slice repet]
    uvascim(i).uvascim.image.reco.scaling_offset = repmat(fixed.uvascim.image.reco.scaling_offset(:,1,1),[uvascim(i).uvascim.image.reco.no_echoes,uvascim(i).uvascim.image.reco.no_slices,uvascim(i).uvascim.image.reco.no_expts]);
    % phaselabel = [echoes slice repet]
    uvascim(i).uvascim.image.reco.phaselabel = repmat(fixed.uvascim.image.reco.phaselabel(:,1,1),[uvascim(i).uvascim.image.reco.no_echoes,uvascim(i).uvascim.image.reco.no_slices,uvascim(i).uvascim.image.reco.no_expts]);
    % label = [echoes slice repet]
    uvascim(i).uvascim.image.reco.label = repmat(fixed.uvascim.image.reco.label(:,1,1),[uvascim(i).uvascim.image.reco.no_echoes,uvascim(i).uvascim.image.reco.no_slices,uvascim(i).uvascim.image.reco.no_expts]);
    
    
    ParamConfig=sprintf('##$SPM_coreg_write algorithum=\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n##$Name of each parameters=\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n##$Parameters set up%s\n',...
        'Ref day', 'Ref Scan', 'Interpolation', 'Warpping', 'Masking', 'Filename Prefix',  add_parameters{1,:});
    uvascim(i).uvascim.image.reco.paramQuantif = ParamConfig;
    
    % Update now the acq field: based on the fixed structure but with info
    % form the original file (the moved structure)
    uvascim(i).uvascim.image.acq= fixed.uvascim.image.acq;
    uvascim(i).uvascim.image.acq.ppl_name = moved(i).uvascim.image.acq.ppl_name;
    uvascim(i).uvascim.image.acq.flip_angle = moved(i).uvascim.image.acq.flip_angle;
    uvascim(i).uvascim.image.acq.echotime = moved(i).uvascim.image.acq.echotime;
    uvascim(i).uvascim.image.acq.champ_magnetique = moved(i).uvascim.image.acq.champ_magnetique;
    uvascim(i).uvascim.image.acq.study = moved(i).uvascim.image.acq.study;
    uvascim(i).uvascim.image.acq.subject = moved(i).uvascim.image.acq.subject;
    uvascim(i).uvascim.image.acq.tr = moved(i).uvascim.image.acq.tr;
    
    % load nii file(s), reorient them and store them in the uvascim structure
    V =spm_vol(fullfile([ PATHSTR,filesep, ['rr' NAME_scan_to_coreg{i}], '.nii']));
    tmp = spm_read_vols(V);
    if numel(size(squeeze(moved(i).uvascim.image.reco.data))) < 5
        if  size(moved(i).uvascim.image.reco.data,3) >  size(moved(i).uvascim.image.reco.data,5)
            uvascim(i).uvascim.image.reco.data = permute(tmp, [1 2 4 3]);
            uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [2 1 3 4]);
        else
            uvascim(i).uvascim.image.reco.data = permute(tmp, [1 2 5 3 4]);
            uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [2 1 3 4 5]);
        end
    else
        % 4D nii file to 5D matlab matrix
        uvascim(i).uvascim.image.reco.data = reshape(tmp, [size(tmp,1), size(tmp,2),size(tmp,3),...
            size(moved(i).uvascim.image.reco.data,3), size(moved(i).uvascim.image.reco.data,5)]);
        uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [1 2 4 3 5]);
        uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [2 1 3 4 5]);
    end
    % end
    
    uvascim(i).uvascim.image.reco.data = flip(uvascim(i).uvascim.image.reco.data, 1);
    uvascim(i).uvascim.image.reco.data = flip(uvascim(i).uvascim.image.reco.data, 2);
    
end

realign_and_coreg_maps = 1;
 
   


