function [moved_maps, uvascim]= Math_module_SPM_realign(Scan_to_realign, add_parameters)

%% use the Coregister function of SPM12 to reslice scans

% load the  scan to realign (if not empty)
for i=1:numel(Scan_to_realign)
    Scan_to_realign_filename = Scan_to_realign{i};
    fid=fopen(Scan_to_realign_filename ,'r');
    if fid>0
        fclose(fid);
        moved = load(Scan_to_realign_filename);
        [PATHSTR, NAME_scan_to_realign, ~] = fileparts(Scan_to_realign_filename);
        NAME_scan_to_realign = [NAME_scan_to_realign '-for_spm'];
        if numel(size(squeeze(moved.uvascim.image.reco.data))) < 5
            NHdr=CreateNiiHdr_from_uvascim(Scan_to_realign_filename, moved, NAME_scan_to_realign);
            if  size(moved.uvascim.image.reco.data,3) >  size(moved.uvascim.image.reco.data,5)
                % 4D multi echo
                data_4d = permute(moved.uvascim.image.reco.data, [1 2 4 3]);
            else
                % 4D multi expt
                data_4d = squeeze(moved.uvascim.image.reco.data(:,:,1,:,:));
            end
            WriteNii(fullfile(PATHSTR, [NAME_scan_to_realign '.nii']),NHdr, data_4d)
            
        else
            % create 4D nii from 5D scan by adding the repetition after the echo
            data_4d = permute(moved.uvascim.image.reco.data, [1 2 4 3 5]);
            data_4d = reshape(data_4d, [size(data_4d, 1), size(data_4d, 2), size(data_4d, 3), size(data_4d, 4)*size(data_4d, 5)]);
            %
            NHdr=CreateNiiHdr_from_uvascim(Scan_to_realign_filename, moved, NAME_scan_to_realign);
            WriteNii(fullfile(PATHSTR, [NAME_scan_to_realign '.nii']),NHdr, data_4d)
            
        end
    else
        warning_text = sprintf('##$ Can not perform the registration because there is\n##$ Somthing wrong with the reference data\n##$Ref scan=%s\n##$',...
            fixed_map_filename);
        msgbox(warning_text, 'Registration map warning') ;
        moved_maps = [];
        return
    end
end

%remove dyn during the bolus 
if strcmp(add_parameters{1,end}, 'Yes')
 % do not realign images (dynamics) during the bolus
    mean_signal=squeeze(mean(mean(mean(data_4d))));
    mean_beg = mean(mean_signal(3:8)); sd_beg = std(mean_signal(3:8));
    mean_end = mean(mean_signal(end-5:end)); sd_end = std(mean_signal(end-5:end));
    
    allscans = 1:size(data_4d,4);
    peak =  find(mean_signal==min(mean_signal));
    % detect the bolus 
    % 1) scan before the peak with a signal < the mean_beg-2*sd_beg
    % 2) scan after the peak with a signal > the mean_end-12*sd_end
    ignorescans = find(mean_signal(1:peak) < mean_beg-2*sd_beg, 1 ) : (size(data_4d,4) -(size(mean_signal(peak:end),1) - find(mean_signal(peak:end) > mean_end-12*sd_end, 1 )));
    keepscans = setdiff(allscans,ignorescans);
    data_4d_old = data_4d;
    data_4d = data_4d(:,:,:,keepscans);
end
%% perform co-reg alogrithum (using SPM 12)


if ~isempty(Scan_to_realign)
    for i = 1:numel(Scan_to_realign)
        [PATHSTR, NAME_scan_to_realign, ~] = fileparts(Scan_to_realign{i});
        NAME_scan_to_realign = [NAME_scan_to_realign '-for_spm'];
        if moved.uvascim.image.reco.no_slices > 1
            for xx = 1:size(data_4d,4)
                if ~exist('Scan_to_realign_nii_filename', 'var')
                    Scan_to_realign_nii_filename{1,:} = [fullfile(PATHSTR, NAME_scan_to_realign) '.nii,1'];
                else
                    Scan_to_realign_nii_filename{size(Scan_to_realign_nii_filename,1)+1,:} =  strcat([fullfile(PATHSTR, NAME_scan_to_realign) '.nii,'], num2str(xx));
                end
            end
        else
            for xx = 1:size(data_4d,3)
                if ~exist('Scan_to_realign_nii_filename', 'var')
                    Scan_to_realign_nii_filename{1,:} = [fullfile(PATHSTR, NAME_scan_to_realign) '.nii,1'];
                else
                    Scan_to_realign_nii_filename{size(Scan_to_realign_nii_filename,1)+1,:} =  strcat([fullfile(PATHSTR, NAME_scan_to_realign) '.nii,'], num2str(xx));
                end
            end
            
        end
    end
    matlabbatch{1}.spm.spatial.realign.estwrite.data  = {Scan_to_realign_nii_filename};
end

%define options
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality =str2double(add_parameters(1));
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = str2double(add_parameters(2));
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = str2double(add_parameters(3));
switch add_parameters{1,4}
    case  'Register to first'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; 
    case 'Register to mean'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
end

% Type of interpolation
switch add_parameters{1,5}
    case 'Trilinear'
       matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 3;
    case '4th Degree B-Spline'
       matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 7;
end
%  Type of Warpping
switch add_parameters{1,6}
    case 'No Wrap'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    case 'Warp X'
       matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap  = [1 1 1];
end
% EstiOption: Weight
switch add_parameters{1,7}
    case '0 files'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    otherwise
         matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = add_parameters(1,7);
end

% 'ResliceOption: Res Image'
%  Mask?
switch add_parameters{1,8}
    case 'All files (1..n)'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0];
    case 'Images 2..n,'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [1 0];
    case 'All Images + Mean Image'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    case  'Mean Images Only'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];   
end

%'ResliceOption: Interpol'
switch add_parameters{1,9}
    case 'Nearest neighbour'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 0;
    case 'Trilinear'
       matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 3;
    case '4th Degree B-Spline'
       matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 7;
end

switch add_parameters{1,10}
    case 'No Wrap'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    case 'Warp X'
       matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap  = [1 1 1];
end

switch add_parameters{1,11}
    case 'Dont mask images'
       matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0;
    case 'Mask images'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
end
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = add_parameters{1,size(add_parameters,2)-1};

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
        [PATHSTR, filesep, NAME_scan_to_realign(1:end-8) '_SPM_realign.ps']); 
end

close(SPMinter)
close(SPMgraph)

% load nii file(s), reorient them and store them in the uvascim structure
if size(Scan_to_realign,1) == 1
    %     if numel(size(squeeze(moved.uvascim.image.reco.data))) < 5
    %         sprintf('##$ Not coded yet');
    %     else
    if strcmp(add_parameters{1,end}, 'Yes')
        tmp = data_4d_old;
        for i = 1:numel(keepscans)
            V =spm_vol(strcat([ PATHSTR,filesep, [add_parameters{1,size(add_parameters,2)-1} NAME_scan_to_realign], '.nii,'], num2str(i)));
            tmp(:,:,:,keepscans(i)) = spm_read_vols(V);
        end
        for i = 1:numel(ignorescans)
            V =spm_vol(strcat([ PATHSTR,filesep, [NAME_scan_to_realign], '.nii,'], num2str(ignorescans(i))));
            tmp(:,:,:,ignorescans(i)) = spm_read_vols(V);
        end
    else
        for i = 1:size(Scan_to_realign_nii_filename,1)
            V =spm_vol(strcat([ PATHSTR,filesep, [add_parameters{1,size(add_parameters,2)-1} NAME_scan_to_realign], '.nii,'], num2str(i)));
            tmp(:,:,:,i) = spm_read_vols(V);
        end
    end
    uvascim.uvascim = moved.uvascim;
    
    if numel(size(squeeze(uvascim.uvascim.image.reco.data))) < 5
        if  size(uvascim.uvascim.image.reco.data,3) >  size(uvascim.uvascim.image.reco.data,5)
            % 4D multi echo
            uvascim.uvascim.image.reco.data = permute(tmp, [1 2 4 3]);
            uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4]);
        else
            % 4D multi expt
            uvascim.uvascim.image.reco.data = permute(tmp, [1 2 5 3 4]);
            uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4 5]);
        end
    else
        % 4D nii file to 5D matlab matrix
        uvascim.uvascim.image.reco.data = reshape(tmp, [size(tmp,1), size(tmp,2),size(tmp,3),...
            size(uvascim.uvascim.image.reco.data,3), size(uvascim.uvascim.image.reco.data,5)]);
        uvascim.uvascim.image.reco.data= permute(uvascim.uvascim.image.reco.data, [1 2 4 3 5]);
        uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4 5]);
    end
    
    uvascim.uvascim.image.reco.data= flip(uvascim.uvascim.image.reco.data, 1);
    uvascim.uvascim.image.reco.data = flip(uvascim.uvascim.image.reco.data, 2);
    
    %
%     uvascim.uvascim.image.reco.data =reshape(tmp, [size(tmp, 1), size(tmp, 2), size(moved.uvascim.image.reco.data, 4), size(moved.uvascim.image.reco.data, 3), size(moved.uvascim.image.reco.data, 5)]);
%     
%     uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [1 2 4 3 5]);
%     uvascim.uvascim.image.reco.data = permute(uvascim.uvascim.image.reco.data, [2 1 3 4 5]);
%     
%     uvascim.uvascim.image.reco.data= flip(uvascim.uvascim.image.reco.data, 1);
%     uvascim.uvascim.image.reco.data = flip(uvascim.uvascim.image.reco.data, 2);
    %     end
    
    
else
    for i = 1:size(Scan_to_realign_nii_filename,1)
        V =spm_vol(Scan_to_realign_nii_filename{i,1});
        %     V =spm_vol(fullfile([ PATHSTR,filesep, [add_parameters{1,size(add_parameters,2)} NAME_scan_to_realign], '.nii']));
        tmp = spm_read_vols(V);
        uvascim(i) = moved;
        if numel(size(squeeze(uvascim(i).uvascim.image.reco.data))) < 5
            if  size(uvascim(i).uvascim.image.reco.data,3) >  size(uvascim(i).uvascim.image.reco.data,5)
                uvascim(i).uvascim.image.reco.data = permute(tmp, [1 2 4 3]);
                uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [2 1 3 4]);
            else
                uvascim(i).uvascim.image.reco.data = permute(tmp, [1 2 5 3 4]);
                
                uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [2 1 3 4 5]);
            end
        else
            % 4D nii file to 5D matlab matrix
            uvascim(i).uvascim.image.reco.data = reshape(tmp, [size(tmp,1), size(tmp,2),size(tmp,3),...
                size(uvascim(i).uvascim.image.reco.data,3), size(uvascim(i).uvascim.image.reco.data,5)]);
            uvascim(i).uvascim.image.reco.data= permute(uvascim(i).uvascim.image.reco.data, [1 2 4 3 5]);
            uvascim(i).uvascim.image.reco.data = permute(uvascim(i).uvascim.image.reco.data, [2 1 3 4 5]);
        end
        uvascim(i).uvascim.image.reco.data= flip(uvascim(i).uvascim.image.reco.data, 1);
        uvascim(i).uvascim.image.reco.data = flip(uvascim(i).uvascim.image.reco.data, 2);
    end
end
moved_maps = 1;

   


