function [files_in,files_out,opt] = Module_DSC(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
% define every option needed to run this module
%     %   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','AllName'};
    module_option(:,4)   = {'output_filename_ext_CBV','CBV'};
    module_option(:,5)   = {'output_filename_ext_CBF','CBF'};
    module_option(:,6)   = {'output_filename_ext_MTT','MTT'};
    module_option(:,7)   = {'output_filename_ext_TMAX','TMAX'};
    module_option(:,8)   = {'output_filename_ext_TTP','TTP'};
    module_option(:,9)   = {'output_filename_ext_T0','T0'};
    module_option(:,10)   = {'RemoveBackground','No'};
    module_option(:,11)   = {'Realign','No'};
    module_option(:,12)   = {'RefInput', 1};
    module_option(:,13)   = {'InputToReshape',1};
    module_option(:,14)   = {'Table_in', table()};
    module_option(:,15)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    opt.NameOutFiles = {opt.Module_settings.output_filename_ext_CBV, opt.Module_settings.output_filename_ext_CBF, opt.Module_settings.output_filename_ext_MTT, opt.Module_settings.output_filename_ext_TMAX, opt.Module_settings.output_filename_ext_TTP, opt.Module_settings.output_filename_ext_T0};
% 
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Select one Dynamic Susceptibility Contrast (DSC) scan as input','1Scan','','',{'SequenceName'},'Mandatory',''};
    user_parameter(:,2)   = {'Parameters','','','','','',''};
    user_parameter(:,3)   = {'   .Remove background','cell',{'No','Yes'},'RemoveBackground','','',{'Remove the background of the perf image before computation. Method by Jean-Albert Lotterie (CHU Toulouse):', newline, ...
        'The following algorithm works on the mean volume of the 4D-input. It works as follow: ', char(13), ...
        ' - Define the 8 voxels at the 8 corners of the volume as the zone 1 (the background).', ...
        ' - -Define the others voxels as the zone 2.', ...
        ' - Define m1 and m2 the means of those 2 zones.', ...
        ' - Compute m the arithmetic mean of m1 and m2.', ...
        ' - Iterative loop: ', ...
        '     - Define the voxels with a signal intensity lower than m as the new zone 1, with a mean m1.', ...
        '     - The other voxels define the zone 2, with the mean m2.', ...
        'The loop stops when floor(m), where m is still the arithmetic mean of m1 and m2, is stable.', ...
        'All the values of the zone 1 voxels are then set to NaN for each volume of the 4D-input.', ...
        }};
    user_parameter(:,4)   = {'   .Realign 4th dimension','cell',{'No','Yes'},'Realign','','',''};
    user_parameter(:,5)   = {'   .Select 1 ROI for the AIF (if empty, the AIF will be autoamtically defined)','1ROI','','',{'SequenceName'}, 'Optional',...
        {'By deflaut, the Arterial Input function (AIF) is automatically determined.'
        'User can manually select one Region-of-interest (ROI). This ROI will be used as AIF'}};
    user_parameter(:,6)   = {'   .Output filename extension CBV','char','CBV','output_filename_ext_CBV','','',''};
    user_parameter(:,7)   = {'   .Output filename extension CBF','char', 'CBF','output_filename_ext_CBF','','',''};
    user_parameter(:,8)   = {'   .Output filename extension MTT','char','MTT','output_filename_ext_MTT','','',''};
    user_parameter(:,9)   = {'   .Output filename extension TMAX','char','TMAX','output_filename_ext_TMAX','','',''};
    user_parameter(:,10)   = {'   .Output filename extension TTP','char','TTP','output_filename_ext_TTP','','',''};
    user_parameter(:,11)   = {'   .Output filename extension T0','char','T0','output_filename_ext_T0','','',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
%     
%%
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%

opt.NameOutFiles = {opt.output_filename_ext_CBV, opt.output_filename_ext_CBF, opt.output_filename_ext_MTT, opt.output_filename_ext_TMAX, opt.output_filename_ext_TTP, opt.output_filename_ext_T0};

if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in(opt.Table_in.Type == 'Scan',:);
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), opt.NameOutFiles{i}]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
end



%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Susceptibility:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end


%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end

[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nifti_header_data = spm_vol(files_in.In1{1});
N = read_volume(nifti_header_data, nifti_header_data, 0);
if size(N,4) == 1
   error([files_in.In1{1} ' is not a 4d image']);
    %return
end
% data have to be converted to double before computing the perfusion map
if isa(class(N), 'double') ==0
    origin_class =  class(N);
    N = double(N);
end

%% Realign if needed
if strcmp(opt.Realign, 'Yes')
    data_4d = squeeze(N);
    
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

    [~, file, ~] = fileparts(files_in.In1{1});
    TmpFile = [opt.folder_out, filesep, 'TMP_FILE_Module_DSC_', file,'.nii'];
    OutTmpFile = [opt.folder_out, filesep, 'rTMP_FILE_Module_DSC_', file, '.nii'];
    copyfile(files_in.In1{1},  TmpFile);
    
    Scan_to_realign_nii_filename = cell(length(keepscans), 1);
    for xx = 1:numel(keepscans)
        Scan_to_realign_nii_filename{xx} =  fullfile([TmpFile ',' num2str(keepscans(xx))]);
    end
     matlabbatch{1}.spm.spatial.realign.estwrite.data  = {Scan_to_realign_nii_filename};
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
    jobs = repmat(matlabbatch, 1, 1);
    inputs = cell(0, 1);
    for crun = 1:1
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
    clear matlabbatch 
    
    N2_header = spm_vol(OutTmpFile);
    N2 = read_volume(N2_header, N2_header, 0);
    
    N3 = N;
    N3(:,:,:,keepscans) = N2(:,:,:,keepscans);
    % the matrix of the DSC is updated using the resut of the realign
    N = N3;
    % replace NaN by 0, otherwise the option RemoveBackground crashes
    N(isnan(N)) = 0;
    % then every data created by SPM are deleted
    delete(TmpFile);
    delete(OutTmpFile);
    delete(strrep(TmpFile, '.nii', '.mat'))
    delete(strrep(OutTmpFile, ['rTMP_FILE_Module_DSC_', file, '.nii'], ['meanTMP_FILE_Module_DSC_', file, '.nii']))
    delete(strrep(OutTmpFile, ['rTMP_FILE_Module_DSC_', file, '.nii'], ['rp_TMP_FILE_Module_DSC_', file, '.txt']))
end


%% Remove Background

if strcmp(opt.RemoveBackground, 'Yes')
    MeanVol = mean(N,4);
    % m1 mean on the 8 voxels on the 8 corners of the volume. Those voxels
    % belongs to the background. It's the first zone.
    Mask1 = false(size(MeanVol));
    Mask1(1,1,1) = 1;
    Mask1(end,1,1) = 1;
    Mask1(1,end,1) = 1;
    Mask1(end,end,1) = 1;
    Mask1(1,1,end) = 1;
    Mask1(end,1,end) = 1;
    Mask1(1,end,end) = 1;
    Mask1(end,end,end) = 1;
    m1 = mean(MeanVol(Mask1(:)));
    %The second zone is the complementary of the first one. Its mean is m2.
    Mask2 = ~Mask1;
    m2 = mean(MeanVol(Mask2(:)));
    
%% Avec Seuil :
%     Seuil = opt.Seuil;
%     
%     ArithmMean = mean([m1,m2]);
%     OldArithmMean = ArithmMean - Seuil;
%     
%     while abs(ArithmMean - OldArithmMean) >= Seuil
%         Mask1 = MeanVol<ArithmMean;
%         m1 = mean(MeanVol(Mask1(:)));
%         Mask2 = ~Mask1;
%         m2 = mean(MeanVol(Mask2(:)));
%         OldArithmMean = ArithmMean;
%         ArithmMean = mean([m1,m2]);
%     end

%% Sans Seuil :
    
    ArithmMean = mean([m1,m2]);
    OldArithmMean = 0;
    
    while floor(ArithmMean) == OldArithmMean
        Mask1 = MeanVol<ArithmMean;
        m1 = mean(MeanVol(Mask1(:)));
        Mask2 = ~Mask1;
        m2 = mean(MeanVol(Mask2(:)));
        OldArithmMean = floor(ArithmMean);
        ArithmMean = mean([m1,m2]);
    end
    
    Mask = MeanVol >= ArithmMean;
    ExcluededVoxels = ~Mask;
    B = repmat(ExcluededVoxels, [1, 1, 1, size(N,4)]);
    N(B) = NaN;
    
    
end



%% defin the Arterial input function.
% could be manual using a ROI or automatique
if isfield(files_in, 'In2')
    nifti_header_ROI = spm_vol(files_in.In2{1});
    ROI = read_volume(nifti_header_ROI, nifti_header_data, 0);
    ROI = double(ROI);
    aif = N.*repmat(ROI, [1 1 1 size(N, 4)]);
    aif(aif == 0)=nan;
    aif = squeeze(nanmean(nanmean(nanmean(aif))))';
else
    meansignal = mean(squeeze(N),4);
    volume_mask = meansignal>max(N(:))*0.01;
    [aif,scores] = extraction_aif_volume(squeeze(N),volume_mask);
end

if exist('scores', 'var') &&  sum(cell2num(scores(:,5))) > 10
    error('No computation because the AIF is not good enough');
else
    [path, name, ext] = fileparts(files_in.In1{1});
    jsonfile = [path, '/', name, '.json'];
    J = ReadJson(jsonfile);
    % We are now using the information about the unit of the repetition
    % time and the echo time in order to automatically convert it to second
        [~,~,~,TMAX,TTP,T0, CBV,CBF,MTT,~,~,~,~] = deconvolution_perfusion_gui(aif,squeeze(N),...
            J.RepetitionTime.value(1) * convert_units(J.RepetitionTime.units, 'sec'),J.EchoTime.value* convert_units(J.RepetitionTime.units, 'sec'));
      % the conversions were hard coded before
%     [~,~,~,TMAX,TTP,T0, CBV,CBF,MTT,~,~,~,~] = deconvolution_perfusion_gui(aif,squeeze(N),J.RepetitionTime.value(1)*10^(-3),J.EchoTime.value*10^(-3));
    if strcmp(opt.RemoveBackground, 'Yes')
        CBV(ExcluededVoxels) = NaN;
        CBF(ExcluededVoxels) = NaN;
        MTT(ExcluededVoxels) = NaN;
        TMAX(ExcluededVoxels) = NaN;
        TTP(ExcluededVoxels) = NaN;
        T0(ExcluededVoxels) = NaN;
    end
    mapsVar = {CBV, CBF, MTT, TMAX, TTP, T0};

end


  
  
%% Reshape to the input scan size
%FilteredImages = reshape(FilteredImages, Size);


info = niftiinfo(files_in.In1{1});


for i=1:length(mapsVar)
  
    
    % reorient the scan
    mapsVar_reoriented{i} = write_volume(mapsVar{i}, nifti_header_data, 0);
    
    % update the header
     info2 = info;
    info2.Filename = files_out.In1{i};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(mapsVar_reoriented{i});
    info2.PixelDimensions = info.PixelDimensions(1:length(size(mapsVar{i})));
    info2.ImageSize = size(mapsVar_reoriented{i});
    info2.MultiplicativeScaling = 1;
    %info2.Description = [info.Description, 'Modified by Susceptibility Module'];
    
    niftiwrite(mapsVar_reoriented{i}, files_out.In1{i}, info2)
    
    % update and save the json 
    J2 = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

    [path, name, ~] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J2, jsonfile)
end
















