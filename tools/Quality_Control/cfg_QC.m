function QC = cfg_QC(varargin)

% ---------------------------------------------------------------------
% nii nii
% ---------------------------------------------------------------------
nii         = cfg_files;
nii.tag     = 'nifti';
nii.name    = 'Image file';
nii.ufilter = '\.nii(\.gz)?';
nii.help = {'Select one or multiple nifti resliced files images (same space, same matrix size)'};
nii.num     = [0 Inf];
% ---------------------------------------------------------------------
% nii nii
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Mask file';
mask.val     = {{''}};
mask.num     = [0 1];
% ---------------------------------------------------------------------
% imtool3D imtool3D
% ---------------------------------------------------------------------
imtool3D         = cfg_exbranch;
imtool3D.tag     = 'imtool3D';
imtool3D.name    = 'imtool3D: 3D Nifti viewer with ROI tool';
imtool3D.val     = {nii mask};
imtool3D.prog = @(job) imtool3D_nii_3planes(job.nifti,job.mask{1});
imtool3D.help = {'NIFTI viewer with ROI tools and statistics'};
% ---------------------------------------------------------------------
% dir html report directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'html report directory';
dir.filter = {'dir'};
dir.num     = [1 1];
dir.help    = {'If directory already exists. qc will append a new line to the report table.'};
% ---------------------------------------------------------------------
% subject Subject Name
% ---------------------------------------------------------------------
inputdate         = cfg_entry;
inputdate.tag     = 'inputdate';
inputdate.name    = 'Date';
inputdate.val     = {datestr(date,'yyyy-mm-dd')};
inputdate.strtype = 's';
inputdate.num     = [0 inf];
inputdate.help    = {'Used in QC table'};
% ---------------------------------------------------------------------
% subject Subject Name
% ---------------------------------------------------------------------
subject         = cfg_entry;
subject.tag     = 'subject';
subject.name    = 'Subject Name';
subject.val     = {'sub01'};
subject.strtype = 's';
subject.num     = [0 inf];
subject.help    = {'Used in QC table'};
% ---------------------------------------------------------------------
% nii nii
% ---------------------------------------------------------------------
contrast         = cfg_entry;
contrast.tag     = 'contrast';
contrast.name    = 'contrast type';
contrast.val     = {'T1'};
contrast.strtype = 's';
contrast.num     = [0 inf];
contrast.help    = {'Used in QC table'};
% ---------------------------------------------------------------------
% nii nii
% ---------------------------------------------------------------------
command         = cfg_entry;
command.tag     = 'command';
command.name    = 'Command';
command.val     = {'denoising'};
command.strtype = 's';
command.num     = [0 inf];
command.help    = {'Used in QC table'};

% ---------------------------------------------------------------------
% nii nii
% ---------------------------------------------------------------------
slices         = cfg_entry;
slices.tag     = 'slices';
slices.name    = 'Slices';
slices.val     = {0};
slices.strtype = 'w';
slices.num     = [0 inf];
slices.help    = {'slices to select. 0: auto mode (6 regularly spaced slices)'};

% ---------------------------------------------------------------------
% QC QC
% ---------------------------------------------------------------------
htmlqc         = cfg_exbranch;
htmlqc.tag     = 'html';
htmlqc.name    = 'HTML report (SCT team)';
htmlqc.val     = {nii mask dir inputdate subject contrast, command, slices};
htmlqc.prog = @(job) qc_write(job.dir,job.inputdate,job.subject,job.nifti,job.mask,job.contrast,job.command, job.slices);
htmlqc.help = {'interactive html report developed by the SCT team'};

% ---------------------------------------------------------------------
% QC QC
% ---------------------------------------------------------------------
QC         = cfg_choice;
QC.tag     = 'QC';
QC.name    = 'Quality Control';
QC.values  = {htmlqc imtool3D};