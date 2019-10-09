function cfg_BIDS = cfg_cfg_BIDS

% ---------------------------------------------------------------------
% parent Parent Directory
% ---------------------------------------------------------------------
parent         = cfg_files;
parent.tag     = 'parent';
parent.name    = 'Parent Directory';
parent.help    = {'Directory where the new directory will be created.'};
parent.filter = {'dir'};
parent.ufilter = '.*';
parent.num     = [1 1];
% ---------------------------------------------------------------------
% scan scan number
% ---------------------------------------------------------------------
bidssessionnum         = cfg_entry;
bidssessionnum.tag     = 'bids_sesnum';
bidssessionnum.name    = 'Number';
bidssessionnum.val     = {1};
bidssessionnum.help    = {'Scanning session to process. 1 = first session of first subject'};
bidssessionnum.strtype = 'n';
bidssessionnum.num     = [1 1];

% ---------------------------------------------------------------------
% scan scan number
% ---------------------------------------------------------------------
bidssessionrel         = cfg_entry;
bidssessionrel.tag     = 'bids_sesrel';
bidssessionrel.name    = 'sub-Name/ses-Session';
bidssessionrel.val     = {'sub-Name/ses-Session'};
bidssessionrel.help    = {'Scanning session to process. Enter a string of this type: sub-Name/ses-Session'};
bidssessionrel.strtype = 's';
bidssessionrel.num     = [1 inf];

% ---------------------------------------------------------------------
% scan scan number
% ---------------------------------------------------------------------
bidssession         = cfg_choice;
bidssession.tag     = 'bids_ses_type';
bidssession.name    = 'BIDS session';
bidssession.values  = {bidssessionnum bidssessionrel};
bidssession.val     = {bidssessionnum};
bidssession.help    = {'Scanning session to process. 1 = first session of first subject',...
                       'Preview ''Parent directory'' to list all sessions'};

% ---------------------------------------------------------------------
% name New Directory Name
% ---------------------------------------------------------------------
derivativesname         = cfg_entry;
derivativesname.tag     = 'derivativesname';
derivativesname.name    = 'New Directory Name (BIDS/derivatives/NewDir)';
derivativesname.help    = {'Name for the new directory.',...
                'Directory is created in the ''derivatives'' folder of your Bids directory',...
                'A virtual output will be accessible in your dependency tree',...
                'Leave empty if you do not want an output directory'};
derivativesname.strtype = 's';
derivativesname.val     = {''};
derivativesname.num     = [0  Inf];
% ---------------------------------------------------------------------
% session_N1 First Session
% ---------------------------------------------------------------------
reffirstsub        = cfg_entry;
reffirstsub.tag    = 'sub_1';
reffirstsub.name   = 'First Subject (first session)';
reffirstsub.val    = {true};
% ---------------------------------------------------------------------
% session_N1 First Session
% ---------------------------------------------------------------------
reffirstsession        = cfg_entry;
reffirstsession.tag    = 'session_N1';
reffirstsession.name   = 'First Session (in each subject)';
reffirstsession.val    = {true};
% ---------------------------------------------------------------------
% session_Nminus1 Session N-1
% ---------------------------------------------------------------------
refprevioussession        = cfg_entry;
refprevioussession.tag    = 'session_Nminus1';
refprevioussession.name   = 'Session N-1';
refprevioussession.val    = {true};
% ---------------------------------------------------------------------
% bidsreference Reference scan
% ---------------------------------------------------------------------
noref         = cfg_const;
noref.tag     = 'noref';
noref.name    = 'No reference';
noref.val     = {true};
% ---------------------------------------------------------------------
% bidsreference Reference scan
% ---------------------------------------------------------------------
bidsreference         = cfg_choice;
bidsreference.tag     = 'bids_ref_type';
bidsreference.name    = 'Reference type';
bidsreference.values  = {noref refprevioussession reffirstsession reffirstsub};
bidsreference.val     = {noref};
bidsreference.help    = {'Reference scan in your dependency tree. 1 = first session of first subject',...
                       'Preview ''Parent directory'' to list all sessions'};
% ---------------------------------------------------------------------
% cfg_parsebids Parse BIDS Directory
% ---------------------------------------------------------------------
parentbids            = parent;
parentbids.name       = 'BIDS directory';
parentbids.preview    = @(parent) cfg_preview_parsebids(parent);
cfg_parsebids         = cfg_exbranch;
cfg_parsebids.tag     = 'cfg_parsebids';
cfg_parsebids.name    = 'Parse BIDS Directory';
cfg_parsebids.val     = {parentbids, bidssession derivativesname bidsreference};
cfg_parsebids.help    = {'Parse a BIDS directory.'};
cfg_parsebids.prog    = @cfg_run_parsebids;
cfg_parsebids.vout    = @cfg_vout_parsebids;
% ---------------------------------------------------------------------
% dicom Dicom Directory
% ---------------------------------------------------------------------
dicom         = cfg_files;
dicom.tag     = 'dicom';
dicom.name    = 'dicom Directory';
dicom.help    = {'Directory that contains dicoms.'};
dicom.filter = {'dir'};
dicom.ufilter = '.*';
dicom.val    = {''};
dicom.num     = [1 1];
% ---------------------------------------------------------------------
% outputdir output Directory
% ---------------------------------------------------------------------
outputdir         = cfg_files;
outputdir.tag     = 'outputdir';
outputdir.name    = 'output Directory';
outputdir.help    = {'Directory where the new BIDS directory will be created.'};
outputdir.filter = {'dir'};
outputdir.val = {''};
outputdir.ufilter = '.*';
outputdir.num     = [1 1];
% ---------------------------------------------------------------------
% cfg_dicm2bids DICOM to BIDS
% ---------------------------------------------------------------------
cfg_dicm2bids         = cfg_exbranch;
cfg_dicm2bids.tag     = 'cfg_dicm2bids';
cfg_dicm2bids.name    = 'DICOM to BIDS';
cfg_dicm2bids.val     = {dicom, outputdir};
cfg_dicm2bids.help    = {'Convert a DICOM directory to BIDS.'};
cfg_dicm2bids.prog    = @(job) dicm2bids(getcell1(job.dicom),getcell1(job.outputdir));
% ---------------------------------------------------------------------
% participant Participant
% ---------------------------------------------------------------------
bidsapptypeparticipant         = cfg_const;
bidsapptypeparticipant.tag     = 'participant';
bidsapptypeparticipant.name    = 'Participant';
bidsapptypeparticipant.val     = {true};
% ---------------------------------------------------------------------
% participant Participant
% ---------------------------------------------------------------------
bidsapptypegroup         = cfg_const;
bidsapptypegroup.tag     = 'group';
bidsapptypegroup.name    = 'Group';
bidsapptypegroup.val     = {true};
% ---------------------------------------------------------------------
% dockerimg docker image
% ---------------------------------------------------------------------
bidsapptype         = cfg_choice;
bidsapptype.tag     = 'bidsapptype';
bidsapptype.name    = 'Type of analysis';
bidsapptype.val     = {bidsapptypeparticipant};
bidsapptype.values  = {bidsapptypeparticipant, bidsapptypegroup};
% ---------------------------------------------------------------------
% dockerimg docker image
% ---------------------------------------------------------------------
dockerimg         = cfg_entry;
dockerimg.tag     = 'dockerimg';
dockerimg.name    = 'docker image';
dockerimg.val     = {'bids/example'};
dockerimg.help    = {'List of available apps on http://bids-apps.neuroimaging.io/apps/'};
dockerimg.strtype = 's';
dockerimg.num     = [1  Inf];
dockerimg.preview = @(val) msgboxdocker(val);
% ---------------------------------------------------------------------
% bidsappOpts docker image
% ---------------------------------------------------------------------
bidsappOpts         = cfg_entry;
bidsappOpts.tag     = 'bidsappOpts';
bidsappOpts.name    = 'BIDS App optional arguments';
bidsappOpts.val     = {''};
bidsappOpts.strtype = 's';
bidsappOpts.num     = [1  Inf];
% ---------------------------------------------------------------------
% cfg_parsebids Parse BIDS Directory
% ---------------------------------------------------------------------
cfg_bidsapp         = cfg_exbranch;
cfg_bidsapp.tag     = 'cfg_bidsapp';
cfg_bidsapp.name    = 'run BIDS app';
cfg_bidsapp.val     = {parentbids, outputdir, bidsapptype, dockerimg bidsappOpts};
cfg_bidsapp.help    = {sprintf('Convert a DICOM directory to BIDS.\nUse Boutiques for a more flexible module')};
cfg_bidsapp.prog    = @(job) system_docker(lower(strrep(job.dockerimg,'BIDS-Apps','bids')),['defaultEntrypoint %i1 %i2 ' getcell1(fieldnames(job.bidsapptype)) ' ' getcell1(job.bidsappOpts)],...
                               job.parent{1},... %i1
                               getcell1(job.outputdir,getcell1(getcell1(regexp(lower(strrep(job.dockerimg,'BIDS-Apps','bids')),'bids/(\w+):?','tokens'))))); %i2
% ---------------------------------------------------------------------
% cfg_BIDS BIDS
% ---------------------------------------------------------------------
cfg_BIDS         = cfg_choice;
cfg_BIDS.tag     = 'cfg_BIDS';
cfg_BIDS.name    = 'BIDS';
cfg_BIDS.help    = {'This toolbox contains BIDS functions. http://bids.neuroimaging.io/'};
cfg_BIDS.values  = {cfg_dicm2bids cfg_parsebids cfg_bidsapp};

function msgboxdocker(val)
val = lower(strrep(val,'BIDS-Apps','bids'));
[~,stdout] = system_docker(val,'defaultEntrypoint -h');
h = msgbox(stdout,['help on ' val]);
set(findall(h,'Type','Text'),'FontName','FixedWidth');
extent =  get(findall(h,'Type','Text'),'Extent');
Pos = get(h,'Position');
Pos(4) = max(Pos(4),extent(4) + extent(2));
Pos(3) = max(Pos(3),extent(3) + extent(1));
set(h,'Position',Pos)




function dicm2bids(src,dest)
if isempty(src) || isempty(dest)
    setpref('dicm2nii_gui_para', 'rstFmt',3) % set bids as defaults
    setpref('dicm2nii_gui_para', 'save_json',true) % export jsons
    dicm2nii
else
    dicm2nii(src,dest,'bids')
end

function out = getcell1(in,default)
if isempty(in)
    if nargin>1
        out = default;
    else
    out = [];
    end
elseif iscell(in)
    out = in{1};
else
    out = in;
end
    
    
