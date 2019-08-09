function out = cfg_run_parsebids(job)

% Make a directory and return its path in out.dir{1}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________

% Tanguy Duval
persistent BIDSfolder
alreadyparsed=0;
if ~isempty(BIDSfolder)
    alreadyparsed = strcmp({BIDSfolder.parent},job.parent{1});
    if any(alreadyparsed) && (now-BIDSfolder(find(alreadyparsed,1)).time)<3e-4
        BIDS = BIDSfolder(alreadyparsed).BIDS;
    end
end
if ~exist('bids_parser','file')
    addpath(fullfile(spm('dir'),'external','bids_tools_matlab'))
    addpath(genpath(fullfile(spm('dir'),'external','bids_tools_matlab','External')))
end

% CALL BIDS_PARSER
if ~exist('BIDS','var')
    BIDS = bids_parser(job.parent{1});
    if ~any(alreadyparsed)
        alreadyparsed = length(BIDSfolder)+1;
    else
        alreadyparsed = find(alreadyparsed,1);
    end
    BIDSfolder(alreadyparsed).parent=job.parent{1};
    BIDSfolder(alreadyparsed).BIDS=BIDS;
    BIDSfolder(alreadyparsed).time=now;
end
% TEST IF SESSION SELECTED
if ~isfield(job,'bids_ses_type')
    out = BIDS;
    return;
end
% TEST IF BIDS
if isempty(BIDS.subjects)
    out = [];
    return
end
% Number or name
if isfield(job.bids_ses_type,'bids_sesnum')
    if isa(job.bids_ses_type.bids_sesnum,'cfg_dep'), job.bids_ses_type.bids_sesnum=1; end
    job.bids_ses = job.bids_ses_type.bids_sesnum;
elseif isfield(job.bids_ses_type,'bids_sesrel')
    if isa(job.bids_ses_type.bids_sesrel,'cfg_dep'), job.bids_ses_type.bids_sesrel=['sub-' BIDS.subjects(1).name '/ses-' BIDS.subjects(1).session]; end
    bids_sesrel = strsplit(job.bids_ses_type.bids_sesrel,{'/' '\'});
    if length(bids_sesrel)==1
        job.bids_ses = strcmp({BIDS.subjects.name},strrep(bids_sesrel{1},'sub-',''));
        bids_sesrel{2} = 'undefined';
    else
        job.bids_ses = strcmp({BIDS.subjects.name},strrep(bids_sesrel{1},'sub-','')) & strcmp({BIDS.subjects.session},strrep(bids_sesrel{2},'ses-',''));
    end
    if any(job.bids_ses)
        job.bids_ses = find(job.bids_ses,1);
    else
        warning(['subject ' bids_sesrel{1} ' with session ' bids_sesrel{2} ' not found']);
        out = [];
        return
    end
end

% SELECT CURRENT SUBJECT
session = min(length(BIDS.subjects),job.bids_ses);
SCAN = BIDS.subjects(session);
disp(['Subject = ' SCAN.name])
disp(['Session = ' SCAN.session])

% Get Path
out.bidsdir         = {BIDS.path};
out.sub             = SCAN.name;
out.ses             = SCAN.session;
if strcmp(job.derivativesname,'<UNDEFINED>'), job.derivativesname = ''; end
out.bidsderivatives = fullfile(out.bidsdir,'derivatives',job.derivativesname,['sub-' SCAN.name],['ses-' SCAN.session]);
out.bidssession = fullfile(out.bidsdir,['sub-' SCAN.name],['ses-' SCAN.session]);
if isfield(SCAN,'tsv') && ~isempty(SCAN.tsv)
    ftsv            = fieldnames(SCAN.tsv);
    for fff = 1:length(ftsv)
        out.(['tsv_' ftsv{fff}]) = SCAN.tsv.(ftsv{fff});
    end
end
if strcmp(SCAN.session,'none')
    out.relpath         = fullfile(['sub-' SCAN.name]);
else
    out.relpath         = fullfile(['sub-' SCAN.name],['ses-' SCAN.session]);
end

if ~isempty(job.derivativesname) && ~exist(out.bidsderivatives{1},'dir')
    mkdir(out.bidsderivatives{1})
end

outses = parseScan(SCAN);
switch char(fieldnames(job.bids_ref_type))
    case 'session_Nminus1'
        if session ==1
            idx_ref = session;
        elseif ~strcmp(BIDS.subjects(session-1).name,BIDS.subjects(session).name)
            idx_ref = session;
        else
            idx_ref = session-1;
        end
    case 'session_N1'
        idx_ref = find(strcmp({BIDS.subjects.name},BIDS.subjects(session).name),1);
    case 'sub_1'
        idx_ref = 1;
    case 'noref'
        idx_ref = 0;
end

if idx_ref
    outref = parseScan(BIDS.subjects(idx_ref));
    disp(['Reference Subject = ' BIDS.subjects(idx_ref).name])
    disp(['Reference Session = ' BIDS.subjects(idx_ref).session])
end

% add files to out
ff = fieldnames(outses);
for iff = 1:length(ff)
    out.(ff{iff}) = outses.(ff{iff});
end

if idx_ref
    ff = fieldnames(outref);
    for iff = 1:length(ff)
        out.(['Reference_' ff{iff}]) = outref.([ff{iff}]);
    end
end




function out = parseScan(SCAN)

% ADD OUTPUT DIRECTORIES

% LOOP OVER MODALITIES AND EXTRACT PATIENT DATA FILENAME AND METADATA
% list = {'anat','T1w',...
%         'anat','T2w',...
%         'anat', 'FLAIR',...
%         'func','bold',...
%         'dwi','dwi'...
%         };
    
mods = setdiff(fieldnames(SCAN),{'name','path','session','tsv'});
list = {};
for imods = 1:length(mods)
    for ifile = 1:length(SCAN.(mods{imods}))
        list{end+1} = mods{imods};
        if isfield(SCAN.(mods{imods}),'modality') && ~isempty(SCAN.(mods{imods})(ifile).modality)
            list{end+1} = SCAN.(mods{imods})(ifile).modality;
        else
            SCAN.(mods{imods})(ifile).modality = 'unknown';
            list{end+1} = regexprep(SCAN.(mods{imods})(ifile).filename,'\.nii(\.gz)?','');
            SCAN.(mods{imods})(ifile).modality = [];
        end
    end
end

for ii=1:2:length(list)
        MODALITY = SCAN.(list{ii})(strcmp({SCAN.(list{ii}).modality},list{ii+1}) | strcmp(regexprep({SCAN.(list{ii}).filename},'\.nii(\.gz)?',''),list{ii+1}));
        if ~isempty(MODALITY)
            % Add nifti
            nii_fname = fullfile(SCAN.path,strrep(list{ii},'image',''),MODALITY(1).filename);
            tag = genvarname([list{ii} '_' list{ii+1}]);
            out.(tag) = {nii_fname};
            
            % Special treatment for dmri
            if strcmp(list{ii},'dwi') && ~isempty(strfind(nii_fname,'_dwi.nii'))
                out.dwi_bvec = {strrep(strrep(nii_fname,'.gz',''),'.nii','.bvec')};
                out.dwi_bval = {strrep(strrep(nii_fname,'.gz',''),'.nii','.bval')};
            end
            
            % Add metadata
            if isfield(MODALITY(1),'meta') && ~isempty(MODALITY(1).meta)
                tag = [tag '_meta'];
                out.(tag) = MODALITY(1).meta;
            end
        end
end
