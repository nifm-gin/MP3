function dep = cfg_vout_parsebids(job)

% Define virtual outputs for cfg_run_parsebids. The directory name can either
% be assigned to a cfg_files directory input or to a evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________

% Tanguy Duval
h_cfg_ui = findobj('Type','figure','Tag','cfg_ui');
t = findobj(h_cfg_ui,'Type','uitoolbar');
if isempty(findobj(t,'Tag','TbFileRunBIDS'))
    try
        icodoublearrow = load('icodoublearrow.mat');
    uipushtool(t,'TooltipString','Run Batch for all subjects',...
        'Tag','TbFileRunBIDS',...
        'ClickedCallback','cfg_ui_loop(guidata(gcbo))',...
        'CData',icodoublearrow.icon);
    end
end

dep               = cfg_dep;
dep(1)            = cfg_dep;
dep(1).sname      = sprintf('BIDS directory path');
dep(1).src_output = substruct('.','bidsdir');
dep(1).tgt_spec   = cfg_findspec({{'filter','dir', 'strtype','e'}});
%
% BIDS = bids_parser(uigetdir([],'BIDS folder on which dependencies will be built'));
% types = setdiff(fieldnames(BIDS.subjects(1)),{'name','path','session'});

dep(2)            = cfg_dep;
dep(2).sname      = fullfile('BIDS/sub-name/ses-session/');
dep(2).src_output = substruct('.','bidssession');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});

dep(3)            = cfg_dep;
dep(3).sname      = sprintf('Subject Name');
dep(3).src_output = substruct('.','sub');
dep(3).tgt_spec   = cfg_findspec({{'class','cfg_entry', 'strtype','s'}});

dep(4)            = cfg_dep;
dep(4).sname      = sprintf('Session Name');
dep(4).src_output = substruct('.','ses');
dep(4).tgt_spec   = cfg_findspec({{'class','cfg_entry', 'strtype','s'}});

dep(5)            = cfg_dep;
dep(5).sname      = sprintf('sub-name/ses-session');
dep(5).src_output = substruct('.','relpath');
dep(5).tgt_spec   = cfg_findspec({{'class','cfg_entry', 'strtype','s'}});

if strcmp(job.parent,'<UNDEFINED>')
    job.parent = {[mfilename('fullpath') '_template']};  
    job.bids_ses = 1;
    job.name   = {'matlabbatch_test'};
end
if strcmp(job.bids_ses_type,'<UNDEFINED>')
    job.bids_ses = 1;
end
out = cfg_run_parsebids(job);

if ~isempty(job.derivativesname)
    dep(6)            = cfg_dep;
    dep(6).sname      = [fullfile('BIDS/derivatives', job.derivativesname, 'sub-name','ses-session') filesep];
    dep(6).src_output = substruct('.','bidsderivatives');
    dep(6).tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});
end

if isempty(out), return; end
for ff = setdiff(fieldnames(out)',{'bidsdir','bidsderivatives','sub','ses','relpath'})
    dep(end+1)            = cfg_dep;
    dep(end).sname      = strrep(strrep(ff{1},'_meta',' metadata'),'_',': ');
    dep(end).src_output = substruct('.',ff{1});
    if strcmp(ff{1}(end-4:end),'_meta') || strcmp(ff{1}(1:min(end,4)),'tsv_')
        dep(end).tgt_spec   = cfg_findspec({{'class','cfg_entry', 'strtype','e'}});
    else
        dep(end).tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});
    end
end