function varargout = cfg_run_call_system(cmd, varargin)
% A generic interface to call any system command through the batch system
% and make its output arguments available as dependencies.
% varargout = cfg_run_call_matlab(cmd, varargin)
% where cmd is one of
% 'run'      - out = cfg_run_call_system('run', job)
%              Run the command, and return the specified output arguments
% 'vout'     - dep = cfg_run_call_system('vout', job)
%              Return dependencies as specified via the output cfg_repeat.
% 'check'    - str = cfg_run_call_system('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = cfg_run_call_system('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              cfg_run_call_system('defaults', key, newval)
%              Set the specified field in the internal def struct to a new
%              value.
% Application specific code needs to be inserted at the following places:
% 'run'      - main switch statement: code to compute the results, based on
%              a filled job
% 'vout'     - main switch statement: code to compute cfg_dep array, based
%              on a job structure that has all leafs, but not necessarily
%              any values filled in
% 'check'    - create and populate switch subcmd switchyard
% 'defaults' - modify initialisation of defaults in subfunction local_defs
% Callbacks can be constructed using anonymous function handles like this:
% 'run'      - @(job)cfg_run_call_system('run', job)
% 'vout'     - @(job)cfg_run_call_system('vout', job)
% 'check'    - @(job)cfg_run_call_system('check', 'subcmd', job)
% 'defaults' - @(val)cfg_run_call_system('defaults', 'defstr', val{:})
%              Note the list expansion val{:} - this is used to emulate a
%              varargin call in this function handle.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2018 Toulouse Neuroimaging Center

% Tanguy Duval

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1},true);
            % convert job structure to simple in and out struct
            in     = cell(size(job.inputs));

            % Test if docker is correctly installed
            if isfield(job.usedocker,'dockerimg')
                [status, res] = system('docker');
                if status
                    if ispc
                        warndlg(sprintf('%s\n ... Docker is not up and running...\nPlease install and run docker.',res))
                    else
                        warndlg(sprintf('%s\n ... Docker is not up and running...\nPlease install and run docker, then launch matlab from a terminal by using the following command:\n "%s"',res,fullfile(matlabroot,'bin','matlab')))
                    end
                end
            end

            
            for k = 1:numel(in)
                intmp = job.inputs{k}.(char(fieldnames(job.inputs{k})));
                intmp = intmp.(char(setdiff(fieldnames(intmp),'help')));
                
                switch char(fieldnames(job.inputs{k}))
                    case 'evaluatedbranch'
                        in{k} = cellfun(@(X) char(X),intmp,'uni',0);
                    case 'anyfilebranch'
                        in{k} = intmp;
                    case 'stringbranch'
                        in{k} = {char(intmp)};
                    otherwise
                        if iscell(intmp)
                            in{k} = cellfun(@(X) char(X),intmp,'uni',0);
                        else
                            in{k} = {char(intmp)};
                        end
                end
            end
            
            % Check if all outputs already exists
            out = job.outputs;
            if ~isempty(out)
                alreadyexist = true;
                for io = 1:length(job.outputs)
                    % replace token in outputs
                    tokens = regexp(job.outputs{io}.outputs.string,'(%i\d)','match');
                    for itok = 1:length(tokens)
                        job.outputs{io}.outputs.string = strrep(job.outputs{io}.outputs.string,tokens{itok},in{str2double(tokens{itok}(3:end))}{1});
                    end
                    tokens = regexp(job.outputs{io}.outputs.directory{1},'(%i\d)','match');
                    for itok = 1:length(tokens)
                        job.outputs{io}.outputs.directory = strrep(job.outputs{io}.outputs.directory,tokens{itok},in{str2double(tokens{itok}(3:end))}{1});
                    end
                    
                    out{io} = {fullfile(job.outputs{io}.outputs.directory{1},job.outputs{io}.outputs.string)};
                    alreadyexist = alreadyexist & exist(out{io}{1},'file');
                end
                if alreadyexist
                    disp(['<strong>output file already exists, assuming that the processing was already done... skipping</strong>'])
                    disp(['Delete output file to restart this job = ' out{1}{1}])
                end
            else
                alreadyexist = false;
            end
            
            % DOCKERIFY
            mountdir = '';
            for k = 1:numel(in)
                type = fieldnames(job.inputs{k});
                if strcmp(type{1},'directorybranch')
                    mountdir = [mountdir '-v "' in{k}{1} ':/i' num2str(k) '" '];
                    dockerinfname{k} = ['/i' num2str(k) '/'];
                else
                    mountdir = [mountdir '-v "' fileparts(in{k}{1}) ':/i' num2str(k) '" '];
                    dockerinfname{k} = strrep(in{k}{1},[fileparts(in{k}{1}) filesep],['/i' num2str(k) '/']);
                    
                    % special case if more than one file choosen in a field
                    for kl = 2:length(in{k})
                        mountdir = [mountdir '-v "' fileparts(in{k}{kl}) ':/i' num2str(k) '" '];
                        dockerinfname{k} = [dockerinfname{k} ' ' strrep(in{k}{kl},[fileparts(in{k}{kl}) filesep],['/i' num2str(k) '/'])];
                    end
                    
                end
            end
            for k = 1:numel(job.outputs)
                mountdir = [mountdir '-v "' fullfile(job.outputs{k}.outputs.directory{1},fileparts(job.outputs{k}.outputs.string)) ':/o' num2str(k) '" '];
                ofullfname = fullfile(job.outputs{k}.outputs.directory{1},job.outputs{k}.outputs.string);
                dockeroutfname{k} = strrep(ofullfname,[fileparts(ofullfname) filesep],['/o' num2str(k) '/']);
            end
            
            
            if ~alreadyexist
                % Replace token i%d and o%d by filenames
                cmd = job.cmd;
                for ii=1:length(in)
                    if isfield(job.usedocker,'dockerimg')
                        cmd = strrep(cmd,sprintf('%%i%d',ii),dockerinfname{ii});
                        cmd = strrep(cmd,[' ' sprintf('i%d',ii)],[' ' dockerinfname{ii} ' ']);
                    else
                        cmd = strrep(cmd,sprintf('%%i%d',ii),in{ii}{1});
                        cmd = strrep(cmd,[' ' sprintf('i%d',ii)],[' ' in{ii}{1} ' ']);
                    end
                end
                for ii=1:length(job.outputs)
                    if isfield(job.usedocker,'dockerimg')
                        cmd = strrep(cmd,sprintf('%%o%d',ii),dockeroutfname{ii});
                        cmd = strrep(cmd,[' ' sprintf('o%d',ii)],[' ' dockeroutfname{ii} ' ']);
                    else
                        cmd = strrep(cmd,sprintf('%%i%d',ii),out{ii}{1});
                        cmd = strrep(cmd,[' ' sprintf('o%d',ii)],[' ' out{ii}{1} ' ']);
                    end
                end
                
                % RUN SYSTEM COMMAND
                if ~isfield(job.usedocker,'dockerimg')
                    disp(['Running terminal command: ' cmd])
                    [status, stdout]=system(cmd,'-echo');
                    if status, error(sprintf('%s\n\nLocal command:\n%s\n',stdout,cmd)); end
                else % docker
                    cmdcell = strsplit(cmd);
                    if strcmpi(cmdcell{1},'defaultEntrypoint')
                        cmddocker = ['docker run ' mountdir job.usedocker.dockerimg ' ' strjoin(cmdcell(2:end))];
                    else
                        cmddocker = ['docker run --entrypoint ' cmdcell{1} ' ' mountdir job.usedocker.dockerimg ' ' strjoin(cmdcell(2:end))];
                    end
                    disp(['Running terminal command: ' cmddocker])
                    [status, stdout]=system(cmddocker,'-echo');
                    if status, error(sprintf('%s\nLocal command: \n%s\n\nCorresponding Docker command:\n%s\n',stdout,cmd,cmddocker)); end
                end
                
                
                
                % gunzip output for FSL if .nii was used
                for io=1:length(out)
                    if ~exist(out{io}{1},'file') && exist([out{io}{1} '.gz'],'file')
                        gunzip([out{io}{1} '.gz'])
                        delete([out{io}{1} '.gz'])
                    end
                end
                for io=1:length(out)
                    if ~exist(out{io}{1},'file') && strcmp(out{io}{1}(end-2:end),'.gz') && exist(out{io}{1}(1:end-3),'file')
                        gzip(out{io}{1}(1:end-3));
                        delete(out{io}{1}(1:end-3));
                    end
                end
                
            end
            
            if nargout > 0
                outjob.outputs = cell(size(job.outputs));
                outjob.outputs  = out;
                varargout{1} = outjob;
            end
        case 'vout'
            h_cfg_ui = findobj('Type','figure','Tag','cfg_ui');
            t = findobj(h_cfg_ui,'Type','uitoolbar');
            if isempty(findobj(t,'Tag','TbFileSaveBoutiques'))
                try
                    icodoublearrow = load('Boutiqueslogo.mat');
                    uipushtool(t,'TooltipString','Save Boutiques Descriptor',...
                        'Tag','TbFileSaveBoutiques',...
                        'ClickedCallback','cfg_save_call_system(guidata(gcbo))',...
                        'CData',icodoublearrow.icon);
                end
            end

            job = local_getjob(varargin{1},true);
            % initialise empty cfg_dep array
            dep = cfg_dep;
            dep = dep(false);
            % determine outputs, return cfg_dep array in variable dep
            if isfield(job,'outputs')
                for k = 1:numel(job.outputs)
                    dep(k)            = cfg_dep;
                    cmdstring = strsplit(job.cmd);
                    desc = strsplit(job.outputs{k}.outputs.help,':'); desc = char(desc{1});
                    dep(k).sname      = sprintf('%s: output %d - %s - %s',cmdstring{1}, k, desc, char(job.outputs{k}.outputs.string));
                    dep(k).src_output = substruct('.','outputs','{}',{k});
                    dep(k).tgt_spec   = cfg_findspec({{'filter', char(fieldnames(job.outputs{k}))}});
                end
            end
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    % implement checks, return status string in variable str
                    otherwise
                        cfg_message('unknown:check', ...
                            'Unknown check subcmd ''%s''.', subcmd);
                end
                varargout{1} = str;
            else
                cfg_message('ischar:check', 'Subcmd must be a string.');
            end
        case 'save'
            job = local_getjob(varargin{1},false);
            parentfolder = textread('boutiques_path.txt','%s');
            parentfolder = GetFullPath_v2(parentfolder,'auto',fileparts(which('boutiques_path.txt')));
            cfg2boutique(job,fileparts(parentfolder{1}))
            %             cmdcell = strsplit(job.cmd);
            %
            %             [directory, tree] = generatetree(['System.' cmdcell{1}],cfg_cfg_call_system);
            %             if isempty(tree), return; end
            %
            %             % generate cfg_CML_def
            %             jobstr = generatecfgdef(tree,job);
            %             cfg_def_fname = fullfile(directory,['cfg_' tree{1} '_def.m']);
            %             fid = fopen(cfg_def_fname,'wt');
            %             fprintf(fid, 'function %s = cfg_%s_def\n',lower(tree{1}),tree{1});
            %             fprintf(fid, '%s\n',jobstr{:});
            %             fclose(fid);
            %             disp(['files added in ' directory])
            
        otherwise
            cfg_message('unknown:cmd', 'Unknown command ''%s''.', cmd);
    end
else
    cfg_message('ischar:cmd', 'Cmd must be a string.');
end

function job = local_getjob(job,rename)
if ~isstruct(job)
    cfg_message('isstruct:job', 'Job must be a struct.');
end
if nargin>1 && rename
    if isfield(job,'inputs_'), job.inputs=job.inputs_; job = rmfield(job,'inputs_');   end
    if isfield(job,'outputs_'), job.outputs=job.outputs_; job = rmfield(job,'outputs_'); end
    if isfield(job,'usedocker_'), job.usedocker=job.usedocker_; job = rmfield(job,'usedocker_'); end
end

function tree = local_tag2cfgtree(cfg,tag)
if strcmpi(gettag(cfg),tag)
    tree = gettag(cfg);
elseif ~isa(cfg,'cfg_exbranch')
    tags = tagnames(cfg,1);
    for ii=1:length(tags)
        tree = local_tag2cfgtree(cfg.values{ii},tag);
        if ~isempty(tree)
            tree = [gettag(cfg) '.' tree];
            break;
        end
    end
else
    tree = [];
end