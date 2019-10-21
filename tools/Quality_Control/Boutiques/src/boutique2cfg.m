function [cfg, def] = boutique2cfg(json)
% Load single boutiques json

% CFG
cfg = cfg_cfg_call_system;
cfg.name = json.name;
cfg.tag = genvarname(lower(json.name));
if isfield(json,'description')
    try
        cfg.help = strsplit(json.description,newline);
    end
end

% DEF
def.cmd = json.command_line;

% groups (choose the first option as mandatory)
if isfield(json,'groups')
for ig = 1:length(json.groups)
    if iscell(json.groups)
        group = json.groups{ig};
    else
        group = json.groups(ig);
    end
    if isfield(group,'one_is_required') && group.one_is_required
        if isfield(group,'mutually_exclusive') && group.mutually_exclusive % set first as mandatory
%             newchoice        = cfg_choice;
%             newchoice.name   = group.name;
%             newchoice.help   = {group.description};
%             newchoice.tag    = group.id;
            for ic = 1:length(group.members)
                inputic = strcmp(group.members{ic},cellfun(@(x) x.id,json.inputs,'uni',0));
%                 switch json.inputs{inputic}.type
%                     case 'File'
%                         Nval = cellfun(@(x) strcmp(x.tag,'anyfilebranch'),cfg.val{1}.values);
%                     case 'String'
%                          Nval = cellfun(@(x) strcmp(x.tag,'stringbranch'),cfg.val{1}.values);
%                     case 'Number'
%                          Nval = cellfun(@(x) strcmp(x.tag,'evaluatedbranch'),cfg.val{1}.values);
%                 end
%                 tmp = cfg.val{1}.values{Nval};
%                 tmp.name = json.inputs{inputic}.name;
%                 tmp.tag  =  json.inputs{inputic}.id;
%                 newchoice.values{end+1} = tmp;
                if ic==1
                    json.inputs{inputic}.optional = 0;
                end
            end
            
          %  cfg.val{1}.values{end+1} = newchoice;
            
        else % set all as mandatory
            for ic = 1:length(group.members)
                inputic = strcmp(group.members{ic},cellfun(@(x) x.id,json.inputs,'uni',0));
                json.inputs{inputic}.optional = 0;
            end
        end
    end
end
end

% Mandatory inputs
inKey = {};
if isfield(json,'inputs')
def.inputs_ = {};
cfg.val{1}.val = {};
for ii = 1:length(json.inputs)
    if iscell(json.inputs)
        in = json.inputs{ii};
    else
        in = json.inputs(ii);
    end
    if ~in.optional
        inKey{end+1} = in.value_key;
        % help
        switch in.type
            case 'File'
                type = 'anyfilebranch';
                cfg.val{1}.val{end+1} = cfg.val{1}.values{3};
            case 'String'
                if strfind(lower(in.name),'directory')
                    type = 'directorybranch';
                    cfg.val{1}.val{end+1} = cfg.val{1}.values{4};
                else
                    type = 'stringbranch';
                    cfg.val{1}.val{end+1} = cfg.val{1}.values{2};
                end
            case 'Number'
                type = 'evaluatedbranch';
                cfg.val{1}.val{end+1} = cfg.val{1}.values{1};
        end
        if isfield(in,'value_choices')
            in.value_choices = cellfun(@char,in.value_choices,'uni',0);
            in.name = [in.name ': {' strjoin(in.value_choices,'/') '}'];
        end
        def.inputs_{end+1}.(type).help  = in.name;
        cfg.val{1}.val{end}.val{1}.val  = {in.name};
        % default-value
        if isfield(in,'default_value')
            def.inputs_{end}.(type).(strrep(type,'branch','')) = in.default_value;
            cfg.val{1}.val{end}.val{2}.val = {in.default_value};
        else
            def.inputs_{end}.(type).(strrep(type,'branch','')) = '<UNDEFINED>';
        end
        % replace flag [FLAG] --> %i1
        if ~isfield(in,'command_line_flag'), in.command_line_flag = ''; end
        def.cmd = strrep(def.cmd,in.value_key,[in.command_line_flag '%i' num2str(length(def.inputs_))]);
    else
        def.cmd = strrep(def.cmd,in.value_key,'');
    end
end
end

% Outputs
if isfield(json,'output_files')
    def.outputs_ = {};
    for io = 1:length(json.output_files)
        if iscell(json.output_files)
            out = json.output_files{io};
        else
            out = json.output_files(io);
        end

        if out.optional == 0
            def.outputs_{end+1}.outputs.help      = out.name;
            cfg.val{2}.val{end+1} = cfg.val{2}.values{1};
            cfg.val{2}.val{end}.val{1}.val = {out.name};
            if isfield(out,'path_template')
                for ikey = 1:length(inKey)
                    out.path_template = strrep(out.path_template,inKey{ikey},['%i' num2str(ikey)]);
                end
                [path,file,ext]                     = fileparts(out.path_template);
                if isempty(path)
                    file = strrep(file,'\',filesep);
                    file = strrep(file,'/',filesep);
                    [path,file,ext]                     = fileparts(file);
                end
                if isempty(path) || strcmp(path,'.')
                    def.outputs_{end}.outputs.directory = '<UNDEFINED>';
                else
                    def.outputs_{end}.outputs.directory = {path};
                    cfg.val{2}.val{end}.val{2}.val = {{path}};
                end
                if isempty(file)
                    def.outputs_{end}.outputs.string    = '<UNDEFINED>';
                else
                    def.outputs_{end}.outputs.string    = [file ext];
                    cfg.val{2}.val{end}.val{3}.val = {[file ext]};
                end
            else
                def.outputs_{end}.outputs.directory = '<UNDEFINED>';
                def.outputs_{end}.outputs.string = '<UNDEFINED>';
            end
            if isfield(out,'value_key')
                if ~isfield(out,'command_line_flag'), out.command_line_flag = ''; end
                def.cmd = strrep(def.cmd,out.value_key,[out.command_line_flag '%o' num2str(length(def.outputs_))]);
            end
        else
            if isfield(out,'value_key')
                def.cmd = strrep(def.cmd,out.value_key,'');
            end
        end
    end
end

if isfield(json,'container_image') && strcmp(json.container_image.type,'docker')
    def.usedocker_.dockerimg = json.container_image.image;
    cfg.val{4}.val{1}.val = {json.container_image.image};
else
    cfg.val{4}.val{1}.val = {''};
end

cfg.val{3}.val = {def.cmd};
        