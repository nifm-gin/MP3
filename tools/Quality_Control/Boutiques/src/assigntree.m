function cfg = assigntree(cfg, cfg_new, tree)
% Assign a config with tree
% cfg = assigntree(cfg, cfg_new, tree)
%
% Example:
%   cfg = assigntree(cfg, cfg_new, {'Boutiques' 'Ants' 'Registration' 'rigid'})

% generate cfg_CML

jj=2;
iscfg_callsystem=false;
new = isempty(cfg);
if ~new
    subs = {};
    cfg_tmp = cfg;
    for jj=2:length(tree)
        ind = cell2mat(cellfun(@(x) strcmpi(x.tag,tree{jj}), cfg_tmp.values,'uni',0));
        if any(ind)
            subs = [subs, {'.','values','{}',{find(ind,1)}}];
            cfg_tmp = subsref(cfg,substruct(subs{:}));
            if isa(cfg_tmp,'cfg_exbranch')
                iscfg_callsystem = true;
            end
        else
            break
        end
    end
    if ~iscfg_callsystem
        subs = [subs, {'.','values','{}',{length(cfg_tmp.values)+1}}];
        subs = substruct(subs{:});
    end
end

if ~iscfg_callsystem
    
    for ii=length(tree)-1:-1:(jj-1)
        cfg_tmp         = cfg_choice;
        cfg_tmp.tag     = genvarname(lower(tree{ii}));
        cfg_tmp.name    = tree{ii};
        cfg_tmp.values  = {cfg_new};
        cfg_new = cfg_tmp;
    end
    if new
        cfg = cfg_new;
    else
        cfg = subsasgn(cfg,subs,cfg_new.values{1});
    end
end
