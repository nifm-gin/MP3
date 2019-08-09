function [cfg,def] = boutiques2cfg(parentfolder)
% load Boutiques jsons in entire folder

cfg = [];
def = struct;
list = spm_select('FPListRec',parentfolder,'.*\.json$');
list = cellstr(list);
path = cellfun(@fileparts,list,'uni',0);

for ifile = 1:length(list)
    % load boutique json
    tree = strsplit(strrep(path{ifile},[fileparts(parentfolder) filesep],''),filesep);
    json = spm_jsonread(list{ifile});
    % convert boutique to cfg
    if isfield(json,'matlab_version')
        [icfg, idef] = boutiquematlab2cfg(json);
    else
        [icfg, idef] = boutique2cfg(json);
    end
    tree{end+1} = json.name;
    
    % generate full tree
    cfg = assigntree(cfg, icfg, tree);
    
    % compute def
    tree{1} = 'def'; 
    code = gencode(idef,strjoin(cellfun(@genvarname,lower(tree),'uni',0),'.'))';
    for iii = 1:length(code); eval(code{iii}); end

end
