function cfg_save_call_system(handles)
try
    dflag = 0;
    udmodlist = get(handles.modlist, 'userdata');
    cmod = get(handles.modlist, 'value');
    udmodule = get(handles.module, 'userdata');
    citem = 1;
    ciid = {udmodlist.cjob udmodlist.id{cmod} udmodule.id{citem}};
    contents = cfg_ui_util('showitem', ciid, dflag);
    [tag, val] = cfg_util('harvest', ciid{:});
    if ~isempty(contents{10}) &&  strcmp(func2str(contents{10}),'@(job)cfg_run_call_system(''save'',job)')
        feval(contents{10}, val);
    else
        warndlg(['"' contents{1} '" is not a "Boutiques>Generic command"'],'Could not save Boutiques Descriptor')
    end
end


