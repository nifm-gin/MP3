function qc_reload(qcdir)
jsonfname = fullfile(qcdir,'qc_results.json');
% read json
fid = fopen(jsonfname,'r');
txt = '';
while ~feof(fid)
    txt = [txt fgetl(fid)];
end
fclose(fid);
txt = regexp(txt,'{.+}','match');
txt = txt{1};

% append to html replace
fid = fopen(fullfile(fileparts(mfilename('fullpath')),'htmltemplate','index.html'),'r');
f = fread(fid,'*char')';
fclose(fid);

if strcmp(txt(1),'['), txt(1)=[]; txt(end)=[]; end
f = regexprep(f,'var sct_data =.*;</script>',['var sct_data = [' txt '];</script>']);
fid = fopen(fullfile(qcdir,'index.html'),'w');
fprintf(fid,'%s',f);
fclose(fid);

copyfile(fullfile(fileparts(mfilename('fullpath')),'htmltemplate','_assets'),fullfile(qcdir,'_assets'))
