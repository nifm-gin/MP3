function cfg_preview_parsebids(bids_dir)
if iscell(bids_dir), bids_dir = bids_dir{1}; end
job.parent{1} = bids_dir;
BIDS = cfg_run_parsebids(job);
hh = figure('Name',['BIDS directory (' bids_dir ')']);
set(hh,'MenuBar','none')
uitable('Data',[{BIDS.subjects.name}', {BIDS.subjects.session}'],'ColumnName',{'subject', 'session'},'Units', 'Normalized', 'Position',[0, 0, 1, 1])