function BIDS = bids_parser(bids_dir)
% BIDS_PARSER   Parse a BIDS (BRAIN IMAGING DATA STRUCTURE) directory path
%                as defined in http://bids.neuroimaging.io/bids_spec.pdf
%
% BIDS = BIDS_PARSER(bids_dir) reads the BIDS architecture located under bids_dir
% Returns the results in a structure with the fields: 
%       participants          -- infos extracted from participants.tsv
%       data_description      -- infos extracted from dataset_description.json
%       path                  -- copy of directory path bids_dir
%       subjects              -- M-by-1 structure with the following fields:
%         name                -- subject ID
%         session             -- session ID
%         path                -- path of subject/session
%        (anat)               -- N-by-1 structure with anatomical files
%        (type)               -- P-by-1 structure with <type> files
%           filename          -- name of the nifti file
%          (parameter)        -- content of file.json (e.g. FlipAngle)
%
% FILTER BY SUBJECT NAME AND SESSION:
%   SCAN = BIDS.subjects(strcmp({BIDS.subjects.name},'MATHIAS') & strcmp({BIDS.subjects.session},'postsurgery'))
%
% FILTER ONE SUBJECT BY MODALITY:
%   MODALITY = SCAN.anat(strcmp({SCAN.anat.modality},'FLAIR'));
%   filepath = fullfile(SCAN.path,'anat',MODALITY.filename);
%
% FILTER ALL SUBJECTS BY MODALITY:
%   filepath={}; 
%   for iscan = 1:length(BIDS.subjects)
%       filepath{iscan} =  fullfile(BIDS.subjects(iscan).path, 'anat', BIDS.subjects(iscan).anat(1).filename); 
%   end

%% Read study infos
if exist(fullfile(bids_dir,'participants.tsv'),'file')
    try
      BIDS.participants = readtable(fullfile(bids_dir,'participants.tsv'),'delimiter','\t');
    catch err
        disp(err.message)
        warning(['Could not read participants.tsv: ' err.message])
    end
else
    warning(['Could not find ' fullfile(bids_dir,'participants.tsv')]);
end
if exist(fullfile(bids_dir,'dataset_description.json'),'file')
    try
        BIDS.data_description = loadjson(fullfile(bids_dir,'dataset_description.json'));
    catch err
        disp(err.message)
        warning(['Could not read dataset_description.json: ' err.message])
    end
else
    warning(['Could not find ' fullfile(bids_dir,'dataset_description.json')]);
end

BIDS.path = bids_dir;
BIDS.subjects=[];

%% Loop over subjects
subjects = strrep(sct_tools_ls(fullfile(bids_dir,'sub-*')),'sub-','');
if isempty(subjects), warning(['no subjects found in ' bids_dir]); end
for isub = 1:length(subjects)
   subdir = fullfile(bids_dir,['sub-' subjects{isub}]);
   ses_list = sct_tools_ls(fullfile(subdir,'ses-*'));
   if isempty(ses_list), ses_list{1}=[]; end
   %% Loop over sessions
   for sesj=1:length(ses_list)
       % read session info if exists
       try
           sesinfo = readtable(fullfile(bids_dir,['sub-' subjects{isub}],['sub-' subjects{isub} '_sessions.tsv']),'delimiter','\t');
       end

       BIDS.subjects(end+1).name = subjects{isub};
       BIDS.subjects(end).path = fullfile(bids_dir,['sub-' subjects{isub}],ses_list{sesj});
       BIDS.subjects(end).session = strrep(ses_list{sesj},'ses-','');
       if isempty(BIDS.subjects(end).session), BIDS.subjects(end).session='none'; end
       
       %% Loop over nifti files
       sesdir = fullfile(subdir,ses_list{sesj});
       [im_list, im_path] = sct_tools_ls(fullfile(sesdir,'*.nii*'),0,1,2,1);
       for imj = 1:length(im_list)
           if strcmp(im_path{imj}(1:end-1),sesdir)
               type = 'image';
           else
               % type and filename
               [~,type] = fileparts(fileparts(im_path{imj}));
           end
           if isfield(BIDS.subjects(end),type)
             BIDS.subjects(end).(type)(end+1).filename = im_list{imj};
             else
              BIDS.subjects(end).(type).filename = im_list{imj};
           end
           
           % parse filename
           labels = regexp(BIDS.subjects(end).(type)(end).filename,[...
               '^sub-[a-zA-Z0-9]+' ...              % sub-<participant_label>
               '(?<ses>_ses-[a-zA-Z0-9]+)?' ...     % ses-<label>
               '(?<acq>_acq-[a-zA-Z0-9]+)?' ...     % acq-<label>
               '(?<run>_run-[a-zA-Z0-9]+)?' ...     % run-<index>
               '(?<rec>_rec-[a-zA-Z0-9]+)?' ...     % rec-<label>
               '.*_(?<modality>[a-zA-Z0-9]+)?' ...    % <modality>
               '\.nii(\.gz)?$'],'names');           % NIfTI file extension
            
           if isempty(labels)
               continue
           end
           for ff = fieldnames(labels)'
               BIDS.subjects(end).(type)(end).(ff{1}) = strrep(strrep(labels.(ff{1}),'_',''),[ff{1} '-'],'');
           end
           
           % read acquisition properties
           infofile = regexprep(fullfile(im_path{imj},im_list{imj}),'\.nii(\.gz)?','.json');
           if exist(infofile,'file')
               props = loadjson(infofile);
               for ff = fieldnames(props)'
                   BIDS.subjects(end).(type)(end).meta.(ff{1}) = props.(ff{1});
               end
           else
               disp(['no json file associated with ' fullfile(im_path{imj},im_list{imj})])
           end
       end
       
       try
           sesinfosesj = strcmp(strrep(sesinfo.session_id,'ses-',''),BIDS.subjects(end).session);
           for iinfo = setdiff(sesinfo.Properties.VariableNames,'session_id')
               BIDS.subjects(end).(iinfo{1}) = sesinfo.(iinfo{1}){sesinfosesj};
           end
       end

   end
    
end



end

function obj = readtable(file, varargin)

  if nargin >= 1
      [~, n, e] = fileparts(file);
      sep = ',';
      heads = 0;
  end

  if nargin >=2
    % check for delimiter
    DelimiterIndex = find( strcmp( varargin, 'delimiter' ) == 1);
    if isempty(DelimiterIndex)
      sep = ',';
    else
      sep = varargin{DelimiterIndex + 1};
    end

    % check for headers
    HeaderIndex = find( strcmp( varargin, 'HeaderLines' ) == 1);
    if isempty(HeaderIndex)
      heads = 0;
    else
      heads = varargin{HeaderIndex + 1};
    end

  end

  % open file
  if strcmp(e, '.xlsx')
    ret = xlsread(file);
  else
    ret = read_file(file, sep, heads);
  end

  if heads == 0
    try
    for ic = 1:size(ret,2)
      obj.(ret{1,ic}) = ret(2:end,ic);
    end
    catch
      obj = cell2table(ret, 'VariableNames', {n});
    end
  else
      obj = cell2table(ret, 'VariableNames', 'Var');
  end

end

function ret = read_file(file, sep, heads)
  % open file
  f = fopen(file);
  if f < 3
    error('Unable to open file')
    return
  end


  % skip headers
  if heads > 0
    for n = 1:heads
      fgetl(f);
    end
  end

  % count separators
  mark = ftell(f);
  tmp = fgetl(f);
  num_sep = numel(strfind(tmp,sprintf(sep)));
  fseek(f, mark, 'bof');

  % read values
  ret = fread (f, 'char=>char').';
  fclose(f);
  
  % check end of line
  if tmp(end)~=sprintf(sep)
	ret = regexprep (ret,'\r\n','\n');
    ret = regexprep (ret,[sep '\n'],'\n');
    ret = regexprep (ret,'\n',[sep '\n']);
    num_sep = num_sep + 1;
  end
  % remove all newlines
  ret=regexprep (ret,'(\n)+','');

  % parsing values
  ret = regexp(ret,sep,'split');

  % format output
  % delete empty last field
  if mod(size(ret,2),num_sep)~=0
    % yes, because we split after each separator
    if numel(ret{1,end}) == 0
      ret = ret(1,1:end-1);
    end

    if mod(size(ret,2),num_sep)==0
      ret = reshape (ret,num_sep,[])';
    end
  else
    ret = reshape (ret,num_sep,[])';
  end


end
