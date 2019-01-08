function system_docker(dockerimage, cmd, varargin)
% Run a system function on a docker container
% system_docker(dockerimage, cmd, varargin)
%
% Example:
%   dockerimage = 'bids/mrtrix3_connectome';   
%   sourcei1 = 'path/source.nii';    
%   output = 'path/ref.nii';    
%   system_docker(dockerimage, 'fslmaths %i1 -Tmean %i2', sourcei1, output)
%
% Tanguy Duval, 2018, Toulouse Neuroimaging Center

in = varargin;
nopath = cellfun(@isempty,cellfun(@fileparts,in,'uni',0));
in(nopath) = cellfun(@(x) fullfile(pwd,x),in(nopath),'uni',0);

% DOCKERIFY
mountdir = '';
for k = 1:numel(in)
    mountdir = [mountdir '-v ' fileparts(in{k}) ':/i' num2str(k) ' '];
    dockerinfname{k} = strrep(in{k},[fileparts(in{k}) filesep],['/i' num2str(k) '/']);
end

% Replace token i%d by filenames
for ii=1:length(in)
    if ~isempty(dockerimage)
        cmd = strrep(cmd,sprintf('%%i%d',ii),dockerinfname{ii});
        cmd = strrep(cmd,[' ' sprintf('i%d',ii)],[' ' dockerinfname{ii} ' ']);
    else
        cmd = strrep(cmd,sprintf('%%i%d',ii),dockerinfname{ii});
        cmd = strrep(cmd,[' ' sprintf('i%d',ii)],[' ' in{ii} ' ']);
    end
end

% RUN SYSTEM COMMAND
if isempty(dockerimage)
    disp(['Running terminal command: ' cmd])
    [status, stdout]=system(cmd,'-echo');
    if status, error(sprintf('%s\n\n%s',cmd,stdout)); end
else % docker
    cmdcell = strsplit(cmd);
    cmddocker = ['docker run --entrypoint ' cmdcell{1} ' ' mountdir dockerimage ' ' strjoin(cmdcell(2:end))];
    disp(['Running terminal command: ' cmddocker])
    [status, stdout]=system(cmddocker,'-echo');
    if status, error(sprintf('%s\n\nRun on docker:\n%s\n\n%s',cmd,cmddocker,stdout)); end
end