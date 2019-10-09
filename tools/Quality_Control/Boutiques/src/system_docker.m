function varargout = system_docker(dockerimage, cmd, varargin)
% Run a system function on a docker container
% system_docker(dockerimage, cmd, varargin)
%
% Example:
%   Test if docker is running correctly
%     system_docker('hello-world','/hello')
%
%   Get ANTS Registration help:
%     system_docker('bids/mrtrix3_connectome', 'antsRegistration --help')
%
%   RUN fslmaths:
%     dockerimage = 'bids/mrtrix3_connectome';   
%     sourcei1 = 'path/source.nii';    
%     output = 'path/out.nii';    
%     system_docker(dockerimage, 'fslmaths %i1 -Tmean %i2', sourcei1, output)
%
% Tanguy Duval, 2018, Toulouse Neuroimaging Center

in = varargin;
nopath = cellfun(@isempty,cellfun(@fileparts,in,'uni',0));
in(nopath) = cellfun(@(x) fullfile(pwd,x),in(nopath),'uni',0);

% Test if docker is correctly installed
if ~isempty(dockerimage)
    [status, res] = system('docker');
    if status
        warndlg(sprintf('%s\n ... Docker is not up and running...\nPlease install and run docker, then launch matlab from a terminal by using the following command:\n "%s"',res,fullfile(matlabroot,'bin','matlab')))
    end
end

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
    if nargout==0 && status, error(sprintf('%s\n\n%s',cmd,stdout)); end
else % docker
    cmdcell = strsplit(cmd);
    if strcmpi(cmdcell{1},'defaultEntrypoint')
        cmddocker = ['docker run ' mountdir dockerimage ' ' strjoin(cmdcell(2:end))];
    else
        cmddocker = ['docker run --entrypoint ' cmdcell{1} ' ' mountdir dockerimage ' ' strjoin(cmdcell(2:end))];
    end
    disp(['Running terminal command: ' cmddocker])
    [status, stdout]=system(cmddocker,'-echo');
    if nargout==0 && status, error(sprintf('%s\n\nRun on docker:\n%s\n\n%s',cmd,cmddocker,stdout)); end
end

varargout{1} = status;
varargout{2} = stdout;
