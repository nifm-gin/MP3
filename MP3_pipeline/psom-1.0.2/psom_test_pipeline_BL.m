function [pipeline,opt] = psom_test_pipeline_BL(path_demo,opt)

psom_gb_vars

if nargin<1||isempty(path_demo)
    local_path_demo = gb_psom_path_demo;    
else
    local_path_demo = path_demo;
end

if ~strcmp(local_path_demo(end),filesep)
    local_path_demo = [local_path_demo filesep];
end

% Set up the options to run the pipeline
opt.path_logs = [local_path_demo 'logs' filesep];  % where to store the log files

if ~isfield(opt,'flag_test')
    flag_test = false;
else
    flag_test = opt.flag_test;
    opt = rmfield(opt,'flag_test');
end

if ~isfield(opt,'flag_pause')
    opt.flag_pause = true;
end

% Job "sample" :    No input, generate a random vector a
command = 'a = randn([opt.nb_samps 1]); save(files_out,''a'')';
pipeline.sample.command      = command;
pipeline.sample.files_out    = [local_path_demo 'sample.mat'];
pipeline.sample.opt.nb_samps = 10;

% % Job "quadratic" : Compute a.^2 and save the results
% command = 'load(files_in); b = a.^2; save(files_out,''b'')';
% pipeline.quadratic.command   = command;
% pipeline.quadratic.files_in  = pipeline.sample.files_out;
% pipeline.quadratic.files_out = [local_path_demo 'quadratic.mat']; 

% Job "quadratic" : Compute a.^2 and save the results
command = 'load(files_in); b = psom_sum(a); save(files_out,''b'')';
pipeline.quadratic.command   = command;
pipeline.quadratic.files_in  = pipeline.sample.files_out;
pipeline.quadratic.files_out = [local_path_demo 'quadratic.mat'];
% 
% Adding a job "cubic" : Compute a.^3 and save the results
command = 'load(files_in); c = a.^3; save(files_out,''c'')';
pipeline.cubic.command       = command;
pipeline.cubic.files_in      = pipeline.sample.files_out;
pipeline.cubic.files_out     = [local_path_demo 'cubic.mat']; 

% Adding a job "sum" : Compute a.^2+a.^3 and save the results
command = 'load(files_in{1}); load(files_in{2}); d = b+c, save(files_out,''d'')';
pipeline.sum.command       = command;
pipeline.sum.files_in{1}   = pipeline.quadratic.files_out;
pipeline.sum.files_in{2}   = pipeline.cubic.files_out;
pipeline.sum.files_out     = [local_path_demo 'sum.mat'];

if flag_test
    return
end
if opt.flag_pause
    psom_visu_dependencies(pipeline);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the demo folder %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msg = 'The demo is about to remove the content of the following folder and save the demo results there:';
msg2 = local_path_demo;
msg3 = 'Press CTRL-C to stop here or any key to continue.';
stars = repmat('*',[1 max([length(msg),length(msg2),length(msg3)])]);
fprintf('\n%s\n%s\n%s\n%s\n%s\n\n',stars,msg,msg2,msg3,stars);
if opt.flag_pause
    pause
end

if exist(local_path_demo,'dir')
    rmdir(local_path_demo,'s');
end
psom_mkdir(local_path_demo);

%%%%%%%%%%%%%%%%%%%%
%% Run a pipeline %%
%%%%%%%%%%%%%%%%%%%%
% msg   = 'The demo is about to execute the toy pipeline.';
% msg2  = 'Press CTRL-C to stop here or any key to continue.';
% stars = repmat('*',[1 max(length(msg),length(msg2))]);
% fprintf('\n%s\n%s\n%s\n%s\n\n',stars,msg,msg2,stars);
% if opt.flag_pause
%     pause
% end

% The following line is running the pipeline manager on the toy pipeline
psom_run_pipeline(pipeline,opt);
