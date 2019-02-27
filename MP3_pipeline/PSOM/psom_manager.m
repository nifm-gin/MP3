function [] = psom_manager(path_logs,time_pipeline)
% Not meant to be used directly. See PSOM_RUN_PIPELINE.
%
% Manage the execution of a pipeline.
% SYNTAX: [] = PSOM_MANAGER(PATH_LOGS,TIME_PIPELINE)
% PATH_LOGS (string) the logs folder.
% TIME_PIPELINE (string) the time at which the pipeline was started.
%
% See licensing information in the code.

% Copyright (c) Pierre Bellec
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2015.
% Maintainer : pierre.bellec@criugm.qc.ca
% Keywords : pipeline
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

psom_gb_vars

%% SYNTAX
if ~exist('path_logs','var')
    error('Syntax: [] = psom_manager(path_logs,opt)')
end

%% Logs folder
if ~strcmp(path_logs(end),filesep)
    path_logs = [ path_logs filesep];
end

%% File names for the pipeline
file_pipeline     = [path_logs 'PIPE.mat'];
file_jobs         = [path_logs 'PIPE_jobs.mat'];
file_status       = [path_logs 'PIPE_status_init.mat'];
file_time         = [path_logs 'PIPE_time.mat'];
file_pipe_running = [path_logs 'PIPE.lock'];
file_heartbeat    = [path_logs 'heartbeat.mat'];
file_kill         = [path_logs 'PIPE.kill'];
file_pipe_end     = [path_logs 'PIPE.end'];
file_news_feed    = [path_logs 'news_feed.csv'];
file_config       = [path_logs 'PIPE_config.mat'];
path_worker       = [path_logs 'worker' filesep];

%% Load configuration
opt = load(file_config);

%% File names for the workers
for num_w = 1:opt.max_queued
    file_worker_heart{num_w} = sprintf('%spsom%i%sheartbeat.mat'  ,path_worker,num_w,filesep);
    file_worker_job{num_w}   = sprintf('%spsom%i%snew_jobs.mat'   ,path_worker,num_w,filesep);
    file_worker_ready{num_w} = sprintf('%spsom%i%snew_jobs.ready' ,path_worker,num_w,filesep);
    file_worker_reset{num_w} = sprintf('%spsom%i%sworker.reset'   ,path_worker,num_w,filesep);
    file_worker_end{num_w}   = sprintf('%spsom%i%sworker.end'     ,path_worker,num_w,filesep);
end
    
%% Check that the time of the pipeline matches the record
%% This check is done to ensure a new pipeline has not been started
%% after the manager was submitted for execution.
logs_time = load(file_time);
if ~strcmp(time_pipeline,logs_time.time_pipeline)
    fprintf('The time of the pipeline does not match the logs. I am quitting.')
    exit
end

%% Start heartbeat
if exist('OCTAVE_VERSION','builtin')    
    main_pid = getpid;
else
    main_pid = feature('getpid');
end
cmd = sprintf('psom_heartbeat(''%s'',''%s'',%i)',file_heartbeat,file_kill,main_pid);
if strcmp(gb_psom_language,'octave')
    instr_heartbeat = sprintf('"%s" %s "addpath(''%s''), %s,exit"',gb_psom_command_octave,gb_psom_opt_matlab,gb_psom_path_psom,cmd);
else 
    instr_heartbeat = sprintf('"%s" %s "addpath(''%s''), %s,exit"',gb_psom_command_matlab,gb_psom_opt_matlab,gb_psom_path_psom,cmd);
end 
system([instr_heartbeat '&']);
    
%% Check for the existence of the pipeline
if ~exist(file_jobs,'file') % Does the pipeline exist ?
    error('Could not find the pipeline file %s. Please use psom_run_pipeline instead of psom_manager directly.',file_jobs);
end

% a try/catch block is used to clean temporary file if the user is
% interrupting the pipeline of if an error occurs
try    
    
    %% Open the news feed file
    hf_news = fopen(file_news_feed,'w');
       
    %% Print general info about the pipeline
    msg_line1 = sprintf('Pipeline started on %s',datestr(clock));
    msg_line2 = sprintf('user: %s, host: %s, system: %s',gb_psom_user,gb_psom_localhost,gb_psom_OS);
    stars = repmat('*',[1 max(length(msg_line1),length(msg_line2))]);
    fprintf('\n%s\n%s\n%s\n',msg_line1,msg_line2,stars);
    
    %% Load the pipeline
    load(file_pipeline,'list_jobs','graph_deps');
    status = load(file_status);
    pipeline = load(file_jobs);
    nb_jobs = length(list_jobs);
    
    %% Initialize the mask of finished jobs
    mask_finished = false([length(list_jobs) 1]);
    for num_j = 1:length(list_jobs)
        mask_finished(num_j) = strcmp(status.(list_jobs{num_j}),'finished');
        if mask_finished(num_j)
            sub_add_line_log(hf_news,sprintf('%s , finished\n',list_jobs{num_j}),false);
        end
    end
    graph_deps(mask_finished,:) = 0;
    
    %% spread the news
    list_finished = list_jobs(mask_finished);
    for ff = 1:length(list_finished)
        sub_add_line_log(hf_news,sprintf('%s , %s\n',list_finished{ff},'finished'),false);
    end
    mask_deps = max(graph_deps,[],1)>0;
    mask_deps = mask_deps(:);
    nb_finished = sum(mask_finished); % The number of finished jobs                         
        
    %% Initialize the to-do list
    mask_todo = false([length(list_jobs) 1]);
    for num_j = 1:length(list_jobs)
        mask_todo(num_j) = strcmp(status.(list_jobs{num_j}),'none');
    end
    nb_todo = sum(mask_todo);    % The number of jobs to do
    
    %% Initialize miscallenaous variables
    psom_plan     = zeros(nb_jobs,1);          % a summary of which worker is running which job
    mask_running  = false(nb_jobs,1);          % A binary mask of running jobs
    mask_failed   = false(nb_jobs,1);          % A binary mask of failed jobs
    nb_failed     = 0;                         % The number of failed jobs
    nb_running    = 0;                         % The number of running jobs
    nb_checks     = 0;                         % The number of checks before printing a point
    worker_reset  = false(opt.max_queued,1);   % A binary list of workers that have been reset
    worker_ready  = false(opt.max_queued,1);   % A binary list of workers that are ready to receive jobs
    worker_active = false(opt.max_queued,1);   % A binary list of workers that are active
    nb_sch_worker = zeros(opt.max_queued,1);   % A list of the number of jobs scheduled for execution per worker   
    flag_point    = false; % A flag to indicate if a . was verbosed last
    
    %% Find the longest job name
    lmax = 0;
    for num_j = 1:length(list_jobs)
        lmax = max(lmax,length(list_jobs{num_j}));
    end
    
    %% Start submitting jobs
    flag_loop = true;
    while flag_loop

        %% Check for workers that have been reset
        for num_w = 1:opt.max_queued
            worker_reset(num_w) = psom_exist(file_worker_reset{num_w});
            if worker_reset(num_w)
                if opt.flag_verbose >= 1
                    fprintf('%s Worker %i has been reset.\n',datestr(clock),num_w);
                end
                psom_clean(file_worker_reset{num_w},struct('flag_verbose',false));
                nb_running = nb_running - sum(mask_running(psom_plan==num_w));
                mask_running(psom_plan==num_w) = false;
                mask_todo(psom_plan==num_w) = true;
                nb_todo = nb_todo + sum(psom_plan==num_w);
                nb_sch_worker(num_w) = Inf;
                psom_plan(psom_plan==num_w) = 0;
            end 
        end
        
        %% Check the state of workers
        %% and read the news
        flag_nothing_happened = true;
        
        %% Look for events
        events = sub_news(path_worker,psom_plan(mask_running),list_jobs(mask_running));
        if ~isempty(events)
            flag_nothing_happened = false;
            if flag_point 
                fprintf('\n')
            end
            flag_point = false;                    
        end
        
        %% Check worker status
        for num_w = 1:opt.max_queued
            worker_active(num_w) = ~psom_exist(file_worker_end{num_w})&&psom_exist(file_worker_heart{num_w}); 
            worker_ready(num_w) = worker_active(num_w)&&~psom_exist(file_worker_ready{num_w});
            if worker_active(num_w) 
                if (nb_sch_worker(num_w)==Inf)
                    nb_sch_worker(num_w) = 0;
                end
            else
                nb_sch_worker(num_w) = Inf; % The worker is not ready. Mark infinite number of jobs running to ensure no new submission will occur.
            end    
        end    
        
        %% Some verbose for the events
        for num_e = 1:size(events,1)
            %% Update status
            name_job = events{num_e,1};
            mask_job = strcmp(list_jobs,name_job);
            switch events{num_e,2}
                case 'failed'
                    nb_sch_worker(psom_plan(mask_job)) = nb_sch_worker(psom_plan(mask_job))-1;
                    nb_running = nb_running-1;
                    nb_failed = nb_failed+1;
                    mask_running(mask_job) = false;
                    mask_todo(mask_job) = false;
                    mask_failed(mask_job) = true;
                    % Remove the children of the failed job from the to-do list
                    mask_child = sub_find_children(mask_job',graph_deps);
                    mask_todo(mask_child) = false; 
                    psom_plan(mask_job) = 0;
                    msg = sprintf('%s %s%s failed    ',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));
                case 'finished'
                    nb_sch_worker(psom_plan(mask_job)) = nb_sch_worker(psom_plan(mask_job))-1;
                    nb_running = nb_running-1;
                    nb_finished = nb_finished+1;
                    mask_running(mask_job) = false;
                    mask_todo(mask_job) = false;
                    mask_finished(mask_job) = true;
                    graph_deps(mask_job,:) = false;
                    psom_plan(mask_job) = 0;
                    msg = sprintf('%s %s%s finished  ',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));
                    otherwise 
                        error('%s is an unkown status',events{num_e,2})
            end
            %% Add to the news feed
            sub_add_line_log(hf_news,sprintf('%s , %s\n',name_job,events{num_e,2}),false);
            if opt.flag_verbose && ~isempty(msg)
                fprintf('%s (%i run | %i fail | %i done | %i left)\n',msg,nb_running,nb_failed,nb_finished,nb_todo);
            end
        end
             
        %% Flush into the news  feed
        if strcmp(gb_psom_language,'octave')
            fflush(hf_news);
        end  
        
        %% Update the dependency mask
        if ~flag_nothing_happened
            mask_deps = max(graph_deps,[],1)>0;
            mask_deps = mask_deps(:);
        end
        
        %% Time to (try to) submit jobs !!
        list_num_to_run = find(mask_todo&~mask_deps);
        mask_new_submit = false(opt.max_queued,1);
        nb_to_submit = length(list_num_to_run);
        tag = [];
        curr_job = 0;
        pipe_sub = cell(opt.max_queued,1);
        nb_sch_worker_r = nb_sch_worker;
        nb_sch_worker_r(~worker_ready) = Inf;
        while (min(nb_sch_worker_r)<=opt.max_buffer)&&(curr_job<length(list_num_to_run))
            curr_job = curr_job+1;
            [val,ind] = min(nb_sch_worker_r);
            name_job = list_jobs{list_num_to_run(curr_job)};
            pipe_sub{ind}.(name_job) = pipeline.(name_job);
            mask_new_submit(ind) = true;
            nb_sch_worker(ind) = nb_sch_worker(ind)+1;
            nb_sch_worker_r(ind) = nb_sch_worker_r(ind)+1;
            mask_todo(list_num_to_run(curr_job)) = false;
            psom_plan(list_num_to_run(curr_job)) = ind;
            mask_running(list_num_to_run(curr_job)) = true;
            nb_running = nb_running+1;
            nb_todo = nb_todo-1;
            if opt.flag_verbose
                msg = sprintf('%s %s%s submitted ',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));
                fprintf('%s (%i run | %i fail | %i done | %i left)\n',msg,nb_running,nb_failed,nb_finished,nb_todo);
            end
        end
                        
        if opt.flag_verbose >= 3
            if ~any(nb_sch_worker~=Inf)
                max_sch_worker = 0;
                min_sch_worker = 0;
            else
                max_sch_worker = max(nb_sch_worker(nb_sch_worker~=Inf));
                min_sch_worker = min(nb_sch_worker(nb_sch_worker~=Inf));
            end
            fprintf('%i/%i jobs submitted, # jobs per worker: %i to %i. %i workers unavailable.\n',curr_job,nb_to_submit,min_sch_worker,max_sch_worker,sum(nb_sch_worker==Inf))
        end
        
        %% Mark new submissions as ready to process
        for num_w = 1:opt.max_queued
            if mask_new_submit(num_w)
                flag_nothing_happened = false;
                tmp = pipe_sub{num_w};
                save(file_worker_job{num_w},'-struct','tmp');
                save(file_worker_ready{num_w},'tag');
            end
        end
        
        %% End workers that do not have jobs left to do
        mask_idle = worker_active&(nb_sch_worker==0);
        if (sum(mask_todo) < sum(worker_active)) && any(mask_idle)
            list_idle = find(mask_idle);
            list_idle = list_idle(1:min(length(list_idle),sum(worker_active)-sum(mask_todo)));
            for num_w = list_idle(:)'
                tmp = [];
                save(file_worker_end{num_w},'tmp');
                if opt.flag_verbose>=1
                    fprintf('%s Stopping idle worker %i (not enough jobs left to do).\n',datestr(clock),num_w);
                end
            end
            if opt.flag_verbose>=2
                fprintf('%i active worker(s) left, %i job(s) left to submit.\n',sum(worker_active)-length(list_idle),sum(mask_todo));
            end
        end
        
        %% Sleep if nothing happened
        flag_loop = (any(mask_todo) || any(mask_running)) && psom_exist(file_pipe_running) && ~psom_exist(file_pipe_end);
        if flag_nothing_happened && flag_loop
            pause(opt.time_between_checks)
         
            if (nb_checks >= opt.nb_checks_per_point)
                nb_checks = 0;
                if opt.flag_verbose
                    fprintf('.');
                    flag_point = true;
                end
                nb_checks = nb_checks+1;
            else
                nb_checks = nb_checks+1;
            end
        end
    end % While there are jobs to do
    
catch
    
    errmsg = lasterror;        
    fprintf('\n\n******************\nSomething went bad ... the pipeline has FAILED !\nThe last error message occured was :\n%s\n',errmsg.message);
    if isfield(errmsg,'stack')
        for num_e = 1:length(errmsg.stack)
            fprintf('File %s at line %i\n',errmsg.stack(num_e).file,errmsg.stack(num_e).line);
        end
    end
    if exist('file_pipe_running','var')
        if exist(file_pipe_running,'file')
            delete(file_pipe_running); % remove the 'running' tag
        end
    end
    
    %% Close the log file
    fclose(hf_news);
    status_pipe = 1;
    return
end

%% Print general info about the pipeline
msg_line1 = sprintf('Pipeline terminated on %s',datestr(now));
stars = repmat('*',[1 length(msg_line1)]);
if opt.flag_verbose
    fprintf('\n%s\n%s\n',stars,msg_line1);
end

%% Report if the lock file was manually removed
if exist('file_pipe_running','var')&&~psom_exist(file_pipe_running)
    fprintf('The pipeline manager was interrupted because the .lock file was manually deleted.\n');
end

%% Report if the pipeline was ended by the deamon 
if psom_exist(file_pipe_end)
    fprintf('The pipeline manager was interrupted because there are no active workers available.\n');
end

%% Print a list of failed jobs
list_num_failed = find(mask_failed);
list_num_failed = list_num_failed(:)';
flag_any_fail = ~isempty(list_num_failed);

if flag_any_fail
    if length(list_num_failed) == 1
        fprintf('1 job has failed.\n',length(list_num_failed));
    else
        fprintf('%i jobs have failed.\n',length(list_num_failed));
    end
end

%% Print a list of jobs that could not be processed
list_num_none = find(~mask_finished&~mask_failed);
list_num_none = list_num_none(:)';
if ~isempty(list_num_none)
    if length(list_num_none) == 1
        fprintf('1 job could not be processed due to failures or interruption.\n');
    else
        fprintf('%i jobs could not be processed due to failures or interruption.\n', length(list_num_none));
    end
end

%% Suggest using psom_pipeline_visu
if flag_any_fail
    fprintf('Use psom_pipeline_visu to access logs, e.g.:\n\n   psom_pipeline_visu(''%s'',''log'',''%s'');\n\n',path_logs,list_jobs{list_num_failed(1)});
    fprintf('You can also get the list of failed jobs using the command:\n\n   psom_pipeline_visu(''%s'',''failed'');\n',path_logs);
end

%% Give a final one-line summary of the processing
if ~flag_any_fail&&isempty(list_num_none)
    fprintf('All jobs have been successfully completed.\n');
end

%% Close the log file
fclose(hf_news);
status_pipe = double(flag_any_fail);

%% Terminate the pipeline
if exist('file_pipe_running','var')&&psom_exist(file_pipe_running)
    psom_clean(file_pipe_running,struct('flag_verbose',false)); % remove the PIPE.lock file
end

%%%%%%%%%%%%%%%%%%
%% subfunctions %%
%%%%%%%%%%%%%%%%%%

%% Find the children of a job
function mask_child = sub_find_children(mask,graph_deps)
% GRAPH_DEPS(J,K) == 1 if and only if JOB K depends on JOB J. GRAPH_DEPS =
% 0 otherwise. This (ugly but reasonably fast) recursive code will work
% only if the directed graph defined by GRAPH_DEPS is acyclic.
% MASK_CHILD(NUM_J) == 1 if the job NUM_J is a children of one of the job
% in the boolean mask MASK and the job is in MASK_TODO.
% This last restriction is used to speed up computation.

if max(double(mask))>0
    mask_child = max(graph_deps(mask,:),[],1)>0;    
    mask_child_strict = mask_child & ~mask;
else
    mask_child = false(size(mask));
end

if any(mask_child)
    mask_child = mask_child | sub_find_children(mask_child_strict,graph_deps);
end

%% Get news
function events = sub_news(path_w,num_worker,list_job_worker);
mask_completed = false(length(list_job_worker),1);
events = cell(length(list_job_worker),2);
if isempty(list_job_worker)
    return
end
events(:,1) = list_job_worker;
for jj = 1:length(list_job_worker)
    job_name = list_job_worker{jj};
    file_finished = [path_w sprintf('psom%i',num_worker(jj)) filesep job_name '.finished'];
    file_failed   = [path_w sprintf('psom%i',num_worker(jj)) filesep job_name '.failed'];
    flag_finished = psom_exist(file_finished,false);
    flag_failed = ~flag_finished&&psom_exist(file_failed,false);
    mask_completed(jj) = flag_failed || flag_finished;
    if flag_finished
        events{jj,2} = 'finished';
    elseif flag_failed
        events{jj,2} = 'failed';
    end
end
events = events(mask_completed,:);

function [] = sub_add_line_log(file_write,str_write,flag_verbose);

if flag_verbose
    fprintf('%s',str_write)
end
fprintf(file_write,'%s',str_write);
