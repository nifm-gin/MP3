function status_pipe = psom_deamon(path_logs)
% Not meant to be used directly. See PSOM_RUN_PIPELINE.
%
% Start workers, a pipeline manager and a garbage collector
% SYNTAX: [] = PSOM_DEAMON(PATH_LOGS)
% PATH_LOGS (string) the logs folder.
% See licensing information in the code.

% Copyright (c) Pierre Bellec, 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up default values for inputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SYNTAX
if ~exist('path_logs','var')
    error('Syntax: [] = psom_deamon(path_logs)')
end

%% fixing the format of path_logs
if ~strcmp(path_logs(end),filesep)
    path_logs = [path_logs filesep];
end

%% Constants
time_death = 300; % Time before a worker is considered dead

%% Pipeline file names
file_pipeline     = [path_logs 'PIPE.mat'];
file_pipe_running = [path_logs 'PIPE.lock'];
file_pipe_heart   = [path_logs 'heartbeat.mat'];
file_kill         = [path_logs 'PIPE.kill'];
file_end          = [path_logs 'PIPE.end'];
file_time         = [path_logs 'PIPE_time.mat'];
path_tmp          = [path_logs 'tmp' filesep];
path_worker       = [path_logs 'worker' filesep];
path_garbage      = [path_logs 'garbage' filesep];
file_conf         = [path_logs 'PIPE_config.mat'];

%% Load opt
opt = load(file_conf);

%% Worker file names
for num_w = 1:opt.max_queued
    name_worker{num_w} = sprintf('psom%i',num_w);
    file_worker_heart{num_w} = [path_worker name_worker{num_w} filesep 'heartbeat.mat'];
    file_worker_reset{num_w} = [path_worker name_worker{num_w} filesep 'worker.reset'];
    file_worker_end{num_w}   = [path_worker name_worker{num_w} filesep 'worker.end'];
end
name_worker{opt.max_queued+1} = 'psom_manager';
file_worker_heart{opt.max_queued+1} = [path_logs 'heartbeat.mat'];
name_worker{opt.max_queued+2} = 'psom_garbage';
file_worker_heart{opt.max_queued+2} = [path_garbage 'heartbeat.mat'];
psom_mkdir(path_tmp);

%% Check for the existence of the pipeline
if ~exist(file_pipeline,'file') % Does the pipeline exist ?
    error('Could not find the pipeline file %s. Please use psom_run_pipeline instead of psom_deamon directly.',file_pipeline);
end

%% Save the time
time_pipeline = datestr(clock);
save(file_time,'time_pipeline');

% a try/catch block is used to clean temporary file if the user is
% interrupting the pipeline of if an error occurs
try    
       
    %% Print general info about the pipeline
    if opt.flag_verbose
        fprintf('Deamon started on %s\n',datestr(clock));
    end
    
    %% Track refresh times for workers
    % (#workers + manager + garbage collector) x 6 (clock info) x 2
    % the first table is to record the last documented active time for the heartbeat
    % the second table is to record the time elapsed since a new heartbeat was detected
    tab_refresh(:,:,1) = -ones(opt.max_queued+2,6);
    tab_refresh(:,:,2) = repmat(clock,[opt.max_queued+2 1]);
    
    %% Initialize miscallenaous variables
    nb_resub    = 0;                   % Number of resubmission               
    nb_checks   = 0;                   % Number of checks to print a points
    nb_points   = 0;                   % Number of printed points
    nb_chars_logs = 0;                 % Number of characters printed from the pipeline history                
    flag_pipe_finished = false;        % Is the pipeline finished?
    flag_started = false([opt.max_queued+2 1]); % Has the worker ever started? two last entries are for the PM and the GC
    flag_alive   = false([opt.max_queued+2 1]); % Is the worker alive? two last entries are for the PM and the GC
    flag_wait    = false([opt.max_queued+2 1]); % Are we waiting for the worker to start? two last entries are for the PM and the GC
    flag_end     = false([opt.max_queued+2 1]); % did the manager request for the worker to end?
    flag_dead    = false([opt.max_queued+2 1]); % did the worker die with no opportunity to restart?
    
    %% Create logs folder for each worker
    path_worker_w = cell(opt.max_queued,1);
    for num_w = 1:opt.max_queued
        path_worker_w{num_w} = sprintf('%spsom%i%s',path_worker,num_w,filesep);
        if psom_exist(path_worker_w{num_w})
            psom_clean(path_worker_w{num_w},struct('flag_verbose',false));
        end
        psom_mkdir(path_worker_w{num_w});
    end
    
    %% General options to submit scripts
    opt_script.path_search    = file_pipeline;
    opt_script.mode           = opt.mode;
    opt_script.init_matlab    = opt.init_matlab;
    opt_script.flag_debug     = opt.flag_verbose >= 2;        
    opt_script.shell_options  = opt.shell_options;
    opt_script.command_matlab = opt.command_matlab;
    opt_script.qsub_options   = opt.qsub_options;
    opt_script.name_job       = ''; % to be specified
    
    %% Options for submission of the pipeline manager
    opt_logs_pipe.txt    = [path_logs 'PIPE_history.txt'];
    opt_logs_pipe.eqsub  = [path_logs 'PIPE.eqsub'];
    opt_logs_pipe.oqsub  = [path_logs 'PIPE.oqsub'];
    opt_logs_pipe.failed = [path_logs 'PIPE.failed'];
    opt_logs_pipe.exit   = [path_logs 'PIPE.exit'];   
    opt_pipe = opt_script;
    opt_pipe.mode = opt.mode_pipeline_manager;
    opt_pipe.name_job = 'psom_manager';   
    cmd_pipe = sprintf('psom_manager(''%s'',''%s'');',path_logs,time_pipeline);    
    if ispc % this is windows
        script_pipe = [path_tmp filesep 'psom_manager.bat'];
    else
        script_pipe = [path_tmp filesep 'psom_manager.sh'];
    end
    
    %% Options for submission of the garbage collector
    opt_logs_garb.txt    = [path_garbage 'garbage_history.txt'];
    opt_logs_garb.eqsub  = [path_garbage 'garbage.eqsub'];
    opt_logs_garb.oqsub  = [path_garbage 'garbage.oqsub'];
    opt_logs_garb.failed = [path_garbage 'garbage.failed'];
    opt_logs_garb.exit   = [path_garbage 'garbage.exit'];   
    opt_garb = opt_script;
    opt_garb.mode = opt.mode_garbage;
    opt_garb.name_job = 'psom_garbage';   
    cmd_garb = sprintf('psom_garbage(''%s'',''%s'',false);',path_logs,time_pipeline);    
    if ispc % this is windows
        script_garb = [path_tmp filesep 'psom_garbage.bat'];
    else
        script_garb = [path_tmp filesep 'psom_garbage.sh'];
    end
    
    %% Options for submission of the workers
    for num_w = 1:opt.max_queued
        opt_logs_worker(num_w).txt    = sprintf('%spsom%i%sworker.log',path_worker,num_w,filesep);
        opt_logs_worker(num_w).eqsub  = sprintf('%spsom%i%sworker.eqsub',path_worker,num_w,filesep);
        opt_logs_worker(num_w).oqsub  = sprintf('%spsom%i%sworker.oqsub',path_worker,num_w,filesep);
        opt_logs_worker(num_w).failed = sprintf('%spsom%i%sworker.failed',path_worker,num_w,filesep);
        opt_logs_worker(num_w).exit   = sprintf('%spsom%i%sworker.exit',path_worker,num_w,filesep);   
        opt_worker(num_w) = opt_script;
        opt_worker(num_w).name_job = name_worker{num_w};
        cmd_worker{num_w} = sprintf('psom_worker(''%s'',''%s'',%i,''%s'');',path_worker_w{num_w},path_logs,num_w,time_pipeline);
        if ispc % this is windows
            script_worker{num_w} = [path_tmp filesep opt_worker(num_w).name_job '.bat'];
        else
            script_worker{num_w} = [path_tmp filesep opt_worker(num_w).name_job '.sh'];
        end
    end
    
    %% Start submitting jobs
    while ~flag_pipe_finished 
    
        %% Check the heartbeats
        for num_w = 1:(opt.max_queued+2)
        
            %% Check for the presence of the heartbeat
            flag_heartbeat = psom_exist(file_worker_heart{num_w});
            
            if ~flag_heartbeat
                if opt.flag_verbose == 3
                    fprintf('No heartbeat for process %s\n',name_worker{num_w})
                end
            else
                if any(tab_refresh(num_w,:,1)<0)
                    % this is the first time an active time is collected
                    % simply update tab_refresh
                    tab_refresh(num_w,:,1) = 0;
                    flag_started(num_w) = true;
                    if opt.flag_verbose == 2
                        fprintf('First time heartbeat for process %s\n',name_worker{num_w})
                    end
                else
                    try
                        refresh_time = load(file_worker_heart{num_w});
                        test_change = etime(refresh_time.curr_time,tab_refresh(num_w,:,1))>1;
                    catch
                        % The heartbeat is unreadable
                        % Assume this is a race condition
                        % Consider no heartbeat was detected
                        if opt.flag_verbose == 2
                            fprintf('There was a problem reading the heartbeat of process %s.\n',name_worker{num_w})
                        end
                        test_change = false;
                    end

                    if test_change
                        % I heard a heartbeat!    
                        tab_refresh(num_w,:,1) = refresh_time.curr_time;
                        tab_refresh(num_w,:,2) = clock;
                        flag_alive(num_w) = true;
                        flag_wait(num_w) = false;
                        if opt.flag_verbose >= 3
                            fprintf('I heard a heartbeat for process %s\n',name_worker{num_w})
                        end
                    else 
                        % I did not hear a heartbeat
                        % how long has it been?
                        elapsed_time = etime(clock,tab_refresh(num_w,:,2));
                        if opt.flag_verbose >= 3
                            fprintf('No heartbeat in %1.2fs for process %s\n',elapsed_time,name_worker{num_w})
                        end
                        
                        % Check if the manager has ended the worker
                        if num_w<=opt.max_queued
                            test_end = psom_exist(file_worker_end{num_w});
                            if test_end&&~flag_end(num_w)
                                if (opt.flag_verbose>=2)
                                    fprintf('%s Worker %s terminated by the manager.\n',datestr(clock),name_worker{num_w});
                                end
                                flag_alive(num_w) = false;
                                flag_wait(num_w) = false;
                            end 
                            flag_end(num_w) = test_end;
                            
                        else
                            flag_end(num_w) = false;
                        end
                        
                        if (elapsed_time > time_death)&&~flag_end(num_w)
                            if opt.flag_verbose
                                fprintf('%s No heartbeat for process %s, counted as dead.\n',datestr(clock),name_worker{num_w});
                            end 
                            % huho it's been too long without a heartbeat, he's dead Jim
                            flag_alive(num_w) = false;
                            flag_wait(num_w) = false;
                        else
                            % Not that long, just go on
                            flag_alive(num_w) = true;
                            flag_wait(num_w) = false;
                        end
                    end
                end
            end
        end
        
        %% Now start workers
        list_worker = [opt.max_queued+[1 2] 1:opt.max_queued]; % first start the manager, then the collector, then workers
        for num_w = list_worker
            if ~flag_end(num_w)&&~flag_wait(num_w)&&~flag_alive(num_w)
             
                if flag_started(num_w)
                    if num_w <= opt.max_queued
                        %% Cleaning the worker folder
                        psom_clean([path_worker name_worker{num_w} filesep],struct('flag_verbose',false));
                        psom_mkdir([path_worker name_worker{num_w} filesep]);
                    elseif num_w == opt.max_queued+1
                        psom_clean_logs(path_logs);
                    elseif num_w == opt.max_queued+2
                        if psom_exist(path_garbage)
                            psom_clean(path_garbage,struct('flag_verbose',false));
                        end
                        psom_mkdir(path_garbage);
                    end
                end
                
                if (nb_resub<opt.nb_resub)||~flag_started(num_w)
                    if num_w <= opt.max_queued
                        [flag_failed,msg] = psom_run_script(cmd_worker{num_w},script_worker{num_w},opt_worker(num_w),opt_logs_worker(num_w));
                    elseif num_w == opt.max_queued+1
                        [flag_failed,msg] = psom_run_script(cmd_pipe,script_pipe,opt_pipe,opt_logs_pipe);
                    elseif num_w == opt.max_queued+2
                        [flag_failed,msg] = psom_run_script(cmd_garb,script_garb,opt_garb,opt_logs_garb);
                    end
                    if ~flag_failed
                        flag_wait(num_w) = true;
                        tab_refresh(num_w,:,1) = -1;
                    else
                        flag_alive(num_w) = false;
                        flag_wait(num_w) = false;
                    end
                end
                
                if flag_started(num_w)
                    if (nb_resub<opt.nb_resub)
                        nb_resub = nb_resub+1;
                        pref_verb = sprintf('Restarting (%i/%i)',nb_resub,opt.nb_resub);
                        if num_w <= opt.max_queued
                            time_reset = clock;
                            save(file_worker_reset{num_w},'time_reset')
                        end
                    elseif ~flag_dead(num_w)
                        pref_verb = sprintf('Marked as dead:');
                        if num_w <= opt.max_queued
                            time_reset = clock;
                            save(file_worker_reset{num_w},'time_reset')
                        end
                        flag_dead(num_w) = true;
                    else
                        pref_verb = '';
                    end
                    
                else
                    pref_verb = 'Starting';
                end
                if (opt.flag_verbose)&&~isempty(pref_verb)
                    if num_w <= opt.max_queued
                        fprintf('%s %s worker number %i...\n',datestr(clock),pref_verb,num_w)
                    elseif num_w == opt.max_queued+1
                        fprintf('%s %s the pipeline manager...\n',datestr(clock),pref_verb);
                    elseif num_w == opt.max_queued+2
                        fprintf('%s %s the garbage collector...\n',datestr(clock),pref_verb);
                    end
                end
                if flag_failed
                    flag_started(num_w) = true;
                end
            end
        end
        if opt.flag_verbose&&strcmp(opt.mode_pipeline_manager,'session')
            nb_chars_logs = psom_pipeline_visu(path_logs,'monitor',nb_chars_logs);      
        end
        
        %% If all workers are dead and no resubmissions are allowed, kill the pipeline
        if psom_exist(file_pipe_running)&&~any(flag_wait(1:opt.max_queued))&&~any(flag_alive(1:opt.max_queued))&&(nb_resub>=opt.nb_resub)&&~psom_exist(file_end)
            if opt.flag_verbose
                 fprintf('%s No active workers left. I am interrupting the pipeline.\n',datestr(now))
            end
            time = clock;
            save(file_end,'time');
        end
        pause(opt.time_between_checks)
        flag_pipe_finished = ~psom_exist(file_pipe_running);
    end

catch
    
    errmsg = lasterror;        
    fprintf('\n\n******************\nSomething went bad ... the PSOM deamon has FAILED !\nThe last error message occured was :\n%s\n',errmsg.message);
    if isfield(errmsg,'stack')
        for num_e = 1:length(errmsg.stack)
            fprintf('File %s at line %i\n',errmsg.stack(num_e).file,errmsg.stack(num_e).line);
        end
    end
    status_pipe = 1;
    return
end

%% Print general info about the pipeline
if opt.flag_verbose&&strcmp(opt.mode_pipeline_manager,'session')
    nb_chars_logs = psom_pipeline_visu(path_logs,'monitor',nb_chars_logs);
end

if psom_exist(file_pipe_running)
    if opt.flag_verbose
        fprintf('The pipeline has not completed, but the max number of allowed submissions was reached or the deamon crashed.\n');
    end
    if ~flag_alive(opt.max_queued+1)
        psom_clean(file_pipe_running,struct('flag_verbose',false))
    end
    status_pipe = 1;
else
    if opt.flag_verbose
        fprintf('Deamon terminated on %s\n',datestr(now));
    end
    status_pipe = 0;
end
