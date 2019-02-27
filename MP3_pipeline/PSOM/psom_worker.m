function status_pipe = psom_worker(path_worker,path_logs,num_worker,time_pipeline)
% Execute jobs.
%
% status = psom_worker( path_worker , path_logs , num_worker )
%
% PATH_WORKER (string) The name of a path where all logs will be saved.
% PATH_LOGS (string)
% NUM_WORKER (integer, default 1)
% TIME_PIPELINE
%
% See licensing information in the code.

% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008-2010.
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2010-2015.
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
if ~exist('path_worker','var')
    error('SYNTAX: status_pipe = psom_worker(path_worker,flag). Type ''help psom_worker'' for more info.')
end

if ~ischar(path_worker)
    error('PATH_WORKER should be a string (name of a path)')
end

if ~strcmp(path_worker(end),filesep)
    path_worker = [path_worker filesep];
end

if ~strcmp(path_logs(end),filesep)
    path_logs = [path_logs filesep];
end

if nargin < 3
    error('Please specify NUM_WORKER')
end

if nargin < 4
    error('Please specify TIME_PIPELINE')
end
     
%% Create folder for worker
if ~psom_exist(path_worker)
    psom_mkdir(path_worker);
end  

%% Generating file names
file_heartbeat = [path_worker filesep 'heartbeat.mat'];
file_kill      = [path_worker filesep 'worker.kill'];
file_end       = [path_worker filesep 'worker.end'];
file_lock      = [path_logs filesep 'PIPE.lock'];
file_time         = [path_logs 'PIPE_time.mat'];

%% Check that the time of the pipeline matches the record
%% This check is done to ensure a new pipeline has not been started
%% since the manager was started
logs_time = load(file_time);
if ~strcmp(time_pipeline,logs_time.time_pipeline)
    fprintf('The time of the pipeline does not match the logs. I am quitting.')
    exit
end

%% Clean-up old submissions
list_ready = dir([path_worker '*.ready']);
list_ready = { list_ready.name };
psom_clean(list_ready);
if ~isempty(list_ready)
    for num_r = 1:length(list_ready)
        [tmp,base_spawn] = fileparts(list_ready{num_r});
        file_spawn = [path_worker base_spawn '.mat'];
        if psom_exist(file_spawn)
            psom_clean(file_spawn);
        end
    end
end

%% Start a heartbeat
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

% a try/catch block is used to crash gracefully if the user is
% interrupting the pipeline of if an error occurs
try    
    %% Initialize and start the execution loop
    test_loop = true;
    num_job = 0;
    flag_any_fail = false;
    time_scheduled = struct();
    time_running = struct();
    list_jobs = {};
    pipeline = struct;
    flag_end = false;
    while test_loop

        %% Check for new spawns
        list_ready = dir([path_worker '*.ready']);
        list_ready = { list_ready.name };
        if ~isempty(list_ready)
            for num_r = 1:length(list_ready)
                [tmp,base_spawn] = fileparts(list_ready{num_r});
                file_spawn = [path_worker base_spawn '.mat'];
                if ~psom_exist(file_spawn)
                    warning('I could not find %s for spawning',file_spawn)
                    continue
                end
                spawn = load(file_spawn);
                list_new_jobs = fieldnames(spawn);
                % Record scheduling time
                for nn = 1:length(list_new_jobs)
                    name_job = list_new_jobs{nn};
                    time_scheduled.(name_job) = clock;
                end
                list_jobs = [ list_jobs ; list_new_jobs ];
                pipeline = psom_merge_pipeline(pipeline,spawn);
                psom_clean({file_spawn,[path_worker list_ready{num_r}]});
            end
        end
        
        %% If there are jobs to run
        if num_job < length(list_jobs)
            num_job = num_job + 1;
            name_job = list_jobs{num_job};
            
            %% Record running time
            time_running.(name_job) = clock;
                
            %% Execute the job in a "shelled" environment
            flag_failed = psom_run_job(pipeline.(name_job),path_worker,name_job);    
            
            %% Update the news feed
            time_job = clock;
            if flag_failed
                save([path_worker name_job '.failed'],'time_job');
            else
                save([path_worker name_job '.finished'],'time_job');
            end
            
            %% Update profile info
            file_prof_job = [path_worker name_job '_profile.mat'];
            new_prof = struct();
            new_prof.time_scheduled = time_scheduled.(name_job);
            new_prof.time_running = time_running.(name_job);
            new_prof.worker = num_worker;
            save(file_prof_job,'-append','-struct','new_prof');
        end 
        
        test_loop = psom_exist(file_lock)&&(~flag_end||(num_job<length(list_jobs)));
        flag_end = psom_exist(file_end);
        if flag_end
            fprintf('%s - The manager has requested to end the work asap!\n',datestr(clock));
        end
        
        if (num_job == length(list_jobs))&&test_loop
            pause(0.1); 
        end
    end % While there are jobs to do
    
    %% Return a 1 status if any job has failed
    status_pipe = double(flag_any_fail);
    
catch
    
    errmsg = lasterror;        
    fprintf('\n\n******************\nSomething went bad ... the pipeline has FAILED !\nThe last error message occured was :\n%s\n',errmsg.message);
    if isfield(errmsg,'stack')
        for num_e = 1:length(errmsg.stack)
            fprintf('File %s at line %i\n',errmsg.stack(num_e).file,errmsg.stack(num_e).line);
        end
    end
   
    status_pipe = 1;
    return
end

%% SUBFUNCTIONS

%% Read a text file
function str_txt = sub_read_txt(file_name)

hf = fopen(file_name,'r');
if hf == -1
    str_txt = '';
else
    str_txt = fread(hf,Inf,'uint8=>char')';
    fclose(hf);    
end
