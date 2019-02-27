function flag_failed = psom_run_job(job,path_logs,name_job)
% Run a PSOM job. 
%
% SYNTAX: flag_failed = psom_run_job(job,path_logs,name_job)
%
% JOB (string) JOB is a structure defining a PSOM job (with
%   COMMAND, FILES_IN, FILES_OUT, FILES_CLEAN, OPT fields). 
% PATH_LOGS (string, default '') the name of a folder for the logs.
% NAME_JOB (string, default 'manual') a name for the job
% FLAG_FAILED (boolean) FLAG_FAILED is true if the job has failed. This 
%   happens if the command of the job generated an error, or if one of 
%   the output files of the job was not successfully generated.
%
% NOTE 1: When running a job, this function will create a global variable named
%   "gb_psom_name_job". This can be accessed by the command executed by the
%   job. This may be useful for example to build unique temporary file names.
% 
% NOTE 2: If PATH_LOGS is specified, this function will generate a log file and 
% a profile file (with extensions .log and _profile.mat, respectively).
%
% See licensing information in the code.

% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008-2010.
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2010-2015.
% Maintainer : pierre.bellec@criugm.qc.ca
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

global gb_psom_name_job
psom_gb_vars
seed = psom_set_rand_seed();

%% Does not allow figure display and reset to default value at end 
%% of function call
fig_is_visible = get(0, 'defaultFigureVisible');
c = onCleanup(@() set(0, 'defaultFigureVisible', fig_is_visible));
set(0,'defaultFigureVisible','off');

%% name of the job
if nargin < 3
    name_job = 'manual';
elseif ~ischar(name_job)
        error('NAME_JOB should be a string.')
end
gb_psom_name_job = name_job;

%% Logs path
if nargin<2
    path_logs = '';
end

if ~ischar(path_logs)
    error('PATH_LOGS should be a string')
end

if ~isempty(path_logs)&&~strcmp(path_logs(end),filesep)
    path_logs = [path_logs filesep];
end

%% Generate file names
if ~isempty(path_logs)
    file_profile = [path_logs filesep name_job '_profile.mat'];
    file_log = [path_logs filesep name_job '.log'];
end 

%% Upload job info
job = psom_struct_defaults( job , ...
   { 'files_in' , 'files_out' , 'files_clean' , 'command','opt' , 'dep' , 'ispipeline' }, ...
   { {}         , {}          , {}            , NaN      , {}   , {}    , false        });

%% The job starts now !
try 
    %% Print general info about the job
    start_time = clock;
    if ~isempty(path_logs)
        if psom_exist(file_log)
            delete(file_log);
        end
        diary(file_log)
    end
    msg = sprintf('Log of the (%s) job : %s\nStarted on %s\nUser: %s\nhost : %s\nsystem : %s',gb_psom_language,name_job,datestr(clock),gb_psom_user,gb_psom_localhost,gb_psom_OS);
    stars = repmat('*',[1 30]);
    fprintf('\n%s\n%s\n%s\n',stars,msg,stars);

    flag_failed = false;
   
    try
        sub_eval(job.command,job.files_in,job.files_out,job.files_clean,job.opt)
        end_time = clock;
    catch
        end_time = clock;
        flag_failed = true;
        errmsg = lasterror;
        fprintf('\n\nSomething went bad ... the job has FAILED !\nThe last error message occured was :\n%s\n',errmsg.message);
        if isfield(errmsg,'stack')
            for num_e = 1:length(errmsg.stack)
                fprintf('File %s at line %i\n',errmsg.stack(num_e).file,errmsg.stack(num_e).line);
            end
        end
    end
    
    %% Checking outputs
    list_files = psom_files2cell(job.files_out);

    if ~isempty(list_files)
        fprintf('\n\n')
    end
    
    for num_f = 1:length(list_files)
        if ~psom_exist(list_files{num_f})
            if size(list_files{num_f},1)>1
                fprintf('The output file or directory %s (...) has not been generated!\n',list_files{num_f}(1,:));
            else
                fprintf('The output file or directory %s has not been generated!\n',list_files{num_f});
            end
            flag_failed = true;
        else
            if size(list_files{num_f},1)>1
                fprintf('The output file or directory %s (...) was successfully generated!\n',list_files{num_f}(1,:));
            else
                fprintf('The output file or directory %s was successfully generated!\n',list_files{num_f});
            end
        end
    end                   

    %% Verbose an epilogue
    elapsed_time = etime(end_time,start_time);
    if flag_failed
        msg1 = sprintf('%s : The job has FAILED',datestr(clock));
    else
        msg1 = sprintf('%s : The job was successfully completed',datestr(clock));
    end
    msg2 = sprintf('Total time used to process the job : %1.2f sec.',elapsed_time);
    stars = repmat('*',[1 max(size(msg1,2),size(msg2,2))]);
    fprintf('\n%s\n%s\n%s\n%s\n',stars,msg1,msg2,stars);
    
    %% Create a tag file for output status
    if ~isempty(path_logs)        
        %% Create a profile
        save(file_profile,'start_time','end_time','elapsed_time','seed');
        diary off
    end
catch
    end_time = clock;
    flag_failed = true;
    errmsg = lasterror;
    fprintf('\n\n%s\nSomething went bad ... the job has FAILED !\nThe last error message occured was :\n%s\n',stars,errmsg.message);
    if isfield(errmsg,'stack')
        for num_e = 1:length(errmsg.stack)
            fprintf('File %s at line %i\n',errmsg.stack(num_e).file,errmsg.stack(num_e).line);
        end
    end
    if ~isempty(path_logs)
        end_time = clock;
        elapsed_time = etime(end_time,start_time);
        save(file_profile,'start_time','end_time','elapsed_time','seed');
        diary off
    end
end

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

function [] = sub_eval(command,files_in,files_out,files_clean,opt)

eval(command)

