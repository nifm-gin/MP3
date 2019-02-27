function [pipel,opt_pipe] = psom_test_hist(path_test,opt)
% A test pipeline with chains of jobs accumulating histograms
%
% [pipe,opt_pipe] = psom_test_hist(path_test,opt)
%
% Parallel computation of the histogram of the square root of the absolute tanh 
% of a normal distribution.
%
% PATH_TEST (string, default current path) where to run the test.
% OPT (structure) any option passed to PSOM will do. In addition the 
%   following options are available:
%   TIME (scalar, default 3) the minimum time (in sec) for the job. The job will 
%       keep looping the estimation until this time has passed. 
%   MINMAX (vector 2x1, default [-5 5]) the min/max of the histogram.
%   SIZEBIN (scalar, default 0.001) the size of the bins of the histogram.
%   NB_SAMP (integer, default 10^7) the number of samples to build one histogram.
%   NB_JOBS (integer, default 100) the number of jobs.
%   FLAG_TEST (boolean, default false) if FLAG_TEST is on, the pipeline
%     is generated but not executed.
%   PATH_LOGS is forced to [path_test filesep 'logs']
% PIPE (structure) the pipeline.
% OPT_PIPE (structure) the options to run the pipeline.
%
% Copyright (c) Pierre Bellec, 
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2015.
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information the LICENSE file.
% Keywords : pipeline, PSOM, test

%% Set up default options
pipel = struct;

if nargin < 2
    opt = struct;
end

list_opt = { 'time' , 'nb_samp' , 'sizebin' , 'nb_jobs' , 'flag_test' };
list_def = { 3      , 10^7      , 0.001     , 100       , false       };
opt = psom_struct_defaults(opt,list_opt,list_def,false);

if (nargin < 1)||isempty(path_test)
    path_test = pwd;
end
if ~strcmp(path_test(end),filesep)
    path_test = [path_test filesep];
end

opt.path_logs = [path_test 'logs'];

%% The options for PSOM
opt_pipe = rmfield(opt,list_opt);

%% Options for job
optj.nb_samp = opt.nb_samp;
optj.sizebin = opt.sizebin;
optj.time = opt.time;

%% Build the pipeline
for jj = 1:opt.nb_jobs
    job_name = sprintf('samp%i',jj);   
    pipel.(job_name).opt = optj;    
    pipel.(job_name).command = sprintf([ ...
                ' edges = 0:opt.sizebin:1;' ...
                ' N = zeros(length(edges),1);' ...
                ' t1 = tic;' ...
                ' nb_samp_tot = 0;' ...
                ' while toc(t1)<opt.time;' ...
                ' data = sqrt(abs(tanh(randn(opt.nb_samp,1))));' ...
                ' N = N + histc(data,edges);' ...
                ' nb_samp_tot = nb_samp_tot + opt.nb_samp;' ...
                ' end;' ...
                ' N = N/nb_samp_tot;' ...
                ' save(files_out,''N'',''nb_samp_tot'');' ...
                ]);      
    samp_name{jj} = sprintf('%ssamp%i.mat',path_test,jj);   
    pipel.(job_name).files_out = samp_name{jj};
end

% now normalize the histograms
pipel.hist.files_in = samp_name;
pipel.hist.files_clean = pipel.hist.files_in;
pipel.hist.files_out = sprintf('%shistogram_gaussian.mat',path_test);
pipel.hist.opt = optj;
pipel.hist.opt.nb_jobs = opt.nb_jobs;
pipel.hist.command = sprintf([ ...
    ' edges = 0:opt.sizebin:1;' ...
    ' N = zeros(length(edges),1);' ...
    ' nb_samp_tot = 0;' ...
    ' for num_c = 1:length(files_in);' ...
    '  data = load(files_in{num_c});' ...
    '  N = N + data.N;' ...
    '  nb_samp_tot = nb_samp_tot + data.nb_samp_tot;' ...
    ' end;' ...
    ' H = N / (opt.sizebin * opt.nb_jobs );' ...
    ' save(files_out,''H'',''edges'',''nb_samp_tot'');' ...
    ' psom_clean(files_in);' ...
    ]);
        
%% Run the pipeline
if ~opt.flag_test
    psom_run_pipeline(pipel,opt_pipe);
end