function [pipe,opt_pipe] = psom_test_pipe_tutorial(path_test,opt)
% The test pipeline from the PSOM tutorial 
% http://psom.simexp-lab.org/how_to_use_psom.html
%
% [pipe,opt_pipe] = psom_test_simple_pipe(path_test,opt)
%
% PATH_TEST (string, default current path) where to run the test.
% OPT (structure) any option passed to PSOM will do. In addition the 
%   FLAG_TEST (boolean, default false) if FLAG_TEST is on, the pipeline
%     is generated but not executed.
% PIPE (structure) the pipeline.
% OPT_PIPE (structure) the options to run the pipeline.
%
% Any additional option passed to OPT will be passed in OPT_PIPE. 
%
% Copyright (c) Pierre Bellec, 
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2015-2017.
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information the LICENSE file.
% Keywords : pipeline, PSOM, test

%% Set up default options
pipel = struct;

if nargin < 2
    opt = struct;
end

list_opt = { 'flag_test' };
list_def = { false       };
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

%% Build the pipeline

% Job "sample" :    No input, generate a random vector a
command = 'a = randn([opt.nb_samps 1]); save(files_out,''a'')';
pipe.sample.command      = command;
pipe.sample.files_out    = [path_test 'sample.mat'];
pipe.sample.opt.nb_samps = 10;

% Job "quadratic" : Compute a.^2 and save the results
command = 'load(files_in); b = a.^2; save(files_out,''b'')';
pipe.quadratic.command   = command;
pipe.quadratic.files_in  = pipe.sample.files_out;
pipe.quadratic.files_out = [path_test 'quadratic.mat']; 

% Adding a job "cubic" : Compute a.^3 and save the results
command = 'load(files_in); c = a.^3; save(files_out,''c'')';
pipe.cubic.command       = command;
pipe.cubic.files_in      = pipe.sample.files_out;
pipe.cubic.files_out     = [path_test 'cubic.mat']; 

% Adding a job "sum" : Compute a.^2+a.^3 and save the results
command = 'load(files_in{1}); load(files_in{2}); d = b+c, save(files_out,''d'')';
pipe.sum.command       = command;
pipe.sum.files_in{1}   = pipe.quadratic.files_out;
pipe.sum.files_in{2}   = pipe.cubic.files_out;
pipe.sum.files_out     = [path_test 'sum.mat'];

%% Adding a job to clean up the outputs of the job "sample"
pipe.cleanup.command     = 'delete(files_clean)';
pipe.cleanup.files_clean = pipe.sample.files_out;
        
%% Run the pipe
if ~opt.flag_test
    psom_run_pipeline(pipe,opt_pipe);
end