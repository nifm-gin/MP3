
clear all

% T2map.command = '[files_in,files_out,opt] = brick_T2map(char(handles.new_module.files_in),handles.new_module.files_out,handles.new_module.opt)';
% T2map.files_in{1} = char(fullfilename(handles, 1, '.nii'));
% T2map.files_in{2} = char(fullfilename(handles, 1, '.json'));
% %         handles.new_module.files_in_index = handles.pipeline.T2.files_in_index;
% T2map.files_out.filename = '';
% T2map.opt.threshold = 5;
% T2map.opt.flag_test =1;
% T2map.opt.output_filename_ext = '_T2map';
% T2map.SequenceName = 'T2map';


% [files_in,files_out,opt] = brick_T2map(handles.new_module.files_in,handles.new_module.files_out,handles.new_module.opt);

database = load('database.mat');
database = database.database;
%% Subject 1
% files_in.subject1.fmri{1} = '/demo_niak/func_motor_subject1.mnc';
% files_in.subject1.fmri{2} = '/demo_niak/func_rest_subject1.mnc';
% files_in.subject1.transf = '/demo_niak/transf_subject1.xfm';

pipeline = struct();
% Get the list of subjects from files_in
list_subject = database(database.SequenceName == 'T2map_MSME',:);
for num_s = 1:size(list_subject,1)
      %% smooth brick
    job_in{1} = ([char(list_subject.Path(num_s)) char(list_subject.Filename(num_s)) '.nii']);
    job_in{2} = ([char(list_subject.Path(num_s)) char(list_subject.Filename(num_s)) '.json']);
   
     job_out.filename = '';
    % Force a specific folder organization for outputs
    opt.smooth.folder_out = [pwd filesep 'logs'];
    opt.smooth.Type = 'gaussian';
    opt.smooth.Sigma = 1;
    opt.smooth.Hsize = 3;
    opt.smooth.flag_test =1;
    opt.smooth.output_filename_ext = '_smoothed';
    % Give a name to the jobs
    job_name = ['smooth_' char(list_subject.Patient(num_s)) '_' char(list_subject.Tp(num_s))];
    % The name of the employed brick
    brick = 'brick_smooth';
    % Add the job to the pipeline
    pipeline =  psom_add_job(pipeline, job_name,brick,job_in,job_out,opt.smooth);
    
    
    % The outputs of this brick are just intermediate outputs :
    % clean these up as soon as possible
%     pipeline = psom_add_clean(pipeline, [job_name '_clean'],pipeline.(job_name).files_out);

    
    %% T2map brick
    % Plug the T2map input files of the subjects in the job
    job_in{1} = pipeline.(job_name).files_out.filename{1} ;
    job_in{2} = pipeline.(job_name).files_out.filename{2} ;
    % Use the default output name
    job_out.filename = '';
    % Force a specific folder organization for outputs
    opt.T2map.folder_out = [pwd filesep 'logs'];
    opt.T2map.threshold = 5;
    opt.T2map.flag_test =1;
    opt.T2map.output_filename_ext = '_T2map';
    % Give a name to the jobs
    job_name = ['T2map_' char(list_subject.Patient(num_s)) '_' char(list_subject.Tp(num_s))];
    % The name of the employed brick
    brick = 'brick_T2map';
    % Add the job to the pipeline
    pipeline =  psom_add_job(pipeline, job_name,brick,job_in,job_out,opt.T2map);
    
    
  
end
opt_pipeline.path_logs = [pwd filesep 'logs' filesep];  % where to store the log files
psom_run_pipeline(pipeline,opt_pipeline);


for num_s = 1:size(list_subject,1)
     %% T2map brick
    % Plug the T2map input files of the subjects in the job
    job_in{1} =([char(list_subject.Path(num_s)) char(list_subject.Filename(num_s)) '.nii']);
    job_in{2} = ([char(list_subject.Path(num_s)) char(list_subject.Filename(num_s)) '.json']);
    % Use the default output name
    job_out.filename = '';
    % Force a specific folder organization for outputs
    opt.T2map.folder_out = '';
    opt.T2map.threshold = 5;
    opt.T2map.flag_test =1;
    opt.T2map.output_filename_ext = '_T2map';
    % Give a name to the jobs
    job_name = ['T2map_' char(list_subject.Patient(num_s)) '_' char(list_subject.Tp(num_s))];
    % The name of the employed brick
    brick = 'brick_T2map';
    % Add the job to the pipeline
    pipeline =  psom_add_job(pipeline, job_name,brick,job_in,job_out,opt.T2map);
    % The outputs of this brick are just intermediate outputs :
    % clean these up as soon as possible
%     pipeline = psom_add_clean(pipeline, [job_name '_clean'],pipeline.(job_name).files_out);
end

    psom_visu_dependencies(pipeline);


% test the brick
% opt.T2map.flag_test =0;
% [files_in,files_out,opt] = brick_T2map(job_in,job_out,opt.T2map);
% opt.smooth.flag_test  =0;
% [files_in,files_out,opt] = brick_smooth(job_in,job_out,opt.smooth);
% Set up the options to run the pipeline
opt.path_logs = [pwd filesep 'logs' filesep];  % where to store the log files
psom_run_pipeline(pipeline,opt);
