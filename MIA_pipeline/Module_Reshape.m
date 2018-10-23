function [files_in,files_out,opt] = Module_Reshape(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

%     %   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values    
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'output_filename_ext','_Reshape'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'Axes','4'};
    module_option(:,6)   = {'IndexVector', '3 4'};
    module_option(:,7)   = {'RefInput',1};
    module_option(:,8)   = {'InputToReshape',1};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)   = {'Table_out', table()};
    
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
% 
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','', '', {''}};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Reshape','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Axes','cell',{'3','4','5'},'Axes','', '',''};
    user_parameter(:,6)   = {'   .IndexVector','char','3 4','IndexVector','', '',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
%%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%


if isempty(files_out)
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end








%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Resphape:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, ext_nii, files_in);  
% end
% 
% %% Building default output names
% if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
%     opt.folder_out = path_nii;    
% end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get data from a specific axe      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = niftiread(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
[path, name, ext] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];

J = ReadJson(jsonfile);

Informations = whos('N');
axes = str2double(opt.Axes);
index = str2double(opt.IndexVector);

if axes < 3
    error('Input axes must > 2');
elseif (length(index) ~=2) && (axes == 3)
    error('For third axes give index as [min,max]');
end


if axes == 3

    NewIm = N(:,:,index(1):index(2));  
    
    OutputMat = info.Transform.T';
    Movement = [0;0;index(1)-1;1];
    OutputMat(:,4) = InputMat * Movement;
    info.Transform.T =  OutputMat';
    
    
elseif axes == 4
 
    NewIm= zeros( [ size(N,1),size(N,2),size(N,3),length(index)] , Informations.class); 
    for i = 1: length(index)
        NewIm(:,:,:,i) = N(:,:,:,index(i));  
    end        
 
    
elseif axes == 5
   
    NewIm = zeros([size(N,1),size(N,2),size(N,3),size(N,4),length(index)], Informations.class);
    
    for i = 1:  length(index)
        NewIm(:,:,:,:,i) = N(:,:,:,:,index(i));     
    end
        
        

    
end


info.ImageSize = size(NewIm);

info.PixelDimensions = info.PixelDimensions(1:length(size(NewIm)));

info.ImageSize = size(NewIm);


info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
%info2.Description = [info.Description, 'Modified by Smoothing Module'];

niftiwrite(NewIm, files_out.In1{1}, info2)

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)
