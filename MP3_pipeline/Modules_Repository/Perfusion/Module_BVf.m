function [files_in,files_out,opt] = Module_BVf(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
  
     % define every option needed to run this module
      %%   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'trash_below',0};
    module_option(:,4)   = {'trash_after',Inf};
    module_option(:,5)   = {'output_name','BVf'};
    module_option(:,6)   = {'B0',4.7};
    module_option(:,7)   = {'Gamma',2.67502};
    module_option(:,8)   = {'deltaxi_CA',0.1867};
    module_option(:,9)   = {'OutputSequenceName','AllName'};
    module_option(:,10)  = {'RefInput',1};
    module_option(:,11)  = {'InputToReshape',1};
    module_option(:,12)  = {'Table_in', table()};
    module_option(:,13)  = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
      user_parameter(:,1)   = {'Description','Text','','','', '',...
        {'Description of the module :',...
        'this module compute a BVf map from a deltaR2* (see module deltaR2*)',...
        'Output map in %', ...
        'If you use this function, please refere to this article : Tropes et al. MRM 2015'}};
user_parameter(:,2)   = {'Select a DeltaR2* scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
user_parameter(:,3)   = {'Parameters','','','','', '',''};
user_parameter(:,4)   = {'   .Output filename','char','BVf','output_name','','',...
    {'Specify the string to be added to the filename input.'
    'Default filename extension is ''BVf''.'}'};
user_parameter(:,5)   = {'   .B0','numeric', '','B0','','',...
    {'Please entre the magnetic field used'}'};
user_parameter(:,6)   = {'   .Gamma (10^8)','numeric','','Gamma','','',...
    {''}'};
user_parameter(:,7)   = {'   .deltaxi HBC (10^-6) (cgs units)','numeric', '','deltaxi_CA','','',...
    {''}'};


VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'BVf'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_name};


if isempty(files_out)
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_name));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_name]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_BVf:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs


% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii');  
% end
% 
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


%% load input Nii file
N = niftiread(files_in.In1{1});
% load nifti info
info = niftiinfo(files_in.In1{1});

%% load input JSON file
%J = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));



B0=opt.B0;
gamma=opt.Gamma; gamma = gamma*10^8;
deltaxi_CA=opt.deltaxi_CA; deltaxi_CA = deltaxi_CA*10^-6;

deltaomega=2*pi*gamma*deltaxi_CA*B0;


% save imformation
BVf_map=N;

% Empty memory
clear data

for i = 1:size(N, 3)
    BVf_map(:,:,i,1)=(3/(2*deltaomega)*BVf_map(:,:,i,1))*100000;
end


OutputImages = BVf_map;
OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 50) = NaN;
% OutputImages(isnan(OutputImages)) = -1;

% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);
%info2.Description = [info.Description, 'Modified by T2star_map Module'];

% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);

% % so far copy the .json file of the first input
% copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))


%% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)




% 
