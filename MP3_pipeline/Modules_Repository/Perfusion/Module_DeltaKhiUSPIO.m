function [files_in,files_out,opt] = Module_DeltaKhiUSPIO(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
  
     % define every option needed to run this module
      %%   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values    
    module_option(:,1)   = {'output_name','DeltaKhiUSPIO'};
    module_option(:,2)   = {'B0',4.7};
    module_option(:,3)   = {'Gamma',2.67502};
    module_option(:,4)   = {'BVf',4};
    module_option(:,5)   = {'update','Yes'};
    
    module_option(:,6)   = {'OutputSequenceName','AllName'};
    module_option(:,7)   = {'RefInput',1};
    module_option(:,8)   = {'InputToReshape',1};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)  = {'Table_out', table()};
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
        {'Blabla'}};
    user_parameter(:,2)   = {'Select a DeltaR2* scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select an ROI of the brain','1ROI','','',{'SequenceName'}, 'Mandatory',''};
    
    user_parameter(:,4)   = {'Parameters','','','','', '',''};
    user_parameter(:,5)   = {'   .Output filename','char','DeltaKhiUSPIO','output_name','','',...
        {'Specify the string to be added to the filename input.'
        'Default filename extension is ''DeltaKhiUSPIO''.'}'};
    user_parameter(:,6)   = {'   .B0','numeric', '','B0','','',...
        {'Please entre the magnetic field used'}'};
    user_parameter(:,7)   = {'   .Gamma (10^8)','numeric','','Gamma','','',...
        {''}'};
    user_parameter(:,8)   = {'   .Mean BVf (in %)','numeric', '','BVf','','',...
        {'This is the mean BVf in the brain'}'};
    user_parameter(:,9)   = {'   .Update register ?','cell', {'Yes','No'},'update','','',...
        {''}'};


    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_name};


if isempty(files_out)
    opt.Table_out = opt.Table_in(1,:);
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


%% Load
input(1).nifti_header = spm_vol(files_in.In1{1});
if isfield(files_in, 'In2')
    input(2).nifti_header = spm_vol(files_in.In2{1});
end

ref_scan = 1;

Scan = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
Types{1} = class(Scan);

if isfield(files_in, 'In2')
    ROI = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
    Types{2} = class(ROI);
    % if there is a missmatch in the data type we have to convert data type
    % before doing any operation except the multiplication to a scalar  
    CommonType = FindCommonDatatype(Types);       
    Scan = cast(Scan, CommonType);
    ROI = cast(ROI, CommonType);
    
end

if length(size(Scan)) > length(size(ROI))
    DeltaR2star = Scan .* repmat(ROI, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]);    
else
    DeltaR2star =  repmat(Scan, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]) .* ROI;    
end

meanDeltaR2star	= 1e3* nanmean(DeltaR2star(DeltaR2star ~= 0));
DeltaKhiUSPIO	= 3/(4*pi*(opt.Gamma*10^8)*opt.B0)  * (meanDeltaR2star/(1e-2*opt.BVf));


if strcmp(opt.update, 'Yes')
    
    folder  = [opt.folder_out '/../Reports'];
    file    = 'DeltaKhiUSPIO.xls';
    if ~exist(folder,'dir')
        mkdir(folder);
    end
    if ~exist([folder '/' file],'file')
        T   = table({char(opt.Table_in.Patient(1))},{char(opt.Table_in.Tp(1))},DeltaKhiUSPIO,'VariableNames',{'Patient','Tp','DeltaKhiUSPIO'}); 
    else 
        T   = sortrows(readtable([folder '/' file]));
        
        if any(strcmp(T.Patient,char(opt.Table_in.Patient(1))) & strcmp(T.Tp,char(opt.Table_in.Tp(1))))
            flag = find(strcmp(T.Patient,char(opt.Table_in.Patient(1))) & strcmp(T.Tp,char(opt.Table_in.Tp(1))) == 1);
        else
            flag = size(T,1)+1;
        end
        
        T.Patient(flag) = {char(opt.Table_in.Patient(1))};
        T.Tp(flag)      = {char(opt.Table_in.Tp(1))};
        T.DeltaKhiUSPIO(flag) = DeltaKhiUSPIO;
    end
    
    writetable(T, [folder '/' file]);
end


% Create output image
OutputImages    = nan(size(Scan));
for s = 1:size(ROI,3)
    OutputImages(:,:,s) = DeltaKhiUSPIO * ROI(:,:,s);
end

% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.In1{1});
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);

% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);


%% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

% 
