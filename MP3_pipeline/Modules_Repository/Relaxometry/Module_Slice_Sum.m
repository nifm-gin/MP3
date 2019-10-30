function [files_in,files_out,opt] = Module_MGEFIDSE_Sum(files_in,files_out,opt)


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
    module_option(:,3)   = {'output_filename_ext','Summed'};
    module_option(:,4)   = {'OutputSequenceName','AllName'};
    module_option(:,5)   = {'first_slice', 3};
    module_option(:,6)   = {'number_of_slices',4};
    module_option(:,7)   = {'echo_to_trash','end'};
    module_option(:,8)   = {'OutputResolution','128'};
    module_option(:,9)   = {'RefInput',1};
    module_option(:,10)   = {'InputToReshape',1};
    module_option(:,11)  = {'Table_in', table()};
    module_option(:,12)  = {'Table_out', table()};
    module_option(:,13)  = {'NbToAdd', 1};
    module_option(:,14)  = {'MakeAverage', 'No'};
    opt.Module_settings  = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
     user_parameter(:,1)   = {'Description','Text','','','','',...
        {'Module description :',...
        'This module compute a MGE2D scan from an MGE3D scan. The idea is to decrease the spatial resolution of the MGE scan, by adding the signal of each voxel, in ordre to decrease the macroscopic inhomogeneity'}'};
user_parameter(:,2)   = {'Select a scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
user_parameter(:,3)   = {'Parameters','','','','', '',''};
user_parameter(:,4)   = {'   .Output filename extension','char','MGE2Dfrom3D','output_filename_ext','','',...
    {'Specify the string to be added to the filename input.'
    'Default filename extension is ''MGE2Dfrom3D''.'}'};
user_parameter(:,5)   = {'   .First slice','numeric','','first_slice','', '',''};
user_parameter(:,6)   = {'   .Number of slices','numeric','','number_of_slices','', '',''};
user_parameter(:,7)   = {'   .Echo to trash','char','','echo_to_trash','', '',''};
user_parameter(:,8)   = {'   .Output Resolution','cell',{'64', '128', '256'},'OutputResolution','', '',''};
user_parameter(:,6)   = {'   .Number of slices to add','numeric','','NbToAdd','', '',''};
user_parameter(:,8)   = {'   .Compute average','cell',{'No', 'Yes'},'MakeAverage','', '',''};


VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'MGE2Dfrom3D'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename_ext};




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
    error('Module_MGE2Dfrom3D:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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


input(1).nifti_header = spm_vol(files_in.In1{1});

J = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
% MGE3D = read_volume(input(1).nifti_header, input(1).nifti_header, 0);
% 
% 
% output.nifti_header = input(1).nifti_header;
Volume = N;
% Empty memory
clear data

first_slice =opt.first_slice;
NbToAdd = opt.NbToAdd;
nbr_of_slice = NbToAdd+1;%opt.number_of_slices;
echo_to_trash = opt.echo_to_trash;

MakeAverage = strcmp(opt.MakeAverage, 'Yes');


% Some echoes need to be remove --> human data
if ~strcmp(echo_to_trash, 'end')
    echo_to_trash = str2double(echo_to_trash);
    Volume = Volume(:,:,:,1:echo_to_trash);
    J.EchoTime.value = J.EchoTime.value(1:echo_to_trash);
end


[~, ~, depths3d, echos3d]=size(Volume);
%nbr_of_final_slice = fix((depths3d-first_slice+1)/nbr_of_slice);
nbr_of_final_slice=1;
if nbr_of_final_slice == 0
    return
end
% temp2 = zeros([size(Volume, 1) size(Volume, 2) nbr_of_final_slice size(Volume, 4)]);
S1 = Volume(:,:, first_slice, :);
for i = 1:NbToAdd
    for k=1:echos3d
%         temp=Volume(:,:, first_slice+(i*nbr_of_slice)-nbr_of_slice:first_slice+(i*nbr_of_slice)-1, k);
        temp=S1(:,:,1,k)+ Volume(:,:, first_slice+i,k);
%         temp=sum(abs(temp),3);
        S1(:,:,1,k)=temp;
    end
end
if MakeAverage
    S1 = S1 ./nbr_of_slice;
end

% x = 1:rows3d;
% y = 1:cols3d;
%MGE2Dfrom3D = NaN(size(MGE3D));
if size(N,1) ~= size(N,2)
    error('The input image hasn''t a squared resolution. This case is not implemented')
end
facteur = size(N,1) / str2double(opt.OutputResolution);
if floor(facteur)~=facteur
    error('The initial resolution is not a multiple of the wanted output one.')
end

mask = ones(facteur);
Test = zeros([size(S1,1)/facteur,size(S1,2)/facteur,size(S1,3), size(S1,4)]);
Test2 = zeros(size(S1));
for i=1:size(S1,3)
    for j=1:size(S1,4)
        Test2(:,:,i,j) = conv2(S1(:,:,i,j), mask, 'same');
        Test(:,:,i,j) = Test2(1:facteur:end, 1:facteur:end,i,j);
    end
end
% 

MGE2Dfrom3D = Test;


J.SliceThickness.value = J.SliceThickness.value*nbr_of_slice;
J.ImageInAcquisition.value = nbr_of_final_slice*length(J.EchoTime.value);



OutputImages = MGE2Dfrom3D;
OutputImages(isinf(OutputImages)) = -1;
OutputImages_reoriented = OutputImages;

% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
% info2.PixelDimensions = zeros(1,4);
% info2.PixelDimensions(1:2) = info.PixelDimensions(1:2);
if size(OutputImages_reoriented, 3) > 1
    info2.PixelDimensions(3) = info.PixelDimensions(3)*nbr_of_slice;
end
info2.Transform.T(1,1) = info2.Transform.T(1,1)*facteur;
info2.Transform.T(2,2) = info2.Transform.T(2,2)*facteur;
info2.Transform.T(4,3) = info2.Transform.T(4,3)+first_slice*info.Transform.T(3,3);
info2.Transform.T(3,3) = info2.Transform.T(3,3)*nbr_of_slice;
info2.ImageSize = size(OutputImages_reoriented);
info2.Description = 'Module Slice Sum';

% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% % so far copy the .json file of the first input
% copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))


%% Json Processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

% 

