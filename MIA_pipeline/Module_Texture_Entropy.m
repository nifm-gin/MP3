function [files_in,files_out,opt] = Module_Texture_Entropy(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename_ext','_Entropy'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'Type','gaussian'};
    module_option(:,6)   = {'Patch_size',9};
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
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
    'This module return a Local entropy of grayscale 2D-3D image.  '
    'The matlab function used is the entropyflit'
    ''
    'User can change the size of the patch used to calculate the Entropy'
    '    --> the default value is 9x9 voxels (Odd value only : 5, 7, 9 ,...)'
    }'};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Entropy','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Patch size','numeric',9,'Patch_size','', '',''};
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
    error('Texture_Entropy:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

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

nifti_header = spm_vol(files_in.In1{1});
Im = read_volume(nifti_header, nifti_header, 0, 'axial');

Im(Im == 0) = NaN;
%clim = [1040.30956306678 1111.37427864798];
%clim = [1050 1100];

clim = [nanmin(Im(:)) nanmax(Im(:))];

origCData = (Im-clim(1))/diff(clim);
clim = [ 0 1];
OutputImage = nan(size(Im));
for i=1:size(Im,3)
    OutputImage(:,:,i) = imadjust(origCData(:,:,i), clim, [0, 1]); % imadjust(im(:,:,i),[nanmin(im(:)) nanmax(im(:))]); 
end

if mod(opt.Patch_size,2)
    Entropy_map = entropyfilt(OutputImage, true(opt.Patch_size));
     %% test other maps
    %Entropy_map =rangefilt(OutputImage, true(opt.Patch_size));
   %Entropy_map = stdfilt(OutputImage, true(opt.Patch_size));
else
    Entropy_map = entropyfilt(OutputImage, true(opt.Patch_size)+1);
    %% test other maps
    %Entropy_map =rangefilt(OutputImage, true(opt.Patch_size)+1);
    %Entropy_map =stdfilt(OutputImage, true(opt.Patch_size)+1);

end

% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)

if ~exist('OutputImages_reoriented', 'var')
    Entropy_map_reoriented = write_volume(Entropy_map, nifti_header, 'axial');
end



% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.In1{1});


info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(Entropy_map_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(Entropy_map_reoriented)));
info2.ImageSize = size(Entropy_map_reoriented);
info2.DisplayIntensityRange = [nanmin(Entropy_map_reoriented(:))  nanmax(Entropy_map_reoriented(:))];
%info2.Description = [info.Description, 'Modified by T2map Module'];


% save the new .nii file
niftiwrite(Entropy_map_reoriented, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))

