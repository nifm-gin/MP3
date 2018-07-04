function [files_in,files_out,opt] = Module_Texture_matlab(files_in,files_out,opt)

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
    
    module_option(:,3)   = {'Local_entropy', 'Yes'};
    module_option(:,4)   = {'Output_filename_ext_Entropy','_Entropy'};
    module_option(:,5)   = {'Local_Entropy_Patch_size',9};
    
    module_option(:,6)   = {'Range_Value', 'Yes'};
    module_option(:,7)   = {'Output_filename_ext_Range','_Range'};
    module_option(:,8)   = {'Range_Value_Patch_size',3};
    
    module_option(:,9)   = {'Standard_Deviation', 'Yes'};
    module_option(:,10)   = {'Output_filename_ext_Std','_Std'};
    module_option(:,11)   = {'Standard_Deviation_Patch_size',3};
    
    module_option(:,12)   = {'OutputSequenceName','Extension'};
    module_option(:,13)   = {'RefInput',1};
    module_option(:,14)   = {'InputToReshape',1};
    module_option(:,15)   = {'Table_in', table()};
    module_option(:,16)   = {'Table_out', table()};
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
    'Serveral texture analysis from a grayscale 2D-3D image can be computed within a patch'
    ' --> Local entropy (matlab function called entropyfilt'
    '           default patch is 9x9 voxels (Odd value only : 5, 7, 9 ,...)'
    ' --> range value (maximum value - minimum value; matlab function called; rangefilt'
    '           default patch is 3x3 voxels (Odd value only : 3, 5, 7,...)'
    ' --> standard deviation value (matlab function called stdfilt'
    '           default patch is 3x3 voxels (Odd value only : 3, 5, 7,...)'
    ''
    'User can change the size of the patch used to calculate each map'
  
    }'};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    
    user_parameter(:,4)   = {'   .Local entropy map','cell',{'No', 'Yes'},'Local_entropy','','',...
        'J = ENTROPYFILT(I) returns the array J, where each output pixel contains the entropy value of the 9-by-9 neighborhood around the corresponding pixel in the input image I. I can have any dimension'};  
    user_parameter(:,5)   = {'   .Output filename extension','char','_Range','Output_filename_ext_Entropy','', '',''};
    user_parameter(:,6)   = {'   .Patch size','numeric',9,'Local_Entropy_Patch_size','', '',''};  
    
    user_parameter(:,7)   = {'   .Rang Value map','cell',{'No', 'Yes'},'Range_Value','','',...
        'J = RANGEFILT(I) returns the array J, where each output pixel contains the range value (maximum value - minimum value) of the 3-by-3 neighborhood around the corresponding pixel in the input image I'};    
    user_parameter(:,8)   = {'   .Output filename extension','char','_Std','Output_filename_ext_Range','', '',''};
    user_parameter(:,9)   = {'   .Patch size','numeric',3,'Range_Value_Patch_size','', '',''};
    
    user_parameter(:,10)   = {'   .Standard Deviation map','cell',{'No', 'Yes'},'Standard_Deviation','','',...
        'J = STDFILT(I) returns the array J, where each output pixel contains the standard deviation value of the 3-by-3 neighborhood around the corresponding pixel in the input image I'};
    user_parameter(:,11)   = {'   .Output filename extension','char','_Std','Output_filename_ext_Std','', '',''};
    user_parameter(:,12)   = {'   .Patch size','numeric',3,'Standard_Deviation_Patch_size','', '',''};
    
    
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
%opt.NameOutFiles = {opt.Output_filename_ext_Entropy, opt.Output_filename_ext_Range, opt.Output_filename_ext_Std};
opt.NameOutFiles = {};
if strcmp(opt.Local_entropy,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Entropy];
end
if strcmp(opt.Range_Value,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Range];
end
if strcmp(opt.Standard_Deviation,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Std];
end


if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in;
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), opt.NameOutFiles{i}]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
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
rescaled_Im = nan(size(Im));
for i=1:size(Im,3)
    rescaled_Im(:,:,i) = imadjust(origCData(:,:,i), clim, [0, 1]); % imadjust(im(:,:,i),[nanmin(im(:)) nanmax(im(:))]);
end

info = niftiinfo(files_in.In1{1});

if strcmp(opt.Local_entropy,'Yes')
    if mod(opt.Local_Entropy_Patch_size,2)
        Entropy_map = entropyfilt(rescaled_Im, true(opt.Local_Entropy_Patch_size));
    else
        Entropy_map = entropyfilt(rescaled_Im, true(opt.Local_Entropy_Patch_size)+1);
    end
     % transform the OutputImages matrix in order to match to the nii header of the
    % first input (rotation/translation)
  
    map_reoriented = write_volume(Entropy_map, nifti_header, 'axial');

    info2 = info;
    info2.Filename = files_out.In1{1};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(map_reoriented);
    info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
    info2.ImageSize = size(map_reoriented);
    info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];

    niftiwrite(map_reoriented, files_out.In1{1}, info2)
    
   % so far copy the .json file of the first input
    copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
    
end
if strcmp(opt.Range_Value,'Yes')
    if mod(opt.Range_Value_Patch_size,2)
        Range_map =rangefilt(rescaled_Im, true(opt.Range_Value_Patch_size));
    else
        Range_map =rangefilt(rescaled_Im, true(opt.Range_Value_Patch_size)+1);
    end
    map_reoriented = write_volume(Range_map, nifti_header, 'axial');
    
    info2 = info;
    info2.Filename = files_out.In1{2};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(map_reoriented);
    info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
    info2.ImageSize = size(map_reoriented);
    info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];

    niftiwrite(map_reoriented, files_out.In1{2}, info2)
    
   % so far copy the .json file of the first input
    copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{2}, '.nii', '.json'))
    
end
if strcmp(opt.Standard_Deviation,'Yes')
    if mod(opt.Standard_Deviation_Patch_size,2)
        Std_map = stdfilt(rescaled_Im, true(opt.Standard_Deviation_Patch_size));
    else
        Std_map =stdfilt(rescaled_Im, true(opt.Standard_Deviation_Patch_size)+1);
    end
     map_reoriented = write_volume(Std_map, nifti_header, 'axial');
    
    info2 = info;
    info2.Filename = files_out.In1{3};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(map_reoriented);
    info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
    info2.ImageSize = size(map_reoriented);
    info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];

    niftiwrite(map_reoriented, files_out.In1{3}, info2)
    
   % so far copy the .json file of the first input
    copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{3}, '.nii', '.json'))
end

%info = niftiinfo(files_in.In1{1});
% for i=1:length(opt.NameOutFiles)
%     
%     % transform the OutputImages matrix in order to match to the nii header of the
%     % first input (rotation/translation)
%     
%     if ~exist('OutputImages_reoriented', 'var')
%         map_reoriented = write_volume(eval(opt.NameOutFiles{:,i}), nifti_header, 'axial');
%     end
% 
%     info2 = info;
%     info2.Filename = files_out.In1{i};
%     info2.Filemoddate = char(datetime('now'));
%     info2.Datatype = class(map_reoriented);
%     info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
%     info2.ImageSize = size(map_reoriented);
%     info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];
% 
%     niftiwrite(map_reoriented, files_out.In1{i}, info2)
%     
%    % so far copy the .json file of the first input
%     copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{i}, '.nii', '.json'))
% end




