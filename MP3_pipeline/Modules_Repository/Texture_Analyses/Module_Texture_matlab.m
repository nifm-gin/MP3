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
    
    % size of the patch
    module_option(:,3)   = {'Patch_size',9};
    
    %% Compute the local mean accros the patch
    module_option(:,4)   = {'Local_Mean', 'Yes'};
    module_option(:,5)   = {'Output_filename_ext_Mean','_LocalMean'};
    
    %% Compute the local Skew accros the patch
    module_option(:,6)   = {'Local_Skewness', 'Yes'};
    module_option(:,7)   = {'Output_filename_ext_Skew','_LocalSkew'};
    
    %% Compute the local Kurtosis accros the patch
    module_option(:,8)   = {'Local_Kurtosis', 'Yes'};
    module_option(:,9)   = {'Output_filename_ext_Kurtosis','_LocalKurtosis'};
    
    %% Compute the Local entropy accros the patch
    module_option(:,10)   = {'Local_Entropy', 'Yes'};
    module_option(:,11)   = {'Output_filename_ext_Entropy','_Entropy'};
    
    %% Compute the Rang accros the patch
    module_option(:,12)   = {'Local_Range', 'Yes'};
    module_option(:,13)   = {'Output_filename_ext_Range','_LocalRange'};
    
    %% Compute the Standard Deviation accros the patch
    module_option(:,14)   = {'Local_SD', 'Yes'};
    module_option(:,15)   = {'Output_filename_ext_Std','_Std'};
    
    %% General's option of the module
    module_option(:,16)   = {'OutputSequenceName','Extension'};
    module_option(:,17)   = {'RefInput',1};
    module_option(:,18)   = {'InputToReshape',1};
    module_option(:,19)   = {'Table_in', table()};
    module_option(:,20)   = {'Table_out', table()};
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
        '           default patch is 9x9 voxels (Odd value only : 5, 7, 9 ,...)'
        ' --> Local Mean: mean across the patch'
        ' --> Local Skewness: skewness of the patch (matlab function called skewness)'
        ' --> Local Kurtosis: Kurtosis of the patch  (matlab function called Kurtosis)'
        ' --> Local entropy: Entropy of hte patch (matlab function called entropyfilt)'
        ' --> range value: range of the patch (maximum value - minimum value; matlab function called rangefilt)'
        ' --> standard deviation value (matlab function called stdfilt)'
        ''
        }'};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Patch size','numeric',9,'Patch_size','', '',...
        'Default patch is 9x9 voxels (Odd value only : 5, 7, 9 ,...)'};
    
    user_parameter(:,5)   = {'   .Local Mean map','cell',{'No', 'Yes'},'Local_Mean','','',...
        'Mean value across the patch'};
    user_parameter(:,6)   = {'   .Output filename extension','char','_LocalMean','Output_filename_ext_Mean','', '',''};   
    
    user_parameter(:,7)   = {'   .Local Skewness map','cell',{'No', 'Yes'},'Local_Skewness','','',...
        'Local Skewness of the patch'};
    user_parameter(:,8)   = {'   .Output filename extension','char','_LocalSkew','Output_filename_ext_Skew','', '',''}; 
    
    user_parameter(:,9)   = {'   .Local Kurtosis map','cell',{'No', 'Yes'},'Local_Kurtosis','','',...
        'Compute the local Kurtosis of the patch'};
    user_parameter(:,10)   = {'   .Output filename extension','char','_LocalKurtosis','Output_filename_ext_Kurtosis','', '',''};  
   
    user_parameter(:,11)   = {'   .Local entropy map','cell',{'No', 'Yes'},'Local_Entropy','','',...
        'J = ENTROPYFILT(I) returns the array J, where each output pixel contains the entropy value of the patch neighborhood around the corresponding pixel in the input image I. I can have any dimension'};
    user_parameter(:,12)   = {'   .Output filename extension','char','_LocalEntropy','Output_filename_ext_Entropy','', '',''};
  
    user_parameter(:,13)   = {'   .Rang Value map','cell',{'No', 'Yes'},'Local_Range','','',...
        'J = RANGEFILT(I) returns the array J, where each output pixel contains the range value (maximum value - minimum value) of the 3-by-3 neighborhood around the corresponding pixel in the input image I'};
    user_parameter(:,14)   = {'   .Output filename extension','char','_Std','Output_filename_ext_Range','', '',''};
    
    user_parameter(:,15)   = {'   .Standard Deviation map','cell',{'No', 'Yes'},'Local_SD','','',...
        'J = STDFILT(I) returns the array J, where each output pixel contains the standard deviation value of the 3-by-3 neighborhood around the corresponding pixel in the input image I'};
    user_parameter(:,16)   = {'   .Output filename extension','char','_LocalSD','Output_filename_ext_Std','', '',''};
    
    
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
if strcmp(opt.Local_Mean,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Mean];
end
if strcmp(opt.Local_Skewness,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Skew];
end
if strcmp(opt.Local_Kurtosis,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Kurtosis];
end
if strcmp(opt.Local_Entropy,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Entropy];
end
if strcmp(opt.Local_Range,'Yes')
    opt.NameOutFiles = [opt.NameOutFiles, opt.Output_filename_ext_Range];
end
if strcmp(opt.Local_SD,'Yes')
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
% data have to be converted to double before computing any
% maps
if isa(class(Im), 'double') ==0
    origin_class =  class(Im);
    Im = double(Im);
end

% data are nomalized and centered
clim = [nanmin(Im(:)) nanmax(Im(:))];
origCData = (Im-clim(1))/diff(clim);

info = niftiinfo(files_in.In1{1});

maps_to_generate = [istrue(opt.Local_Mean) istrue(opt.Local_Skewness)...
    istrue(opt.Local_Kurtosis) istrue(opt.Local_Entropy) istrue(opt.Local_Range) istrue(opt.Local_SD)];

if mod(opt.Patch_size,2)
    patch_2D =  true(opt.Patch_size);
    patch_in_vector = [opt.Patch_size, opt.Patch_size];
else
    patch_2D =  true(opt.Patch_size+1);
    patch_in_vector = [opt.Patch_size+1, opt.Patch_size+1];
end
files_out_nbr = 0;
for i=1:numel(maps_to_generate)
    if maps_to_generate(i)
        switch i
            % Local_Mean
            case 1
                for j= 1:size(origCData,3)
                    New_texture_map(:,:,j) = colfilt(origCData(:,:,j),patch_in_vector,'sliding',@mean);
                end
                
                % Local_Skewness
            case 2
                for j= 1:size(origCData,3)
                    New_texture_map(:,:,j) = colfilt(origCData(:,:,j),patch_in_vector,'sliding',@skewness);
                end

                % Local_Kurosis
            case 3
                for j= 1:size(origCData,3)
                    New_texture_map(:,:,j) = colfilt(origCData(:,:,j),patch_in_vector,'sliding',@kurtosis);
                end
                
                % Local_Entropy
            case 4
                New_texture_map = entropyfilt(origCData, patch_2D);
                
                % Local_Range
            case 5
                for j= 1:size(origCData,3)
                    New_texture_map(:,:,j) = colfilt(origCData(:,:,j),patch_in_vector,'sliding',@range);
                end
                
                % Local_SD
            case 6
                New_texture_map = stdfilt(origCData, patch_2D);   
        end
        
        files_out_nbr = files_out_nbr+1;
        %transform the OutputImages matrix in order to match to the nii header of the
        % first input (rotation/translation)
        
        map_reoriented = write_volume(New_texture_map, nifti_header, 'axial');
        
        info2 = info;
        info2.Filename = files_out.In1{files_out_nbr};
        info2.Filemoddate = char(datetime('now'));
        info2.Datatype = class(map_reoriented);
        info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
        info2.ImageSize = size(map_reoriented);
        info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];
        
        niftiwrite(map_reoriented, files_out.In1{files_out_nbr}, info2)
        
        % so far copy the .json file of the first input
       % copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{files_out_nbr}, '.nii', '.json'))

        if exist(strrep(files_in.In1{1}, '.nii', '.json'), 'file') == 2
            J = ReadJson(strrep(files_in.In1{1}, '.nii', '.json'));
            J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
            
            [path, name, ~] = fileparts(files_out.In1{files_out_nbr});
            jsonfile = [path, '/', name, '.json'];
            WriteJson(J, jsonfile)
        end
    end
end




%
% if strcmp(opt.Local_Entropy,'Yes')
%     if mod(opt.Patch_size,2)
%         Entropy_map = entropyfilt(rescaled_Im, true(opt.Patch_size));
%     else
%         Entropy_map = entropyfilt(rescaled_Im, true(opt.Patch_size)+1);
%     end
%     % transform the OutputImages matrix in order to match to the nii header of the
%     % first input (rotation/translation)
%     
%     map_reoriented = write_volume(Entropy_map, nifti_header, 'axial');
%     
%     info2 = info;
%     info2.Filename = files_out.In1{1};
%     info2.Filemoddate = char(datetime('now'));
%     info2.Datatype = class(map_reoriented);
%     info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
%     info2.ImageSize = size(map_reoriented);
%     info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];
%     
%     niftiwrite(map_reoriented, files_out.In1{1}, info2)
%     
%     % so far copy the .json file of the first input
%     copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
%     
% end
% if strcmp(opt.Local_Range,'Yes')
%     if mod(opt.Patch_size,2)
%         Range_map =rangefilt(rescaled_Im, true(opt.Patch_size));
%     else
%         Range_map =rangefilt(rescaled_Im, true(opt.Patch_size)+1);
%     end
%     map_reoriented = write_volume(Range_map, nifti_header, 'axial');
%     
%     info2 = info;
%     info2.Filename = files_out.In1{2};
%     info2.Filemoddate = char(datetime('now'));
%     info2.Datatype = class(map_reoriented);
%     info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
%     info2.ImageSize = size(map_reoriented);
%     info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];
%     
%     niftiwrite(map_reoriented, files_out.In1{2}, info2)
%     
%     % so far copy the .json file of the first input
%     copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{2}, '.nii', '.json'))
%     
% end
% if strcmp(opt.Local_SD,'Yes')
%     if mod(opt.Patch_size,2)
%         Std_map = stdfilt(rescaled_Im, true(opt.Patch_size));
%     else
%         Std_map =stdfilt(rescaled_Im, true(opt.Patch_size)+1);
%     end
%     map_reoriented = write_volume(Std_map, nifti_header, 'axial');
%     
%     info2 = info;
%     info2.Filename = files_out.In1{3};
%     info2.Filemoddate = char(datetime('now'));
%     info2.Datatype = class(map_reoriented);
%     info2.PixelDimensions = info.PixelDimensions(1:length(size(map_reoriented)));
%     info2.ImageSize = size(map_reoriented);
%     info2.DisplayIntensityRange = [nanmin(map_reoriented(:))  nanmax(map_reoriented(:))];
%     
%     niftiwrite(map_reoriented, files_out.In1{3}, info2)
%     
%     % so far copy the .json file of the first input
%     copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{3}, '.nii', '.json'))
% end
