function [files_in,files_out,opt] = Module_R2Prim(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the module's parameters with default values
if isempty(opt)
    %
    %     %%   % define every option needed to run this module
    %     % --> module_option(1,:) = field names
    %     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','AllName'};
    module_option(:,4)   = {'RefInput',1};
    module_option(:,5)   = {'InputToReshape',1};
    module_option(:,6)   = {'Table_in', table()};
    module_option(:,7)   = {'Table_out', table()};
    module_option(:,8)   = {'output_filename_ext','R2Prim'};
    module_option(:,9)   = {'Output_orientation','First input'};
    
    
    
    
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
    user_parameter(:,1)   = {'Description','Text','','','', '',{'Description of the module',...
        'this module compute a R2prim map from a T2map and from a T2corrected3D (see module MGE2Dfrom3D) ',...
        'If you use this function, please refere to this article : Christen et al. NMR in bimed 2011'}};
    
    user_parameter(:,2)   = {'Select the T2_map scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the T2*corr3D scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','R2Prim','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''R2Prim''.'}'};
    user_parameter(:,5)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        }'};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    %%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
    
end
%%%%%%%%

if isempty(files_out)
    opt.Table_out = opt.Table_in(opt.RefInput,:);
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), '_', opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_R2Prim:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%% load input Nii file
input(1).nifti_header = spm_vol(files_in.In1{1});
input(2).nifti_header = spm_vol(files_in.In2{1});

JT2Map = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
JT2StarCorr3D = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));


if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
T2Map = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
T2StarCorr3D = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
%T2Map = niftiread(files_in.In1{1});
%T2StarCorr3D = niftiread(files_in.In2{1});

% check data compatibility (slice thickness and slice number)
% if JT2StarCorr3D.SliceThickness.value ~= JT2Map.SliceThickness.value
%     warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ a slice thickness missmatch between\n##$T2*map=%s\n##$ and \n##$Tmap=%s',...
%         T2star_filename,T2map_filename);
%     msgbox(warning_text, 'R2prim map warning') ;
%     R2prim_map = [];
%     return
% end


%R2prim_map = T2StarCorr3D;
% if size(T2Map,3) ~= size(T2StarCorr3D,3)
%     for i = 1:size(T2StarCorr3D, 3)
%         for j = 1:size(T2Map, 3)
%             R2prim_slice_nbr = R2prim_slice_nbr+1;
%             % Compute the CMRO2 map each slice with the same offset
%             R2prim_map(:,:,R2prim_slice_nbr,1)=1./(RAW_T2star.data(:,:,i,1)*10^-3)-1./(RAW_T2map.data(:,:,j,1)*10^-3);
%         end
%     end
%     if R2prim_slice_nbr == 0
%         warning_text = sprintf('##$ Can not calculate the R2prim map because there is\n##$ no slice offset match between\n##$T2*map=%s\n##$ and \n##$T2map=%s',...
%             files_in.In2{1},files_in.In1{1});
%         msgbox(warning_text, 'R2prim map warning') ;
%         return
%     end
% else
%     final_res = max([size(T2Map,1), size(T2Map, 2), size(T2StarCorr3D, 1), size(T2StarCorr3D, 2)]);
%     if size(T2StarCorr3D,1) ~= final_res || size(T2StarCorr3D,2) ~= final_res
%         for i = 1:size(T2StarCorr3D, 3)
%             tmp(:,:,i,1)= imresize(T2StarCorr3D(:,:,i,1),[final_res final_res],'bilinear');
%         end
%         T2StarCorr3D = tmp;
%     end
%     if size(T2Map,1) ~= final_res || size(T2Map,2) ~= final_res
%         for i = 1:size(T2Map,3)
%             tmp(:,:,i,1)= imresize(T2Map(:,:,i,1),[final_res final_res],'bilinear');
%         end
%         T2Map = tmp;
%     end
%     clear tmp
%     
%     % Compute the R2prim map for all slices
%     R2prim_map=1./(T2StarCorr3D*10^-3)-1./(T2Map*10^-3);
% end


R2prim_map=1./(T2StarCorr3D*10^-3)-1./(T2Map*10^-3);




%replace negative value of R2prim by the interpolation of their neibourgh 
R2prim_map_m = convn(R2prim_map,ones(3,3,3)./26,'same');
R2prim_map(R2prim_map(:)<=0) =R2prim_map_m(R2prim_map(:)<=0);


OutputImages = R2prim_map;
OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 25) = NaN;
%OutputImages(isnan(OutputImages)) = NaN;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header, 'Axial');
end


% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(1)]){1});

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
info2.ImageSize = size(OutputImages_reoriented);

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


