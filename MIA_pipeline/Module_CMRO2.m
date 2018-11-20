function [files_in,files_out,opt] = Module_CMRO2(files_in,files_out,opt)

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
    module_option(:,5)   = {'output_filename','CMRO2'};
    module_option(:,6)   = {'Ca',20.43};
    module_option(:,7)  = {'OutputSequenceName','AllName'};
    module_option(:,8)  = {'RefInput',1};
    module_option(:,9)  = {'InputToReshape',1};
    module_option(:,10)  = {'Table_in', table()};
    module_option(:,11)  = {'Table_out', table()};
    module_option(:,12)   = {'Output_orientation','First input'};
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
     user_parameter(:,1)   = {'Description','Text','','','','',...
        {'This module compute a CMRO2 (oxygen consumption) using a SO2 and a CBF map (see SO2 and CBF modules for more information)'}'};
    user_parameter(:,2)   = {'Select a CBF scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select a SO2 scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'Parameters','','','','', '',''};
    user_parameter(:,5)   = {'   .Output filename','char','CMRO2','output_filename','','','Specify the name of the output file.'};
    user_parameter(:,6)   = {'   .Constant (Ca)','numeric', '','Ca','','',{''}'};
    user_parameter(:,7)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        }'};

VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'CMRO2'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename};

if isempty(files_out)
    opt.Table_out = opt.Table_in(opt.RefInput,:);
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_CMRO2:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

if isempty(files_out)
   files_out.In1 = {cat(2,opt.folder_out,filesep,name_nii,'_',opt.output_filename,ext_nii)};
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


input(1).nifti_header = spm_vol(files_in.In1{1});
input(2).nifti_header = spm_vol(files_in.In2{1});
J_CBF = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
J_SO2 = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));

if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
CBF = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
SO2 = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);


if size(CBF, 4)>1
    CBF = CBF(:,:,:,1);
end


% data = CBF.reco.data.*(1/SO2map.reco.data).*Ca;

% formule Julien
% data = Ca*(CBF./100).*(OEF./100);

Ca = opt.Ca; %14.7*1.39;  % 14.7g/dl d'Hb      1g Hb transporte 1.39mL d'O2



% check data compatibility (slice thickness and slice number)
% if J_CBF.SliceThickness.value ~= J_SO2.SliceThickness.value
%     warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ a slice thickness missmatch between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
%         files_in.In2{1},files_in.In1{1});
%     msgbox(warning_text, 'CMRO2 map warning') ;
%     return
% end


if size(CBF,3) ~= size(SO2,3)
    for i = 1:size(CBF, 3)
        for j = 1:size(SO2, 3)
                CMRO2_slice_nbr = CMRO2_slice_nbr+1;
                % Compute the CMRO2 map each slice with the same offset
                %                 data(:,:,:,CMRO2_slice_nbr)=Ca.*(imresize(squeeze(CBF.reco.data(:,:,:,i)), 0.5)./100).*((100-SO2map.reco.data(:,:,:,j))./100);
                data(:,:,CMRO2_slice_nbr,:)=Ca.*(CBF(:,:,i,:)./100).*((100-SO2(:,:,j,:))./100);

        end
    end
    if CMRO2_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the CMRO2 map because there is\n##$ no slice offset match between\n##$SO2map=%s\n##$ and \n##$CBFmap=%s',...
            files_in.In2{1},files_in.In1{1});
        msgbox(warning_text, 'CMRO2 map warning') ;
        return
    end
else
    % Compute the CMRO2 map for all slices
    data=Ca.*(CBF./100).*((100-SO2)./100);
end

CMRO2 = data;


OutputImages = CMRO2;
% OutputImages(OutputImages < 0) = -1;
% OutputImages(OutputImages > 1000) = -1;
% OutputImages(isnan(OutputImages)) = -1;


OutputImages(isinf(OutputImages)) = -1;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
end


info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});
% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
info2.ImageSize = size(OutputImages_reoriented);
%info2.Description = [info.Description, 'Modified by T2star_map Module'];

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
