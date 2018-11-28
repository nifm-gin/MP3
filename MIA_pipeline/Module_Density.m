function [files_in,files_out,opt] = Module_Density(files_in,files_out,opt)

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
    module_option(:,8)   = {'output_filename_ext','Density'};
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
          user_parameter(:,1)   = {'Description','Text','','','', '',...
        {'Description of the module :',...
        'this module compute a Vessel Density map from a deltaR2* (see module deltaR2*) and a deltaR2 (see module deltaR2) ',...
        'If you use this function, please refere to this article : Tropes et al. MRM 2015'}};
    
    user_parameter(:,2)   = {'Select the DeltaR2 scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the DeltaR2Star scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','Density','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''Density''.'}'};
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
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Density:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

J1 = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
J2 = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));


if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
DeltaR2 = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
DeltaR2Star = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);

%DeltaR2 = niftiread(files_in.In1{1});
%DeltaR2Star = niftiread(files_in.In2{1});




% 
% if size(DeltaR2Star, 3) ~= size(DeltaR2, 3)
%     for i = 1:size(DeltaR2, 3)
%         for j = 1:size(DeltaR2Star, 3)
% %            if abs(deltaR2.reco.fov_offsets(3,1,i) - deltaR2star.reco.fov_offsets(3,1,j)) < 1e-5
%                 Density_map_slice_nbr = Density_map_slice_nbr+1;
%                 %ref lemasson et al 2013
%                 data(:,:,Density_map_slice_nbr,:)=329000.*((DeltaR2(:,:,i,:).^3)./(DeltaR2Star(:,:,j,:).^2));
%                 
% %            end
%         end
%     end
%     if Density_map_slice_nbr == 0
%         warning_text = sprintf('##$ Can not calculate the Density map because there is\n##$ no slice offset match between\n##$DeltaR2=%s\n##$ and \n##$DeltaR2*=%s',...
%             filename_SO2map,filename_CBF);
%         msgbox(warning_text, 'Density map warning') ;
%         return
%     end
%     Density_map.reco.no_slices=Density_map_slice_nbr;
% else
%     % Compute theDensity_map for all slices
%     data=329000.*((DeltaR2.^3)./(DeltaR2Star.^2));
% end
data=329000.*((DeltaR2.^3)./(DeltaR2Star.^2));


Density_map = data;


% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)

OutputImages = Density_map;
OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 10000) =NaN;
OutputImages(isinf(OutputImages)) = NaN;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
else 
    OutputImages_reoriented = OutputImages;
end





% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});


info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
info2.ImageSize = size(OutputImages_reoriented);
%info2.Description = [info.Description, 'Modified by T2map Module'];


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
