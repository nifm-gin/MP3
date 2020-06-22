function [files_in,files_out,opt] = Module_VSI(files_in,files_out,opt)

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
    module_option(:,4)   = {'B0',4.7};
    module_option(:,5)   = {'Gamma',2.67502};
    module_option(:,6)   = {'deltaxi_CA',0.1867};
    module_option(:,7)   = {'RefInput',1};
    module_option(:,8)   = {'InputToReshape',1};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)   = {'Table_out', table()};
    module_option(:,11)   = {'output_filename_ext','VSI'};
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
    user_parameter(:,1)   = {'Description','Text','','','', '',...
        {'Description of the module :',...
        'this module compute a Vessel Sice Index map (VSI) from a deltaR2* (see module deltaR2*), a deltaR2 (see module deltaR2) and de diffusion map (ADC)',...
        'If you use this function, please refere to this article : Tropes et al. MRM 2015'}};
    
    user_parameter(:,2)   = {'Select the DeltaR2 scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the DeltaR2* scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'Select the ADC scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,5)   = {'   .Output filename','char','VSI','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''VSI''.'}'};
    user_parameter(:,6)   = {'   .Output orientation','cell',{'First input', 'Second input', 'Third Input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        '--> Output orientation = Third input'
        }'};
    user_parameter(:,7)   = {'   .B0','numeric', '','B0','','',...
    {'Please entre the magnetic field used'}'};
    user_parameter(:,8)   = {'   .Gamma (10^8)','numeric','','Gamma','','',...
    {''}'};
    user_parameter(:,9)   = {'   .deltaxi HBC (10^-6) (cgs units)','numeric', '','deltaxi_CA','','',...
    {''}'};
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
    error('Module_VSI:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
input(3).nifti_header = spm_vol(files_in.In3{1});

J_Delta = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
J_DeltaStar = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));
J_ADC = spm_jsonread(strrep(files_in.In3{1}, '.nii', '.json'));


if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
elseif strcmp(opt.Output_orientation, 'Second input')
    ref_scan = 2;
else
    ref_scan = 3;
end
% all data have to be in double before any calculs
DeltaR2 = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
if isa(class(DeltaR2), 'double') ==0
    DeltaR2 = double(DeltaR2);
end
DeltaR2Star = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);
if isa(class(DeltaR2Star), 'double') ==0
    DeltaR2Star = double(DeltaR2Star);
end
ADC = read_volume(input(3).nifti_header, input(ref_scan).nifti_header, 0);
if isa(class(DeltaR2Star), 'double') ==0
    DeltaR2Star = double(DeltaR2Star);
end

B0=opt.B0;
gamma=opt.Gamma; gamma = gamma*10^8;
deltaxi=opt.deltaxi_CA; deltaxi = deltaxi*10^-6;


temp_thickness = [J_DeltaStar.SliceThickness J_Delta.SliceThickness  J_ADC.SliceThickness];
temp_slice_nbr = [size(DeltaR2Star,3) size(DeltaR2,3)  size(ADC,3)];
temp_resolution = [size(DeltaR2Star,1)*size(DeltaR2Star,2) size(DeltaR2,1)*size(DeltaR2,2)  size(ADC,1)*size(ADC,2)];

% check data compatibility (slice thickness and slice number)
% if  length(find(temp_thickness == deltaR2star.reco.thickness)) ~= numel(temp_thickness)
%     warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ a slice thickness missmatch between\n##$deltaR2star=%s\n##$ and \n##$deltaR2=%s\n##$ and\n##$ADC=%s\n',...
%         deltaR2star_filename,deltaR2_filename,ADC_filename);
%     msgbox(warning_text, 'VSI map warning') ;
%      VSImap = [];
%     return
% end
% if length(find(temp_resolution == deltaR2star.reco.no_samples)) ~= numel(temp_resolution)
%      warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ a resolution missmatch between\n##$deltaR2star=%s\n##$ and \n##$deltaR2=%s\n##$ and\n##$ADC=%s\n',...
%        deltaR2star_filename,deltaR2_filename,ADC_filename);
%     msgbox(warning_text, 'VSI map warning') ;
%      VSImap = [];
%     return
% end

VSImap = DeltaR2;

if length(find(temp_slice_nbr == size(DeltaR2Star,3))) ~= numel(temp_slice_nbr)
    VSImap_slice_nbr = 0;
    for i = 1:size(DeltaR2Star, 3)
        for j = 1:size(DeltaR2, 3)
                for x = 1:size(ADC, 3)
                        VSImap_slice_nbr = VSImap_slice_nbr+1;
                        % Compute the VSI map each slice with the same offset
                        ratio = DeltaR2Star(:,:,i,1) ./ DeltaR2(:,:,j,1); % ratio : no unit
                        index_ratiofinite= find(~isfinite(ratio));
                        index_rationan= find(isnan(ratio));
                        if ~isempty(index_ratiofinite)
                            ratio(index_ratiofinite) = 0;
                        end
                        if ~isempty(index_rationan)
                            ratio(index_rationan)    = 0;
                        end
                        index_neg= find(ratio<0);
                        if ~isempty(index_neg)
                            ratio(index_neg)    = 0;
                        end
                        ratio = ratio.^1.5;
                          slice_data = 1.77^(-1.5) * sqrt(squeeze(ADC(:,:,x,1)) ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
                       %%%% Rui data
%                        test = imresize(squeeze(ADC.reco.data(:,:,1,x)), [128 128]);
%                        index_inferieurzero=find(test < 0);
%                        test(index_inferieurzero) = 0;
%                        slice_data = 1.77^(-1.5) * sqrt(test ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
                        %%%
                        
                        index_vsifinite=find(~isfinite(slice_data));
                        slice_data(index_vsifinite)= nan;
                        index_vsinan=find(isnan(slice_data));
                        slice_data(index_vsinan)= nan;
                        index_infzero=find(slice_data < 0);
                        slice_data(index_infzero)= nan;
                        VSImap(:,:,VSImap_slice_nbr,1) = slice_data;
                        
                end

        end
    end
    if VSImap_slice_nbr == 0
        warning_text = sprintf('##$ Can not calculate the VSI map because there is\n##$ no slice offset match between\n##$deltaR2star=%s\n##$ and \n##$deltaR2=%s\n##$ and\n##$ADC=%s\n',...
            deltaR2star_filename,deltaR2_filename,ADC_filename);
        msgbox(warning_text, 'VSI map warning') ;
        return
    end
else
    % calculate VSI map for all slices
    ratio = DeltaR2Star ./ DeltaR2; % ratio : no unit
    index_ratiofinite= find(~isfinite(ratio));
    index_rationan= find(isnan(ratio));
    if ~isempty(index_ratiofinite)
        ratio(index_ratiofinite) = NaN;
    end
    if ~isempty(index_rationan)
        ratio(index_rationan)    = NaN;
    end
     index_neg= find(ratio<0);
    if ~isempty(index_neg)
        ratio(index_neg)    = NaN;
    end
    ratio = ratio.^1.5;
    VSImap(:,:,:,1) = 1.77^(-1.5) * sqrt(squeeze(ADC(:,:,:,1)) ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
    % VSImap(:,:,:,1) = 1.77^(-1.5) * sqrt(squeeze(ADC) ./ (gamma * deltaxi * B0)) .* squeeze(ratio);
    
%     index_vsifinite=find(~isfinite(VSImap));
%     VSImap(index_vsifinite)= 0;
%     index_vsinan=find(isnan(VSImap));
%     VSImap(index_vsinan)= 0;
%     index_infzero=find(VSImap < 0);
%     VSImap(index_infzero)= 0;
    
end


OutputImages = VSImap;
OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 50) = NaN;
%OutputImages(isnan(OutputImages)) = NaN;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
end



% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
%info2.Description = [info.Description, 'Modified by the DeltaR2 Module'];

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


