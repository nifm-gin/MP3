function [files_in,files_out,opt] = Module_DeltaR2Star(files_in,files_out,opt)

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
    module_option(:,4)   = {'RefInput',2};
    module_option(:,5)   = {'InputToReshape',2};
    module_option(:,6)   = {'Table_in', table()};
    module_option(:,7)   = {'Table_out', table()};
    module_option(:,8)   = {'output_filename_ext','DeltaR2Star'};
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
        {'Description of the module : this module compute a deltaR2* map from two T2* maps acquired before and after a UPSIO injection',...
        'Output map in ms-1', ...
        'If you use this function, please refere to this article : Tropes et al. MRM 2015'}};
    user_parameter(:,2)   = {'Select the T2_star_Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the T2_star_Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','DeltaR2Star','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''DeltaR2Star''.'}'};
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
    error('Module_DeltaR2Star:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
input_pre = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
input_post = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);


if size(input_pre, 3) ~= size(input_post, 3)
    deltaR2star_slice_nbr = 0;
    
    for i = 1:size(input_pre, 3)
        for j = 1:size(input_post, 3)
            deltaR2star_slice_nbr = deltaR2star_slice_nbr+1;
            r2savant=input_pre(:,:,i,1);
            r2sapres=input_post(:,:,j,1);
            index_avant=find(r2savant<10e-2);% looking for the T2*pre < 100�s
            warning off					%To avoid division by zero warning
            r2savant=1./r2savant;       % R2*pre (ms-1)
            warning on
            r2savant(index_avant)=0;
            index_apres=find(r2sapres<10e-2);% olooking for the T2*post  < 100�s
            warning off					%To avoid division by zero warning
            r2sapres=1./r2sapres;	% R2*post (ms-1)
            warning on %#ok<*WNON>
            r2sapres(index_apres)=0;
            DeltaR2StarMap(:,:,deltaR2star_slice_nbr,1)=r2sapres-r2savant;
            index_nan=find(input_pre(:,:,j,1)==NaN); %#ok<FNAN>
            DeltaR2StarMap(index_nan)=0;
        end
    end
else
    %Initialisation des variables
    r2savant=input_pre(:,:,:,1);
    r2sapres=input_post(:,:,:,1);
    DeltaR2StarMap=zeros(size(input_pre));     %param�tres (DeltaR2*)
    index_avant=find(r2savant<10e-2);% on cherche les T2* < 100�s
    warning off					%To avoid division by zero warning
    r2savant=1./r2savant;	% R2*avant en ms-1
    warning on
    r2savant(index_avant)=0;
    index_apres=find(r2sapres<10e-2);% on cherche les T2* < 100�s
    warning off					%To avoid division by zero warning
    r2sapres=1./r2sapres;	% R2*apres en ms-1
    warning on %#ok<*WNON>
    r2sapres(index_apres)=0;
    DeltaR2StarMap(:,:,:,1)=r2sapres-r2savant;
    index_nan=find(input_pre==NaN); %#ok<FNAN>
    DeltaR2StarMap(index_nan)=0;
end

OutputImages = DeltaR2StarMap;
% OutputImages(OutputImages < 0) = -1;
% OutputImages(OutputImages > 5000) = -1;
% OutputImages(isnan(OutputImages)) = -1;
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




% 


