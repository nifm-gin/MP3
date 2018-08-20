function [files_in,files_out,opt] = Module_Susceptibility(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
% define every option needed to run this module
%     %   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','AllName'};
    module_option(:,4)   = {'output_filename_ext_CBV','CBV'};
    module_option(:,5)   = {'output_filename_ext_CBF','CBF'};
    module_option(:,6)   = {'output_filename_ext_MTT','MTT'};
    module_option(:,7)   = {'output_filename_ext_TMAX','TMAX'};
    module_option(:,8)   = {'output_filename_ext_TTP','TTP'};
    module_option(:,9)   = {'output_filename_ext_T0','T0'}; 
    module_option(:,10)   = {'RefInput', 1};
    module_option(:,11)   = {'InputToReshape',1};
    module_option(:,12)   = {'Table_in', table()};
    module_option(:,13)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    opt.NameOutFiles = {opt.Module_settings.output_filename_ext_CBV, opt.Module_settings.output_filename_ext_CBF, opt.Module_settings.output_filename_ext_MTT, opt.Module_settings.output_filename_ext_TMAX, opt.Module_settings.output_filename_ext_TTP, opt.Module_settings.output_filename_ext_T0};
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
    user_parameter(:,1)   = {'Select one PERF scan as input','1Scan','','',{'SequenceName'},'Mandatory',''};
    user_parameter(:,2)   = {'Parameters','','','','','',''};
    user_parameter(:,3)   = {'   .Output filename extension CBV','char','CBV','output_filename_ext_CBV','','',''};
    user_parameter(:,4)   = {'   .Output filename extension CBF','char', 'CBF','output_filename_ext_CBF','','',''};
    user_parameter(:,5)   = {'   .Output filename extension MTT','char','MTT','output_filename_ext_MTT','','',''};
    user_parameter(:,6)   = {'   .Output filename extension TMAX','char','TMAX','output_filename_ext_TMAX','','',''};
    user_parameter(:,7)   = {'   .Output filename extension TTP','char','TTP','output_filename_ext_TTP','','',''};
    user_parameter(:,8)   = {'   .Output filename extension T0','char','T0','output_filename_ext_T0','','',''};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
%     
%%
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%

opt.NameOutFiles = {opt.output_filename_ext_CBV, opt.output_filename_ext_CBF, opt.output_filename_ext_MTT, opt.output_filename_ext_TMAX, opt.output_filename_ext_TTP, opt.output_filename_ext_T0};




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
    error('Module_Susceptibility:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end


% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, ext_nii, files_in);  
% end
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



N = niftiread(files_in.In1{1});
if size(N,4) == 1
   error([files_in.In1{1} ' is not a 4d image']);
    %return
end

N = double(N);
info = niftiinfo(files_in.In1{1});
[path, name, ext] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
J = jsondecode(raw);


%% Processing
meansignal = mean(squeeze(N),4);
volume_mask = meansignal>max(N(:))*0.01;
[aif,scores] = extraction_aif_volume(squeeze(N),volume_mask);

if sum(cell2num(scores(:,5))) > 10
    error('No computation because the AIF is not good enough');
else
    [~,~,~,TMAX,TTP,T0, CBV,CBF,MTT,~,~,~,~] = deconvolution_perfusion_gui(aif,squeeze(N),J.RepetitionTime.value(1)*10^(-3),J.EchoTime.value*10^(-3));
    %maps = {'CBV','CBF','MTT','TMAX','TTP','T0'};
    mapsVar = {CBV, CBF, MTT, TMAX, TTP, T0};

end


%% Reshape to the input scan size
%FilteredImages = reshape(FilteredImages, Size);


for i=1:length(mapsVar)
    info2 = info;
    info2.Filename = files_out.In1{i};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(mapsVar{i});
    info2.PixelDimensions = info.PixelDimensions(1:length(size(mapsVar{i})));
    info2.ImageSize = size(mapsVar{i});
    %info2.Description = [info.Description, 'Modified by Susceptibility Module'];
    niftiwrite(mapsVar{i}, files_out.In1{i}, info2)
    

    JMod = jsonencode(J);
    [path, name, ext] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    fidmod = fopen(jsonfile, 'w');
    fwrite(fidmod, JMod, 'uint8');
    fclose(fidmod);
end



