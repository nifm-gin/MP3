function [files_in,files_out,opt] = Module_Brain_Mask_PCNN3D(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename','Brain'};
    module_option(:,4)   = {'OutputSequenceName','AllName'};
    module_option(:,5)   = {'Min_BrainSize', 1200};
    module_option(:,6)   = {'Max_BrainSize', 4400};
    module_option(:,7)   = {'Radius',5};
    module_option(:,8)   = {'RefInput',1};
    module_option(:,9)   = {'InputToReshape',1};
    module_option(:,10)   = {'Table_in', table()};
    module_option(:,11)   = {'Table_out', table()};
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
        'Brain Extraction Tool from PCNN3D tool'
        ''
        'This module needs a brain anatomical image as input'
        'Output --> a brain mask saved as ROI'
        ''
        'Tips :'
        '      1. Set proper brain volume range is most critical for getting good result.'
        '         If your data has only partial brain coverage, reduce the max volume.'
        '         Also, if the data has resolution magnified by 10x, the volume range should be magnified by 1000x.'
        ''
         '     2. Effect of "Radius of the structural element" can be seen below. Note: larger radius, slower the processing.'     
        ''
        'Ref:'
        'N Chou, J Wu, J Bai, A Qiu, KH Chuang. Robust automatic rodent brain extraction using 3D Pulse-coupled NeuralNetworks (PCNN). IEEE Trans Imag Proc 20(9):2554-64, 2011.'
        }'};
    
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory','Please select an anatomical image'};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename','char','Brain','output_filename','', '',''};
    user_parameter(:,5)   = {'   .Min Brain size','numeric', 1200,'Min_BrainSize','', '',...
        'Min BrainSize: [min max] brain volume sizes (in mm3) for estimating optimal iteration; default values [100 550] is for mouse brain. For adult rat brain, one may try [1200 4400]'};
    user_parameter(:,6)   = {'   .Max Brain size','numeric', 4400,'Max_BrainSize','', '',...
        'Max BrainSize: [min max] brain volume sizes (in mm3) for estimating optimal iteration; default values [100 550] is for mouse brain. For adult rat brain, one may try [1200 4400]'};
     user_parameter(:,7)   = {'   .Radius (in pixel)?','numeric',5,'Radius','', '',...
         'radius of the structural element (in pixel); For example for a 256x256 field-of-view 5 is a good number'};
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
    opt.Table_out.Type = categorical(cellstr('ROI'));
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
    error('Smoothing:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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


nifti_info = niftiinfo(files_in.In1{1});
I = niftiread(files_in.In1{1});
voxel_volume = prod(nifti_info.raw.pixdim(2:4));
voxel_surface = prod(nifti_info.raw.pixdim(2:3));

I(I(:) < 0) = nan; % PCNN3D does not support negative values

maxiter = 200;
OptMask = 0;
[I_border, Plot, optG] = PCNN3D(I, opt.Radius, nifti_info.raw.pixdim(2:4), [opt.Min_BrainSize opt.Max_BrainSize], maxiter, OptMask);
% Automatically close the figure which contains the plot
close(gcf)
% Compute new optimal 
Plot(Plot < opt.Min_BrainSize) = NaN;
[~,New_optG] = min(gradient(Plot));
% dispo for info
diff_optG = New_optG -  optG;
disp(strcat('Automatic Optimum versus calculated one:', num2str(optG), 'vs ' ,num2str(New_optG), '=', num2str(diff_optG)))
mtx=size(I);
Mask=zeros(mtx);
for n=1:mtx(3)
    Mask(:,:,n)=I_border{New_optG}{n}; % the "optimal" segmentation is always to high
end
Mask = single(Mask);

info = niftiinfo(files_in.In1{1});
info2 = nifti_info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(Mask);

niftiwrite(Mask, files_out.In1{1}, info2)



% %% Json processing
% [path, name, ~] = fileparts(files_in.In1{1});
% jsonfile = [path, '/', name, '.json'];
% J = ReadJson(jsonfile);
% 
% J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
% 
% [path, name, ~] = fileparts(files_out.In1{1});
% jsonfile = [path, '/', name, '.json'];
% WriteJson(J, jsonfile)
