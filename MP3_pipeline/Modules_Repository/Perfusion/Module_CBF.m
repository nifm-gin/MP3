function [files_in,files_out,opt] = Module_CBF(files_in,files_out,opt)
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
    module_option(:,4)   = {'InvEffValue','Default'};
    module_option(:,5)   = {'QuantifMethod','Buxton1998'};
    module_option(:,6)   = {'Mean_Dyn','Mean'};
    module_option(:,7)   = {'OutExcl','Yes'};
    module_option(:,8)   = {'OutPlot','Yes'};
    module_option(:,9)   = {'M0Value','Default'};
    module_option(:,10)  = {'T1Value','Default'};
    module_option(:,11)  = {'TransitValue','Default'};
    module_option(:,12)  = {'RefInput',1};
    module_option(:,13)  = {'InputToReshape',1};
    module_option(:,14)  = {'Table_in', table()};
    module_option(:,15)  = {'Table_out', table()};
    module_option(:,16)  = {'output_filename_ext','CBF'};
    module_option(:,17)  = {'Output_orientation','First input'};
    

%                 additional_info_name = {'ASL_InvEff', 'InvEff ROI', 'InvEff value', ...
%                 'Method quantif',  'Mask?', 'Output: Mean or Dyn?', 'Outlier: Exclude?',  'Outlier: Plot?',...
%                 'M0 map', 'M0 value', 'T1 map', 'T1 value', 'Transit map', 'Transit value'};
%             
%             
%                         additional_info_data = {ASL_InvEffState, 'None', 'Default',...
%                 'Buxton1998', 'No', 'Mean', 'Yes', 'Yes',...
%                 'None', 'Default', T1mapState, 'Default','None', 'Default',};
%             additional_info_format = {['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], ['None' handles.MIA_data.VOIs(1:end-2)], 'numeric',...
%                 {'Buxton1998' 'ASL-diff' 'ASL-Hadamard-Bolus' 'HadamardFitTransitTime' 'Alsop1996'  'Parkes2002'}, ['None' handles.MIA_data.VOIs(1:end-2)],  {'Mean' 'Dyn'},...
%                 {'Yes' 'No'}, {'Yes' 'No'},...
%                 ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], 'numeric', ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], 'numeric',...
%                 ['None' 'NOT COMPUTED YET' handles.CPM_parameters_list], 'numeric'};
    
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
    user_parameter(:,1)   = {'Description','Text','','','', '','Description of the module'}  ;
    user_parameter(:,2)   = {'Select the pCasl scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the ASL_InvEff scan','1Scan','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,4)   = {'   .Output filename','char','CBF','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''CBF''.'}'};
    user_parameter(:,5)   = {'   .Output orientation','cell',{'First input', 'Second input', 'Third Input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        '--> Output orientation = Third input'
        }'};
    user_parameter(:,6)   = {'   .InvEff ROI','1ROI','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,7)   = {'   .Inversion efficiency value','numeric', 'Default','InvEffValue','','',{''}'};
    user_parameter(:,8)   = {'   .Quantification Method','cell', {'Buxton1998' 'ASL-diff' 'ASL-Hadamard-Bolus' 'HadamardFitTransitTime' 'Alsop1996'  'Parkes2002'},'QuantifMethod','','',...
    {''}'};
    user_parameter(:,9)   = {'   .Output: Mean or Dyn?','cell',{'Mean' 'Dyn'},'Mean_Dyn','','',...
    {''}'};
    user_parameter(:,10)   = {'   .Outlier: Exclude?','cell', {'Yes' 'No'},'OutExcl','','',...
    {''}'};
    user_parameter(:,11)   = {'   .Outlier: Plot?','cell', {'Yes' 'No'},'OutPlot','','',{''}'};
    user_parameter(:,12)   = {'   .M0 map','1Scan', '','',{'SequenceName'},'Optional',{''}'};
    user_parameter(:,13)   = {'   .M0 value','numeric', 'Default','M0Value','','',{''}'};
    user_parameter(:,14)   = {'   .T1 map','1Scan', '','',{'SequenceName'},'Optional',{''}'};
    user_parameter(:,15)   = {'   .T1 value','numeric', 'Default','T1Value','','',{''}'};
    user_parameter(:,16)   = {'   .Transit map','1Scan', '','',{'SequenceName'},'Optional',{''}'};
    user_parameter(:,17)   = {'   .Transit value','numeric', 'Default','TransitValue','','',{''}'};
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
    error('ModuleCBF:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

% if strcmp(opt.Output_orientation, 'First input')   
%     ref_scan = 1;
% else
%     ref_scan = 2;
% end
ref_scan = 1 ;

input(1).nifti_header = spm_vol(files_in.In1{1});
pCasl = read_volume(input(1).nifti_header, input(1).nifti_header, 0);
J_pCasl = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));

if isfield(files_in, 'In2')
    input(2).nifti_header = spm_vol(files_in.In2{1});
    ASL_InvEff = read_volume(input(2).nifti_header, input(2).nifti_header, 0);
    J_ASL_InvEff = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));
end



InvEff_ROI = [];
M0Map = [];
T1Map = [];
TransitMap = [];

if isfield(files_in, 'In3')
    input(3).nifti_header = spm_vol(files_in.In3{1});
    InvEff_ROI = read_volume(input(3).nifti_header, input(2).nifti_header, 0);
end
if isfield(files_in, 'In4')
    input(4).nifti_header = spm_vol(files_in.In4{1});
    M0Map = read_volume(input(4).nifti_header, input(ref_scan).nifti_header, 0);
    J_M0Map = spm_jsonread(strrep(files_in.In4{1}, '.nii', '.json'));
end

if isfield(files_in, 'In5')
    input(5).nifti_header = spm_vol(files_in.In5{1});
    T1Map = read_volume(input(5).nifti_header, input(ref_scan).nifti_header, 0);
    J_T1Map = spm_jsonread(strrep(files_in.In5{1}, '.nii', '.json'));
end

if isfield(files_in, 'In6')
    input(6).nifti_header = spm_vol(files_in.In6{1});
    TransitMap = read_volume(input(6).nifti_header, input(ref_scan).nifti_header, 0);
    J_TransitMap = spm_jsonread(strrep(files_in.In6{1}, '.nii', '.json'));
end

if ~isnan(str2double(opt.InvEffValue))
    InvEffValue = str2double(opt.InvEffValue);
else
    InvEffValue = opt.InvEffValue;
end
if ~isnan(str2double(opt.TransitValue)) 
    TTvalue = str2double(opt.TransitValue);
else
    TTvalue = opt.TransitValue;
end
if ~isnan(str2double(opt.M0Value))
    M0value = str2double(opt.M0Value);
else
    M0value = opt.M0Value;
end
if ~isnan(str2double(opt.T1Value))
    T1value = str2double(opt.T1Value);
else
    T1value = opt.T1Value;
end

T1Sang = 'Default';
if strcmp(opt.Mean_Dyn, 'Mean')          % Mean of CBF accros expt at exit of program
    MeanYn = 1;
else
    MeanYn = 0;
end

% generate a CBF from an ASL scan
% this code comes from the the ASL module coded by C. Debacker
% Adaptation Paravision 6 et gestion ASL 3D : L. Hirschler

progression_Hadamard = 0;
%% Retrieve parameters
lambda = 0.9;                                            % brain-blood partition coefficient of water (ml/g) from biblio
quantification_method = opt.QuantifMethod;   % Choice of quantification method
if strcmp(quantification_method, 'ASL-Hadamard-Bolus')
    alpha = 'Default';
else
    if isempty(InvEff_ROI)      % inversion efficiency
        if strcmp(opt.InvEffValue, 'Default')
            alpha = 'Default';
        else
            alpha = InvEffValue;
        end
    else
        alpha = InvEff_ROI(:,:,1);
    end
end


if contains(J_pCasl.MethodDiffusion.value,'casl', 'IgnoreCase', true) && ~contains(J_pCasl.MethodDiffusion.value,'pcasl', 'IgnoreCase', true)
    ASLprefix = 'C';
else
    ASLprefix = 'pC';
end






% %% Alpha value : inversion efficiency % On utilise La carte Inv_eff et la ROI pour sortir une unique valeur (un pourcentage)
% % Load Alpha scan if necessary
% AlphaScanMethod = 'Not loaded';
% AlphaMethodNeed = 'ASL_InvEff';
% [AlphaScan,AlphaScan_status] = extract_map(ASL_InvEff_filename,finalMap,'map');
% 
% if AlphaScan_status && ~isempty(AlphaScan)
%     AlphaScanMethod = scan_acqp('##$QuantifMethod=',AlphaScan.reco.paramQuantif,0);
% end
% % Test if alpha measure is possible
% % Test if we have Voi and Scan for alpha
% if ~isempty( ASL_InvEff_filename ) && ~AlphaVoiYn && strcmp(extract_param(add_parameters,3), 'Default')
%     warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ You must enter a valid VOI name for Alpha value \n##$');
%     msgbox(warning_text, 'CBF map warning') ;
%     ASL_maps = [];
%     return
%     % Test if we have Voi and Scan for alpha
% elseif (~strcmp(AlphaScanMethod,AlphaMethodNeed) && ~isempty( ASL_InvEff_filename )) || (isempty( ASL_InvEff_filename ) && AlphaVoiYn )
%     warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ Scan for alpha must be ASL_InvEff type\n##$AlphaScan=%s\n##$',...
%         ASL_InvEff_filename);
%     msgbox(warning_text, 'CBF map warning') ;
%     ASL_maps = [];
%     return
% elseif ~AlphaVoiYn && ischar(alpha) && ~strcmp('Default',alpha)
%     warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ Voi don t exist with this pattern\n##$AlphaScan=%s\n##$',...
%         ASL_InvEff_filename);
%     msgbox(warning_text, 'CBF map warning') ;
%     ASL_maps = [];
%     return
% end
% 
% % Find good slice of Voi for Alpha
% if AlphaVoiYn
%     AlphaScanOS = AlphaScan.acq.fov_offsets(:,3);
%     NumVoiA = [];
%     for itVoi=1:numel(VoiAlpha)
%         VoiAlphaOS = VoiAlpha(itVoi).fov_offsets(3);
%         if VoiAlphaOS==AlphaScanOS
%             NumVoiA = itVoi;
%         end
%     end
%     
%     if isempty(NumVoiA)
%         warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ Voi offset and scan offset different \n##$AlphaScan=%s\n##$',...
%             ASL_InvEff_filename);
%         msgbox(warning_text, 'CBF map warning') ;
%         ASL_maps = [];
%         return
%     end
%     
%     VoiAlpha = VoiAlpha(NumVoiA);
%     if exist('VoiAlpha2','var')
%         VoiAlpha2 = VoiAlpha2(NumVoiA);
%     end
% end

% Alpha = value specified by user (nothing to do)
% Alpha = value from ASL_GEFC scan in a VOI specified by AlphaValue
% if( ~isempty(ASL_InvEff_filename) && AlphaVoiYn )
%     if isempty(VoiAlpha)
%         VoiAlpha.value = true(size(AlphaScan.reco.data));
%     end
%     factorResize = size(AlphaScan.reco.data,1) / size(VoiAlpha.value,1);
%     VoiAlphaTemp = logical(imresize(VoiAlpha.value,factorResize,'nearest'));
%     if exist('VoiAlpha2','var')
%         VoiAlphaTemp = VoiAlphaTemp + logical(imresize(VoiAlpha2.value,factorResize,'nearest'));
%     end
%     AlphaScanTemp = squeeze(AlphaScan.reco.data(:,:,1,1,1));
%     alpha = nanmean(AlphaScanTemp(logical(VoiAlphaTemp))) / 100;                 % /100 because ASL_GEFC give a value of alpha in percent
% end


if exist('ASL_InvEff', 'var')
    if isempty(InvEff_ROI)
        InvEff_ROI = true(size(ASL_InvEff));
    end
    alpha = nanmean(ASL_InvEff(logical(InvEff_ROI))) / 100; 
end







freq = J_pCasl.ImagingFrequency.value*1e6; % Mz
Rap_Gyr = 2.68e8;
B0 = 2*pi*freq/Rap_Gyr;
B0 = round(B0*10)/10; %Arrondi au 10eme

%% Set default value according to B0 magnetic field of scanner
% B0 = ASL.acq.champ_magnetique;
if  B0<3                          % 3T
    alphaPC = 0.88;  % not calculated
    alphaC = 0.81;  % not calculated
    T1sangD = 1660;  % Lu 2004
    T1valueD = 1100; % Stanisz 2005, gray matter
elseif  B0>=3 && B0<4.7                          % 3T
    alphaPC = 0.88;  % not calculated
    alphaC = 0.81;  % not calculated
    T1sangD = 1660;  % Lu 2004
    T1valueD = 1200; % Stanisz 2005, gray matter
elseif B0 >= 4.7 && B0<7            % 4.7T
    alphaPC = 0.88;
    alphaC = 0.81;
    T1sangD = 1880;  % Dobre 2007
    T1valueD = 1500;
elseif B0 >= 7 && B0<9.4            % 7T
    alphaPC = 0.88;
    alphaC = 0.84;
    T1sangD = 2230;  % Dobre 2007
    T1valueD = 1600;
elseif B0 >= 9.4 && B0<11.7         % 9.4T
    alphaPC = 0.88;  % not calculated
    alphaC = 0.8;  % not calculated
    T1sangD = 2430;  % NASRALLAH 2012
    T1valueD = 1650;
elseif B0 >= 11.7 && B0<17          % 11.7T
    alphaPC = 0.78;
    alphaC = 0.72;
    T1sangD = 2800;  % Lin 2012
    T1valueD = 1700;
elseif B0 >= 17              % 17.2T PARTIE A OPTIMISER
    alphaPC = 0.65;
    alphaC = 0.72;
    T1sangD = 3500;
    T1valueD = 2230;
else % default value for 4.7T
    alphaPC = 0.88;
    alphaC = 0.81;
    T1sangD = 1880;
    T1valueD = 1500;
end
if strcmp(ASLprefix,'pC') && strcmp(alpha,'Default')
    alpha = alphaPC;
elseif strcmp(ASLprefix,'C') && strcmp(alpha,'Default')
    alpha = alphaC;
end
if strcmp(T1Sang,'Default'), T1sang = T1sangD; end
if strcmp(T1value,'Default'), T1value = T1valueD; end

% %% Apply script of exclusion of outlier
% if size(pCasl,5) < 2 && strcmp(opt.OutExcl, 'Yes')
%     warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ The data haven''t several repetition but you ask to repair outlier, may be remove outlier detection \n##$ASL_Scan=%s\n##$',...
%         file_in.In2{1});
%     msgbox(warning_text, 'CBF map warning') ;
% elseif strcmp(opt.OutExcl, 'Yes')
%     alphaT = 1e-6;
%     interpolation = false;
%     plotDiagnostics = strcmp(opt.OutPlot, 'Yes');
%     detrend = false;
%     MeanMap = true;
%     
%     dASL = pCasl;
%     if strcmp(quantification_method, 'ASL-Hadamard-Bolus')
%         repair_outlier_slices({ASL_filename},{VoiMask_filename},alphaT,interpolation,plotDiagnostics,detrend,MeanMap);
%     else
%         dASL.reco.data = ASL.reco.data(:,:,2,:,:) - ASL.reco.data(:,:,1,:,:);
%         output = repair_outlier_slices({ASL_filename,dASL},{VoiMask_filename},alphaT,interpolation,plotDiagnostics,detrend,MeanMap);
%         if MeanYn
%             ASL = output{2};
%         else
%             ASL = output{1};
%         end
%     end
% end

% AcqDim = scan_acqp('##$PVM_SpatDimEnum=',ASL.texte,2);
%AcqDim = scan_acqp('##$PVM_SpatDimEnum=',ASL.texte,0);
AcqDim = J_pCasl.SpatDimEnum.value;

%% Initialisation parameter of ASL

% retrieve parameter values
%LabelTime=scan_acqp('ASL_LabelTime=',ASL.texte,1);
%PostLabelTime=scan_acqp('ASL_PostLabelTime=',ASL.texte,1);
LabelTime = J_pCasl.LabelTime.value;
PostLabelTime = J_pCasl.PostLabelTime.value;
% Calcul du temps interslice
if strcmp(AcqDim, '3D') == 1
    interSliceTime=0;
else
    interSliceTime = InterSliceTimeCaslEpi_Adapted(J_pCasl);
end
% Ordre d'acquisition des slices
if strcmp(AcqDim, '3D') == 1 || (size(pCasl,3)==1)
    PVM_ObjOrderList = 0;
else
    %ObjOrderListName = sprintf('##$PVM_ObjOrderList=( %d )',size(pCasl,3));
    %PVM_ObjOrderList=scan_acqp(ObjOrderListName,ASL.texte,1)';
    PVM_ObjOrderList = J_pCasl.ObjOrderList.value; % Est ce le bon JSON ...?
end
UnitCoeff = 60*100*1000;    % Conversion coefficient (ml blood/min/100ml tissue) en (ml blood/ms/ml tissue)

%% Initialized variables
if strcmp(quantification_method, 'ASL-Hadamard-Bolus')
    repetitions = size(pCasl,5);
    Hadamard_matrix = scan_acqp('CASL_HadamardMatrix=',ASL.texte,1);
    dataASL=zeros(size(pCasl,1),size(pCasl,2),size(pCasl,3),size(Hadamard_matrix,2),repetitions);
else
    dataASL=zeros(size(pCasl,1),size(pCasl,2),size(pCasl,3),1,size(pCasl,5));
end


%% Begin of calculation
switch quantification_method
    case 'ASL-diff'
        
        method = scan_acqp('##$Method=',ASL.texte,0);
        dataASL = ASL.reco.data(:,:,2,:,:) - ASL.reco.data(:,:,1,:,:);
        clipASL{1} = [min(dataASL(:)) max(dataASL(:)) 1];
        
        if strcmp(M0value,'Yes')
            dataASL = dataASL ./ ASL.reco.data(:,:,2,:,:)*100;
            clipASL{1} = [-10 30 1];
        end
        
        if AlphaVoiYn && ~isempty(strfind(method,'PhiSw'))
            pCASL_PR = scan_acqp('##$pCASL_PulseRate=',ASL.texte,1);
            pCASL_PhaseStart=scan_acqp('##$pCASL_PhaseStart=',ASL.texte,1);
            pCASL_PhaseEnd=scan_acqp('##$pCASL_PhaseEnd=',ASL.texte,1);
            pCASL_PhaseInc=scan_acqp('##$pCASL_PhaseInc=',ASL.texte,1);
            pCASL_Phase = pCASL_PhaseStart:pCASL_PhaseInc:pCASL_PhaseEnd;
            for itSlice = 1 : size(ASL.reco.data,4)
                factorResize = size(ASL.reco.data,1) / size(VoiAlpha(itSlice).value,1);
                VoiAlphaTemp = imresize(VoiAlpha(itSlice).value,factorResize,'nearest');
                dataPhiSw = squeeze(dataASL(:,:,1,itSlice,:));
                VoiAlphaTemp = reshape(VoiAlphaTemp,[1,size(VoiAlphaTemp,1)*size(VoiAlphaTemp,2)]);
                dataPhiSw = reshape(dataPhiSw,[size(dataPhiSw,1)*size(dataPhiSw,2),size(dataPhiSw,3)]);
                AlphaPhiSw = nanmean(dataPhiSw(VoiAlphaTemp,:),1);
                [Amax,indAmax]=max(AlphaPhiSw);
                hPS = figure;
                plot(pCASL_Phase,AlphaPhiSw)
                maxHz = pCASL_Phase(indAmax);
                maxDeg = maxHz * pCASL_PR*1e-6 * 360 + 360;
                title(sprintf('pCASL phase sweep \n Max = %0.1f %% for phi = %0.0f Hz or %0.1f deg',Amax,maxHz,maxDeg))
                
                patternFile = 'Image_Analyses_data';
                IndxPath = strfind(ASL_filename,patternFile);
                SavePath = ASL_filename(1:IndxPath+numel(patternFile));
                FileNameMat = ASL_filename(IndxPath+numel(patternFile)+1:end-4);
                
                PhiSwFigName = fullfile(SavePath,sprintf('pCASL_Phase_Sweep_file-%s.fig',FileNameMat));
                saveas(hPS, PhiSwFigName)
            end
        end
%         ASL_maps.reco.texte = 'ASL-diff';
    case 'ASL-Hadamard-Bolus'
        Hadamard_matrix = scan_acqp('CASL_HadamardMatrix=',ASL.texte,1);
        Subboli_Duration = scan_acqp('CASL_SubboliTime=',ASL.texte,1);
        TotalLabeling_Duration = scan_acqp('CASL_LabelTime=',ASL.texte,1);
        PostLabeling_Duration = scan_acqp('CASL_PostLabelTime=',ASL.texte,1);
        Effective_PLD = zeros(1,size(Hadamard_matrix,2));
        for subbolus = 1:size(Hadamard_matrix,2)
            Effective_PLD(subbolus) = TotalLabeling_Duration - sum(Subboli_Duration(1:subbolus)) + PostLabeling_Duration;
            for pCASL_image=1:size(Hadamard_matrix,1)
                coefficient = Hadamard_matrix(pCASL_image,subbolus);
                current_image = ASL.reco.data(:,:,pCASL_image,:,:);
                dataASL(:,:,subbolus,:,:) = dataASL(:,:,subbolus,:,:) - coefficient.* current_image ;
            end
        end
        
        dataASL= dataASL/(size(Hadamard_matrix,1)/2);
        clipASL{1} = [min(dataASL(:)) max(dataASL(:)) 1];
%         ASL_maps.reco.texte = 'ASL-Hadamard-Bolus';
    case 'Alsop1996'
        clipASL{1} = [0 500 1];
        for it_sl=1:size(ASL.reco.data,4)    % boucle sur les slices
            counter = 0;
            for it_expt=1:size(ASL.reco.data,5)    % boucle sur les experiments
                dM = squeeze(( ASL.reco.data(:,:,2,it_sl,it_expt) - ASL.reco.data(:,:,1,it_sl,it_expt) ));
                % T1 map with FAIR
                num_echo_t1_map = 1;
                T1app = T1map.reco.data(:,:,num_echo_t1_map,it_sl,1);
                if strcmp(AcqDim, '3D') == 1 % en 3D toutes les coupes sont acquises simultanemment
                    NumSlice = 0;
                else
                    % T1 apparent : map
                    NumSlice = find(PVM_ObjOrderList==it_sl-1)-1;
                end
                SliceDec = exp(NumSlice*interSliceTime./T1app);
                % M0 map approximate to control image
                M0 = ASL.reco.data(:,:,2,it_sl,it_expt);
                ObjOrderListNameM0 = sprintf('##$PVM_ObjOrderList=( %d )',size(ASL.reco.data,4));
                PVM_ObjOrderListM0=scan_acqp(ObjOrderListNameM0,ASL.texte,1)';
                InvTimeM0 = LabelTime + PostLabelTime;
                NumSliceM0 = find(PVM_ObjOrderListM0==it_sl-1)-1;
                SliceDecM0 = ((1-exp(-(InvTimeM0+NumSliceM0*interSliceTime)./T1app))).^(-1);
                % CBF in ml/100g/min
                f = dM .* SliceDec .* lambda ./ ( 2 .* alpha .* M0 .* SliceDecM0 .* T1app );
                dataASL(:,:,1,it_sl,it_expt) = f * UnitCoeff;             % passage de (ml/g/ms) en (ml/100g/min);
                counter = counter+1;
                omit(it_sl, it_expt) = 1;
            end
        end
%         ASL_maps.reco.texte = 'CBF';
        
    case {'Buxton1998','Parkes2002','HadamardFitTransitTime'}
        
        %clipASL{1} = [0 500 1];
        % boucle sur les slices
        for it_sl=1:size(pCasl,3)
            
            %% M0 map
            % M0 = value given by user
            if( isempty(M0Map) && max(~isnan(M0value))==1 && ~ischar(M0value) )
                M0 = repmat(M0value, [size(pCasl, 1) size(pCasl, 2)]);
                PVM_ObjOrderListM0 = 0;
                InvTimeM0 = Inf;
                interSliceTimeM0 = 0;
                M0SatYN=true;
                % M0 = map of mean control acquisition
            elseif isempty(M0Map) && strcmp(M0value, 'Default')
                M0 = nanmean(pCasl(:,:,it_sl,2,:),5);
                PVM_ObjOrderListM0=PVM_ObjOrderList;
                InvTimeM0 = LabelTime + PostLabelTime;
                interSliceTimeM0 = interSliceTime;
                %M0SatYN=char2logical(scan_acqp('##$PVM_FovSatOnOff=',ASL.texte,0));
                M0SatYN = strcmp(J_pCasl.FovSatOnOff.value{1}, 'On');
                %M0_TR = scan_acqp('##$PVM_RepetitionTime=',ASL.texte,1);
                M0_TR = J_pCasl.RepetitionTime.value;
                % M0 = map from a given scan (see which condition)
            elseif( ~isempty(M0Map) )
                %methodM0 = scan_acqp('##$Method=',M0map.texte,0);
                methodM0 = J_M0Map.MethodDiffusion.value{1};
%                 ObjOrderListNameM0 = sprintf('##$PVM_ObjOrderList=( %d )',size(M0map.reco.data,4));
%                 PVM_ObjOrderListM0=scan_acqp(ObjOrderListNameM0,M0map.texte,1)';
                PVM_ObjOrderListM0 = J_M0Map.ObjOrderList.value;
                %AcqDimM0 = scan_acqp('##$PVM_SpatDimEnum=',M0map.texte,2);
                AcqDimM0 = J_M0Map.SpatDimEnum.value;
                if(size(M0map,3)==1) || strcmp(AcqDimM0, '3D') == 1
                    PVM_ObjOrderListM0 = 0;
                end
                % M0 = map from FAIR scan (see which condition)
                if(~isempty(regexpi(methodM0,'\w*fair\w*')))
                    NaASL = ASL.acq.no_averages;
                    NaFair = M0map.acq.no_averages;
                    % M0 = map from fit of inversion recovery scan (FAIR)
                    if(regexpi(M0map.reco.echo_label{1},'fit_M0\w*'))
                        M0 = M0map(:,:,it_sl,1,1)*NaASL/NaFair;
                        InvTimeM0 = Inf;
                        interSliceTimeM0 = 0; % already taken into account in T1 fitting
                        M0SatYN=true;
                        % M0 = map from last echo (inversion time) of a FAIR scan
                    elseif(isempty(M0map.reco.echo_label{1}))
                        % inter slice time
                        interSliceTimeM0 = InterSliceTimeFairEpi_Adapted(M0map);
                        M0 = M0map(:,:,it_sl,1,end)*NaASL/NaFair;
                        %FairTIR_Arr=(scan_acqp('##$FairTIR_Arr=',M0map.texte,1))';
                        FairTIR_Arr = J_M0Map.FairTIR_Arr.value;
                        InvTimeM0 = FairTIR_Arr(end);
                        M0SatYN=true;
                    end
                    % M0 = map from control image of a CASL or pCASL scan
                elseif(~isempty(regexpi(methodM0,'\w*casl\w*')))
                    M0 = mean(M0map.reco.data(:,:,2,it_sl,:),5);
                    %LabelTimeM0=scan_acqp('ASL_LabelTime=',M0map.texte,1);
                    LabelTimeM0 = J_M0Map.LabelTime.value;
                    %PostLabelTimeM0=scan_acqp('ASL_PostLabelTime=',M0map.texte,1);
                    PostLabelTimeM0 = J_M0Map.PostLabelTime.value;
                    InvTimeM0 = LabelTimeM0 + PostLabelTimeM0;
                    %AcqDimM0 = scan_acqp('##$PVM_SpatDimEnum=',M0map.texte,2);
                    AcqDimM0 = J_M0Map.SpatDimEnum.value;
                    if strcmp(AcqDimM0, '3D') == 1
                        interSliceTimeM0=0;
                    else
                        interSliceTimeM0 = InterSliceTimeCaslEpi(M0map);
                    end
                    %M0SatYN=char2logical(scan_acqp('##$PVM_FovSatOnOff=',M0map.texte,0));
                    M0SatYN = J_M0Map.FovSatOnOff.value;
                    %M0_TR = scan_acqp('##$PVM_RepetitionTime=',M0map.texte,1);
                    M0_TR = J_M0Map.RepetitionTime.value;
                    % M0 = map from first echo and first experiment of selected scan
                else
                    M0 = M0map.reco.data(:,:,1,it_sl,1);
                    InvTimeM0 = Inf;
                    interSliceTimeM0 = 0;
                    M0SatYN=true;
                end
            elseif isempty( M0Map ) && ~isnumeric(M0value)
                warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ M0map and M0value are empty \n##$');
                msgbox(warning_text, 'CBF map warning') ;
                ASL_maps = [];
                return
            end
            
            %% T1 map
            % T1 apparent = constante
            if( isempty(T1Map) )
                T1app = repmat(T1value, [size(pCasl, 1) size(pCasl, 2)]);
                % T1 apparent = map
            elseif( ~isempty(T1Map))
%                 if ~isfield(J_T1Map,'echo_label')
%                      J_T1Map.echo_label.value= '';
%                 end
                T1app = T1Map(:,:,it_sl);
                T1app(isnan(T1app)) = T1value;
%                 if(strcmp(J_T1Map.echo_label.value(1),'CBF')==1)
%                     num_echo_t1_map = 3;
%                     T1app = T1Map(:,:,it_sl, num_echo_t1_map,1);
%                 elseif(strcmp(J_T1Map.echo_label.value(1),'fit_T1')==1)
%                     num_echo_t1_map = 1;
%                     T1app = T1Map(:,:,it_sl,num_echo_t1_map,1);
%                 else
%                     num_echo_t1_map = 1;
%                     T1app = T1Map(:,:,it_sl,num_echo_t1_map,1);
%                     %                     num_echo_t1_map = 5;
%                     %                     T1app = T1map.reco.data(:,:,1,it_sl,num_echo_t1_map);
%                 end
%                 
            end
            
            %% Transit delay map
            if( isempty(TransitMap) && ~ischar(TTvalue) )
                deltaTemp = TTvalue;
                delta = deltaTemp * ones(size(pCasl,1),size(pCasl,2));
            elseif( isempty(TransitMap) && strcmp(TTvalue,'Default') )
                deltaTemp = PostLabelTime;
                delta = deltaTemp * ones(size(pCasl,1),size(pCasl,2));
            elseif( ~isempty(TransitMap) )
                delta = TransitMap(:,:,it_sl,1,1);
            end
            
            %% Slice order and slice timing correction
            if(size(PVM_ObjOrderListM0)==size(PVM_ObjOrderList))
                NumSliceM0 = find(PVM_ObjOrderListM0==it_sl-1)-1;
            end
            if(size(PVM_ObjOrderListM0)==1)
                NumSliceM0 = 0;
            end
%             if size(PVM_ObjOrderListM0) == 1 
%                 NumSliceM0=0;
%             end

            %% Post labeling time correction
            if M0SatYN
                SliceDecM0 = (1-exp(-(InvTimeM0+NumSliceM0*interSliceTimeM0)./T1app)).^(-1);
            else
                SliceDecM0 = (1-exp(-M0_TR./T1app)).^(-1);
            end
            % if Buxton or Parkes
            if strcmp(quantification_method, 'HadamardFitTransitTime')~=1
                %% Calul of CBF
                NumSlice = find(PVM_ObjOrderList==it_sl-1)-1;
                if size(PVM_ObjOrderList) == 1 % cas 1 slice ou 3D
                    NumSlice=0;
                end                
                PostLabelTimeASL = PostLabelTime + NumSlice*interSliceTime;
                tau = LabelTime;        % (ms) labeling time
                omega = PostLabelTimeASL;  % (ms) transit delay
                qss1 = 1 - exp(-(tau + omega - delta)./T1app);      % terms related to the approach to the steady state factors
                qss2 = 1 - exp(-tau ./T1app);                       % terms related to the approach to the steady state factors
                
                % boucle sur les experiments
                for it_expt=1:size(pCasl,5)    % boucle sur les experiments
                    
                    dM = squeeze(( pCasl(:,:,it_sl,2,it_expt) - pCasl(:,:,it_sl,1,it_expt) ));
                    
                    dM1 = zeros(size(dM));
                    dM2 = zeros(size(dM));
                    
                    % single compartment model from Buxton 1998
                    switch quantification_method
                        case 'Buxton1998'
                            dM1(delta > omega) = dM(delta > omega);
                            dM2(delta <= omega) = dM(delta <= omega);
                            f1 = dM1 ./ ( 2 .* M0.* SliceDecM0 .* T1app .* alpha .* exp(-delta ./ T1sang) .* qss1 );
                            f2 = dM2 ./ ( 2 .* M0.* SliceDecM0 .* T1app .* alpha .* exp(-delta ./ T1sang) .* exp(-(omega-delta) ./ T1app) .* qss2 );
                            
                            % two compartment model from Parkes 2002
                        case 'Parkes2002'
                            
                            %% Parameters for 2 compartment model
                            % Expected values of parameters
                            PSExp = 155;                % (ml/100g/min) takagi Stroke 1987 and Bos Stroke 2012
                            %                         deltaExpCx = 230;         % cf. Thomas JCBFM 2006
                            %                         deltaExpSt = 275;         % cf. Thomas JCBFM 2006
                            %                         delta = deltaExpCx * ones(size(ASL.reco.data,1),size(ASL.reco.data,2));
                            %                         vbCx = 0.031;         	% (ml blood/ml tissue) entre 0.03 et 0.05 ???; cx:3.1%, str:2.8% selon valable 2008 ; str: 2.4% selon Beaumont 2009
                            %                         vbSt = 0.028;         	% (ml blood/ml tissue) entre 0.03 et 0.05 ???; cx:3.1%, str:2.8% selon valable 2008 ; str: 2.4% selon Beaumont 2009
                            %                         vb = vbCx;
                            vb = 0.03;
                            vwb = 0.7;                  % (ml water/ml blood)
                            lambda = 0.9;               % (ml water(ml tissue)-1/ml water(ml blood)-1)
                            % Other parameters
                            PS = PSExp / UnitCoeff;     % (ml water/ms/ml tissue)
                            Solution = 'slow';          % 'fast' or 'slow'
                            T1b = T1sang;               % (ms)
                            tL = LabelTime;             % (ms) labeling time
                            tW = PostLabelTime;         % (ms) postlabeling time
                            t = tL + tW;                % (ms)
                            D = 1/T1b;                  % (ms-1)
                            tA = delta;                 % (ms) transit delay
                            tp = t - tA;                % (ms)
                            vbw = vb * vwb;             % (ml water/ml tissue)
                            ma0 = M0;                   % (a.u.)
                            T1e = T1app;                % (ms)
                            
                            %% Elements of equations
                            A = PS / vbw;               % (ms-1)
                            C = 1./T1e;                 % (ms-1)
                            if(strcmp(Solution,'slow'))
                                J = A + D;              % (ms-1)
                            elseif(strcmp(Solution,'fast'))
                                J = A + D + f / vb;     % (ms-1)
                            end
                            
                            %% 2 compartment model equation
                            % dM1 = 2 * f * ma0 * alpha * exp(-D * tA) * ( (1-exp(-J*tp))/J + A * ( (J-C+C*exp(-J*tp)-J*exp(-C*tp))/(J*C*(J-C)) ) );                                  % (tA < t < tA + tL)
                            % dM2 = 2 * f * ma0 * alpha * exp(-D * tA) * ( ( 1/J + A/(J*(C-J)) ) * ( (exp(J*tL)-1) * exp(-J*tp) ) - A / (C*(C-J)) * exp(-C*tp) * (exp(C*tL)-1) );     % (t > tA + tL)
                            dM1(tA < t <= tA + tL) = dM(tA < t <= tA + tL);
                            dM2(t > tA + tL) = dM(t > tA + tL);
                            f1 = dM1 ./ ( 2 * ma0.* SliceDecM0 .* alpha .* exp(-D .* tA) .* ( (1-exp(-J.*tp))./J + A .* ( (J-C+C.*exp(-J.*tp)-J.*exp(-C.*tp))./(J.*C.*(J-C)) ) ) );
                            f2 = dM2 ./ ( 2 * ma0.* SliceDecM0 .* alpha .* exp(-D .* tA) .* ( ( 1./J + A./(J.*(C-J)) ) .* ( (exp(J.*tL)-1) .* exp(-J.*tp) ) - A ./ (C.*(C-J)) .* exp(-C.*tp) .* (exp(C.*tL)-1) ) );
                            
                    end
                    f = f1 + f2;
                    dataASL(:,:,it_sl,1,it_expt) = f * lambda * UnitCoeff;
                end
                
                
%                 ASL_maps.echo_label(1,it_sl) = {'Buxton1998'};
%                 ASL_maps.unit(1,it_sl) = {'(ml/100g/min)'};
%                 ASL_maps.reco.texte = 'CBF';
            elseif strcmp(quantification_method, 'HadamardFitTransitTime')==1
                if it_sl == 1
                    fit_transit_result = NaN(size(ASL.reco.data,1)*size(ASL.reco.data,2),size(ASL.reco.data,4));
                    %fit_transit_err= NaN(size(data_in_vector,1),size(data_in_vector,3));
                    
                    fit_CBF_result = NaN(size(ASL.reco.data,1)*size(ASL.reco.data,2),size(ASL.reco.data,4));
                    %fit_CBF_err= NaN(size(data_in_vector,1),size(data_in_vector,3));
                    
                    fit_T1app_result = NaN(size(ASL.reco.data,1)*size(ASL.reco.data,2),size(ASL.reco.data,4));
                    %fit_T1app_err= NaN(size(data_in_vector,1),size(data_in_vector,3));
                end
                
                [fit_CBF_result(:,it_sl), fit_transit_result(:,it_sl),fit_T1app_result(:,it_sl)] =ASL_HadamardFit(it_sl,ASL,VoiMask_filename,M0,SliceDecM0,T1app,PVM_ObjOrderList,...
                    interSliceTime,alpha,T1sang,lambda,progression_Hadamard,fit_CBF_result,fit_transit_result,fit_T1app_result);

%                 [fit_CBF_result(:,it_sl), fit_transit_result(:,it_sl),fit_T1app_result(:,it_sl)] =parametric_map_ASL_dyn(it_sl,ASL,VoiMask_filename,PVM_ObjOrderList,...
%                      interSliceTime,alpha,T1sang,lambda,fit_CBF_result,fit_transit_result,fit_T1app_result);
            end
        end
        
end

% Mise en forme donnees
if strcmp(quantification_method, 'HadamardFitTransitTime') ~= 1
    if MeanYn
        output_data{1} = nanmean(dataASL(:,:,:,:,1:end),5);
    else
         output_data{1} = dataASL;
    end
else
    
    output_data{1} = reshape(fit_CBF_result, [size(pCasl,1),size(pCasl,2),size(pCasl,3),1,1]) * 6000; % from ml/g/sec/ to  ml/100g/min
    output_data{2} = reshape(fit_transit_result, [size(pCasl,1),size(pCasl,2),size(pCasl,3),1,1]) * 1000; % sec to ms 
    output_data{3} =reshape(fit_T1app_result, [size(pCasl,1),size(pCasl,2),size(pCasl,3),1,1]) *1000; % sec to ms   
    
    clipASL{1} = [0 500 1];
    clipASL{2} = [0 700 1];
    clipASL{3} = [0 3500 1];
end



% for i = 1:numel(output_data)
%     ASL_maps{i}=Template;
%     if size(output_data{i},4) == 1
%         ASL_maps{i}.reco.data(:,:,1,:)=output_data{i};     
%     else
%         ASL_maps{i}.reco.data = output_data{i}; 
%     end
%     ASL_maps{i}.reco.no_echoes = size(ASL_maps{i}.reco.data,3);
%     ASL_maps{i}.reco.no_expts   = size(ASL_maps{i}.reco.data,5);
%     
%     ASL_maps{i}.reco.displayedecho=ASL.reco.no_echoes;
%     for m_expt=1:ASL_maps{i}.reco.no_expts,
%         for m_slice=1:ASL_maps{i}.reco.no_slices,
%             for m_echo=1:ASL_maps{i}.reco.no_echoes
%                 ASL_maps{i}.reco.fov_offsets(:,m_echo,m_slice,m_expt) = ASL.reco.fov_offsets(:,1,m_slice,m_expt);
%                 ASL_maps{i}.reco.fov_orientation(:,m_echo,m_slice,m_expt) = ASL.reco.fov_orientation(:,1,m_slice,m_expt);
%                 ASL_maps{i}.reco.label(m_echo,m_slice,m_expt) = ASL.reco.label(1,m_slice,m_expt);
%                 ASL_maps{i}.reco.phaselabel(m_echo,m_slice,m_expt) = ASL.reco.phaselabel(1,m_slice,m_expt);
%                 ASL_maps{i}.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = ASL.reco.fov_phase_orientation(1,m_slice,m_expt);
%                 ASL_maps{i}.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
%                 ASL_maps{i}.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
%             end
%         end
%     end
% 
%     ASL_maps{i}.reco.globalmin=min(ASL_maps{i}.reco.data(:));
%     ASL_maps{i}.reco.globalmax=max(ASL_maps{i}.reco.data(:)); 
% %     ASL_maps{i}.echo_label(1,it_sl) = {'Buxton1998'};
%     ASL_maps{i}.unit = {'(ml/100g/min)'};
%     ASL_maps{i}.reco.texte = 'CBF';
% 
%     ASL_maps{i}.reco=orderfields(ASL_maps{i}.reco);
%     
%     ASL_maps{i}.clip=clipASL{i};
%     
% end







OutputImages = output_data{1};
OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 10000) = NaN;
%OutputImages(isnan(OutputImages)) = NaN;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
end

%% Initial function 


%% 




% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
% info2.ImageSize = ones(1, length(info2.ImageSize));
% info2.ImageSize(1:length(size(OutputImages_reoriented))) = size(OutputImages_reoriented);
info2.ImageSize = size(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));


% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% % so far copy the .json file of the first input
% copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% % 


%% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);

J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

end

