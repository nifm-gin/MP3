function ASL_maps = parametric_map_CBF_ASL(InvEff_ROI_filename, VoiMask_filename, ASL_InvEff_filename, T1map_filename, M0map_filename, TTmap_filename, ASL_filename, add_parameters)

% generate a CBF from an ASL scan
% this code comes from the the ASL module coded by C. Debacker
% Adaptation Paravision 6 et gestion ASL 3D : L. Hirschler

progression_Hadamard = 0;
%% Retrieve parameters
lambda = 0.9;                                            % brain-blood partition coefficient of water (ml/g) from biblio
quantification_method = extract_param(add_parameters,4);   % Choice of quantification method
if strcmp(quantification_method, 'ASL-Hadamard-Bolus')
    alpha = 'Default';
    ASL_InvEff_filename = '';
else
    if strcmp(extract_param(add_parameters,2), 'None')      % inversion efficiency
        if strcmp(extract_param(add_parameters,3), 'Default')
            alpha = 'Default';
        else
            alpha = extract_param(add_parameters,3);
        end
    else
        alpha = extract_param(add_parameters,2);
         if isempty(ASL_InvEff_filename)
                warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ No ASL_InvEff file found \n##$ For %s', ASL_filename);
                msgbox(warning_text, 'CBF map warning') ;
             ASL_maps = [];
             return
         end
    end
end
T1sang = 'Default';               % T1 du sang (ms)

M0value =  extract_param(add_parameters,10);                % M0 value used (theoriticaly M0 of blood)
T1value = extract_param(add_parameters,12);                 % Apparent tissue T1
TTvalue = extract_param(add_parameters,14);                 % Transit time between labeling plane and imaging plane
ExcludeOutlier = char2logical(add_parameters{:}(7));        % Run script for remove outlier image of ASL scan
if strcmp(extract_param(add_parameters,6), 'Mean')          % Mean of CBF accros expt at exit of program
    MeanYn = 1;
else
    MeanYn = 0;
end


%% Load maps and VOIs
finalMap = 'CBF';
ListMap = {'ASL','T1map','M0map','TTmap'};

% Function : filename, map return by this file, type of file (VOI or map)
[ASL, ~] = extract_map(ASL_filename,finalMap,'map');               % load the ASL file
[T1map, T1map_status] = extract_map(T1map_filename,finalMap,'map');         % load the T1map file
[M0map, M0map_status] = extract_map(M0map_filename,finalMap,'map');         % load the M0map file
[TTmap, TTmap_status] = extract_map(TTmap_filename,finalMap,'map');         % load the TTmap file
[VoiAlpha, AlphaVoiYn] = extract_map(InvEff_ROI_filename,finalMap,'VOI');     % load ROI for alpha if necessary

if isempty(M0map)
    VoiM0 = [];
    M0map_filename = '';
end
if isempty(T1map)
    VoiT1 = [];
end

%% Compare compatibility between two map : give name of variable
Status = [];
for it = 1:numel(ListMap)-1
    [StatusTemp,ValidSlicesTemp] = check_compatibility(eval(ListMap{1}),eval(ListMap{it+1}));
    
    if StatusTemp == 0
        ASL_maps = [];
        return
    end
    eval(sprintf('%s.reco.data = %s.reco.data(:,:,:,%s==1,:);',ListMap{1},ListMap{1},'ValidSlicesTemp{1}'));
    eval(sprintf('%s.acq.no_slices = %i;',ListMap{1},sum(ValidSlicesTemp{1})));
    if eval(sprintf('~isempty(%s_filename)',ListMap{it+1}))
        eval(sprintf('%s.reco.data = %s.reco.data(:,:,:,%s==1,:);',ListMap{it+1},ListMap{it+1},'ValidSlicesTemp{2}'));
        eval(sprintf('%s.acq.no_slices = %i;',ListMap{it+1},sum(ValidSlicesTemp{2})));
    end
    Status = cat(1,Status,StatusTemp);
end

% Exit of program if an error is detected
if( isempty(ASL) || ~T1map_status || ~M0map_status || ~TTmap_status || ~AlphaVoiYn || min(Status)==0 )
    eval(sprintf('%s = []',finalMap));
    return
end

% Use of VOI : true or false
if isempty(InvEff_ROI_filename)
    AlphaVoiYn = false;
end

%% Prefix for CASL or pCASL sequence
if( ~isempty(regexpi(ASL.acq.ppl_name,'\w*casl\w*')) && isempty(regexpi(ASL.acq.ppl_name,'\w*pcasl\w*')) )
    ASLprefix = 'C';
else
    ASLprefix = 'pC';
end

%% Alpha value : inversion efficiency
% Load Alpha scan if necessary
AlphaScanMethod = 'Not loaded';
AlphaMethodNeed = 'ASL_InvEff';
[AlphaScan,AlphaScan_status] = extract_map(ASL_InvEff_filename,finalMap,'map');

if AlphaScan_status && ~isempty(AlphaScan)
    AlphaScanMethod = scan_acqp('##$QuantifMethod=',AlphaScan.reco.paramQuantif,0);
end
% Test if alpha measure is possible
% Test if we have Voi and Scan for alpha
if ~isempty( ASL_InvEff_filename ) && ~AlphaVoiYn && strcmp(extract_param(add_parameters,3), 'Default')
    warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ You must enter a valid VOI name for Alpha value \n##$');
    msgbox(warning_text, 'CBF map warning') ;
    ASL_maps = [];
    return
    % Test if we have Voi and Scan for alpha
elseif (~strcmp(AlphaScanMethod,AlphaMethodNeed) && ~isempty( ASL_InvEff_filename )) || (isempty( ASL_InvEff_filename ) && AlphaVoiYn )
    warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ Scan for alpha must be ASL_InvEff type\n##$AlphaScan=%s\n##$',...
        ASL_InvEff_filename);
    msgbox(warning_text, 'CBF map warning') ;
    ASL_maps = [];
    return
elseif ~AlphaVoiYn && ischar(alpha) && ~strcmp('Default',alpha)
    warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ Voi don t exist with this pattern\n##$AlphaScan=%s\n##$',...
        ASL_InvEff_filename);
    msgbox(warning_text, 'CBF map warning') ;
    ASL_maps = [];
    return
end

% Find good slice of Voi for Alpha
if AlphaVoiYn
    AlphaScanOS = AlphaScan.acq.fov_offsets(:,3);
    NumVoiA = [];
    for itVoi=1:numel(VoiAlpha)
        VoiAlphaOS = VoiAlpha(itVoi).fov_offsets(3);
        if VoiAlphaOS==AlphaScanOS
            NumVoiA = itVoi;
        end
    end
    
    if isempty(NumVoiA)
        warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ Voi offset and scan offset different \n##$AlphaScan=%s\n##$',...
            ASL_InvEff_filename);
        msgbox(warning_text, 'CBF map warning') ;
        ASL_maps = [];
        return
    end
    
    VoiAlpha = VoiAlpha(NumVoiA);
    if exist('VoiAlpha2','var')
        VoiAlpha2 = VoiAlpha2(NumVoiA);
    end
end

% Alpha = value specified by user (nothing to do)
% Alpha = value from ASL_GEFC scan in a VOI specified by AlphaValue
if( ~isempty(ASL_InvEff_filename) && AlphaVoiYn )
    if isempty(VoiAlpha)
        VoiAlpha.value = true(size(AlphaScan.reco.data));
    end
    factorResize = size(AlphaScan.reco.data,1) / size(VoiAlpha.value,1);
    VoiAlphaTemp = logical(imresize(VoiAlpha.value,factorResize,'nearest'));
    if exist('VoiAlpha2','var')
        VoiAlphaTemp = VoiAlphaTemp + logical(imresize(VoiAlpha2.value,factorResize,'nearest'));
    end
    AlphaScanTemp = squeeze(AlphaScan.reco.data(:,:,1,1,1));
    alpha = nanmean(AlphaScanTemp(logical(VoiAlphaTemp))) / 100;                 % /100 because ASL_GEFC give a value of alpha in percent
end

%% Set default value according to B0 magnetic field of scanner
B0 = ASL.acq.champ_magnetique;
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
if strcmp(T1sang,'Default'), T1sang = T1sangD; end
if strcmp(T1value,'Default'), T1value = T1valueD; end

%% Apply script of exclusion of outlier
if size(ASL.reco.data,5) < 2 && ExcludeOutlier
    warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ The data haven''t several repetition but you ask to repair outlier, may be remove outlier detection \n##$ASL_Scan=%s\n##$',...
        ASL_filename);
    msgbox(warning_text, 'CBF map warning') ;
    ASL_maps = [];
    return
elseif ExcludeOutlier
    alphaT = 1e-6;
    interpolation = false;
    plotDiagnostics = char2logical(add_parameters{:}{8});
    detrend = false;
    MeanMap = true;
    
    dASL = ASL;
    if strcmp(quantification_method, 'ASL-Hadamard-Bolus')
        repair_outlier_slices({ASL_filename},{VoiMask_filename},alphaT,interpolation,plotDiagnostics,detrend,MeanMap);
    else
        dASL.reco.data = ASL.reco.data(:,:,2,:,:) - ASL.reco.data(:,:,1,:,:);
        output = repair_outlier_slices({ASL_filename,dASL},{VoiMask_filename},alphaT,interpolation,plotDiagnostics,detrend,MeanMap);
        if MeanYn
            ASL = output{2};
        else
            ASL = output{1};
        end
    end
end

% AcqDim = scan_acqp('##$PVM_SpatDimEnum=',ASL.texte,2);
AcqDim = scan_acqp('##$PVM_SpatDimEnum=',ASL.texte,0);


%% Initialisation parameter of ASL
% retrieve parameter values
LabelTime=scan_acqp('ASL_LabelTime=',ASL.texte,1);
PostLabelTime=scan_acqp('ASL_PostLabelTime=',ASL.texte,1);
% Calcul du temps interslice
if strcmp(AcqDim, '3D') == 1
    interSliceTime=0;
else
    interSliceTime = InterSliceTimeCaslEpi(ASL);
end
% Ordre d'acquisition des slices
if strcmp(AcqDim, '3D') == 1 || (size(ASL.reco.data,4)==1)
    PVM_ObjOrderList = 0;
else
    ObjOrderListName = sprintf('##$PVM_ObjOrderList=( %d )',size(ASL.reco.data,4));
    PVM_ObjOrderList=scan_acqp(ObjOrderListName,ASL.texte,1)';
end
UnitCoeff = 60*100*1000;    % Conversion coefficient (ml blood/min/100ml tissue) en (ml blood/ms/ml tissue)

%% Initialized variables
if strcmp(quantification_method, 'ASL-Hadamard-Bolus')
    repetitions = size(ASL.reco.data,5);
    Hadamard_matrix = scan_acqp('CASL_HadamardMatrix=',ASL.texte,1);
    dataASL=zeros(size(ASL.reco.data,1),size(ASL.reco.data,2),size(Hadamard_matrix,2),size(ASL.reco.data,4),repetitions);
else
    dataASL=zeros(size(ASL.reco.data,1),size(ASL.reco.data,2),1,size(ASL.reco.data,4),size(ASL.reco.data,5));
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
        
        clipASL{1} = [0 500 1];
        % boucle sur les slices
        for it_sl=1:size(ASL.reco.data,4)
            
            %% M0 map
            % M0 = value given by user
            if( isempty(M0map_filename) && max(~isnan(M0value))==1 && ~ischar(M0value) )
                M0 = repmat(M0value, [size(ASL.reco.data, 1) size(ASL.reco.data, 2)]);
                PVM_ObjOrderListM0 = 0;
                InvTimeM0 = Inf;
                interSliceTimeM0 = 0;
                M0SatYN=true;
                % M0 = map of mean control acquisition
            elseif isempty(M0map_filename) && strcmp(M0value, 'Default')
                M0 = nanmean(ASL.reco.data(:,:,2,it_sl,:),5);
                PVM_ObjOrderListM0=PVM_ObjOrderList;
                InvTimeM0 = LabelTime + PostLabelTime;
                interSliceTimeM0 = interSliceTime;
                M0SatYN=char2logical(scan_acqp('##$PVM_FovSatOnOff=',ASL.texte,0));
                M0_TR = scan_acqp('##$PVM_RepetitionTime=',ASL.texte,1);
                % M0 = map from a given scan (see which condition)
            elseif( ~isempty(M0map_filename) )
                methodM0 = scan_acqp('##$Method=',M0map.texte,0);
                ObjOrderListNameM0 = sprintf('##$PVM_ObjOrderList=( %d )',size(M0map.reco.data,4));
                PVM_ObjOrderListM0=scan_acqp(ObjOrderListNameM0,M0map.texte,1)';
                AcqDimM0 = scan_acqp('##$PVM_SpatDimEnum=',M0map.texte,2);
                if(size(M0map.reco.data,4)==1) || strcmp(AcqDimM0, '3D') == 1
                    PVM_ObjOrderListM0 = 0;
                end
                % M0 = map from FAIR scan (see which condition)
                if(~isempty(regexpi(methodM0,'\w*fair\w*')))
                    NaASL = ASL.acq.no_averages;
                    NaFair = M0map.acq.no_averages;
                    % M0 = map from fit of inversion recovery scan (FAIR)
                    if(regexpi(M0map.reco.echo_label{1},'fit_M0\w*'))
                        M0 = M0map.reco.data(:,:,1,it_sl,1)*NaASL/NaFair;
                        InvTimeM0 = Inf;
                        interSliceTimeM0 = 0; % already taken into account in T1 fitting
                        M0SatYN=true;
                        % M0 = map from last echo (inversion time) of a FAIR scan
                    elseif(isempty(M0map.reco.echo_label{1}))
                        % inter slice time
                        interSliceTimeM0 = InterSliceTimeFairEpi(M0map);
                        M0 = M0map.reco.data(:,:,1,it_sl,end)*NaASL/NaFair;
                        FairTIR_Arr=(scan_acqp('##$FairTIR_Arr=',M0map.texte,1))';
                        InvTimeM0 = FairTIR_Arr(end);
                        M0SatYN=true;
                    end
                    % M0 = map from control image of a CASL or pCASL scan
                elseif(~isempty(regexpi(methodM0,'\w*casl\w*')))
                    M0 = mean(M0map.reco.data(:,:,2,it_sl,:),5);
                    LabelTimeM0=scan_acqp('ASL_LabelTime=',M0map.texte,1);
                    PostLabelTimeM0=scan_acqp('ASL_PostLabelTime=',M0map.texte,1);
                    InvTimeM0 = LabelTimeM0 + PostLabelTimeM0;
                    AcqDimM0 = scan_acqp('##$PVM_SpatDimEnum=',M0map.texte,2);
                    if strcmp(AcqDimM0, '3D') == 1
                        interSliceTimeM0=0;
                    else
                        interSliceTimeM0 = InterSliceTimeCaslEpi(M0map);
                    end
                    M0SatYN=char2logical(scan_acqp('##$PVM_FovSatOnOff=',M0map.texte,0));
                    M0_TR = scan_acqp('##$PVM_RepetitionTime=',M0map.texte,1);
                    % M0 = map from first echo and first experiment of selected scan
                else
                    M0 = M0map.reco.data(:,:,1,it_sl,1);
                    InvTimeM0 = Inf;
                    interSliceTimeM0 = 0;
                    M0SatYN=true;
                end
            elseif isempty( M0map_filename ) && ~isnumeric(M0value)
                warning_text = sprintf('##$ Can not calculate the CBF map because there is\n##$ Somthing wrong with the data\n##$ M0map and M0value are empty \n##$');
                msgbox(warning_text, 'CBF map warning') ;
                ASL_maps = [];
                return
            end
            
            %% T1 map
            % T1 apparent = constante
            if( isempty(T1map_filename) )
                T1app = repmat(T1value, [size(ASL.reco.data, 1) size(ASL.reco.data, 2)]);
                % T1 apparent = map
            elseif( ~isempty(T1map_filename) )
                if ~isfield(T1map.reco,'echo_label')
                     T1map.reco.echo_label= '';
                end
                if(strcmp(T1map.reco.echo_label(1),'CBF')==1)
                    num_echo_t1_map = 3;
                    T1app = T1map.reco.data(:,:,num_echo_t1_map,it_sl,1);
                elseif(strcmp(T1map.reco.echo_label(1),'fit_T1')==1)
                    num_echo_t1_map = 1;
                    T1app = T1map.reco.data(:,:,num_echo_t1_map,it_sl,1);
                else
                    num_echo_t1_map = 1;
                    T1app = T1map.reco.data(:,:,num_echo_t1_map,it_sl,1);
                    %                     num_echo_t1_map = 5;
                    %                     T1app = T1map.reco.data(:,:,1,it_sl,num_echo_t1_map);
                end
                
            end
            
            %% Transit delay map
            if( isempty(TTmap_filename) && ~ischar(TTvalue) )
                deltaTemp = TTvalue;
                delta = deltaTemp * ones(size(ASL.reco.data,1),size(ASL.reco.data,2));
            elseif( isempty(TTmap_filename) && strcmp(TTvalue,'Default') )
                deltaTemp = PostLabelTime;
                delta = deltaTemp * ones(size(ASL.reco.data,1),size(ASL.reco.data,2));
            elseif( ~isempty(TTmap_filename) )
                delta = TTmap.reco.data(:,:,1,it_sl,1);
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
                for it_expt=1:size(ASL.reco.data,5)    % boucle sur les experiments
                    
                    dM = squeeze(( ASL.reco.data(:,:,2,it_sl,it_expt) - ASL.reco.data(:,:,1,it_sl,it_expt) ));
                    
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
                    dataASL(:,:,1,it_sl,it_expt) = f * lambda * UnitCoeff;
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
    
    output_data{1} = reshape(fit_CBF_result, [size(ASL.reco.data,1),size(ASL.reco.data,2),1,size(ASL.reco.data,4),1]) * 6000; % from ml/g/sec/ to  ml/100g/min
    output_data{2} = reshape(fit_transit_result, [size(ASL.reco.data,1),size(ASL.reco.data,2),1,size(ASL.reco.data,4),1]) * 1000; % sec to ms 
    output_data{3} =reshape(fit_T1app_result, [size(ASL.reco.data,1),size(ASL.reco.data,2),1,size(ASL.reco.data,4),1]) *1000; % sec to ms   
    
    clipASL{1} = [0 500 1];
    clipASL{2} = [0 700 1];
    clipASL{3} = [0 3500 1];
end

Template.acq=ASL.acq;
Template.filename=ASL.filename;
Template.texte=ASL.texte;
if isfield(ASL.reco, 'iminfos')
    
    Raw_data_info = ASL.reco.iminfos{:};
else
    Raw_data_info = ASL.reco.paramQuantif ;
    
end
ParamConfig=sprintf('##$QuantifMethod=%s\n##$t1_map=%i\n##$t1_value=%i\n##$M0_map=%i\n##$M0_value=%s\n##$TT_map=%i\n##$TT_value=%s\n##$AlphaScan=%i\n##$Alpha_value=%s\n##$Scan=%d\n##$T1sang=%d\n##$lambda=%d\n##$T1map=%s\n##$M0map=%s\n##$TTmap=%s\n##$%sASL=%s\n##$%sASL scan info\n%s\n##$Outlier=%s\n##END=',...
    quantification_method,...
    ~isempty(T1map_filename),...
    T1value,...
    ~isempty(M0map_filename),...
    num2str(M0value),...
    ~isempty(TTmap_filename),...
    num2str(TTvalue),...
    ~isempty(ASL_InvEff_filename),...
    num2str(alpha),...
    ASL.scan_number,...
    T1sang,lambda,...
    T1map_filename, M0map_filename, TTmap_filename, ASLprefix, ASL_filename,...
    ASLprefix, Raw_data_info,...
    num2str(ExcludeOutlier));


if ~isempty(T1map_filename)
    ParamConfig=sprintf('%s\n\n##$T1map scan info %s',ParamConfig, T1map.reco.paramQuantif);
end
%% define a template 
Template.reco.paramQuantif = ParamConfig;
Template.reco.date = date;
Template.reco_number = ASL.reco_number;
Template.scan_number = ASL.scan_number;
Template.reco.no_slices = ASL.reco.no_slices;
if isfield(ASL, 'reco_number')
    Template.reco_number=ASL.reco_number;
end
if isfield(ASL, 'scan_number')
    Template.scan_number=ASL.scan_number;
end
Template.reco.displayedslice=ASL.reco.no_slices;
Template.reco.displayedexpt=1;
Template.reco.thickness=ASL.acq.thickness;
Template.reco.no_views=ASL.reco.no_views;
Template.reco.no_samples=ASL.reco.no_samples;
Template.reco.fov=ASL.reco.fov;
Template.reco.angulation=ASL.reco.angulation;
Template.reco.bitpix=ASL.reco.bitpix;


for i = 1:numel(output_data)
    ASL_maps{i}=Template;
    if size(output_data{i},4) == 1
        ASL_maps{i}.reco.data(:,:,1,:)=output_data{i};     
    else
        ASL_maps{i}.reco.data = output_data{i}; 
    end
    ASL_maps{i}.reco.no_echoes = size(ASL_maps{i}.reco.data,3);
    ASL_maps{i}.reco.no_expts   = size(ASL_maps{i}.reco.data,5);
    
    ASL_maps{i}.reco.displayedecho=ASL.reco.no_echoes;
    for m_expt=1:ASL_maps{i}.reco.no_expts,
        for m_slice=1:ASL_maps{i}.reco.no_slices,
            for m_echo=1:ASL_maps{i}.reco.no_echoes
                ASL_maps{i}.reco.fov_offsets(:,m_echo,m_slice,m_expt) = ASL.reco.fov_offsets(:,1,m_slice,m_expt);
                ASL_maps{i}.reco.fov_orientation(:,m_echo,m_slice,m_expt) = ASL.reco.fov_orientation(:,1,m_slice,m_expt);
                ASL_maps{i}.reco.label(m_echo,m_slice,m_expt) = ASL.reco.label(1,m_slice,m_expt);
                ASL_maps{i}.reco.phaselabel(m_echo,m_slice,m_expt) = ASL.reco.phaselabel(1,m_slice,m_expt);
                ASL_maps{i}.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = ASL.reco.fov_phase_orientation(1,m_slice,m_expt);
                ASL_maps{i}.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
                ASL_maps{i}.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
            end
        end
    end

    ASL_maps{i}.reco.globalmin=min(ASL_maps{i}.reco.data(:));
    ASL_maps{i}.reco.globalmax=max(ASL_maps{i}.reco.data(:)); 
%     ASL_maps{i}.echo_label(1,it_sl) = {'Buxton1998'};
    ASL_maps{i}.unit = {'(ml/100g/min)'};
    ASL_maps{i}.reco.texte = 'CBF';

    ASL_maps{i}.reco=orderfields(ASL_maps{i}.reco);
    
    ASL_maps{i}.clip=clipASL{i};
    
end

%% Function for extract good value of parameters
    function Param = extract_param(add_parameters,Number)
        ParamD = str2double(add_parameters{:}(Number));
        ParamS = add_parameters{:}(Number);
        if ~strcmp(ParamS{1},'NaN') && isnan(ParamD)
            Param = ParamS{1};
        else
            Param = ParamD;
        end
        
    end

end
