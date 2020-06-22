function [files_in,files_out,opt] = Module_vMRF_Rest_1D(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the module's parameters with default values

if isempty(opt)
    
    %   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'dictionary_folder_filename',  'Dictionary Folder'};
    module_option(:,2)   = {'combUsed', 'Post/Pre'};
    module_option(:,3)   = {'indivNorm', 'No'};
    module_option(:,4)   = {'finalNorm', 'Yes'};
    module_option(:,5)   = {'prefix',           'MRF_Rest_'};
    module_option(:,6)   = {'method',           'DBM'};
    module_option(:,7)   = {'filtered',         'No'};
    
    module_option(:,8)   = {'RefInput',         1};
    module_option(:,9)   = {'InputToReshape',   1};
    module_option(:,10)   = {'Table_in',         table()};
    module_option(:,11)   = {'Table_out',        table()};
    module_option(:,12)   = {'folder',           table()};
    module_option(:,13)  = {'OutputSequenceName','AllName'};
    module_option(:,14)  = {'Params',           'Vf'};
    module_option(:,15)  = {'K',                50};
    module_option(:,16)  = {'Lw',               0};
    module_option(:,17)  = {'cstrS',            'd'};
    module_option(:,18)  = {'cstrG',            'd'};
    module_option(:,19)  = {'RelErr', 0.05};
    module_option(:,20)   = {'RestType', 'T2'};
    module_option(:,21)  = {'mkScoreMap', 'Yes'};
    module_option(:,22)  = {'removed',          0};
    module_option(:,23)  = {'augment',          60};
    opt.Module_settings  = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    
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
        'This module allows to perform vascular MRF with different signals (MGEFIDSE, MSME), combined as different patterns (ratio or concatenations)'
        'A quantitative map of either ADC or T2 is fed to the module to help the dictionary matching'
        ''
        'The dictionary-based matching (DBM) method is based on the paper: Ma, Dan, et al. "Magnetic resonance fingerprinting." Nature (2013)'
        ['This specific vascular MRF method is based on : Christen T. et al. "MR vascular fingerprinting: ' ...
            'A new approach to compute cerebral blood volume, mean vessel radius, and oxygenation maps in the human brain" NeuroImage (2014)']
        'The dictionary-based learning (DBL) method is based on the paper: Boux, Fabien, et al. [work in progress]'
        ''
        'Prerequisite:'
        '      - Put your dictionary files (pre and post simulated scans and/or MSME) in the ''data/dictionaries'' folder'
        '      - (or) Put your ''DICO.mat'' dictionary file ratio between the post and pre simulated scans in the ''data/dictionaries'' folder'
        '      - (or performing the regression method) Put your ''MODEL.mat'' model file'
        ''
        'The dictionaries are designed with the Mrvox simulation tool'
        }'};
    
    user_parameter(:,2)   = {'Select the MGEFIDSE Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the MGEFIDSE Post scan','1Scan','','',{'SequenceName'}, '',''};
    user_parameter(:,4)   = {'Select the MSME scan','1Scan','','',{'SequenceName'}, '',''};
    user_parameter(:,5)   = {'Select the sequence combination to use', 'cell', {'Post/Pre','Pre-Post','MSME-Post/Pre','MSME-Pre-Post'}, 'combUsed', '','Mandatory', 'Combination of signals to use'};
    user_parameter(:,6)   = {'Normalize each signals', 'cell', {'Yes', 'No'}, 'indivNorm', '','Optional', 'Do you want to normalize (N2) each signal before concatenation? NB: MGEFIDSE are not normalized in case of Ratio'};
    user_parameter(:,7)   = {'Normalize final signal', 'cell', {'Yes', 'No'}, 'finalNorm', '','Optional', 'Do you want to normalize (N2) the final combination before MRF?'};
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    folder_files	= dir(fullfile(s{1}, 'data/dictionaries/'));
    folder_files    = folder_files([folder_files.isdir]);
    opt.Module_settings.folder = fullfile(s{1}, 'data/dictionaries/');
    if isempty(folder_files), folder_files(1).name = ' '; end
    user_parameter(:,8)   = {'Dictionary folder','cell', {folder_files.name}, 'dictionary_folder_filename','','Mandatory',...
        {'Select the folder containing dico files (.json), ratio dico file (.mat) and/or model file (.mat)'}};
    
    user_parameter(:,9)   = {'Prefix','char', '', 'prefix', '', '',...
        {'Choose a prefix for your output maps'}};
    user_parameter(:,10)   = {'Parameters','check', ...
        {'Vf', 'VSI', 'R', 'SO2', 'DH2O', 'B0theta', 'khi', 'Hct', 'T2'},...
        'Params', '', '',...
        {'Select the parameters considered in the model (default ''Vf'')'
        }'};
    user_parameter(:,11)   = {'Remove last echoes?','numeric','', 'removed', '', '',...
        {'WARNING: Only applies to MGEFIDSE signals'}};
    user_parameter(:,12)   = {'Smooth?','cell', {'Yes','No'}, 'filtered', '', '',...
        {'Select ''Yes'' to smooth the signals  (recommanded ''No'')'}};
    user_parameter(:,13)   = {'Make ScoreMap?','cell', {'Yes','No'}, 'mkScoreMap', '', '',...
        {'Select ''Yes'' add score map in outputs  (recommanded ''Yes'')'; 'WARNING: does not work with DBL'}'};
    user_parameter(:,14)   = {'Select the scans to use for restriction','1Scan','','',{'SequenceName'}, '','A T2 map should be in ms, ADC map in mm²/s'};
    user_parameter(:,15)   = {'Indicate the type of restriction','cell', {'T2', 'DH2O'}, 'RestType', '', '',...
        {}'};
    user_parameter(:,16)   = {'Relative error tolerated on prior input value','numeric','','RelErr','','',''};
    user_parameter(:,17)   = {'Method','cell', {'DBM', 'DBL'}, 'method', '', '',...
        { 'Choose:'
        '	- ''DBM'' to use the Dan Ma method'
        '	- ''DBL'' to use the regression method'
        }'};
    user_parameter(:,18)   = {'   .Model settings (if the DBL method is chosen)','Text','','','','',...
        {'Recommanded:'
        'K = 50'
        'Lw = 0'
        'cstr on Sigma = ''d*'''
        'cstr on Gamma = '' '''
        'Dictionary augmentation = 60'
        }'};
    user_parameter(:,19)   = {'       .Number of regions','numeric','','K','','',...
        {'Recommanded: K = 50'
        'If K is -1, an automatic tuning of the parameter is performed using BIC (time-consuming)'
        }'};
    user_parameter(:,20)  = {'       .Number of additional unsupervised parameter','numeric','','Lw','','',...
        {'Recommanded: Kw = 0'
        'If Lw is -1, an automatic tuning of the parameter is performed using BIC (time-consuming)'
        }'};
    user_parameter(:,21)  = {'       .Model constraint on Sigma','cell',{'i*','i','d*','d',' '},'cstrS','','',...
        {'''d'' = diagonal'
        '''i'' = isotropic'
        '''*'' = equal for all K regions'
        }'};
    user_parameter(:,22)  = {'       .Model constraint on Gamma','cell',{'i*','i','d*','d',' '},'cstrG','','',...
        {'''d'' = diagonal'
        '''i'' = isotropic'
        '''*'' = equal for all K regions'
        }'};
    user_parameter(:,23)   = {'       .Dictionary augmentation?','numeric','', 'augment', '', '',...
        {'Recommanded: = 60'
         'This value represents the signal-to-noise ratio (SNR) between dictionary signals and noise added on these signals'
         '0 and Inf values lead to no augmentation by noise addition'
         }};    
    
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)', 'VariableNames', VariableNames);
    %%
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1    = {''};
    files_out.In1   = {''};
    return
    
end
%%%%%%%%

if isempty(files_out)
    
    opt.Params = opt.Params(cell2mat(opt.Params(:,2)),1);
    
    opt.Table_out = opt.Table_in(1,:);
    
    for i = 1:numel(opt.Params)
        opt.Table_out(i,:) = opt.Table_out(1,:);
        opt.Table_out(i,:).Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            opt.Table_out.SequenceName(i) = categorical(cellstr([char(opt.prefix), char(opt.Params{i})]));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            opt.Table_out.SequenceName(i) = categorical(cellstr([char(opt.Table_out.SequenceName(i)), opt.Params{i}]));
        end
        opt.Table_out.Filename(i) = categorical(cellstr([char(opt.Table_out.Patient(i)), '_', char(opt.Table_out.Tp(i)), '_', char(opt.Table_out.SequenceName(i))]));
        opt.Table_out.IsRaw(i) = categorical(cellstr('0'));
        files_out.In1{i} = [char(opt.Table_out.Path(i)), char(opt.Table_out.Filename(i)) '.nii'];
    end
    
    if strcmp(opt.method,'DBL')
        nb = i;
        
        for i = 1:numel(opt.Params)
            opt.Table_out(nb+i,:) = opt.Table_out(1,:);
            
            opt.Table_out.Filename(nb+i) =  categorical(cellstr([char(opt.Table_out.Patient(i)), '_', char(opt.Table_out.Tp(i)), '_', char(opt.Table_out.SequenceName(i)) '_confidence']));
            opt.Table_out.SequenceName(nb+i) = categorical(cellstr([char(opt.prefix), char(opt.Params{i}), '_confidence']));
            
            files_out.In2{i} = [char(opt.Table_out.Path(nb+i)), char(opt.Table_out.Filename(nb+i)) '.nii'];
        end
    end
    
    if strcmp(opt.mkScoreMap, 'Yes')
        scoreRow.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            scoreRow.SequenceName = categorical(cellstr([char(opt.prefix), 'ScoreMap']));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            scoreRow.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), 'ScoreMap']));
        end
        scoreRow.Filename = categorical(cellstr([char(opt.Table_out.Patient(1)), '_', char(opt.Table_out.Tp(1)), '_', char(scoreRow.SequenceName)]));
        scoreRow.IsRaw = categorical(cellstr('0'));
        
        scoreRow.Group = opt.Table_out.Group(1);
        scoreRow.Patient = opt.Table_out.Patient(1);
        scoreRow.Tp = opt.Table_out.Tp(1);
        scoreRow.Type = opt.Table_out.Type(1);
        opt.Table_out = [opt.Table_out; struct2table(scoreRow)];
        
        files_out.In1{end+1} = [char(opt.Table_out.Path(end)), char(opt.Table_out.Filename(end)) '.nii'];
    end
end


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('MRF_MultiSeq:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
opt.finalNorm = strcmp(opt.finalNorm, 'Yes');

% Read json files and create ratio dictionary
opt.dictionary_folder_filename = fullfile(opt.folder, opt.dictionary_folder_filename);

d = dir(opt.dictionary_folder_filename);
if contains(opt.combUsed, 'Pre')
    opt.dictionary_pre_filename     = d(contains({d.name}, 'PRE_'));
    if ~isempty(opt.dictionary_pre_filename)
        opt.dictionary_pre_filename     = opt.dictionary_pre_filename.name;
    end
    opt.dictionary_post_filename    = d(contains({d.name}, 'POST_'));
    if ~isempty(opt.dictionary_post_filename)
        opt.dictionary_post_filename    = opt.dictionary_post_filename.name;
    end
end

if contains(opt.combUsed, 'MSME')
    opt.dictionary_MSME_filename    = d(contains({d.name}, 'MSME'));
    if ~isempty(opt.dictionary_MSME_filename)
        opt.dictionary_MSME_filename    = opt.dictionary_MSME_filename.name;
    end
end

% Name the model using a random ID and check if this model already exist
if strcmp(opt.method, 'DBL')
    model_filename = [opt.dictionary_folder_filename filesep 'MODEL_' num2str(10^5+randi(10^6-10^5,1)) '.mat'];

    list_models = d(contains({d.name}, 'MODEL_'));
    flag_model_exist = exist(model_filename,'file');
    flag_valid = false(size(list_models));
    for i = 1:length(list_models)

        load([opt.dictionary_folder_filename filesep list_models(i).name],'Parameters')
        
        try
            flag_valid(i) = strcmp(opt.cstrG, Parameters.cstr.Gammat);
            flag_valid(i) = flag_valid(i) && strcmp(opt.cstrS, Parameters.cstr.Sigma);
            flag_valid(i) = flag_valid(i) && opt.K == Parameters.K;
            flag_valid(i) = flag_valid(i) && opt.Lw == Parameters.Lw;
            flag_valid(i) = flag_valid(i) && opt.augment == Parameters.data_augmentation;
        catch
            flag_valid(i) = 0;
        end
    end
    flag_model_exist = flag_model_exist || any(flag_valid);
    if flag_model_exist
        l = find(flag_valid == true);
        model_filename = [opt.dictionary_folder_filename filesep list_models(l(1)).name];
    end
end

if (strcmp(opt.method, 'DBL') && ~exist(model_filename,'file')) || strcmp(opt.method, 'DBM')
    
    % Dictionary filename
    dico_filename   = [opt.dictionary_folder_filename filesep 'DICO.mat'];
    
    % If, dico exists, load it, else, create it and save it
    if exist(dico_filename,'file')
        load(dico_filename)
    else
        switch opt.combUsed
            case {'Post/Pre', 'Pre-Post'}
                Pre     = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
                Post    = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_post_filename]);
                
%                 Pre     = ReadJson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
%                 Pre.MRSignals = reshape(Pre.MRSignals.x_ArrayData_(:,1) + 1i*Pre.MRSignals.x_ArrayData_(:,2),...
%                     [Pre.MRSignals.x_ArraySize_(1), Pre.MRSignals.x_ArraySize_(2)]);
%                 
%                 Post     = ReadJson([opt.dictionary_folder_filename filesep opt.dictionary_post_filename]);
%                 Post.MRSignals = reshape(Post.MRSignals.x_ArrayData_(:,1) + 1i*Post.MRSignals.x_ArrayData_(:,2),...
%                     [Post.MRSignals.x_ArraySize_(1), Post.MRSignals.x_ArraySize_(2)]);
                
                Dico.MRSignals{1}       = abs(Pre.MRSignals);
                Dico.MRSignals{2}       = abs(Post.MRSignals);
                Dico.Tacq               = Pre.Sequence.Tacq;
                Dico.Parameters.Par     = Pre.Parameters.Par; % Parameters used to simulate X signals
%                 Dico.Parameters.Par     = ReplaceNaNCell(Pre.Parameters.Par);
                Dico.Parameters.Labels  = Pre.Parameters.Labels;
                clear Pre Post
                save(dico_filename,'Dico')
            case {'MSME-Post/Pre','MSME-Pre-Post'}
                Pre     = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
                Post    = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_post_filename]);
                
%                 Pre     = ReadJson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
%                 Pre.MRSignals = reshape(Pre.MRSignals.x_ArrayData_(:,1) + 1i*Pre.MRSignals.x_ArrayData_(:,2),...
%                     [Pre.MRSignals.x_ArraySize_(1), Pre.MRSignals.x_ArraySize_(2)]);
%                 
%                 Post     = ReadJson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
%                 Post.MRSignals = reshape(Post.MRSignals.x_ArrayData_(:,1) + 1i*Post.MRSignals.x_ArrayData_(:,2),...
%                     [Post.MRSignals.x_ArraySize_(1), Post.MRSignals.x_ArraySize_(2)]);
                
                MSME    = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_MSME_filename]);
                
%                 MSME     = ReadJson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
%                 MSME.MRSignals = reshape(MSME.MRSignals.x_ArrayData_(:,1) + 1i*MSME.MRSignals.x_ArrayData_(:,2),...
%                     [MSME.MRSignals.x_ArraySize_(1), MSME.MRSignals.x_ArraySize_(2)]);
                
                Dico.MRSignals{1}       = abs(Pre.MRSignals);
                Dico.MRSignals{2}       = abs(Post.MRSignals);
                Dico.MRSignals{3}       = abs(MSME.MRSignals);
                Dico.Tacq{1}            = Pre.Sequence.Tacq;
                Dico.Tacq{2}            = MSME.Sequence.Tacq;
                Dico.Parameters.Par     = Pre.Parameters.Par; % Parameters used to simulate X signals
%                 Dico.Parameters.Par     = ReplaceNaNCell(Pre.Parameters.Par);
                Dico.Parameters.Labels  = Pre.Parameters.Labels;
                clear Pre Post MSME
                save(dico_filename,'Dico')
        end
    end
end

% After loading or creating the 'raw' dico, operations are still needed to
% adapt to the combination of sequences

switch opt.combUsed
    case 'Post/Pre'
        % Generate ratio signals from scans (and TODO: ROI if given)
        % TODO: what if In1 is the post and In2 the pre scan
        Xobs            = niftiread(files_in.In2{1}) ./ niftiread(files_in.In1{1});
        json_filename   = split(files_in.In2{1},'.');
        json_filename{2} = '.json';
        Obs             = ReadJson([json_filename{1} json_filename{2}]);
        
        % If necessary smooth observations
        if strcmp(opt.filtered, 'Yes') == 1
            p = 2;
            for x = 1:size(Xobs,1); for y = 1:size(Xobs,2); for z = 1:size(Xobs,3)
                        signal = squeeze(Xobs(x,y,z,:));
                        signal = filter(ones(1, p)/p, 1, [signal(1); signal; signal(end)]);
                        %     	Xobs(x,y,z,:) = conv(signal, gaussian_window, 'valid');
                        Xobs(x,y,z,:) = signal(2:end-1);
                    end; end; end
        end
        
        % If necessary, remove echoes
        if opt.removed > 0
            clear tmp
            for x = 1:size(Xobs,1); for y = 1:size(Xobs,2); for z = 1:size(Xobs,3)
                tmp(x,y,z,:) = Xobs(x,y,z,1:end-opt.removed);
            end; end; end
            Xobs    = tmp;
            Obs.EchoTime.value = Obs.EchoTime.value(1:end-opt.removed);
        end
        
        Xobs        = permute(Xobs, [1 2 4 3]);
        
        % Reformat dico (not needed if MODEL is already computed)
        if strcmp(opt.method, 'DBM') || ~exist(model_filename,'file')
            tmp = nan(size(Dico.MRSignals{1},1), length(Obs.EchoTime.value'));
            if size(Xobs,3) ~= size(Dico.MRSignals{1},2)
                warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
                for i = 1:size(Dico.MRSignals{1},1)
                    tmp(i,:) = interp1(Dico.Tacq(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{2}(i,:)./Dico.MRSignals{1}(i,:), Obs.EchoTime.value'*1e-3);
                end
            else
                tmp = Dico.MRSignals{2}./Dico.MRSignals{1};
            end
            Dico.MRSignals = tmp;
            Dico.Tacq   = Obs.EchoTime.value'*1e-3;
            %remove row containning nan values
            nn = ~any(isnan(Dico.MRSignals),2);
            Dico.MRSignals = Dico.MRSignals(nn,:);
            Dico.Parameters.Par = Dico.Parameters.Par(nn,:);
            %             Tmp{1}      = Dico;
        end
               
    case 'Pre-Post'
        XobsPre             = niftiread(files_in.In1{1});
        XobsPost            = niftiread(files_in.In2{1});
        json_filename       = split(files_in.In2{1},'.');
        json_filename{2}    = '.json';
        Obs                 = ReadJson([json_filename{1} json_filename{2}]);
        
        % If necessary smooth observations
        if strcmp(opt.filtered, 'Yes') == 1
            p = 2;
            for x = 1:size(XobsPre,1); for y = 1:size(XobsPre,2); for z = 1:size(XobsPre,3)
                        % Pre
                        signalPre = squeeze(XobsPre(x,y,z,:));
                        signalPre = filter(ones(1, p)/p, 1, [signalPre(1); signalPre; signalPre(end)]);
                        XobsPre(x,y,z,:) = signalPre(2:end-1);
                        % Post
                        signalPost = squeeze(XobsPost(x,y,z,:));
                        signalPost = filter(ones(1, p)/p, 1, [signalPost(1); signalPost; signalPost(end)]);
                        XobsPost(x,y,z,:) = signalPost(2:end-1);
                    end; end; end
        end
        
        if strcmp(opt.indivNorm, 'Yes') % Normalize nifti signals
            for x = 1:size(XobsPre,1); for y = 1:size(XobsPre,2); for z = 1:size(XobsPre,3)
                        XobsPre(x,y,z,:) = XobsPre(x,y,z,:)./(sqrt(sum(XobsPre(x,y,z,:).^2)));
                        XobsPost(x,y,z,:) = XobsPost(x,y,z,:)./(sqrt(sum(XobsPost(x,y,z,:).^2)));
                    end; end; end
        end
        
        % If necessary, remove echoes
        if opt.removed > 0
            for x = 1:size(XobsPre,1); for y = 1:size(XobsPre,2); for z = 1:size(XobsPre,3)
                clear tmpPre tmpPost
                tmpPre(x,y,z,:) = XobsPre(x,y,z,1:end-opt.removed);
                tmpPost(x,y,z,:) = XobsPost(x,y,z,1:end-opt.removed);
            end; end; end
            XobsPre    = tmpPre;
            XobsPost    = tmpPost;
            Obs.EchoTime.value = Obs.EchoTime.value(1:end-opt.removed);
        end
        
        timeDim         = 4;%find(size(XobsPre)==length(Obs.EchoTime.value));
        Xobs            = cat(timeDim, XobsPre, XobsPost);
        Xobs            = permute(Xobs, [1 2 4 3]);
        XobsPre         = permute(XobsPre, [1 2 4 3]);
        
        % Reformat dico (not needed if MODEL is already computed)
        if strcmp(opt.method, 'DBM') || ~exist(model_filename,'file')
            tmpPre = nan(size(Dico.MRSignals{1},1), length(Obs.EchoTime.value'));
            tmpPost = nan(size(Dico.MRSignals{1},1), length(Obs.EchoTime.value'));
            if size(XobsPre,3) ~= size(Dico.MRSignals{1},2)
                warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
                for i = 1:size(Dico.MRSignals{1},1)
                    tmpPre(i,:) = interp1(Dico.Tacq(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{1}(i,:), Obs.EchoTime.value'*1e-3);
%                     if strcmp(opt.indivNorm, 'Yes') % Normalize dico
%                         tmpPre(i,:) = tmpPre(i,:)./(sqrt(sum(tmpPre(i,:).^2)));
%                     end
                    tmpPost(i,:) = interp1(Dico.Tacq(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{2}(i,:), Obs.EchoTime.value'*1e-3);
%                     if strcmp(opt.indivNorm, 'Yes') % Normalize dico
%                         tmpPost(i,:) = tmpPost(i,:)./(sqrt(sum(tmpPost(i,:).^2)));
%                     end
                end
                if strcmp(opt.indivNorm, 'Yes') % Normalize dico
                    tmpPre = tmpPre./vecnorm(tmpPre,2,2);
                    tmpPost = tmpPost./vecnorm(tmpPost,2,2);
                end
                tmp = cat(2, tmpPre, tmpPost);
            else
                if strcmp(opt.indivNorm, 'Yes') % Normalize dico
                    Dico.MRSignals{1} = Dico.MRSignals{1}./vecnorm(Dico.MRSignals{1},2,2);
                    Dico.MRSignals{2} = Dico.MRSignals{2}./vecnorm(Dico.MRSignals{2},2,2);
                end             
                tmp = cat(2, Dico.MRSignals{1}, Dico.MRSignals{2});
            end
            
            Dico.MRSignals      = tmp;
            clear tmp;
            Dico.Tacq           = Obs.EchoTime.value'*1e-3;
            %remove row containning nan values
            nn                  = ~any(isnan(Dico.MRSignals),2);
            Dico.MRSignals      = Dico.MRSignals(nn,:);
            Dico.Parameters.Par = Dico.Parameters.Par(nn,:);
            %Tmp{1}              = Dico;
        end
        %end
        
    case 'MSME-Pre-Post'
        XobsPre             = niftiread(files_in.In1{1});
        XobsPost            = niftiread(files_in.In2{1});
        XobsMSME            = niftiread(files_in.In3{1});
        json_filename       = split(files_in.In2{1},'.');
        json_filename{2}    = '.json';
        Obs                 = ReadJson([json_filename{1} json_filename{2}]);
        
        % If necessary smooth observations
        if strcmp(opt.filtered, 'Yes') == 1
            p = 2;
            for x = 1:size(XobsPre,1); for y = 1:size(XobsPre,2); for z = 1:size(XobsPre,3)
                        % Pre
                        signalPre = squeeze(XobsPre(x,y,z,:));
                        signalPre = filter(ones(1, p)/p, 1, [signalPre(1); signalPre; signalPre(end)]);
                        XobsPre(x,y,z,:) = signalPre(2:end-1);
                        % Post
                        signalPost = squeeze(XobsPost(x,y,z,:));
                        signalPost = filter(ones(1, p)/p, 1, [signalPost(1); signalPost; signalPost(end)]);
                        XobsPost(x,y,z,:) = signalPost(2:end-1);
                        % MSME
                        signalMSME = squeeze(XobsMSME(x,y,z,:));
                        signalMSME = filter(ones(1, p)/p, 1, [signalMSME(1); signalMSME; signalMSME(end)]);
                        XobsMSME(x,y,z,:) = signalMSME(2:end-1);
                    end; end; end
        end
        if strcmp(opt.indivNorm, 'Yes') % Normalize nifti signals
            for x = 1:size(XobsPre,1); for y = 1:size(XobsPre,2); for z = 1:size(XobsPre,3)
                        XobsPre(x,y,z,:) = XobsPre(x,y,z,:)./(sqrt(sum(XobsPre(x,y,z,:).^2)));
                        XobsPost(x,y,z,:) = XobsPost(x,y,z,:)./(sqrt(sum(XobsPost(x,y,z,:).^2)));
                        XobsMSME(x,y,z,:) = XobsMSME(x,y,z,:)./(sqrt(sum(XobsMSME(x,y,z,:).^2)));
                    end; end; end
        end
        
        % If necessary, remove echoes
        if opt.removed > 0
            for x = 1:size(XobsPre,1); for y = 1:size(XobsPre,2); for z = 1:size(XobsPre,3)
                clear tmpPre tmpPost
                tmpPre(x,y,z,:) = XobsPre(x,y,z,1:end-opt.removed);
                tmpPost(x,y,z,:) = XobsPost(x,y,z,1:end-opt.removed);
%                 tmpMSME(x,y,z,:) = XobsMSME(x,y,z,1:end-opt.removed);
            end; end; end
            XobsPre     = tmpPre;
            XobsPost    = tmpPost;
%             XobsMSME    = tmpMSME;
            Obs.EchoTime.value = Obs.EchoTime.value(1:end-opt.removed);
        end
        
        timeDim         = find(size(XobsPre)==length(Obs.EchoTime.value));
        Xobs            = cat(timeDim, XobsMSME, XobsPre, XobsPost);
        Xobs            = permute(Xobs, [1 2 4 3]);
        XobsPre         = permute(XobsPre, [1 2 4 3]);
        
        % Reformat dico (not needed if MODEL is already computed)
        if strcmp(opt.method, 'DBM') || ~exist(model_filename,'file')
            tmpPre = nan(size(Dico.MRSignals{1},1), length(Obs.EchoTime.value'));
            tmpPost = nan(size(Dico.MRSignals{1},1), length(Obs.EchoTime.value'));
            
            if size(XobsPre,length(size(XobsPre))) ~= size(Dico.MRSignals{1},2)
                warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
                for i = 1:size(Dico.MRSignals{1},1)
                    tmpPre(i,:) = interp1(Dico.Tacq{1}(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{1}(i,:), Obs.EchoTime.value'*1e-3);
                    tmpPost(i,:) = interp1(Dico.Tacq{1}(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{2}(i,:), Obs.EchoTime.value'*1e-3);
                end
            end
            if strcmp(opt.indivNorm, 'Yes') % Normalize dico
                for i = 1:size(Dico.MRSignals{1},1)
                    tmpPre(i,:) = tmpPre(i,:)./(sqrt(sum(tmpPre(i,:).^2)));
                    tmpPost(i,:) = tmpPost(i,:)./(sqrt(sum(tmpPost(i,:).^2)));
                    Dico.MRSignals{3}(i,:) = Dico.MRSignals{3}(i,:)./(sqrt(sum(Dico.MRSignals{3}(i,:).^2)));
                end
            end
            tmp                 = cat(2, Dico.MRSignals{3}, tmpPre, tmpPost);
            Dico.MRSignals      = tmp;
            Dico.Tacq           = Obs.EchoTime.value'*1e-3;
            %remove row containning nan values
            nn                  = ~any(isnan(Dico.MRSignals),2);
            Dico.MRSignals      = Dico.MRSignals(nn,:);
            Dico.Parameters.Par = Dico.Parameters.Par(nn,:);
            %             Tmp{1}              = Dico;
        end
        
    case 'MSME-Post/Pre'
        Xobs                = niftiread(files_in.In2{1}) ./ niftiread(files_in.In1{1});
        XobsMSME            = niftiread(files_in.In3{1});
        json_filename       = split(files_in.In2{1},'.');
        json_filename{2}    = '.json';
        Obs                 = ReadJson([json_filename{1} json_filename{2}]);
        
        % If necessary smooth observations
        if strcmp(opt.filtered, 'Yes') == 1
            p = 2;
            for x = 1:size(Xobs,1); for y = 1:size(Xobs,2); for z = 1:size(Xobs,3)
                        % Ratio
                        signal = squeeze(Xobs(x,y,z,:));
                        signal = filter(ones(1, p)/p, 1, [signal(1); signal; signal(end)]);
                        Xobs(x,y,z,:) = signal(2:end-1);
                        % MSME
                        signalMSME = squeeze(XobsMSME(x,y,z,:));
                        signalMSME = filter(ones(1, p)/p, 1, [signalMSME(1); signalMSME; signalMSME(end)]);
                        XobsMSME(x,y,z,:) = signalMSME(2:end-1);
                    end; end; end
        end
        
        if strcmp(opt.indivNorm, 'Yes') % Normalize nifti signals
            for x = 1:size(Xobs,1); for y = 1:size(Xobs,2); for z = 1:size(Xobs,3)
                        XobsMSME(x,y,z,:) = XobsMSME(x,y,z,:)./(sqrt(sum(XobsMSME(x,y,z,:).^2)));
                    end; end; end
        end
        
         % If necessary, remove echoes
        if opt.removed > 0
            clear tmp
            for x = 1:size(Xobs,1); for y = 1:size(Xobs,2); for z = 1:size(Xobs,3)
                tmp(x,y,z,:) = Xobs(x,y,z,1:end-opt.removed);
            end; end; end
            Xobs    = tmp;
            Obs.EchoTime.value = Obs.EchoTime.value(1:end-opt.removed);
        end
        
        timeDim         = find(size(Xobs)==length(Obs.EchoTime.value));
        sizeLastDim     = size(Xobs,length(size(XobsRatio)));
%         XobsRatio       = Xobs; % create copy for dimension
        Xobs            = cat(timeDim, XobsMSME, Xobs);
        Xobs            = permute(Xobs, [1 2 4 3]);
        
        
        % Reformat dico (not needed if MODEL is already computed)
        if strcmp(opt.method, 'DBM') || ~exist(model_filename,'file')
            tmp = nan(size(Dico.MRSignals{1},1), length(Obs.EchoTime.value'));
            if sizeLastDim ~= size(Dico.MRSignals{1},2)
                warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
                for i = 1:size(Dico.MRSignals{1},1)
                    tmp(i,:) = interp1(Dico.Tacq(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{2}(i,:)./Dico.MRSignals{1}(i,:), Obs.EchoTime.value'*1e-3);
                end
            end
            if strcmp(opt.indivNorm, 'Yes') % Normalize dico
                for i = 1:size(Dico.MRSignals{1},1)
                    Dico.MRSignals{3}(i,:) = Dico.MRSignals{3}(i,:)./(sqrt(sum(Dico.MRSignals{3}(i,:).^2)));
                end
            end
            tmp                 = cat(2, Dico.MRSignals{3}, tmp);
            Dico.MRSignals      = tmp;
            Dico.Tacq           = Obs.EchoTime.value'*1e-3;
            %remove row containning nan values
            nn                  = ~any(isnan(Dico.MRSignals),2);
            Dico.MRSignals      = Dico.MRSignals(nn,:);
            Dico.Parameters.Par = Dico.Parameters.Par(nn,:);
            %             Tmp{1}              = Dico;
        end
        
end %end switch

%% Dico restriction
% Check input restriction maps
% inPar = {'T2', 'DH20'};
% inMaps = [0,0];
% if isfield(files_in, 'In4')
%     inMaps(1)=1;
% end
% if isfield(files_in, 'In5')
%     inMaps(2)=1;
% end
switch opt.RestType
    case 'T2'
        inPar = 'T2';
        Fact = 1e3;
    case 'DH2O'
        inPar = 'DH2O';
        Fact = 1e12;
end

%% Iterate on restriction maps
% for i =1:numel(inMaps)
%     if inMaps(i)
%% Reduce the number of search
colNb               = find(contains(cellfun(@char, Dico.Parameters.Labels, 'UniformOutput', 0), inPar));
%inNumber        = sprintf('In%i', 4);
MapPath             = files_in.In4;
RestMap             = niftiread(MapPath{1});
RestMap             = round(RestMap); % round T2 values to reduce the amount of restrictions to perform
Values              = unique(RestMap); % get unique values
minDicoVal          = min(Dico.Parameters.Par(:,colNb))*Fact; % min of the dico
maxDicoVal          = max(Dico.Parameters.Par(:,colNb))*Fact; % max of the dico
Values(Values*(1+opt.RelErr) < minDicoVal) = nan;
Values(Values*(1-opt.RelErr) > maxDicoVal) = nan;
if nnz(isnan(Values))
    warning('%i voxels will not be evaluated as their %s value falls outside the dictionary range', nnz(isnan(Values)), inPar)
end
%         [nanRow, nanCol, nanSl] = ind2sub(size(RestMap), find(RestMap < minDicoVal & RestMap > maxDicoVal)); % Get coordinates of points that will not fit the dico
%% Dico Restriction
for v=1:numel(Values)
    % Removing dico entries where parameter of interest is out of
    % the tolerated range
    if isnan(Values(v)) %If value at this iteration is nan, don't consider it
        continue
    end
    % if not nan, get coordinates of corresponding voxels
    [row, col, sl] = ind2sub(size(RestMap), find(RestMap == Values(v))); % Get coordinates of voxels considered at this iteration

    % remove dico entries out of tolerated range
    toRemoveInf     = Dico.Parameters.Par(:,colNb)*Fact < Values(v)*(1-opt.RelErr);
    toRemoveSup     = Dico.Parameters.Par(:,colNb)*Fact > Values(v)*(1+opt.RelErr);
    toKeep          = ~(toRemoveInf + toRemoveSup);
    
    if nnz(toKeep) == 0
        warning('Value %i could not be evaluated as restricted dico is empty, %i voxels concerned\n', Values(v), numel(row))
        continue
    end
    
    % Copying the restricted dico
    Tmp{1}.MRSignals = Dico.MRSignals(toKeep, :);
    Tmp{1}.Parameters.Par = Dico.Parameters.Par(toKeep, :);
    Tmp{1}.Parameters.Labels = Dico.Parameters.Labels;
    
    localXobs = zeros(numel(row), size(Xobs, 3));
    for k=1:numel(row)
        localXobs(k,:) = Xobs(row(k), col(k), :, sl(k)); % Keep only the corresponding voxels in the observation
    end
    %% Compute MRF/regression
    switch opt.method
        case 'DBM'
            % TODO: find something nicer than this permute trick
            Estimation  = AnalyzeMRImages_Chunk(localXobs,Tmp,opt.method, [], [], [], opt.finalNorm);
            Map.Y       = permute(Estimation.GridSearch.Y, [1 2 4 3]);
            
        case 'DBL'
            %Compute the learing only one time per dictionar
            if exist(model_filename,'file')
                load(model_filename,'Parameters','labels');
                Estimation = AnalyzeMRImages(Xobs, [], opt.method, Parameters);
                
                Tmp{1}.Parameters.Labels = labels;
            else
                count = 1;
                
                for z = 1:length(Dico.Parameters.Labels)
                    tmp = split(Dico.Parameters.Labels{z},'.',2);
                    if any(strcmp(tmp{end}, opt.Params))
                        params.Par(:,count)     = Tmp{1}.Parameters.Par(:,z);
                        params.Labels{count}    = Tmp{1}.Parameters.Labels{z};
                        count = count +1;
                    end
                end
                Tmp{1}.Parameters = params;
                
                % Parameters of the regression
                clear Parameters
                if opt.K >= 0,  Parameters.K = opt.K; end
                if opt.Lw >= 0, Parameters.Lw = opt.Lw; end
                if strcmp(opt.cstrS,' '), opt.cstrS = ''; end
                if strcmp(opt.cstrG,' '), opt.cstrG = ''; end
                Parameters.cstr.Sigma   = opt.cstrS;
                Parameters.cstr.Gammat  = opt.cstrG;
                Parameters.cstr.Gammaw  = '';
                
                [Estimation, Parameters] = AnalyzeMRImages(Xobs, Tmp, opt.method, Parameters);
                
                labels      = Tmp{1}.Parameters.Labels;
                save(model_filename,'Parameters', 'labels')
            end
            
            Map.Y   = permute(Estimation.Regression.Y, [1 2 4 3]);
            Map.Std	= permute(Estimation.Regression.Cov, [1 2 4 3]).^.5;
    end
    %% Extract maps (and modify unit if necessary)
    count = 1;
    for k = 1:length(Tmp{1}.Parameters.Labels)
        tmp = split(Tmp{1}.Parameters.Labels{k},'.',2);
        if any(strcmp(tmp{end}, opt.Params))
            
            Labels{count} = tmp{end};
            
            switch tmp{end}
                %If the ression method is performed, extract also the confidence maps
                case {'Vf', 'SO2'} %convert to percent
                    for j = 1:numel(row)
                        MapStruct{count}(row(j), col(j), sl(j))     = 100*Map.Y(j,:,:,k);
                        if strcmp(opt.method, 'DBL')
                            StdStruct{count}(row(j), col(j), sl(j)) = 100*Map.Std(j,:,:,k);
                        end
                    end
                case {'VSI', 'R'} % convert m to µm
                    for j = 1:numel(row)
                        MapStruct{count}(row(j), col(j), sl(j))     = 1e6*Map.Y(j,:,:,k);
                        if strcmp(opt.method, 'DBL')
                            StdStruct{count}(row(j), col(j), sl(j))	= 1e6*Map.Std(j,:,:,k);
                        end
                    end
                case 'T2' % convert s to ms
                    for j = 1:numel(row)
                        MapStruct{count}(row(j), col(j), sl(j))     = 1e3*Map.Y(j,:,:,k);
                        if strcmp(opt.method, 'DBL')
                            StdStruct{count}(row(j), col(j), sl(j))	= 1e6*Map.Std(j,:,:,k);
                        end
                    end
                otherwise
                    for j = 1:numel(row)
                        MapStruct{count}(row(j), col(j), sl(j)) 	= Map.Y(j,:,:,k);
                        if strcmp(opt.method, 'DBL')
                            StdStruct{count}(row(j), col(j), sl(j)) = 100*Map.Std(j,:,:,k);
                        end
                    end
            end   % end switch tmp
            count   = count +1;
        end %end if any
    end % end for k in labels
    
    for j = 1:numel(row)
        ScoreMap(row(j), col(j), sl(j)) = Estimation.scoreMap(j); %%%%%%%%%%%%%%%% GERER LA DIMENSION DE MAPSCORE
    end
        
end % end for(Values)
%     end % end if(inMaps)
%
% end % end for(inMaps)

%% Ensure output dimension is correct
if any(size(MapStruct{1}) ~= size(RestMap))
    [s1, s2, s3] = size(RestMap);
    for i = 1:numel(MapStruct)
        MapStruct{i}(s1, s2, s3) = 0;
    end
end

%% Put NaN where out of the dico
for i = 1:numel(MapStruct)
    MapStruct{i}(isnan(Values))= nan;
end

if strcmp(opt.mkScoreMap, 'Yes')
    MapStruct{end+1} = ScoreMap;
    Labels{end+1} = 'ScoreMap';
end

%% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J       = ReadJson(jsonfile);
J       = KeepModuleHistory(J, struct('files_in',files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

% Reoriented and save nifti maps
% nifti_header = spm_vol(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
info.Filemoddate = char(datetime('now'));

for i = 1:length(files_out.In1)
    
    for j = 1:numel(Labels)
        
        if contains(files_out.In1{i},Labels{j})
            [path, name, ~] = fileparts(files_out.In1{i});
            WriteJson(J, [path, '/', name, '.json'])
            
            info.Filename = files_out.In1{i};
            info.ImageSize = size(MapStruct{j});
            info.PixelDimensions = info.PixelDimensions(1:length(size(MapStruct{j})));
            info.Datatype = class(MapStruct{j});
            niftiwrite(MapStruct{j}, files_out.In1{i}, info);
            
            
            if strcmp(opt.method, 'DBL')
                [path, name, ~] = fileparts(files_out.In2{i});
                WriteJson(J, [path, '/', name, '.json'])
                
                info.Filename = files_out.In2{i};
                info.ImageSize = size(StdStruct{j});
                info.PixelDimensions = info.PixelDimensions(1:length(size(StdStruct{j})));
                info.Datatype = class(StdStruct{j});
                niftiwrite(StdStruct{j}, files_out.In2{i}, info);
            end
        end
    end
end

