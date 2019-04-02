function [files_in,files_out,opt] = Module_MGEFIDSE_ratio(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

	%   % define every option needed to run this module
	% --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'dictionary_folder_filename',  'Dictionary Folder'};
    module_option(:,2)   = {'prefix',           'MRF_'};
    module_option(:,3)   = {'method',           'ClassicMRF'};
    module_option(:,4)   = {'filtered',         'No'};
    
    module_option(:,5)   = {'RefInput',         1};
    module_option(:,6)   = {'InputToReshape',   1};
    module_option(:,7)   = {'Table_in',         table()};
    module_option(:,8)   = {'Table_out',        table()};
    module_option(:,9)   = {'folder',           table()};
    module_option(:,10)  = {'OutputSequenceName','AllName'};
    module_option(:,11)  = {'Params',           'Vf'};
    
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
        'The MRF method is based on the paper : Ma, Dan, et al. "Magnetic resonance fingerprinting." Nature (2013)'
        'The regression method is based on the paper : Boux, Fabien, et al. [work in progress]'
        ''
        'Prerequisite:'
        '      - Put your ''PRE_*.json'' and ''POST_*'' dictionary files (pre and post simulated scans) in the ''data/dictionaries'' folder'
        '      - (or) Put your ''DICO.mat'' dictionary file ratio between the post and pre simulated scans in the ''data/dictionaries'' folder'
        '      - (or performing the regression method) Put your ''MODEL.mat'' model file'
        ''
        'The dictionaries are designed with the Mrvox simulation tool'
        }'};
    
    user_parameter(:,2)   = {'Select the MGEFIDSE Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the MGEFIDSE Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    folder_files	= dir(fullfile(s{1}, 'data/dictionaries/'));
    folder_files    = folder_files([folder_files.isdir]);
    opt.Module_settings.folder = fullfile(s{1}, 'data/dictionaries/');
    if isempty(folder_files), folder_files(1).name = ' '; end
    user_parameter(:,4)   = {'   .Dictionary Pre/Post file folder','cell', {folder_files.name}, 'dictionary_folder_filename','','Mandatory',...
        {'Select the folder containing Pre/Post dico files (.json), ratio dico file (.mat) and/or model file (.mat)'}};
    
    user_parameter(:,5)   = {'   .Prefix','char', '', 'prefix', '', '',...
        {'Choose a prefix for your output maps'}};
    user_parameter(:,6)   = {'   .Smooth?','cell', {'Yes','No'}, 'filtered', '', '',...
        {'Select ''Yes'' to smooth the signals  (recommanded ''No'')'}};
    user_parameter(:,7)   = {'   .Method','cell', {'ClassicMRF', 'RegressionMRF'}, 'method', '', '',...
        { 'Choose:'
        '	- ''ClassicMRF'' to use the Dan Ma method'
        '	- ''RegressionMRF'' to use the regression method'
        }'};
    
    user_parameter(:,8)   = {'   .Parameters','check', ...
        {'Vf', 'VSI', 'R', 'SO2', 'DH2O', 'B0theta', 'khi', 'Hct'},...
        'Params', '', '',...
        {'Select the parameters considered in the model (default ''Vf'')'
        }'};
    
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
    
    if strcmp(opt.method,'RegressionMRF')
        nb = i;
        
        for i = 1:numel(opt.Params)
            opt.Table_out(nb+i,:) = opt.Table_out(1,:);

            opt.Table_out.Filename(nb+i) =  categorical(cellstr([char(opt.Table_out.Patient(i)), '_', char(opt.Table_out.Tp(i)), '_', char(opt.Table_out.SequenceName(i)) '_confidence'])); 
            opt.Table_out.SequenceName(nb+i) = categorical(cellstr([char(opt.prefix), char(opt.Params{i}), '_confidence']));

            files_out.In2{i} = [char(opt.Table_out.Path(nb+i)), char(opt.Table_out.Filename(nb+i)) '.nii'];
        end
    end
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


% Read json files and create ratio dictionary
opt.dictionary_folder_filename = fullfile(opt.folder, opt.dictionary_folder_filename);

d   = dir(opt.dictionary_folder_filename);
opt.dictionary_pre_filename     = d(contains({d.name}, 'PRE_'));
if ~isempty(opt.dictionary_pre_filename)
    opt.dictionary_pre_filename     = opt.dictionary_pre_filename.name;
end
opt.dictionary_post_filename    = d(contains({d.name}, 'POST_'));
if ~isempty(opt.dictionary_post_filename)
    opt.dictionary_post_filename    = opt.dictionary_post_filename.name;
end

dico_filename   = [opt.dictionary_folder_filename filesep 'DICO.mat'];
model_filename  = [opt.dictionary_folder_filename filesep 'MODEL.mat'];


if (strcmp(opt.method, 'RegressionMRF') && ~exist(model_filename,'file')) || strcmp(opt.method, 'ClassicMRF')
    
    % If, dico exists, load it, else, create it and save it
    if exist(dico_filename,'file')
        load(dico_filename)
    else
        Pre     = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
        Post    = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_post_filename]);

        Dico.MRSignals      = abs(Post.MRSignals ./ Pre.MRSignals); % Ratio post/pre signals 
        Dico.Tacq           = Pre.Sequence.Tacq;
        Dico.Parameters.Par = Pre.Parameters.Par; % Parameters used to simulate X signals
        Dico.Parameters.Labels = Pre.Parameters.Labels;
        clear Pre Post
        save(dico_filename,'Dico')
    end
end
        

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
Xobs        = permute(Xobs, [1 2 4 3]);


% Reformat dico (not needed if MODEL is already computed)
if strcmp(opt.method, 'ClassicMRF') || ~exist(model_filename,'file')
    tmp = nan(size(Dico.MRSignals,1), length(Obs.EchoTime.value'));
    if size(Xobs,length(size(Xobs))) ~= size(Dico.MRSignals,2)
        warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
        for i = 1:size(Dico.MRSignals,1)
            tmp(i,:) = interp1(Dico.Tacq(1:size(Dico.MRSignals,2)), Dico.MRSignals(i,:), Obs.EchoTime.value'*1e-3);
        end
    end
    Dico.MRSignals = tmp;
    Dico.Tacq   = Obs.EchoTime.value'*1e-3;
    %remove row containning nan values
    nn = ~any(isnan(Dico.MRSignals),2);
    Dico.MRSignals = Dico.MRSignals(nn,:);
    Dico.Parameters.Par = Dico.Parameters.Par(nn,:);
    Tmp{1}      = Dico;
end


% Compute MRF/regression
switch opt.method
    
    case 'ClassicMRF'
        % TODO: find something nicer than this permute trick
        Estimation  = AnalyzeMRImages(Xobs,Tmp,opt.method,[]);
        Map.Y       = permute(Estimation.GridSearch.Y, [1 2 4 3]);
        
    case 'RegressionMRF'
        
        %Compute the learing only one time per dictionar        
        if exist(model_filename,'file')
            load(model_filename,'Parameters','labels');
            Estimation = AnalyzeMRImages(Xobs, [], opt.method, Parameters);
            
            Tmp{1}.Parameters.Labels = labels;
        else
            count = 1;
            
            for i = 1:length(Dico.Parameters.Labels)
                tmp = split(Dico.Parameters.Labels{i},'.',2);
                if any(strcmp(tmp{end}, opt.Params))
                    params.Par(:,count)     = Tmp{1}.Parameters.Par(:,i);
                    params.Labels{count}    = Tmp{1}.Parameters.Labels{i};
                    count = count +1;
                end
            end
            Tmp{1}.Parameters = params;
            
            % Parameters of the regression
            clear Parameters
            Parameters.K = 30;
            Parameters.cstr.Sigma   = 'i'; %'d' can be used
            Parameters.cstr.Gammat  = 'd'; 
            Parameters.cstr.Gammaw  = '';
            Parameters.Lw           = 0;
            
            [Estimation, Parameters] = AnalyzeMRImages(Xobs, Tmp, opt.method, Parameters);
            
            labels      = Tmp{1}.Parameters.Labels; 
            save(model_filename,'Parameters', 'labels')
        end
        
        Map.Y   = permute(Estimation.Regression.Y, [1 2 4 3]);
        Map.Std	= permute(Estimation.Regression.Cov, [1 2 4 3]).^.5;
end


% Extract maps (and modify unit if necessary)
count = 1;
for i = 1:length(Tmp{1}.Parameters.Labels)
    tmp = split(Tmp{1}.Parameters.Labels{i},'.',2);
    if any(strcmp(tmp{end}, opt.Params))
        
        Labels{count} = tmp{end};
        
        
        switch tmp{end}
        %If the ression method is performed, extract also the confidence maps
            case {'Vf', 'SO2'} %convert to percent 
                MapStruct{count}    = 100*Map.Y(:,:,:,i);
                if strcmp(opt.method, 'RegressionMRF')
                    StdStruct{count}    = 100*Map.Std(:,:,:,i);
                end
            case {'VSI', 'R'} % convert m to µm
                MapStruct{count}    = 1e6*Map.Y(:,:,:,i);
                if strcmp(opt.method, 'RegressionMRF')
                    StdStruct{count}	= 1e6*Map.Std(:,:,:,i);
                end
            otherwise
                MapStruct{count} 	= Map.Y(:,:,:,i);
                if strcmp(opt.method, 'RegressionMRF')
                    StdStruct{count}    = 100*Map.Std(:,:,:,i);
                end
        end
        
        count   = count +1;
    end
end

    
% Json processing
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


            if strcmp(opt.method, 'RegressionMRF')
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