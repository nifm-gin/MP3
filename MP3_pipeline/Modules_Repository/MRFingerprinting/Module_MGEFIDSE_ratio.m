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
    module_option(:,11)  = {'output_name_bvf',  'BVf'};
    module_option(:,12)  = {'output_name_vsi',  'VSI'};
    module_option(:,13)  = {'output_name_sto2', 'StO2'};
    
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
        'Based on the Dan Ma paper [TODO: add ref]'
        ''
        'Prerequisite: Put your ''.json'' dictionaries (pre and post scans) in the ''data/dictionaries'' folder'
        }'};
    
    user_parameter(:,2)   = {'Select the MGEFIDSE Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the MGEFIDSE Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    folder_files	= dir(fullfile(s{1}, 'data/dictionaries/'));
    folder_files    = folder_files([folder_files.isdir]);
    opt.Module_settings.folder = fullfile(s{1}, 'data/dictionaries/');
    if isempty(folder_files), folder_files(1).name = ' '; end
    user_parameter(:,4)   = {'   .Dictionary Pre/Post file folder','cell', {folder_files.name}, 'dictionary_folder_filename','','Mandatory',...
        {'Select the folder containing Pre/Post dico files'}};
    
    user_parameter(:,5)   = {'   .Prefix','char', '', 'prefix', '', '',...
        {'Choose a prefix'}};
    user_parameter(:,6)   = {'   .Smooth?','cell', {'Yes','No'}, 'filtered', '', '',...
        {''}};
    user_parameter(:,7)   = {'   .Method','cell', {'ClassicMRF', 'RegressionMRF'}, 'method', '', '',...
        { 'Choose:'
        '	- ''ClassicMRF'' to use Dan Ma method'
        '	- ''RegressionMRF'' to use proposed approach'
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
    opt.Table_out = opt.Table_in(1,:);
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.prefix), char(opt.output_name_bvf)]));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_name_bvf]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    opt.Table_out.IsRaw = categorical(cellstr('0'));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Filename) '.nii'];
    files_out.In1{1} = f_out;
    
    opt.Table_out(2,:) = opt.Table_out(1,:);
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName(2) = categorical(cellstr([char(opt.prefix), char(opt.output_name_vsi)]));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName(2) = categorical(cellstr([char(opt.Table_out.SequenceName(2)), opt.output_name_vsi]));
    end
    opt.Table_out.Filename(2) = categorical(cellstr([char(opt.Table_out.Patient(2)), '_', char(opt.Table_out.Tp(2)), '_', char(opt.Table_out.SequenceName(2))])); 
    f_out = [char(opt.Table_out.Path(2)), char(opt.Table_out.Filename(2)) '.nii'];
    files_out.In1{2} = f_out;
    
    opt.Table_out(3,:) = opt.Table_out(1,:);
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName(3) = categorical(cellstr([char(opt.prefix), char(opt.output_name_sto2)]));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName(3) = categorical(cellstr([char(opt.Table_out.SequenceName(3)), opt.output_name_sto2]));
    end
    opt.Table_out.Filename(3) = categorical(cellstr([char(opt.Table_out.Patient(3)), '_', char(opt.Table_out.Tp(3)), '_', char(opt.Table_out.SequenceName(3))])); 
    f_out = [char(opt.Table_out.Path(3)), char(opt.Table_out.Filename(3)) '.nii'];
    files_out.In1{3} = f_out;
    
    if strcmp(opt.method,'RegressionMRF')
        opt.Table_out(4,:) = opt.Table_out(1,:);
        opt.Table_out(5,:) = opt.Table_out(2,:);
        opt.Table_out(6,:) = opt.Table_out(3,:);
        
        opt.Table_out.Filename(4) =  categorical(cellstr([char(opt.Table_out.Patient(1)), '_', char(opt.Table_out.Tp(1)), '_', char(opt.Table_out.SequenceName(1)) '_confidence'])); 
        opt.Table_out.Filename(5) =  categorical(cellstr([char(opt.Table_out.Patient(2)), '_', char(opt.Table_out.Tp(2)), '_', char(opt.Table_out.SequenceName(2)) '_confidence'])); 
        opt.Table_out.Filename(6) =  categorical(cellstr([char(opt.Table_out.Patient(3)), '_', char(opt.Table_out.Tp(3)), '_', char(opt.Table_out.SequenceName(3)) '_confidence'])); 
        opt.Table_out.SequenceName(4) = categorical(cellstr([char(opt.prefix), char(opt.output_name_bvf), '_confidence']));
        opt.Table_out.SequenceName(5) = categorical(cellstr([char(opt.prefix), char(opt.output_name_vsi), '_confidence']));
        opt.Table_out.SequenceName(6) = categorical(cellstr([char(opt.prefix), char(opt.output_name_sto2), '_confidence']));
            
        files_out.In2{1} = [char(opt.Table_out.Path(2)), char(opt.Table_out.Filename(4)) '.nii'];
        files_out.In2{2} = [char(opt.Table_out.Path(2)), char(opt.Table_out.Filename(5)) '.nii'];
        files_out.In2{3} = [char(opt.Table_out.Path(2)), char(opt.Table_out.Filename(6)) '.nii'];
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


opt.dictionary_folder_filename = fullfile(opt.folder, opt.dictionary_folder_filename);

% Read json files and create ratio dictionary
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

switch opt.method
    
    case 'ClassicMRF'
        
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
        
    case 'RegressionMRF'
        
        if ~exist(model_filename,'file')
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
end
        

% Generate ratio signals from scans (and ROI if given)
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

% ind = 1;
% figure
% plot(Dico.Tacq(1:size(Dico.MRSignals,2)), Dico.MRSignals(ind,:), '.-')
% hold on
% plot(Obs.EchoTime.value'*1e-3, interp1(Dico.Tacq(1:size(Dico.MRSignals,2)), Dico.MRSignals(ind,:), Obs.EchoTime.value'*1e-3), '.-')

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

% Compute MRF
switch opt.method
    
    case 'ClassicMRF'
        % Find something nicer than this permute trick
        Estimation  = AnalyzeMRImages(Xobs,Tmp,opt.method,[]);
        Map.Y       = permute(Estimation.GridSearch.Y, [1 2 4 3]);
        
    case 'RegressionMRF'
        
        %Compute the learing only one time per dictionar        
        if exist(model_filename,'file')
            load(model_filename,'Parameters','labels');
            Estimation = AnalyzeMRImages(Xobs, [], opt.method, Parameters);
            
            Tmp{1}.Parameters.Labels = labels;
        else
            
            for i = 1:length(Dico.Parameters.Labels)
                tmp = split(Dico.Parameters.Labels{i},'.',2);
                switch tmp{end}
                    case 'Vf'
                        params.Par(:,1)     = Tmp{1}.Parameters.Par(:,i);
                        params.Labels{1}    = Tmp{1}.Parameters.Labels{i};
                    case 'VSI'
                        params.Par(:,2)     = Tmp{1}.Parameters.Par(:,i);
                        params.Labels{2}    = Tmp{1}.Parameters.Labels{i};
                    case 'SO2'
                        params.Par(:,3)     = Tmp{1}.Parameters.Par(:,i);
                        params.Labels{3}    = Tmp{1}.Parameters.Labels{i};
                end
            end
            Tmp{1}.Parameters = params;
            
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

% Extract maps
for i = 1:length(Tmp{1}.Parameters.Labels)
    tmp = split(Tmp{1}.Parameters.Labels{i},'.',2);
    switch tmp{end}
        case 'Vf'
            BVf     = 100*Map.Y(:,:,:,i);
        case 'VSI'
            VSI     = 1e6*Map.Y(:,:,:,i);
        case 'SO2'
            StO2    = 100*Map.Y(:,:,:,i);
    end
end
if strcmp(opt.method, 'RegressionMRF')
    for i = 1:length(Tmp{1}.Parameters.Labels)
        tmp = split(Tmp{1}.Parameters.Labels{i},'.',2);
        switch tmp{end}
            case 'Vf'
                BVf_conf    = 100*Map.Std(:,:,:,i);
            case 'VSI'
                VSI_conf    = 1e6*Map.Std(:,:,:,i);
            case 'SO2'
                StO2_conf   = 100*Map.Std(:,:,:,i);
        end
    end
end

    
% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J       = ReadJson(jsonfile);
J       = KeepModuleHistory(J, struct('files_in',files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 


% Reoriented and save nifti maps
nifti_header = spm_vol(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
info.Filemoddate = char(datetime('now'));

for i = 1:length(files_out.In1)
    if contains(files_out.In1{i},'BVf')
            [path, name, ~] = fileparts(files_out.In1{i});
            jsonfile = [path, '/', name, '.json'];
            WriteJson(J, jsonfile)
            
            %BVf     = write_volume(BVf, nifti_header);
            info.Filename = files_out.In1{i};
            info.ImageSize = size(BVf);
            info.PixelDimensions = info.PixelDimensions(1:length(size(BVf)));
            info.Datatype = class(BVf);
            
            niftiwrite(BVf, files_out.In1{i}, info);
    end
    if contains(files_out.In1{i},'VSI')
            [path, name, ~] = fileparts(files_out.In1{i});
            jsonfile = [path, '/', name, '.json'];
            WriteJson(J, jsonfile)
            
            %VSI     = write_volume(VSI, nifti_header);
            info.Filename = files_out.In1{i};
            info.ImageSize = size(VSI);
            info.PixelDimensions = info.PixelDimensions(1:length(size(VSI)));
            info.Datatype = class(VSI);
            
            niftiwrite(VSI, files_out.In1{i}, info);
    end
    if contains(files_out.In1{i},'StO2')
            [path, name, ~] = fileparts(files_out.In1{i});
            jsonfile = [path, '/', name, '.json'];
            WriteJson(J, jsonfile)
            
            %VSI     = write_volume(VSI, nifti_header);
            info.Filename = files_out.In1{i};
            info.ImageSize = size(StO2);
            info.PixelDimensions = info.PixelDimensions(1:length(size(StO2)));
            info.Datatype = class(StO2);
            
            niftiwrite(StO2, files_out.In1{i}, info);
    end
end

if strcmp(opt.method, 'RegressionMRF')
    for i = 1:length(files_out.In2)
        if contains(files_out.In2{i},'BVf')
                [path, name, ~] = fileparts(files_out.In2{i});
                jsonfile = [path, '/', name, '.json'];
                WriteJson(J, jsonfile)

                %BVf     = write_volume(BVf, nifti_header);
                info.Filename = files_out.In2{i};
                info.ImageSize = size(BVf_conf);
                info.PixelDimensions = info.PixelDimensions(1:length(size(BVf_conf)));
                info.Datatype = class(BVf_conf);

                niftiwrite(BVf_conf, files_out.In2{i}, info);
        end
        if contains(files_out.In2{i},'VSI')
                [path, name, ~] = fileparts(files_out.In2{i});
                jsonfile = [path, '/', name, '.json'];
                WriteJson(J, jsonfile)

                %VSI     = write_volume(VSI, nifti_header);
                info.Filename = files_out.In2{i};
                info.ImageSize  = size(VSI_conf);
                info.PixelDimensions = info.PixelDimensions(1:length(size(VSI_conf)));
                info.Datatype = class(VSI_conf);

                niftiwrite(VSI_conf, files_out.In2{i}, info);
        end
        if contains(files_out.In2{i},'StO2')
                [path, name, ~] = fileparts(files_out.In2{i});
                jsonfile = [path, '/', name, '.json'];
                WriteJson(J, jsonfile)

                %VSI     = write_volume(VSI, nifti_header);
                info.Filename = files_out.In2{i};
                info.ImageSize  = size(StO2_conf);
                info.PixelDimensions = info.PixelDimensions(1:length(size(StO2_conf)));
                info.Datatype = class(StO2_conf);

                niftiwrite(StO2_conf, files_out.In2{i}, info);
        end
    end

end
