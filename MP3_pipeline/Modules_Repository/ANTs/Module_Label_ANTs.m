function [files_in,files_out,opt] = Module_Label_ANTs(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

	%   % define every option needed to run this module
	% --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'rois', ''};
    module_option(:,2)   = {'output_extension', '_labels'};
    module_option(:,3)   = {'resume', 'Yes'};
    
    module_option(:,4)   = {'RefInput',         1};
    module_option(:,5)   = {'InputToReshape',   1};
    module_option(:,6)   = {'Table_in',         table()};
    module_option(:,7)   = {'Table_out',        table()};
    module_option(:,8)   = {'OutputSequenceName','Extension'};
    
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
        'ANTs Atlas Matching:'
        '    ANTS was created by:'
        '    Brian B. Avants, Nick Tustison and Gang Song'
        '    Penn Image Computing And Science Laboratory'
        '    University of Pennsylvania'
        '    https://github.com/ANTsX/ANTs'
        }'};
    
    user_parameter(:,2)   = {'Select label map','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select maps as input','XScan','','',{'SequenceName'}, 'Optionnal',''};
    user_parameter(:,4)   = {'Parameters','','','','', '', ''};
    user_parameter(:,5)   = {'   .ROIs','char','','rois','', 'Mandatory',...
        {'Select all ROIs to extract separeted by '','' character (no space)'}};
    user_parameter(:,6)   = {'   .Ouput extension','char','','output_extension','', '',...
        {'Give the output extension'}};
    user_parameter(:,7)   = {'   .Resume data','cell',{'Yes', 'No'},'resume','', '',...
        {''}};
    
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
    
    tmp_Table_out = opt.Table_in(1,:);
    tmp_Table_out.IsRaw = categorical(0);   
    tmp_Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        tmp_Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext_atlas));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        tmp_Table_out.SequenceName = categorical(cellstr([char(opt.Table_in(1,:).SequenceName), opt.output_extension]));
    end
    tmp_Table_out.Filename = categorical(cellstr([char(opt.Table_in(1,:).Patient), '_', char(opt.Table_in(1,:).Tp), '_', char(tmp_Table_out.SequenceName)]));
    tmp_Table_out.Type = categorical(cellstr('Cluster'));
    opt.Table_out = tmp_Table_out;
        
    files_out.In1{1} = [char(tmp_Table_out.Path), char(opt.Table_in.Patient(1)), '_', char(opt.Table_in.Tp(1)), '_', char(tmp_Table_out.SequenceName), '.nii'];

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


nb_maps = length(files_in.In2);


% Create label maps
if strcmp(opt.resume,'No')
    
    label_map = niftiread(files_in.In1{1});
    if isempty(opt.rois)
        labels = unique(label_map(:));
    else
        labels = split(opt.rois, ',');
    end

    new_label_map = nan(size(label_map));
    for l = 1:length(labels)
        new_label_map(label_map == str2double(labels{l})) = 1;
    end


    for nn = 1:nb_maps

        complete_map = niftiread(files_in.In2{nn});
        complete_map = imresize3(complete_map, size(label_map)); %TODO

        labelized_map = complete_map .* new_label_map;


        % save nifti file
        % update the header before saving the new .nii
        info = niftiinfo(files_in.In2{nn});
        info.Filename       = files_out.In1{1};
        info.Filemoddate    = char(datetime('now'));
        info.ImageSize      = size(labelized_map);
        info.PixelDimensions = info.PixelDimensions(1:length(size(labelized_map)));

        niftiwrite(labelized_map, files_out.In1{nn}, info);


        % Json processing
        [path, name, ~] = fileparts(files_in.In1{1});
        jsonfile = [path, '/', name, '.json'];

        J = ReadJson(jsonfile);
        J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

        [path, name, ~] = fileparts(files_out.In1{nn});
        jsonfile = [path, '/', name, '.json'];
        WriteJson(J, jsonfile)    
    end

    
%Create a full structure to resume data
elseif strcmp(opt.resume,'Yes')
    
    Resume = struct();
    
    label_map = niftiread(files_in.In1{1});
    Label   = unique(label_map(:));
    
    for nn = 1:nb_maps
        if exist(files_in.In2{nn},'file')
            complete_map = niftiread(files_in.In2{nn});
            complete_map = imresize3(complete_map, size(label_map)); %TODO

            for l = 1:length(Label)
                h       = reshape(complete_map(label_map == Label(l)), 1,[]);
                Mean(l,1) = nanmean(h);
                Variance(l,1) = nanvar(h);
                Number(l,1) = length(h);
                Hist{l} = h;
            end

            name = char(opt.Table_in(nn+1,:).SequenceName);
            eval(['Resume.' name ' = table(Label, Mean, Variance, Number);']);
            eval(['Complete.' name ' = Hist;']);
        else
            for l = 1:length(Label)
                Mean(l,1) = nan;
                Variance(l,1) = nan;
                Number(l,1) = nan;
                Hist{l} = nan;
            end
            
            eval(['Resume.' name ' = table(Label, Mean, Variance, Number);']);
            eval(['Complete.' name ' = Hist;']);
        end
    end
    
    niftiwrite([], files_out.In1{1});
    [path,filename] = fileparts(files_out.In1{1});
    save([path filesep filename '.mat'], 'Resume', 'Complete')
end



