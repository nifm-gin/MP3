function [files_in,files_out,opt] = Module_DeltaR2(files_in,files_out,opt)

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
    module_option(:,4)   = {'NbEchoUsed','3'};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    module_option(:,9)   = {'output_filename_ext','DeltaR2'};
    module_option(:,10)   = {'Output_orientation','First input'};
    
    
    
    
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    %
    %      additional_info_name = {'Operation', 'Day of the 2nd scan', 'Name of the 2nd scan'};
    %         additional_info_data = {'Addition', 'first', handles.Math_parameters_list{1}};
    %         additional_info_format = {{'Addition', 'Subtraction', 'Multiplication', 'Division', 'Percentage', 'Concentration'},...
    %             ['first', '-1', '+1',  'same', handles.Math_tp_list], handles.Math_parameters_list};
    %         additional_info_editable = [1 1 1];
    
    
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
        {'Description of the module : this module compute a deltaR2 map from 2 MGEFIDSE scans acquired before and after a UPSIO injection',...
        'If you use this function, please refere to this article : Tropes et al. MRM 2015'}}  ;
    user_parameter(:,2)   = {'Select the MGEFIDSE Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the MGEFIDSE Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename','char','DeltaR2','output_filename_ext','','',...
        {'Specify the name of the output scan.'
        'Default filename extension is ''DeltaR2''.'}'};
    user_parameter(:,5)   = {'   .Number of echo used','cell',{'1', '3'},'NbEchoUsed','', '','Please select the number of echo used to compute the deltaR2 map. Is 3 is selected, a mean value bewteen the 3 echoes around the SpinEcho will be computed'};
    user_parameter(:,6)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','','',...
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
    error('Module_DeltaR2:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%input_pre=niftiread(files_in.In1{1});
%input_post=niftiread(files_in.In2{1});

%% WTF
% %%% ? What's the difference between se_echotime and RAW_data_pre.echotime ?
% %%% I don't think those informations are available in the nifti files.
% se_echotime = scan_acqp('##$SpinEchoTime=',data_pre.uvascim.image.texte,1);
% if isnan(se_echotime)
%     se_echotime = scan_acqp('##$PVM_EchoTime=',data_pre.uvascim.image.texte,1);
% end
% se_echo_pos = abs(RAW_data_pre.echotime - se_echotime) <= 2;
% se_echo_pos = find(se_echo_pos == 1);
%%
%nb_echo_used=opt.nb_echos;
se_echotime = J1.SpinEchoTime.value; 
%se_echotime = J1.EchoTime.value;
%se_echo_pos = 1:length(J1.EchoTime.value);
[~,se_echo_pos] = min(abs(J1.EchoTime.value - se_echotime));
%se_echo_pos = abs(J1.EchoTime.value/1000 - se_echotime) <= 1/1000; %Conversion from ms to sec
%se_echo_pos = find(se_echo_pos == 1);
NbEchoUsed = str2double(opt.NbEchoUsed);

    % Compute the deltaR2 map for all slices
deltaR2_map=zeros([size(input_pre,1), size(input_pre, 2), size(input_pre, 3)]);     % DeltaR2 map
for m_slice=1:size(input_pre, 3)
    if NbEchoUsed == 1
        temp_avant=squeeze(input_pre(:,:,m_slice,se_echo_pos));
        temp_apres=squeeze(input_post(:,:,m_slice,se_echo_pos));
    else
        temp_avant=mean(squeeze(input_pre(:,:,m_slice,se_echo_pos-1:se_echo_pos+1)),3);
        temp_apres=mean(squeeze(input_post(:,:,m_slice,se_echo_pos-1:se_echo_pos+1)),3);
    end
    index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
    index_apres=find(temp_apres<(1e-3*mean(temp_apres(:)))); 
    %To avoid division by small numbers
    temp_apres(index_apres)=1;
    warning off %#ok<WNOFF>
    deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1

%         deltaR2_slice_nbr = deltaR2_slice_nbr+1;
%         % Compute the CMRO2 map each slice with the same offset
%         if nb_echo_used == 1
%             if rescale == 0
%                 temp_avant=squeeze(RAW_data_pre.data(:,:,se_echo_pos,i));
%                 temp_apres=squeeze(RAW_data_post.data(:,:,se_echo_pos,j));
%             else
%                 temp_avant=imresize(squeeze(RAW_data_pre.data(:,:,se_echo_pos,i)),[final_res final_res],'bilinear');
%                 temp_apres=imresize(squeeze(RAW_data_post.data(:,:,se_echo_pos,j)),[final_res final_res],'bilinear');
%             end
%         else
%             if rescale == 0
%                 temp_avant=squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,i),3));
%                 temp_apres=squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,j),3));
%             else
%                 temp_avant=imresize(squeeze(mean(RAW_data_pre.data(:,:,se_echo_pos-1:se_echo_pos+1,i),3)),[final_res final_res],'bilinear');
%                 temp_apres=imresize(squeeze(mean(RAW_data_post.data(:,:,se_echo_pos-1:se_echo_pos+1,j),3)),[final_res final_res],'bilinear');
%             end
%         end
%         index_avant=find(temp_avant<(1e-3*mean(temp_avant(:))));
%         index_apres=find(temp_apres<(1e-3*mean(temp_apres(:))));
% 
%         %To avoid division by small numbers
%         temp_apres(index_apres)=1;
%         warning off %#ok<WNOFF>
%         deltaR2=-(1/se_echotime) * log(temp_apres ./ temp_avant);%ms-1
%         warning on %#ok<WNON>
%         deltaR2(index_avant)=0;
%         deltaR2(index_apres)=0;
%         deltaR2_map.reco.data(:,:,1,deltaR2_slice_nbr)=deltaR2;








%         A = -(1./se_echotime)';
%         B = log(temp_apres ./ temp_avant);
%         B = permute(B, [3,2,1]);
%         for i=1:size(B,3)
%             deltaR2(i,:) = A*squeeze(B(:,:,i));
%         end
    %deltaR2= A.*B ;%ms-1
    %deltaR2 = permute(deltaR2, [3,2,1]);
    warning on %#ok<WNON>
    %deltaR2(index_avant)=0;
    %deltaR2(index_apres)=0;
    deltaR2_map(:,:,m_slice,1)=deltaR2;
    %Compute the err map
    %proc.err(:,:,1,m_slice)=tmpim;
end



% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)

OutputImages = deltaR2_map;
%OutputImages(OutputImages < 0) = -1;
OutputImages(OutputImages > 5000) = -1;
OutputImages(isnan(OutputImages)) = -1;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header, 'axial');
else
    OutputImages_reoriented = OutputImages;
end





% save the new files (.nii & .json)
% update the header before saving the new .nii
info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});


info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
info2.ImageSize = size(OutputImages_reoriented);
%info2.Description = [info.Description, 'Modified by T2map Module'];


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