function [files_in,files_out,opt] = Module_MGE2Dfrom3D(files_in,files_out,opt)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
     % define every option needed to run this module
      %%   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'output_filename_ext','MGE2Dfrom3D'};
    module_option(:,4)   = {'OutputSequenceName','AllName'};
    module_option(:,5)   = {'first_slice', 4};
    module_option(:,6)   = {'number_of_slices',4};
    module_option(:,7)   = {'echo_to_trash','end'};
    module_option(:,8)   = {'OutputResolution','128'};
    module_option(:,9)   = {'RefInput',1};
    module_option(:,10)   = {'InputToReshape',1};
    module_option(:,11)  = {'Table_in', table()};
    module_option(:,12)  = {'Table_out', table()};
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
        {'ASL_InvEff Description'}'};
user_parameter(:,2)   = {'Select a MGE-3D scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
user_parameter(:,3)   = {'Parameters','','','','', '',''};
user_parameter(:,4)   = {'   .Output filename extension','char','MGE2Dfrom3D','output_filename_ext','','',...
    {'Specify the string to be added to the filename input.'
    'Default filename extension is ''MGE2Dfrom3D''.'}'};
user_parameter(:,5)   = {'   .First slice','numeric','','first_slice','', '',''};
user_parameter(:,6)   = {'   .Number of slices','numeric','','number_of_slices','', '',''};
user_parameter(:,7)   = {'   .Echo to trash','char','','echo_to_trash','', '',''};
user_parameter(:,8)   = {'   .Output Resolution','cell',{'64', '128', '256'},'OutputResolution','', '',''};



VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'MGE2Dfrom3D'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename_ext};




if isempty(files_out)
    opt.Table_out = opt.Table_in;
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
    error('Module_MGE2Dfrom3D:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs

% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii');  
% end
% 
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


%% load input Nii file
N = niftiread(files_in.In1{1});
% load nifti info
info = niftiinfo(files_in.In1{1});

%% load input JSON file


input(1).nifti_header = spm_vol(files_in.In1{1});

J = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
% MGE3D = read_volume(input(1).nifti_header, input(1).nifti_header, 0);
% 
% 
% output.nifti_header = input(1).nifti_header;
MGE3D = N;
% Empty memory
clear data

first_slice =opt.first_slice;
nbr_of_slice = opt.number_of_slices;
echo_to_trash = opt.echo_to_trash;

% Some echoes need to be remove --> human data
if ~strcmp(echo_to_trash, 'end')
    echo_to_trash = str2double(echo_to_trash);
    MGE3D = MGE3D(:,:,:,1:echo_to_trash);
    J.EchoTime.value = J.EchoTime.value(1:echo_to_trash);
end


[rows3d, cols3d, depths3d, echos3d]=size(MGE3D);
nbr_of_final_slice = fix((depths3d-first_slice+1)/nbr_of_slice);
if nbr_of_final_slice == 0
    return
end
temp2 = zeros([size(MGE3D, 1) size(MGE3D, 2) nbr_of_final_slice size(MGE3D, 4)]);
for i = 1:nbr_of_final_slice
    for k=1:echos3d
        temp=MGE3D(:,:,first_slice+(i*nbr_of_slice)-nbr_of_slice:first_slice+(i*nbr_of_slice)-1,k);
        temp=sum(abs(temp),3);
        temp2(:,:,i,k)=temp;
    end
end

% x = 1:rows3d;
% y = 1:cols3d;
%MGE2Dfrom3D = NaN(size(MGE3D));
if size(N,1) ~= size(N,2)
    error('The input image hasn''t a squared resolution. This case is not implemented')
end
facteur = size(N,1) / str2double(opt.OutputResolution);
if floor(facteur)==facteur
    error('The initial resolution is not a multiple of the wanted output one.')
end

mask = ones(facteur);
Test = zeros([size(temp2,1)/facteur,size(temp2,2)/facteur,size(temp2,3), size(temp2,4)]);
Test2 = zeros(size(temp2));
for i=1:size(temp2,3)
    for j=1:size(temp2,4)
        Test2(:,:,i,j) = conv2(temp2(:,:,i,j), mask, 'same');
        Test(:,:,i,j) = Test2(1:facteur:end, 1:facteur:end,i,j);
    end
end
% 

MGE2Dfrom3D = Test;


J.SliceThickness.value = J.SliceThickness.value*nbr_of_slice;
J.ImageInAcquisition.value = nbr_of_final_slice*length(J.EchoTime.value);



OutputImages = MGE2Dfrom3D;
OutputImages(isinf(OutputImages)) = -1;
OutputImages_reoriented = OutputImages;
% if ~exist('OutputImages_reoriented', 'var')
%     OutputImages_reoriented = write_volume(OutputImages, input(1).nifti_header);
% end

% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
info2.PixelDimensions(3) = info2.PixelDimensions(3)*nbr_of_slice;
info2.Transform.T(1,1) = info2.Transform.T(1,1)*facteur;
info2.Transform.T(2,2) = info2.Transform.T(2,2)*facteur;
info2.Transform.T(4,3) = info2.Transform.T(4,3)+first_slice*info.Transform.T(3,3);
info2.Transform.T(3,3) = info2.Transform.T(3,3)*nbr_of_slice;
info2.ImageSize = size(OutputImages_reoriented);
%info2.Description = [info.Description, 'MGE2Dfrom3D'];

% save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 



%% Initial function
% fid=fopen(MGE3D_map_filename ,'r');
% if fid>0
%     fclose(fid);
%     data = load(MGE3D_map_filename);
%     reco = data.uvascim.image.reco;
% else
%     warning_text = sprintf('##$ Can not calculate the MGE2Dfrom3D map because there is\n##$ Somthing wrong with the data\n##$MGE3D=%s\n##$',...
%         MGE3D_map_filename);
%     msgbox(warning_text, 'MGE2Dfrom3D map warning') ;
%     MGE2Dfrom3D = [];
%     return
% end
% MGE2Dfrom3D = data.uvascim.image;
% 
% % Empty memory
% clear data
% 
% fist_slice =str2double(add_parameters{:}(3));
% nbr_of_slice = str2double(add_parameters{:}(4));
% resolution = str2double(add_parameters{:}(5));
% echo_to_trash = add_parameters{:}(8);
% 
% % Some echoes need to be remove --> human data
% if ~strcmp(echo_to_trash, 'end')
%     reco.data= reco.data(:,:,1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.echo_label = reco.echo_label(1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.echotime = reco.echotime(1:str2double(echo_to_trash));
%     reco.fov_offsets = reco.fov_offsets(:,1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.fov_orientation = reco.fov_orientation(:,1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.fov_phase_orientation = reco.fov_phase_orientation(1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.label = reco.label(1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.no_echoes = str2double(echo_to_trash);
%     MGE2Dfrom3D.reco.phaselabel = reco.phaselabel(1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.scaling_factor = reco.scaling_factor(1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.scaling_offset = reco.scaling_offset(1:str2double(echo_to_trash),:);
%     MGE2Dfrom3D.reco.unit = reco.unit(1:str2double(echo_to_trash),1:str2double(echo_to_trash));
% end
% 
% 
% [rows3d, cols3d, echos3d, depths3d]=size(reco.data);
% % nbr_of_final_slice = round((depths3d-fist_slice+1)/nbr_of_slice)-1;
% nbr_of_final_slice = fix((depths3d-fist_slice+1)/nbr_of_slice);
% if nbr_of_final_slice ==0
%     return
% end
% temp2  =zeros([size(reco.data, 1) size(reco.data, 2) size(reco.data, 3) nbr_of_final_slice]);
% for i = 1:nbr_of_final_slice
%     for k=1:echos3d
%         temp=(reco.data(:,:,k,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
%         temp=sum(abs(temp),4);
%         temp2(:,:,k,i)=temp;
%     end
% end
% 
% xfactor = rows3d/resolution;
% yfactor = cols3d/resolution;
% 
% x = 1:xfactor:reco.no_samples;
% y = 1:yfactor:reco.no_views;
% MGE2Dfrom3D.reco.data = NaN(size(reco.data((x+xfactor-1)/2,(y+yfactor-1)/2,:,1:nbr_of_final_slice)));
% 
% for i = 1:nbr_of_final_slice
%     MGE2Dfrom3D.reco.data((x+xfactor-1)/2,(y+yfactor-1)/2,:,i)=(temp2(x,y,:,i)+...
%         (temp2(x+xfactor-1,y,:,i)+...
%         temp2(x,y+yfactor-1,:,i))+...
%          temp2(x+xfactor-1,y+yfactor-1,:,i));
% end
% if length(MGE2Dfrom3D.reco.fov) == 3
%     MGE2Dfrom3D.reco.fov(3)=[];
% end
% MGE2Dfrom3D.reco.fov_offsets(:,:,nbr_of_final_slice+1:end)=[];
% for i = 1:nbr_of_final_slice
%     temp=squeeze(reco.fov_offsets(:,1,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
%     new_offset=median(temp,2);
%     MGE2Dfrom3D.reco.fov_offsets(:,:,i) = repmat(new_offset, [1 MGE2Dfrom3D.reco.no_echoes]);
% end
% for i = 1:nbr_of_final_slice
%     temp=squeeze(reco.fov_orientation(:,1,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
%     new_offset=median(temp,2);
%     MGE2Dfrom3D.reco.fov_orientation(:,:,i) = repmat(new_offset, [1 MGE2Dfrom3D.reco.no_echoes]);
% end
% for i = 1:nbr_of_final_slice
%     temp=squeeze(reco.fov_phase_orientation(:,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
%     MGE2Dfrom3D.reco.fov_phase_orientation(:,i)=median(temp,2);
% end
%  
% 
% MGE2Dfrom3D.reco.globalmax = max(MGE2Dfrom3D.reco.data(:));
% MGE2Dfrom3D.reco.globalmin = min(MGE2Dfrom3D.reco.data(:));
% MGE2Dfrom3D.reco.iminfos = 'MGE2Dfrom3D';
% MGE2Dfrom3D.reco.no_samples = resolution;
% MGE2Dfrom3D.reco.no_views = resolution;
% MGE2Dfrom3D.reco.phaselabel(:,nbr_of_final_slice+1:end)=[];
% MGE2Dfrom3D.reco.reco_meth = 'MGE2Dfrom3D';
% for i = 1:nbr_of_final_slice
%     temp=squeeze(reco.scaling_factor(:,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
%     MGE2Dfrom3D.reco.scaling_factor(:,i)=median(temp,2);
% end
% for i = 1:nbr_of_final_slice
%     temp=squeeze(reco.scaling_offset(:,fist_slice+(i*nbr_of_slice)-nbr_of_slice:fist_slice+(i*nbr_of_slice)-1));
%     MGE2Dfrom3D.reco.scaling_offset(:,i)=median(temp,2);
% end
% 
% MGE2Dfrom3D.reco.thickness = MGE2Dfrom3D.reco.thickness*nbr_of_slice;
% MGE2Dfrom3D.reco.unit(:,nbr_of_final_slice+1:end)=[];
% MGE2Dfrom3D.reco.no_slices = nbr_of_final_slice;
% 
% MGE2Dfrom3D.acq.thickness = MGE2Dfrom3D.reco.thickness;
% MGE2Dfrom3D.acq.pix_spacing = [MGE2Dfrom3D.reco.fov(1)/MGE2Dfrom3D.reco.no_samples MGE2Dfrom3D.reco.fov(2)/MGE2Dfrom3D.reco.no_views]';
% MGE2Dfrom3D.acq.fov_orientation =   median(squeeze(MGE2Dfrom3D.reco.fov_orientation(1:3,1,:)),2)';
% MGE2Dfrom3D.acq.fov_offsets =  median(squeeze(MGE2Dfrom3D.reco.fov_offsets(1:3,1,:)),2)';
% 
% 
% ParamConfig=sprintf('##$QuantifMethod=Merge slices\n##$First slice used=%s\n##$Number of slice merged=%s\n##$New resolution=%s\n##$Raw scan used=%s\n',...
%     add_parameters{:}{3},...
%     add_parameters{:}{4},...
%     add_parameters{:}{5},...
%     MGE3D_map_filename);
% MGE2Dfrom3D.reco.paramQuantif = ParamConfig;
% MGE2Dfrom3D.reco=orderfields(MGE2Dfrom3D.reco);
% 
% MGE2Dfrom3D.clip = [MGE2Dfrom3D.reco.globalmin MGE2Dfrom3D.reco.globalmax 1];