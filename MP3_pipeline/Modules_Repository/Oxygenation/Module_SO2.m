function [files_in,files_out,opt] = Module_SO2(files_in,files_out,opt)

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
    module_option(:,3)   = {'trash_below',0};
    module_option(:,4)   = {'trash_after',Inf};
    module_option(:,5)   = {'output_filename','SO2'};
    module_option(:,6)   = {'B0',4.7};
    module_option(:,7)   = {'Gamma',2.67502};
    module_option(:,8)   = {'deltaxi_HBC',0.264};
    module_option(:,9)   = {'HCT',0.375};
    module_option(:,10)  = {'UseT1Map','No'};
    module_option(:,11)  = {'HTC_map','None'};
    module_option(:,12)  = {'OutputSequenceName','AllName'};
    module_option(:,13)  = {'RefInput',1};
    module_option(:,14)  = {'InputToReshape',1};
    module_option(:,15)  = {'Table_in', table()};
    module_option(:,16)  = {'Table_out', table()};
    module_option(:,17)   = {'Output_orientation','First input'};
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
    user_parameter(:,1)   = {'Description','Text','','','', '',{'Description of the module',...
        'this module compute a SO2 map from a R2prim and a BVf maps (see R2prim and BVf maps for more information) ',...
        'If you use this function, please refere to this article : Christen et al. NMR in bimed 2011'}};
    user_parameter(:,2)   = {'Select a R2Prim scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select a BVf scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'Parameters','','','','', '',''};
    user_parameter(:,5)   = {'   .Output filename','char','SO2','output_filename','','',...
        {'Specify the string to be added to the filename input.'
        'Default filename extension is ''SO2''.'}'};
    user_parameter(:,6)   = {'   .B0','numeric', '','B0','','',...
        {''}'};
    user_parameter(:,7)   = {'   .Gamma (10^8)','numeric','','Gamma','','',...
        {''}'};
    user_parameter(:,8)   = {'   .deltaxi HBC (10^-6)','numeric', '','deltaxi_HBC','','',...
        {''}'};
    user_parameter(:,9)   = {'   .HCT','numeric','','HCT','', '',''};
    user_parameter(:,10)   = {'   .Using T1map ?','1Scan','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,11)   = {'   .Using HCT_map ?','1Scan','','',{'SequenceName'}, 'Optional',''};
    user_parameter(:,12)   = {'   .Output orientation','cell',{'First input', 'Second input'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = First input'
        '--> Output orientation = Second input'
        }'};

VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'StO2'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename};

if isempty(files_out)
    opt.Table_out = opt.Table_in(opt.RefInput,:);
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
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
    error('Module_SO2:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end


% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii');  
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


input(1).nifti_header = spm_vol(files_in.In1{1});
input(2).nifti_header = spm_vol(files_in.In2{1});
J1 = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
J2 = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));

if strcmp(opt.Output_orientation, 'First input')   
    ref_scan = 1;
else
    ref_scan = 2;
end
R2Prim = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0);
BVf = read_volume(input(2).nifti_header, input(ref_scan).nifti_header, 0);

if isfield(files_in, 'In3') &&  ~isempty(files_in.In3{1})
    input(3).nifti_header = spm_vol(files_in.In3{1});
    J3 = spm_jsonread(strrep(files_in.In3{1}, '.nii', '.json'));
    T1Map = read_volume(input(3).nifti_header, input(ref_scan).nifti_header, 0);
else
    T1Map = [];
end
if isfield(files_in, 'In4') && ~isempty(files_in.In4{1})
    input(4).nifti_header = spm_vol(files_in.In4{1});
    J4 = spm_jsonread(strrep(files_in.In4{1}, '.nii', '.json'));
    HCTMap = read_volume(input(4).nifti_header, input(ref_scan).nifti_header, 0);
else
    HCTMap = [];
end




B0=opt.B0;
gamma=opt.Gamma; gamma = gamma*10^8;
deltaxi_HBC=opt.deltaxi_HBC; deltaxi_HBC = deltaxi_HBC*10^-6;
HCT = opt.HCT;



if isempty(T1Map)
    % Apply the same DeltaKhi AND use HCT variable (using a HCTmap)
    if ~isempty(HCTMap)
        temp_slice_nbr = [size(R2Prim,3) size(BVf,3)  size(HCTMap,3)];
        for a=1:size(HCTMap,3)
            tmp(:,:,a,1) = imresize(HCTMap(:,:,a,1),[size(R2Prim, 1) size(R2Prim,2)],'bilinear');
        end
        HCTmap_resized = tmp;
        SO2map = [];
        if length(find(temp_slice_nbr == size(R2Prim,3))) ~= numel(temp_slice_nbr)
            for i = 1:size(R2Prim, 3)
                for j = 1:size(BVf, 3)
                    for x = 1:size(HCTMap, 3)
                        SO2map_slice_nbr = SO2map_slice_nbr+1;
                        % Compute the SO2map each slice with the same offset
                        SO2map(:,:,SO2map_slice_nbr,1)=(1-(R2Prim(:,:,i,1))./...
                            ((BVf(:,:,j,1)/100)*gamma*4/3*pi*deltaxi_HBC.*(HCTmap_resized(:,:,x,1)./100)*B0))*100;%en %
                    end
                    
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s ##$T1map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
        else
            % calculate SO2 map for all slices
            %             SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
        end
        
        
        % Apply the same DeltaKhi AND use constant HCT value
    else
        % check data compatibility (slice thickness and slice number)
        if J1.SliceThickness.value ~= J2.SliceThickness.value
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a slice thickness missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s',...
                R2prim_filename,BVf_map_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        SO2map = BVf;
        if size(R2Prim,3) ~= size(BVf,3)
            for i = 1:size(R2Prim, 3)
                for j = 1:size(BVf, 3)
                        SO2map_slice_nbr = SO2map_slice_nbr+1;
                        % Compute the SO2map map each slice with the same offset
                        SO2map(:,:,SO2map_slice_nbr,1)=(1-(R2Prim(:,:,i,1))./((BVf(:,:,j,1)/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
        else
            % calculate SO2 map for all slices
            SO2map=(1-(R2Prim)./((BVf/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
        end
    end
else
    % case using the T1map to segment the brain and apply 3 differents
    % DeltaKhi: gray, white and CF   AND use HCT variable (using a HCTmap)
    if ~isempty(HCTMap)
        % check data compatibility (slice thickness and slice number)
        temp_thickness = [J1.SliceThickness.value J2.SliceThickness.value  J3.SliceThickness.value];
        if length(find(temp_thickness == J1.SliceThickness.value)) ~= numel(temp_thickness)
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a slice thickness missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s and \n##$T1map=%s',...
                R2prim_filename,BVf_map_filename, T1map_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end
        temp_slice_nbr = [size(R2Prim,3) size(BVf,3)  size(HCTMap,3)];
        HCTmap_resized = HCTMap;
        
        for a=1:size(HCTMap,3)
            tmp(:,:,a,1) = imresize(HCTMap(:,:,a,1),[size(R2Prim, 1) size(R2Prim,2)],'bilinear');
        end
        HCTMap = tmp;
        SO2map = [];
        if length(find(temp_slice_nbr == size(R2Prim,3))) ~= numel(temp_slice_nbr)
            for i = 1:size(R2Prim, 3)
                for j = 1:size(BVf, 3)
                        for x = 1:size(HCTMap, 3)
                                SO2map_slice_nbr = SO2map_slice_nbr+1;
                                % Compute the SO2map each slice with the same offset
                                
                                mask_white = T1Map(:,:,j,1) <= 1300;
                                mask_gray = logical(T1Map(:,:,j,1) > 1300) & logical(T1Map(:,:,j,1) < 1900);
                                mask_CF = T1Map(:,:,j,1) >= 1900;
                                %
                                % calculate SO2 map for all slices
                                tmp=(1-(R2Prim(:,:,i,1))./((BVf(:,:,j,1)/100)*gamma*4/3*pi*5.36e-07.*(HCTmap_resized.reco.data(:,:,x,1)./100)*B0))*100;%en %
                                tmp(isnan(tmp)) = 0;
                                mask_white2 = mask_white.*tmp;
                                %
                                tmp=(1-(R2Prim(:,:,i,1))./((BVf(:,:,j,1)/100)*gamma*4/3*pi*2.6e-07.*(HCTmap_resized(:,:,x,1)./100)*B0))*100;%en %
                                tmp(isnan(tmp)) = 0;
                                mask_gray2 = mask_gray.*tmp;
                                
                                tmp=(1-(R2Prim(:,:,i,1))./((BVf(:,:,j,1)/100)*gamma*4/3*pi*0.33e-07.*(HCTmap_resized(:,:,x,1)./100)*B0))*100;%en %
                                tmp(isnan(tmp)) = 0;
                                mask_CF2 = mask_CF.*tmp;
                                
                                SO2map(:,:,SO2map_slice_nbr,1)= mask_white2 + mask_gray2 + mask_CF2;
                                
                                
                                %                                 SO2map.reco.data(:,:,1,SO2map_slice_nbr)=(1-(R2prim.reco.data(:,:,1,i))./...
                                %                                     ((BVf_map.reco.data(:,:,1,j)/100)*gamma*4/3*pi*deltaxi_HBC.*(HCTmap_resized.reco.data(:,:,1,x)./100)*B0))*100;%en %
                                % Update the SO2map structure
                        end
                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s ##$T1map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
            %         else
            %             % calculate SO2 map for all slices
            % %             SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %
            %
            %             % classification of the brain in 3 classes (Gray, White and
            %             % Cerebrospinal fluid (CF))based on the T1 value
            %             mask_white = T1map.reco.data <= 1300;
            %             mask_gray = logical(T1map.reco.data > 1300) & logical(T1map.reco.data < 1900);
            %             mask_CF = T1map.reco.data >= 1900;
            %             %
            %             % calculate SO2 map for all slices
            %             tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*5.36e-07*HCT*B0))*100;%en %
            %             tmp(isnan(tmp)) = 0;
            %             mask_white2 = mask_white.*tmp;
            %             %
            %             tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*2.6e-07*HCT*B0))*100;%en %
            %             tmp(isnan(tmp)) = 0;
            %             mask_gray2 = mask_gray.*tmp;
            %
            %             tmp=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi*0.33e-07*HCT*B0))*100;%en %
            %             tmp(isnan(tmp)) = 0;
            %             mask_CF2 = mask_CF.*tmp;
            %
            %             SO2map.reco.data = mask_white2 + mask_gray2 + mask_CF2;
        end
        
        % case using the T1map to segment the brain and apply 3 differents
        % DeltaKhi: gray, white and CF   AND use constant HCT value
    else
        % check data compatibility (slice thickness and slice number)
        temp_thickness = [J1.SliceThickness.value J2.SliceThickness.value  J3.SliceThickness.value];
        temp_slice_nbr = [size(R2Prim,3) size(BVf,3)  size(T1Map,3)];
        if length(find(temp_thickness == J1.SliceThickness.value)) ~= numel(temp_thickness)
            warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ a slice thickness missmatch between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s and \n##$T1map=%s',...
                R2prim_filename,BVf_map_filename, T1map_filename);
            msgbox(warning_text, 'SO2 map warning') ;
            SO2map = [];
            return
        end

        
        SO2map = BVf;
        if length(find(temp_slice_nbr == size(R2Prim,3))) ~= numel(temp_slice_nbr)
            for i = 1:size(R2Prim, 3)
                for j = 1:size(BVf, 3)
                        for x = 1:size(T1Map, 3)
                                SO2map_slice_nbr = SO2map_slice_nbr+1;
                                % Compute the SO2map each slice with the same offset
                                SO2map(:,:,SO2map_slice_nbr,1)=(1-(R2Prim(:,:,i,1))./((BVf(:,:,j,1)/100)*gamma*4/3*pi*deltaxi_HBC*HCT*B0))*100;%en %

                        end

                end
            end
            if SO2map_slice_nbr == 0
                warning_text = sprintf('##$ Can not calculate the SO2 map because there is\n##$ no slice offset match between\n##$R2prim=%s\n##$ and \n##$BVf_map=%s ##$T1map=%s',...
                    R2prim_filename,BVf_map_filename);
                msgbox(warning_text, 'SO2 map warning') ;
                return
            end
        else
            % classification of the brain in 3 classes (Gray, White and
            % Cerebrospinal fluid (CF))based on the T1 value
            mask_white = T1Map <= 1300;
            mask_gray = logical(T1Map > 1300) & logical(T1Map < 1900);
            mask_CF = T1Map >= 1900;
            %
            % calculate SO2 map for all slices
            tmp=(1-(R2Prim)./((BVf/100)*gamma*4/3*pi*5.36e-07*HCT*B0))*100;%en %
            tmp(isnan(tmp)) = 0;
            mask_white2 = mask_white.*tmp;
            %
            tmp=(1-(R2Prim)./((BVf/100)*gamma*4/3*pi*2.6e-07*HCT*B0))*100;%en %
            tmp(isnan(tmp)) = 0;
            mask_gray2 = mask_gray.*tmp;
            
            tmp=(1-(R2Prim)./((BVf/100)*gamma*4/3*pi*0.33e-07*HCT*B0))*100;%en %
            tmp(isnan(tmp)) = 0;
            mask_CF2 = mask_CF.*tmp;
            
            SO2map = mask_white2 + mask_gray2 + mask_CF2;
            
            % test if a relationship between the T1 value and the deltaKhi
            % exist
            %         deltakhi_map = ((T1map.reco.data.*-0.0041)+9.23).*0.0000001;
            %         SO2map.reco.data=(1-(R2prim.reco.data)./((BVf_map.reco.data/100)*gamma*4/3*pi.*deltakhi_map*HCT*B0))*100;%en %
            
        end
    end
end


% complete structure









OutputImages = SO2map;
OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 5000) = NaN;
%OutputImages(isnan(OutputImages)) = NaN;


OutputImages(isinf(OutputImages)) = -1;
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header);
end


info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});
% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages_reoriented);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages_reoriented)));
info2.ImageSize = size(OutputImages_reoriented);
%info2.Description = [info.Description, 'Modified by T2star_map Module'];

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
