function [files_in,files_out,opt] = Module_Export_HTML_Quality_Control(files_in,files_out,opt)

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
    module_option(:,3)   = {'RefInput',1};
    module_option(:,4)   = {'InputToReshape',1};
    module_option(:,5)   = {'Table_in', table()};
    module_option(:,6)   = {'Table_out', table()};
    module_option(:,7)   = {'Output_HTML_Folder','HTML_Quality_Control'};
    module_option(:,8)   = {'Contrast','0.5 0.99'};
    module_option(:,9)   = {'NumSlices', 'Max(12)'};
    module_option(:,10)   = {'AutomaticJobsCreation', 'No'};
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
         {'This module create a Quelity Control report as an HTML page. Max 12 images.'
         'This module uses the toolbox "Spinal Cord Toolbox" (https://sourceforge.net/projects/spinalcordtoolbox/) and its adaptation by Tanguy Duval (TONIC - Toulouse)'}
        };
    user_parameter(:,2)   = {'   .Scan','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the scan that will be displayed'};
    user_parameter(:,3)   = {'   .ROI','1ROIOr1Cluster','','',{'SequenceName'},'Optional',...
         'You can also select an ROI or a cluster that will be displayed on the scan'};
    user_parameter(:,4)   = {'   .Output HTML folder','char','','Output_HTML_Folder','', '',''};
    user_parameter(:,5)   = {'   .Contrast tolerence','char','','Contrast','', '',''};
    user_parameter(:,6)   = {'   .Index of the slices to display','char','','NumSlices','', '',''};

    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional', 'Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)', 'VariableNames', VariableNames);
%%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
  
end
%%%%%%%%


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Export_Values_VoxelByVoxel:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
databScans = opt.Table_in(opt.Table_in.Type == categorical(cellstr('Scan')),:);
databROIs = opt.Table_in(opt.Table_in.Type == categorical(cellstr('ROI')) | opt.Table_in.Type == categorical(cellstr('Cluster')),:);
for i=1:size(databScans,1)
    
    
    
    
    file1 = [char(databScans.Path(i)), char(databScans.Filename(i)), '.nii'];
    file_h = spm_vol(file1);
    databR = databROIs(databROIs.Patient == databScans.Patient(i), :);
    databRo = databR(databR.Tp == databScans.Tp(i), :);
    roi = [char(databRo.Path), char(databRo.Filename), '.nii'];
    if ~isempty(databROIs) && isempty(databRo)
        continue
    elseif ~isempty(databROIs)
        roi_h = spm_vol(roi);
        roi_vol = read_volume(roi_h, roi_h, 0, 'Axial');
        file_vol = read_volume(file_h, roi_h, 0, 'Axial');
    else
        file_vol = read_volume(file_h, file_h, 0, 'Axial');
    end

    Path_html = [opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'HTML', filesep, opt.Output_HTML_Folder, filesep];
    Path_nii = [Path_html, 'TMP_NIFTIS', filesep];
    if ~exist(Path_html, 'dir') || ~exist(Path_nii, 'dir')
        mkdir(Path_nii);
    end
    
    
    file_vol(isnan(file_vol)) = 0;
%     %file_vol = histeq(file_vol, 1024);
%     


     file_vol = (file_vol-min(file_vol(:)));
     file_vol = file_vol./max(file_vol(:));
%     %N2 = (N2_brut-min(N2_brut(:)))./max(N2_brut(:));
%     m = nanmean(file_vol(:));
%     std = nanstd(file_vol(:));
%     N = 1;
%     if m-N*std<0
%         MIN = 0.0001;
%     else
%         MIN = m-N*std;
%     end
%     if m+N*std>1
%         MAX = 0.9999;
%     else
%         MAX = m+N*std;
%     end
%     
    %filevol(file_vol>1.35) = 1;
    %file_vol(file_vol<0.06) = 0;
    
    %file_vol = imadjustn(file_vol);
    file_vol = imadjustn(file_vol, stretchlim(file_vol(:), str2num(opt.Contrast)));%, [MIN; MAX]);%, [0,128]);
    %file_vol = imadjustn(file_vol, [MIN MAX], []);
    %file_vol = histeq(file_vol, 2048);
% %     if ind==0
% %         file_vol = imadjustn(file_vol);
% %         hist = histogram(file_vol, 100000);
% %         Val = hist.Values;
% %         ind = 1;
% %     else
% %         file_vol = histeq(file_vol, Val);
% %     end
%     %N2_adj = imadjustn(N2_brut);
%     
    overlay = '';
    command = '';
    
    if ~isempty(files_in.In2)
        info = niftiinfo(roi);
        B = write_volume(cast(roi_vol, info.Datatype), roi_h, 'Axial');
        niftiwrite(B, [Path_nii, 'ROI.nii'], info);
        overlay = [Path_nii, 'ROI.nii'];
        command = char(databROIs.SequenceName(1));
        C = write_volume(file_vol, roi_h, 'Axial');
    else
        info = niftiinfo(file1);
        C = write_volume(file_vol, file_h, 'Axial');
    end
    info.Datatype = class(C);
    niftiwrite(C, [Path_nii, 'T2.nii'], info);
    qcdir = Path_html;
    inputdate = datestr(datetime('now'));
    subject = [char(databScans.Patient(i)), ' - ', char(databScans.Tp(i))];
    img = [Path_nii, 'T2.nii'];
    contrast = char(databScans.SequenceName(i));
    
    %z = [2 3 4 5 6 7 8 9 10 11 12 13 14];%2:round(size(file_vol,3)/4);
    if strcmp(opt.NumSlices, 'Max(12)')
        if isempty(files_in.In2) % No ROI/Cluster
            NSlicesTot = size(file_vol,3);
            if NSlicesTot<12
                z = 1:NSlicesTot;
            else
                Space = round(NSlicesTot/12);
                Indexes = round(Space/2):Space:NSlicesTot;
                z = unique(Indexes);
            end
        else
            NSlicesTot = size(roi_vol,3);
            if NSlicesTot<12
                z = 1:NSlicesTot;
            else
                Space = round(NSlicesTot/12);
                Indexes = round(Space/2):Space:NSlicesTot;
                z = unique(Indexes);
                if length(z)>12
                    z = z(1:12);
                end
            end
        end
    else
        z = str2num(opt.NumSlices); %#ok<ST2NM>
    end
    qc_write(qcdir,inputdate,subject,{img},{overlay},contrast,command, z)

end

rmdir(Path_nii, 's');
%web([Path_html, 'index.html']);
