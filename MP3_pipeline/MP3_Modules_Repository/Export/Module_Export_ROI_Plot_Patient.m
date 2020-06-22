function [files_in,files_out,opt] = Module_Export_ROI_Plot_Patient(files_in,files_out,opt)

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
    module_option(:,7)   = {'Output_Folder','Plot'};
    module_option(:,8)   = {'Normalisation','No'};
    module_option(:,9)   = {'AutomaticJobsCreation', 'No'};
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
         {'This module plot the mean values of the voxels of selected scans inside a ROI.'}
        };
    user_parameter(:,2)   = {'   .Scan','XScan','','', {'SequenceName'},'Optional',...
         'Please select the scans that will be analysed'};
    user_parameter(:,3)   = {'   .ROI','1ROI','','',{'SequenceName'},'Mandatory',...
         'Please select the ROI that will be analysed'};
    user_parameter(:,4)   = {'   .Output folder','char','','Output_Folder','', '','the name of the folder where will be saved your figure, inside your project folder.'};
    user_parameter(:,4)   = {'   .Normalization ?','cell',{'Yes', 'No'},'Normalisation','', '','It will plot the evolution of  the mean, between 0 and 1 for each scan.'};
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

UPat = unique(opt.Table_in.Patient);
FirstVoxels = cell(1,length(UPat));
LastVoxels = cell(1,length(UPat));
for i=1:length(UPat)
    FirstVox = [];
    LastVox = [];
    DBPat = opt.Table_in(opt.Table_in.Patient == UPat(i),:);
    UTp = unique(DBPat.Tp);
    TP = nan(1,length(UTp));
    UScans = unique(DBPat.SequenceName(DBPat.Type=='Scan',:));
    Means = NaN(length(UTp),length(UScans));
    STDs = NaN(length(UTp),length(UScans));
    for j=1:length(UTp)
        t = char(UTp(j));
        TP(j) = str2double(t(2:end));
        DBPatTp = DBPat(DBPat.Tp == UTp(j),:);
        ROIDB = DBPatTp(DBPatTp.Type == 'ROI',:);       
        if isempty(ROIDB)
            continue
        elseif height(ROIDB)>1
            warning('Plus d''une ROI trouvée, pas de figure générée')
            continue
        end
        ROIFilename = [char(ROIDB.Path), char(ROIDB.Filename), '.nii'];
        ROI_h = spm_vol(ROIFilename);
        ROI_vol = read_volume(ROI_h, ROI_h, 0, 'Axial');
        for k=1:length(UScans)
            DBScan = DBPatTp(DBPatTp.SequenceName == UScans(k),:);       
            if isempty(DBScan)
                continue
            elseif height(DBScan)>1
                warning('Plus d''un scan trouvé, pas de figure générée pour ce scan')
                continue
            end
            ScanFilename = [char(DBScan.Path), char(DBScan.Filename), '.nii'];
            scan_h = spm_vol(ScanFilename);
            scan_vol = read_volume(scan_h, ROI_h, 0, 'Axial');
            for l=1:size(scan_vol,3)
                if ~any(any(scan_vol(:,:,l)))
                    scan_vol(:,:,l) = nan(size(scan_vol(:,:,l)));
                end
            end
            
            Voxels = scan_vol(logical(ROI_vol(:)));
            Means(j, k) = nanmean(Voxels);
            STDs(j, k) = nanstd(Voxels);
            if isempty(FirstVox)
                FirstVox = [FirstVox;Voxels];
            end
            LastVox = Voxels;
        end
    end
    
    FirstVoxels{i} = FirstVox;
    LastVoxels{i} = LastVox;
    
    if any(isnan(TP))
        TP = 1:length(UTp);
    end
    MatTp = repmat(TP.', 1, length(UScans));
    if strcmp(opt.Normalisation, 'Yes')
        NMeans_tmp = Means - min(Means, [],1);
        NMeans = NMeans_tmp./max(NMeans_tmp, [],1);
        NSTDs = STDs./max(NMeans_tmp, [], 1);
        figure
        errorbar(MatTp, NMeans, [], NSTDs, '-*');
        xticks(TP)
        xticklabels(cellstr(UTp))
        warning off
        legend(cellstr(UScans), 'Interpreter', 'None')
        warning on
        title(['Normalized evolution of the mean (+std): ', char(UPat(i))])
        xlabel('Timepoints')
        ylabel('Normalized mean')
    else
        figure
        errorbar(MatTp, Means, [], STDs, '-*');
        xticks(TP)
        xticklabels(cellstr(UTp))
        warning off
        legend(cellstr(UScans), 'Interpreter', 'None')
        warning on
        title(['Evolution of the mean (+std): ', char(UPat(i))])
        xlabel('Timepoints')
        ylabel('Mean')
    end    
end





% 
% 
% for ind = 1:length(FirstVoxels)
%     m = nanmean(FirstVoxels{ind});
%     sig = nanstd(FirstVoxels{ind});
%     T = LastVoxels{ind};
%     T = (T-m)./sig;
%     [Hsw,psw,Wsw] = adtest(T);
%     if Hsw
%         [h, p, stats] = ttest(T);
%         if h
%             ['La médiane de la répartition des voxels de la tumeur au dernier jour sur le scan ', char(UScans), ' est significativement différente ce celle au premier jour, avec p = ' num2str(p)]
%         else
%             ['EH NON, car p = ' num2str(p)]
%         end
%     end
% end
% 
% test = 8;





