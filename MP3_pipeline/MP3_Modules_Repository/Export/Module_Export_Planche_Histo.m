function [files_in,files_out,opt] = Module_Export_Planche_Histo(files_in,files_out,opt)

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
    module_option(:,7)   = {'Output_PNG_Folder','PNG_folder'};
    module_option(:,8)   = {'AutomaticJobsCreation', 'Yes'};
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
         {''}
        };
    user_parameter(:,2)   = {'   .Scan Anatomic','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the anatomical scan that will be displayed'};
    user_parameter(:,3)   = {'   .Scan Multiparametric','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the multiparametrical scan whose the field of view will be displayed on the anatomical scan.'};
    user_parameter(:,4)   = {'   .ROI','1ROIOr1Cluster','','',{'SequenceName'},'Optional',...
         'You can also select an ROI to get zoomed images. Ex: the brain'};
    user_parameter(:,5)   = {'   .Output PNG folder','char','','Output_PNG_Folder','', '','the name of the folder where will be saved your PNGs, inside your project folder.'};
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
fileA = files_in.In1{1};
fileM = files_in.In2{1};
if isfield(files_in, 'In3') && ~isempty(files_in.In3)
    fileROI = files_in.In3{1};
    ROI_h = spm_vol(fileROI);
    A_h = spm_vol(fileA);
    ROI = read_volume(ROI_h, A_h, 0, 'Axial');
else
    A_h = spm_vol(fileA);
    ROI = [];
end

M_h = spm_vol(fileM);

J_A = ReadJson(strrep(files_in.In1{1}, '.nii', '.json'));
%J_M = ReadJson(strrep(files_in.In2{1}, '.nii', '.json'));
AcqDate = J_A.AcquisitionDate.value{1};
Day = strsplit(AcqDate);
Day = Day{1};
%A_s = read_volume(A_h, A_h, 0, 'Saggital');
A_a = read_volume(A_h, A_h, 0, 'Axial');
%A_c = read_volume(A_h, A_h, 0, 'Coronal');

M_init = read_volume(M_h, M_h, 0, 'Axial');
M_Nb_Slices = size(M_init,3);
%M_s = read_volume(M_h, A_h, 0, 'Saggital');
M_a = read_volume(M_h, A_h, 0, 'Axial');
%M_c = read_volume(M_h, A_h, 0, 'Coronal');

A_Size_vox = diag(A_h.mat(1:3,1:3)).';
M_Size_vox = diag(M_h(1).mat(1:3,1:3)).';
Facteur = round(M_Size_vox(3)/A_Size_vox(3));

%%


if ~isempty(ROI)
    Vol_Clean = A_a.*ROI;
    [rows, columns] = find(sum(ROI,3));
    row1 = min(rows);
    row2 = max(rows);
    col1 = min(columns);
    col2 = max(columns);
    Vol_Clean = Vol_Clean(row1:row2, col1:col2,:);
    MidInd_11 = row1 + round(size(Vol_Clean,1)/2);
    MidInd_12 = row1 + round(size(Vol_Clean,1)/2+1);
    MidInd_21 = col1 + round(size(Vol_Clean,2)/2);
    MidInd_22 = col1 + round(size(Vol_Clean,2)/2+1);
    
else
    Vol_Clean = A_a;
    MidInd_11 = round(size(A_a,1)/2);
    MidInd_12 = round(size(A_a,1)/2+1);
    MidInd_21 = round(size(A_a,2)/2);
    MidInd_22 = round(size(A_a,2)/2+1);
end



fig = figure('Name','Planche Histologie','units','normalized','outerposition',[0 0 1 1]);
[ha, pos] = tight_subplot(3,4,[.001 .001],[.001 .001],[.001 .001]);
subplot(3,4,[1,2])
%[ha, pos] = tight_subplot(2,2,[.03 .01],[.1 .03],[.01 .01]);
%axes(ha(1))


[~, fnA, ~] = fileparts(fileA);
[~, fnM, ~] = fileparts(fileM);



image(squeeze(A_a(MidInd_11,:,:)), 'CDataMapping', 'scaled')
Ind_a = [];
for i=1:size(M_a,3)
    if any(M_a(MidInd_11,:,i))
        Ind_a = [Ind_a, i];
    end
end
colormap('Gray')
hold on
axis off
%for i=Ind_a(1):Facteur:Ind_a(end)
for i=Ind_a(1)-Facteur:Facteur:Ind_a(end)+Facteur % Pour tracer 2 rectangles de plus : un avant et un apres
    r = rectangle('Position', [i-0.5,0.5, Facteur, size(A_a,1)], 'EdgeColor', [1-eps 1 1], 'Linewidth', 3);
end
xlabel(['Anat: ', fnA], 'Interpreter', 'none')
title(['Rat ', char(opt.Table_in.Patient(1)),' -- Date : ', Day], 'FontSize', 20, 'Interpreter', 'none');
h=get(gca);
h.XAxis.Label.Visible='on';

%axes(ha(2))
subplot(3,4,[3,4])
image(squeeze(A_a(MidInd_12,:,:)), 'CDataMapping', 'scaled')
colormap('Gray')
hold on
axis off
%for i=Ind_a(1):Facteur:Ind_a(end)
for i=Ind_a(1)-Facteur:Facteur:Ind_a(end)+Facteur % Pour tracer 2 rectangles de plus : un avant et un apres
    r = rectangle('Position', [i-0.5,0.5, Facteur, size(A_a,1)], 'EdgeColor', [1-eps 1 1], 'Linewidth', 3);
end

if size(A_a, 1) == size(A_a,2)
    %title(['Numeros Tranches (HG-HD/BG-BD)', ' : ', num2str(MidInd_11), ' - ', num2str(MidInd_12), ' sur ', num2str(size(A_a,1)),' / ', num2str(MidInd_21), ' - ', num2str(MidInd_22), ' sur ', num2str(size(A_a,2))], 'FontSize', 16);
    title(['Epaisseur : ',num2str(A_Size_vox(3)), ' / ', num2str(M_Size_vox(3)), 'mm',' -- N° tranches', ' : ', num2str(MidInd_11), '&', num2str(MidInd_12), ' sur ', num2str(size(A_a,1))], 'FontSize', 11, 'Interpreter', 'none');
end
xlabel(['MultiPara: ', fnM], 'Interpreter', 'none')
h=get(gca);
h.XAxis.Label.Visible='on';

% axes(ha(3))
% 
% image(squeeze(A_a(:,MidInd_21,:)), 'CDataMapping', 'scaled');
% hold on
% axis off
% %for i=Ind_a(1):Facteur:Ind_a(end)
% for i=Ind_a(1)-Facteur:Facteur:Ind_a(end)+Facteur % Pour tracer 2 rectangles de plus : un avant et un apres
%     r = rectangle('Position', [i-0.5,0.5, Facteur, size(A_a,1)], 'EdgeColor', [1-eps 1 1], 'Linewidth', 3);
% end
% title('Vue Axiale')
% title(['Vue Saggitale - Sujet ', char(opt.Table_in.Patient(1)), ' -  ', Day], 'FontSize', 16, 'Interpreter', 'none');
% 
% 
% axes(ha(4))
% 
% image(squeeze(A_a(:,MidInd_22,:)), 'CDataMapping', 'scaled');
% hold on
% axis off
% %for i=Ind_a(1):Facteur:Ind_a(end)
% for i=Ind_a(1)-Facteur:Facteur:Ind_a(end)+Facteur % Pour tracer 2 rectangles de plus : un avant et un apres
%     r = rectangle('Position', [i-0.5,0.5, Facteur, size(A_a,1)], 'EdgeColor', [1-eps 1 1], 'Linewidth', 3);
% end
% [~, fnA, ~] = fileparts(fileA);
% [~, fnM, ~] = fileparts(fileM);
% title(['A : ', fnA, ' - M : ', fnM], 'Interpreter', 'none', 'FontSize', 16);

% %123
% %321retourné
% %132 de coté
% %231retourné
% %312 de coté
% %213 ax au sol
% Vol = permute(A_a, [2,1,3]);
% %Vol = flip(Vol,2);
% %figure;
% %subplot(2,2,4)
% axes(ha([2,4]))
% subplot(2,2,[2,4])
% colormap('Gray')
% s = slice(1:256, 1:256, 1:31,Vol, MidInd_1 ,MidInd_2, [Ind_c(1), Ind_c(end)]);
% s(1).MeshStyle = 'none';
% s(2).MeshStyle = 'none';
% s(3).MeshStyle = 'none';
% s(4).MeshStyle = 'none';
% hold on
% ax = gca;
% ax.Visible = 'Off';
% %ax.CameraPosition = [700, -1000, -326];
% for i=1:length(Ind_c)+1
%     p = plot3([0 0 256 256 0], [0 256 256 0 0], [Ind_c(1)+i-1, Ind_c(1)+i-1, Ind_c(1)+i-1, Ind_c(1)+i-1, Ind_c(1)+i-1], 'Linewidth', 2);
%     p.Color = [1 1 1];
% end
% 
% p = plot3([0 0], [Ind_c(1), Ind_c(end)+1], [256 256], 'Linewidth', 2);
% p.Color = [1 1 1];
% p = plot3([256 256], [Ind_c(1), Ind_c(end)+1], [256 256], 'Linewidth', 2);
% p.Color = [1 1 1];
% p = plot3([0 0], [Ind_c(1), Ind_c(end)+1], [0 0], 'Linewidth', 2);
% p.Color = [1 1 1];
% p = plot3([256 256], [Ind_c(1), Ind_c(end)+1], [0 0], 'Linewidth', 2);
% p.Color = [1 1 1];

if ~exist([opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG'], 'dir')
    mkdir([opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG']);
end
if ~exist([opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG', filesep, opt.Output_PNG_Folder, filesep], 'dir')
    mkdir([opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG', filesep, opt.Output_PNG_Folder, filesep]);
end

%file1 = [opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG', filesep, opt.Output_PNG_Folder, filesep, [char(opt.Table_in.Patient(1)), '_', char(opt.Table_in.Tp(1)), '_3D.png', ]];
%saveas(gcf, file1)
%close 'Vues Saggitale et Axiale'

% 
% 
% figure;
% Vol = A_a;
% colormap('Gray')
% s = slice(1:size(Vol,2), 1:size(Vol,1), 1:size(Vol,3),Vol, floor(size(Vol,2)/2)+2, floor(size(Vol,1)/2)+2,[Ind_c(1), Ind_c(end)]);
% s(1).MeshStyle = 'none';
% s(2).MeshStyle = 'none';
% s(3).MeshStyle = 'none';
% s(4).MeshStyle = 'none';
% hold on
% ax = gca;
% ax.Visible = 'Off';
% ax.CameraPosition = [400 -500 -300];
% 
% for i=1:length(Ind_c)+1
%     p = plot3([0 0 256 256 0], [0 256 256 0 0], [Ind_c(1)+i-1, Ind_c(1)+i-1, Ind_c(1)+i-1, Ind_c(1)+i-1, Ind_c(1)+i-1]);
%     p.Color = [1 0 0];
% end
% 
% p = plot3([0 0], [Ind_c(1), Ind_c(end)+1], [256 256]);
% p.Color = [1 0 0];
% p = plot3([256 256], [Ind_c(1), Ind_c(end)+1], [256 256]);
% p.Color = [1 0 0];
% p = plot3([0 0], [Ind_c(1), Ind_c(end)+1], [0 0]);
% p.Color = [1 0 0];
% p = plot3([256 256], [Ind_c(1), Ind_c(end)+1], [0 0]);
% p.Color = [1 0 0];



if M_Nb_Slices ~= 5
    warning('Plus de 5 tranches imagées en multiparamétrique. Il faut adapter le subplot')
end

if Facteur == 1 % Même résolution en z sur la carte anat que sur la multiparametrique.
    %figure('Name','Vues Coronales','units','normalized','outerposition',[0 0 1 1])
    %[ha, pos] = tight_subplot(2,4,[.01 .01],[.1 .01],[.01 .01]);
    
    for i=1:length(Ind_a)+3
        %axes(ha(i))
        subplot(3,4,4+i)
        imshow(flip(Vol_Clean(:,:,Ind_a(1)+i-3),2), [])
        if i==1
            title(['Sujet ', char(opt.Table_in.Patient(1)), ' Timepoint ', char(opt.Table_in.Tp(1)), ' - ', num2str(size(Vol_Clean,3)), ' tranches'], 'FontSize', 11, 'Interpreter', 'none');
        else
            title(['[TRANCHE ', num2str(Ind_a(1)+i-3), ']'], 'FontSize', 11);
%         else
%             title(['Tranche ', num2str(Ind_c(1)+i-3)], 'FontWeight', 'normal')
        end
    end
elseif Facteur == 2
    %figure('Name','Vues Coronales','units','normalized','outerposition',[0 0 1 1])
    %[ha, pos] = tight_subplot(2,4,[.01 .01],[.1 .01],[.01 .01]);
    for i=1:(round(length(Ind_a)/Facteur)+3)
        %axes(ha(i))
        subplot(3,4,4+i)
        imshow(flip(nanmean(Vol_Clean(:,:,(Ind_a(1)+(i-3)*2):(Ind_a(1)+i*2-5)),3),2), [])
        if i==1
            title(['Sujet ', char(opt.Table_in.Patient(1)), ' - Timepoint ', char(opt.Table_in.Tp(1)), ' - ', num2str(size(Vol_Clean,3)), ' tranches'], 'FontSize', 11, 'Interpreter', 'none');
        else
            title(['[TRANCHES ', num2str((Ind_a(1)+i*2-6):(Ind_a(1)+i*2-5)), ']'], 'FontSize', 11);
%         else
%             title(['Tranche ', num2str(Ind_c(1)+i-3)], 'FontWeight', 'normal')
        end
    end
else
    warning('Difference de resolution d''un facteur autre que 1 et 2. Cas non pris en charge.')
end

file = [opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG', filesep, opt.Output_PNG_Folder, filesep, [char(opt.Table_in.Patient(1)), '_', char(opt.Table_in.Tp(1)), '_Planche.png', ]];
set(fig, 'PaperOrientation', 'landscape')
saveas(gcf, file)
close 'Planche Histologie'
%print( fig,'-depsc', file)
% file2 = [opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG', filesep, opt.Output_PNG_Folder, filesep, [char(opt.Table_in.Patient(1)), '_', char(opt.Table_in.Tp(1)), '_Planche.png', ]];
% saveas(gcf, file2)
% 
% F1 = imread(file1);
% F2 = imread(file2);
% Image = [F1;F2];
% file3 = [opt.Table_in.Properties.UserData.MP3_data_path, filesep, 'PNG', filesep, opt.Output_PNG_Folder, filesep, [char(opt.Table_in.Patient(1)), '_', char(opt.Table_in.Tp(1)), '_Total.png', ]];
% imwrite(Image, file3)
% close 'Vues Coronales'



%close all





% 
% 
% L = findobj('Name','Vues Coronales');
% copyobj(L, findobj('Name', 'Vues Saggitale et Axiale'));


% % Load saved figures
% c=hgload(file1);
% k=hgload(file2);
% % Prepare subplots
% figure
% h(1)=subplot(1,2,1);
% h(2)=subplot(1,2,2);
% % Paste figures on the subplots
% copyobj(allchild(get(c,'CurrentAxes')),h(1));
% copyobj(allchild(get(k,'CurrentAxes')),h(2));
% % Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')





