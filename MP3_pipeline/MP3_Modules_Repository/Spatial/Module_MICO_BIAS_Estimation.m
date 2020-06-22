function [files_in,files_out,opt] = Module_MICO_BIAS_Estimation(files_in,files_out,opt)

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
    module_option(:,3)   = {'OutputSequenceName','Extension'};
    module_option(:,4)   = {'iterNum_outer',8};
    module_option(:,5)   = {'iterCM',2};
    module_option(:,6)   = {'iter_b',1};
    module_option(:,7)   = {'q',1.5};
    module_option(:,8)   = {'th_bg',0};
    module_option(:,9)   = {'N_regions',3};
    module_option(:,10)   = {'tissueLabel','1 2 3'};
    module_option(:,11)   = {'output_filename_ext_Field','_BiasField'};
    module_option(:,12)   = {'output_filename_ext_Scan','_bc'};
    module_option(:,13)   = {'Segmentation','No'};
    module_option(:,14)   = {'output_filename_ext_seg','_Seg'};
    module_option(:,15)   = {'RefInput',1};
    module_option(:,16)   = {'InputToReshape',1};
    module_option(:,17)   = {'Table_in', table()};
    module_option(:,18)   = {'Table_out', table()};
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
        {
    'MICO algorithm as described in the Chunming Li, John C. Gore and Christos Davatzikos'' paper ''Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation''.'
    'The module is basically an application of the sequence provided by those authors on their webpages.'
    'This sequence consist in a preprocessing (thresholding, cropping ...) and then the application of the MICO algorithm.'
    }'};

    user_parameter(:,2)   = {'Reference Image','XScan','','', {'SequenceName'},'Mandatory',...
         'This is the image on which estimate the bias field map.'};
    user_parameter(:,3)   = {'Parameters','','','','','',''};
    user_parameter(:,4)   = {'   .Number of outer iterations','numeric','','iterNum_outer','','',...
        'Number of outer iterations'};
    user_parameter(:,5)   = {'   .Number of inner iterations for segmentation estimation','numeric','','iterCM','','',...
        'Number of inner iterations for segmentation estimation'};
    user_parameter(:,6)   = {'   .Number of bias estimation iterations','numeric','','iter_b','','',...
        'Number of bias computation iterations'};
    user_parameter(:,7)  = {'   .Fuzzifier','numeric','','q','','',...
        'Fuzzifier'};
    user_parameter(:,8)  = {'   .Threshold to remove background','numeric','','th_bg','','',...
        'Threshold to remove background. My default the threshold = 0 --> no threashold. But you can increase this number such as 5 to remove background voxels'};
    user_parameter(:,9)  = {'   .Number of regions used to perform the segmentation','numeric','','N_regions','','',...
        'Number of regions used to perform the segmentation'};
    user_parameter(:,10)  = {'   .Labels of the resulting segmentation','char','','tissueLabel','','',...
        'Labels of the resulting segmentation'};
    user_parameter(:,11)  = {'   .Bias Corrected Scan Filename extension','char','','output_filename_ext_Scan','','',...
        'Specify the string to be added to the filename of the bias corrected scan. Default extension is ''_bc''.'};
    user_parameter(:,12)  = {'   .Bias Field map Filename extension','char','','output_filename_ext_Field','','',...
        'Specify the string to be added to the filename of the bias field map. Default extension is ''_BiasField''.'};
    user_parameter(:,13)  = {'   .Save the segmentation Map ?','cell',{'No', 'Yes'},'Segmentation','','',...
        'The computation of the bias field map uses a segmentation computation. Save this resulting Cluster ?'};
    user_parameter(:,14)  = {'       .Segmentation Filename extension','char','','output_filename_ext_seg','','',...
        'Specify the string to be added to the filename of the segmentation Cluster. Default extension is ''_Seg''.'};


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


if isempty(files_out)
    for i=1:size(opt.Table_in,1)
        tag1 = opt.Table_in(i,:);
        tag1.IsRaw = categorical(0);   
        tag1.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            tag1.SequenceName = categorical(cellstr(opt.output_filename_ext_Field));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            tag1.SequenceName = categorical(cellstr([char(tag1.SequenceName), opt.output_filename_ext_Field]));
        end
        tag1.Filename = categorical(cellstr([char(tag1.Patient), '_', char(tag1.Tp), '_', char(tag1.SequenceName)]));
        f_out = [char(tag1.Path), char(tag1.Patient), '_', char(tag1.Tp), '_', char(tag1.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        tag2 = tag1;
        tag2.SequenceName = categorical(cellstr([char(opt.Table_in.SequenceName(i)), opt.output_filename_ext_Scan]));
        tag2.Filename = categorical(cellstr([char(tag2.Patient), '_', char(tag2.Tp), '_', char(tag2.SequenceName)]));
        f_out2 = [char(tag2.Path), char(tag2.Patient), '_', char(tag2.Tp), '_', char(tag2.SequenceName), '.nii'];
        files_out.In2{i} = f_out2;
        TMP_Table_out = [tag1; tag2];
        if strcmp(opt.Segmentation, 'Yes')
            tag3 = tag1;
            tag3.Type = categorical(cellstr('Cluster'));
            if strcmp(opt.OutputSequenceName, 'AllName')
                tag3.SequenceName = categorical(cellstr(opt.output_filename_ext_seg));
            elseif strcmp(opt.OutputSequenceName, 'Extension')
                tag3.SequenceName = categorical(cellstr([char(opt.Table_in.SequenceName(i)), opt.output_filename_ext_seg]));
            end
            tag3.Filename = categorical(cellstr([char(tag3.Patient), '_', char(tag3.Tp), '_', char(tag3.SequenceName)]));
            f_out3 = [char(tag3.Path), char(tag3.Patient), '_', char(tag3.Tp), '_', char(tag3.SequenceName), '.nii'];
            files_out.In3{i} = f_out3;
            TMP_Table_out = [TMP_Table_out; tag3];
        end
        opt.Table_out = [opt.Table_out; TMP_Table_out];
    end
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_MICO:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

% FixedImInfo = niftiinfo(files_in.In1{1});
% [path, name, ~] = fileparts(files_in.In1{1});
% FixedImJsonfile = [path, filesep, name, '.json'];
% fid = fopen(FixedImJsonfile, 'r');
% raw = fread(fid, inf, 'uint8=>char');
% fclose(fid);
% %raw = reshape(raw, 1,length(raw));
% FixedImJSON = jsondecode(raw);

iterNum_outer=opt.iterNum_outer;  % outer iteration
iterCM=opt.iterCM;  % inner interation for C and M
Iter_b=opt.iter_b;  % inner iteration for bias
q = opt.q;   % fuzzifier
th_bg = opt.th_bg;  %% threshold for removing background
N_region = opt.N_regions; %% number of tissue types, e.g. WM, GM, CSF
tissueLabel=opt.tissueLabel;

unique_name = tempname;
Fold = [opt.folder_out, unique_name, filesep];

while exist(Fold, 'dir')
    unique_name = tempname;
    Fold = [opt.folder_out, unique_name, filesep];
end
mkdir(Fold);


N_scan = length(files_in.In1);
for nn = 1:N_scan
    str=files_in.In1{nn};
    data = load_untouch_nii(str);
    %hist = data.hdr.hist;
    Img=data.img;
    Img = double(Img);      
    save(strcat(opt.folder_out, filesep, 'temp_info.mat'), 'data');
    clear data;    
    
    [x1, x2, y1, y2, z1, z2] = cropImg(Img, th_bg);   %% crop image
    [DimX1, DimY1, DimZ1]=size(Img);
    x1 = max(1,x1-2); x2 = min(DimX1,x2+2);
    y1 = max(1,y1-2); y2 = min(DimY1,y2+2);
    z1 = max(1,z1-2); z2 = min(DimZ1,z2+2);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Img3D=Img(x1:x2,y1:y2,z1:z2);
    clear Img;
    [DimX, DimY, DimZ] = size(Img3D);    
    ROI = (Img3D>th_bg);
    saveBasisOrder3_3D(ROI, Fold);    
    Img3D = Img3D.*ROI;
    %%%%%%%%%%%%%%%%%%%%%% Initialization
    
    A=max(Img3D(:));
    C= linspace(0.1,0.9,N_region)*A;
    b=ones(size(Img3D));
    
    
    M=rand(DimX, DimY, DimZ,N_region);
    a=sum(M,4);
    for k = 1 : N_region
        M(:,:,:,k)=M(:,:,:,k)./a;
    end
    clear a;
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    totaltime = 0;
    M_old = M; chg=10000;
    save(strcat(opt.folder_out, filesep, 'M_old.mat'), 'M_old');
    clear M_old;
   
    C_old =C;
    for n_iter= 1:iterNum_outer 
        n_iter
        tic
        [M, b, C]=   MICO_3D(Img3D,ROI,M,C,b,Iter_b,iterCM,q, Fold);
      
        totaltime = totaltime + toc
     
        [M, C]=sortMemC(M, C);
%         PC2d=zeros(size(Img3D(:,:,1)));
         PC3d=zeros(DimX1, DimY1, DimZ1);
%         N_slc=70;
%         for k=1:N_region
%             PC2d = PC2d +  tissueLabel(k)*M(:,:,N_slc,k);
%         end
%         pause(0.1);
%         figure(1);
%         subplot(1,2,1);
%         imagesc(Img3D(:,:,N_slc));colormap(gray);
%         title('a 2D slice of 3D image');
%         subplot(1,2,2);
%         imagesc(PC2d);colormap(gray);  
%         title('a 2D slice of 3D segmentation result');
        C_new = C/norm(C);

        chg= max(abs(C_new(:)-C_old(:)));
        C_old = C_new;
        if chg<0.0001    % check convergence 
            break
        end         
    end   
   
  
    %[U, C]=sortMemC(M, C);
    % save the segmentation result if needed
    if strcmp(opt.Segmentation, 'Yes')    
        U=maxMembership(M);
        clear M;
        Membership = zeros(DimX1, DimY1, DimZ1,N_region);
        Membership(x1:x2,y1:y2,z1:z2,:)=U;
        for k=1:N_region
            PC3d(:,:,:)=PC3d(:,:,:)+tissueLabel(k)*Membership(:,:,:,k);
        end
        % save the nii file
        info = niftiinfo(files_in.In1{1}); % use the head of files_in.In1
        info.Filename = files_out.In3{nn};
        info.Filemoddate = char(datetime('now'));
        niftiwrite(single(PC3d), files_out.In3{nn},  info)
        
       % save_nii(pc3d,files_out.In3{nn});
    end
    img_bc = zeros(DimX1, DimY1, DimZ1);
    img_bc(x1:x2,y1:y2,z1:z2)=Img3D./b;
    B_field = zeros(size(img_bc));
    B_field(x1:x2, y1:y2, z1:z2) = b;
   
    % save the B_field nifti file
    info = niftiinfo(files_in.In1{1});
    info.Filename = files_out.In1{1};
    info.Filemoddate = char(datetime('now'));
    niftiwrite(single(B_field), files_out.In1{nn},  info)
    
    % save the B_field Json file
    %% Json Processing
    [path, name, ~] = fileparts(files_in.In1{nn});
    jsonfile = [path, '/', name, '.json'];
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

    [path, name, ~] = fileparts(files_out.In1{nn});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)

    clear b Img3D;
    
     % save the Image_corrected nifti file
    info.Filename = files_out.In1{1};
    info.Filemoddate = char(datetime('now'));
    niftiwrite(single(img_bc), files_out.In2{nn},  info)
    
    % save the Image_corrected Json file
    %% Json Processing
    [path, name, ~] = fileparts(files_in.In1{nn});
    jsonfile = [path, '/', name, '.json'];
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
    [path, name, ~] = fileparts(files_out.In2{nn});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
    
    clear Membership Bias U image_bc;
    
    for basis_index=1:20  % delete basis files
        filename = [Fold, 'basis_',num2str(basis_index),'.mat'];
        delete(filename);
    end  
    delete(strcat(opt.folder_out, filesep, 'M_old.mat'));
    delete(strcat(opt.folder_out, filesep, 'temp_info.mat'));
    
end


% sort the constants c1, c2, ..., and change the order of the membership
% functions accordingly.
function [M_out, C_out]=sortMemC(M, C)

[C_out IDX]=sort(C);

for k = 1 : length(C)
    M_out(:,:,:,k) = M(:,:,:,IDX(k));
end

% This Matlab function binarizes a fuzzy membership function 
function M_out = maxMembership(M)

if size(M,4)==1
    
    N_class=size(M,3);
    ROI = (sum(M,3)>0);
    
    M_out = zeros(size(M));
    
    [e_min,N_min] = max(M,[], 3);   % do not consider 2 minimum in this version
    for kk=1:N_class
        M_out(:,:,kk) = ROI.*(N_min == kk);
    end
    
elseif size(M,4)>1
    N_class=size(M,4);    
    ROI = (sum(M,4)>0);
       
    [e_min,N_min] = max(M,[], 4);   % do not consider 2 minimum in this version
    for kk=1:N_class
        M_out(:,:,:,kk) = ROI.*(N_min == kk);
    end
else
    error('wrong dimension: maxMembership');
end
