function [BVf_Num_struct, stO2_Num_struct,   Vessel_orient_struct, ADC_Num_struct, R_Num_struct, R2_Num_struct]  = parametric_map_BVf_stO2_R(MGEFIDSEpre_filename, MGEFIDSEpost_filename,ADC_filename, Mask_filename, add_parameters)
   % Coucou 
% generate a 3 vascular maps using a numerical fit from a 2 scans (MGEFIDSEpre and MGEFIDSEpost contrast agent injection: 
%  i) BVf-Num (blood volume fraction)
%  ii)stO2-Num (saturation of tissu oxygenation)
%  iii) R-Num (vessel radii)
%
% this code come from Nicolas Pannetier
% additional_parameters correspond to B0 (T) and the TE (ms) used

% B0=str2double(add_parameters{:}(1));
% Mask_name = add_parameters{:}(2);
final_res = add_parameters{:}(3);
Hct = add_parameters{:}(4);
Hct = str2double(Hct{:});
With_mask = add_parameters{:}(5);
Dico_selected = add_parameters{:}(6);
MaxTimPt = str2double(add_parameters{:}(7));


if strcmp(final_res, 'Original')
   rescale = 0;
else
    rescale = 1;
    final_res = str2double(add_parameters{:}(3));
end


% load the MGEFIDSEpre file
fid=fopen(MGEFIDSEpre_filename ,'r');
if fid>0
    fclose(fid);
    MGEFIDSEpre = load(MGEFIDSEpre_filename);
    MGEFIDSEpre = MGEFIDSEpre.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the BVf_stO2_R maps because there is\n##$ Somthing wrong with the data\n##$MGEFIDSEpre*=%s\n##$',...
        MGEFIDSEpre_filename);
    msgbox(warning_text, 'BVf_stO2_R maps warning') ;
    BVf_Num_struct = []; stO2_Num_struct = []; R_Num_struct = []; R2_Num_struct = []; Vessel_orient_struct= []; ADC_Num_struct= [];
    return
end
% load the MGEFIDSEpost file
fid=fopen(MGEFIDSEpost_filename ,'r');
if fid>0
    fclose(fid);
    MGEFIDSEpost = load(MGEFIDSEpost_filename);
    MGEFIDSEpost = MGEFIDSEpost.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the BVf_stO2_R maps because there is\n##$ Somthing wrong with the data\n##$MGEFIDSEpost=%s\n##$',...
        MGEFIDSEpost_filename);
    msgbox(warning_text, 'BVf_stO2_R maps warning') ;
    BVf_Num_struct = []; stO2_Num_struct = []; R_Num_struct = []; R2_Num_struct = [];Vessel_orient_struct= []; ADC_Num_struct= [];
    return
end
% load the ADC file

fid=fopen(ADC_filename ,'r');
if fid>0
    fclose(fid);
    ADC = load(ADC_filename);
    ADC = ADC.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the  BVf_stO2_R maps because there is\n##$ Somthing wrong with the data\n##$ADC=%s\n##$',...
        ADC_filename);
    msgbox(warning_text, 'BVf_stO2_R maps warning') ;
    BVf_Num_struct = []; stO2_Num_struct = []; R_Num_struct = []; R2_Num_struct = [];Vessel_orient_struct= []; ADC_Num_struct= [];
    return
end

% Load mask (if selected and/or exist)

if ~isempty(Mask_filename)
    fid=fopen(Mask_filename ,'r');
    if fid>0
        fclose(fid);
        Mask = load(Mask_filename);
        Mask = Mask.uvascroi;
    else
        Mask = []; 
    end
else
    Mask = [];
end
 
 
temp_thickness = [MGEFIDSEpre.reco.thickness MGEFIDSEpost.reco.thickness  ADC.reco.thickness];
temp_slice_nbr = [MGEFIDSEpre.reco.no_slices MGEFIDSEpost.reco.no_slices  ADC.reco.no_slices];
temp_resolution = [MGEFIDSEpre.reco.no_samples MGEFIDSEpost.reco.no_samples  ADC.reco.no_samples];

% check data compatibility (slice thickness and slice number)
if  length(find(temp_thickness == MGEFIDSEpre.reco.thickness)) ~= numel(temp_thickness)
    warning_text = sprintf('##$ Can not calculate the BVf_stO2_R maps because there is\n##$ a slice thickness missmatch between\n##$MGEFIDSEpre=%s\n##$ and \n##$MGEFIDSEpost=%s\n##$ and\n##$ADC=%s\n',...
        MGEFIDSEpre_filename,MGEFIDSEpost_filename,ADC_filename);
    msgbox(warning_text, 'BVf_stO2_R maps warning') ;
    BVf_Num_struct = []; stO2_Num_struct = []; R_Num_struct = [];
    return
end
if length(find(temp_resolution == MGEFIDSEpre.reco.no_samples)) ~= numel(temp_resolution)
     warning_text = sprintf('##$ Can not calculate the BVf_stO2_R maps because there is\n##$ a resolution missmatch between\n##$MGEFIDSEpre=%s\n##$ and \n##$MGEFIDSEpost=%s\n##$ and\n##$ADC=%s\n',...
       MGEFIDSEpre_filename,MGEFIDSEpost_filename,ADC_filename);
    msgbox(warning_text, 'BVf_stO2_R maps warning') ;
    BVf_Num_struct = []; stO2_Num_struct = []; R_Num_struct = [];
    return
end

% get sequence timing
tout  = MGEFIDSEpre.reco.echotime*1e-3;

% Process dictionary 

% % if exist('Dic','var') ~= 1
% 	Dic = BuiltUpDico(DicoPrePath,DicoPostPath,tout);
if ispc
    
    switch Dico_selected{:}
        case 'TE60ms'
            load('D:\Google Drive\matlab\Dictionaries\dico_dense_se_60ms.mat');
        case {'dico_full_se_60ms', 'dico_full_se_60ms+ADC'}
            load('D:\Google Drive\matlab\Dictionaries\dico_full_se_60ms.mat');
       
    end
else
    switch Dico_selected{:}
        case 'TE60ms'
            load('/Users/blemasso/Google Drive/matlab/Dictionariesdico_dense_se_60ms.mat');
        case {'dico_full_se_60ms', 'dico_full_se_60ms+ADC'}
            load('/Users/blemasso/Google Drive/matlab/Dictionaries/dico_full_se_60ms.mat');
        case 'dico_full_se_60ms_extendeddiff'
            load('/Users/blemasso/Google Drive/matlab/Dictionaries/dico_full_se_60ms_extendeddiff.mat')
    end
end


% START TO BE ADAPTED
PostDataIn = MGEFIDSEpost.reco.data;
PreDataIn  = MGEFIDSEpre.reco.data;
ADCdata = ADC.reco.data;
[X, Y, E, Z] = size(PostDataIn);

% Change resolution (by adding the signal) to increase SNR
if rescale == 1
    xfactor = X/final_res;
    yfactor = Y/final_res;
    
    x = 1:xfactor:MGEFIDSEpost.reco.no_samples;
    y = 1:yfactor:MGEFIDSEpost.reco.no_views;
    PostDataIn_tmp = NaN(size(PostDataIn((x+xfactor-1)/2,(y+yfactor-1)/2,:,1:Z)));
    PreDataIn_tmp =PostDataIn_tmp;
    for i = 1:Z
        % add signal for MGEFIDSEpre
        PreDataIn_tmp((x+xfactor-1)/2,(y+yfactor-1)/2,:,i)=(PreDataIn(x,y,:,i)+...
            (PreDataIn(x+xfactor-1,y,:,i)+...
            PreDataIn(x,y+yfactor-1,:,i))+...
            PreDataIn(x+xfactor-1,y+yfactor-1,:,i));
        
        % add signal for MGEFIDSEpost
        PostDataIn_tmp((x+xfactor-1)/2,(y+yfactor-1)/2,:,i)=(PostDataIn(x,y,:,i)+...
            (PostDataIn(x+xfactor-1,y,:,i)+...
            PostDataIn(x,y+yfactor-1,:,i))+...
            PostDataIn(x+xfactor-1,y+yfactor-1,:,i));
        
        % merge voxels for ADC map
        ADCdata_tmp(:,:,1,i)=imresize(ADCdata(:,:,1,i),[final_res final_res],'bilinear');
    end
    PostDataIn = PostDataIn_tmp;
    PreDataIn= PreDataIn_tmp;
    ADCdata = ADCdata_tmp;
    clear PostDataIn_tmp PreDataIn_tmp ADCdata_tmp
end
    
if ~isempty(Mask)
    if size(Mask(1).value, 1) ~= X || size(Mask(1).value, 1) ~= Y
        for i=1:numel(Mask)
            tmp(i).value = imresize(Mask(i).value, [X, Y]);
        end
        Mask = rmfield(Mask, 'value');
        for i=1:numel(Mask)
            Mask(i).value = tmp(i).value;
        end
    end
    for i = 1:Z
        for j = 1:numel(Mask)
            if abs(MGEFIDSEpost.reco.fov_offsets(3,1,i) - Mask(j).fov_offsets(3)) < 1e-5
                 PostDataIn(:,:,:,i) = PostDataIn(:,:,:,i).*repmat(Mask(j).value,[1 1 E 1]);
                 PreDataIn(:,:,:,i) = PreDataIn(:,:,:,i).*repmat(Mask(j).value,[1 1 E 1]);
                 ADCdata(:,:,:,i) = ADCdata(:,:,:,i).*repmat(Mask(j).value,[1 1 1 1]);             
            end
        end
    end
end
%figure;imshow3D(squeeze(PostDataIn(:,:,1,:)))
% Squeeze temporal point to avoid bias induced by very load SNR for long
% echo time
% MaxTimPt = 26; % 23 semble pas mal
PostDataIn = PostDataIn(:,:,1:MaxTimPt,:);
PreDataIn = PreDataIn(:,:,1:MaxTimPt,:);
Dic.dat = Dic.dat(:,1:MaxTimPt); %#ok<NODEF>

[X, Y, E, Z] = size(PostDataIn);
PostDataIn = permute(PostDataIn,[1 2 4 3]);
PostDataIn = reshape(PostDataIn,[X*Y*Z E]);
PreDataIn = permute(PreDataIn,[1 2 4 3]);
PreDataIn = reshape(PreDataIn,[X*Y*Z E]);
% for i=1:size(PostDataIn,1)
%     if sum(PostDataIn(i,:)) ~= 0
%         PostDataIn(i,:) = smooth(PostDataIn(i,:), 'lowess')';
%         PreDataIn(i,:) = smooth(PreDataIn(i,:), 'lowess')';
%     end
% end

DataIn = PostDataIn./PreDataIn;

%%Compute qBOLD based 
%########################

% Reshape data
% DataIn = permute(DataIn,[1 2 4 3]);
% DataIn = reshape(DataIn,[X*Y*Z E]);
ADCdata = permute(ADCdata,[1 2 4 3]); 
if size(ADCdata,3) == 9
    ADCdata = ADCdata(:,:,3:7);
    disp([ADC_filename, ' : Only the slice 3 4 5 6 7 have been keept']);
end
ADCdata = reshape(ADCdata,[X*Y*Z 1]);
% disp('warning: ADC fixed at 800');
% ADCdata = 800*ones([size(ADCdata,1) 1]);
% Normalized DataIn and Dico by blocks: FID, SE
dtout = diff(tout);
Gap = find(single(dtout) ~= single(dtout(1))); % where the gap between FID and SE is 
Dic.dat(:,1:Gap) = bsxfun(@rdivide,Dic.dat(:,1:Gap),mean(squeeze(Dic.dat(:,1:Gap)),2));
Dic.dat(:,Gap+1:end) = bsxfun(@rdivide,Dic.dat(:,Gap+1:end),mean(squeeze(Dic.dat(:,Gap+1:end)),2));
DataIn(:,1:Gap)  = bsxfun(@rdivide,DataIn(:,1:Gap),mean(squeeze(DataIn(:,1:Gap)),2));
DataIn(:,Gap+1:end)  = bsxfun(@rdivide,DataIn(:,Gap+1:end),mean(squeeze(DataIn(:,Gap+1:end)),2));

%% Optimized with DICO ONLY (much faster, less accurate)
switch Dico_selected{:}
    case 'TE60ms'
        DataOut = qBoldL2DicoCompute(DataIn,ADCdata(:),Dic);
%         DataOut= qBoldL2DicoCompute_rmse(DataIn,ADCdata(:),Dic);

    case {'dico_full_se_60ms', 'dico_full_se_60ms_extendeddiff'}  %'DenseLargeVessel'
       DataOut = qBoldL2DicoCompute2(DataIn,ADCdata(:),Dic);
    case 'dico_full_se_60ms+ADC' 
        DataOut = qBoldL2DicoCompute_withD(DataIn,Dic);
   
    otherwise
        warndlg('This dico is not coded','Warning');
        return
end


%% Reshape Output
switch Dico_selected{:}
    case {'dico_full_se_60ms', 'dico_full_se_60ms_extendeddiff'}  %'DenseLargeVessel'
        DataOut = reshape(DataOut,[X Y Z 5]);
   case 'dico_full_se_60ms+ADC'
       DataOut = reshape(DataOut,[X Y Z 6]);
    otherwise
        DataOut = reshape(DataOut,[X Y Z 4]);
end
% permut slices <--> parameters
DataOut = permute(DataOut,[1 2 4 3]);
switch Dico_selected{:}
    case 'TE=75'
        %Convert Dk value in lSO2 values
        Dk0 = 0.264 * 4 * pi;
        DataOut(:,:,3,:) = (1- DataOut(:,:,3,:)/(Dk0 * Hct))*100;
        % Format radius and Vf
        DataOut(:,:,1,:) = DataOut(:,:,1,:);
        DataOut(:,:,2,:) = DataOut(:,:,2,:)*100;
    otherwise 
        Dk0 = 0.264 * 4 * pi * 1e-6;
        DataOut(:,:,3,:) = (1- DataOut(:,:,3,:)/(Dk0 * Hct))*100;
%         % Format BVf 
        DataOut(:,:,1,:) = DataOut(:,:,1,:)*100;
        % Format radius
        DataOut(:,:,2,:) = DataOut(:,:,2,:)*1e6;

        
% %         % Format indice 
%          DataOut(:,:,1,:) = DataOut(:,:,1,:);
%         % Format ADC-Num
%         DataOut(:,:,2,:) = DataOut(:,:,2,:)*1e12;
end

% Output: - DataOut: N x 4 containing the estimates in this order: R, Vf, Dk, r2
% save imformation
tempstruct=MGEFIDSEpre;
tempstruct.reco.no_echoes = 1;
tempstruct.reco.displayedecho = 1;
tempstruct.reco.displayedslice = 1;
tempstruct.reco.displayedexpt = 1;
tempstruct.reco.no_expts = 1;
tempstruct.reco.echotime = NaN('single');

%Adapt the fov offsets and orientation infos
tmpfov=MGEFIDSEpre.reco.fov_offsets(:,1,:,1);
tmpori=MGEFIDSEpre.reco.fov_orientation(:,1,:,1);
% tmpphaselabel=MGEFIDSEpre.reco.phaselabel(1,:,1);
if rescale == 1
tempstruct.reco.no_samples = final_res;
tempstruct.reco.no_views = final_res;
end
tempstruct.reco.fov_offsets = [];
tempstruct.reco.fov_orientation = [];
tempstruct.reco.scaling_factor = [];
tempstruct.reco.scaling_offset = [];
tempstruct.reco.label = {};
tempstruct.reco.phaselabel = {};

for m_slice=1:MGEFIDSEpre.reco.no_slices,
    tempstruct.reco.fov_offsets(:,:,m_slice)=squeeze(tmpfov(:,m_slice));
    tempstruct.reco.fov_orientation(:,:,m_slice)=squeeze(tmpori(:,m_slice));
    tempstruct.reco.label(1,m_slice)=MGEFIDSEpre.reco.label(1,m_slice,1);
    tempstruct.reco.phaselabel(1,m_slice)=MGEFIDSEpre.reco.phaselabel(1,m_slice,1);
    tempstruct.reco.scaling_factor(1,m_slice) = MGEFIDSEpre.reco.scaling_factor(1,m_slice,1);
    tempstruct.reco.scaling_offset(1,m_slice) = MGEFIDSEpre.reco.scaling_offset(1,m_slice,1);
end
tempstruct.reco.fov_phase_orientation=MGEFIDSEpre.reco.fov_phase_orientation(:,:,1);


%delete useless field
tempstruct.reco = rmfield(tempstruct.reco, 'angAP');
tempstruct.reco = rmfield(tempstruct.reco, 'angFH');
tempstruct.reco = rmfield(tempstruct.reco, 'angRL');
tempstruct.reco = rmfield(tempstruct.reco, 'angulation');
tempstruct.reco = rmfield(tempstruct.reco, 'echotime');
tempstruct.reco = rmfield(tempstruct.reco, 'iminfos');
tempstruct.reco = rmfield(tempstruct.reco, 'reco_meth');

ParamConfig=sprintf('##$QuantifMethod=BVf stO2 and R maps from a numerical analysis of the MGEFIDSE pre and post scans\n##$Magnetic field=%s\n##$TE(ms)=%s\n##$Final resolution=%s\n##$Hct used=%s\n##$Dico used=%s\n##$MaxTE (nbr) used=%s\n##$Raw scans used=\n%s\n%s\n%s',...
        add_parameters{:}{1},...
        add_parameters{:}{2},...
        add_parameters{:}{3},...
        add_parameters{:}{4},...
        add_parameters{:}{6},...
        add_parameters{:}{7},...
        MGEFIDSEpre_filename,MGEFIDSEpost_filename,ADC_filename);
tempstruct.reco.paramQuantif = ParamConfig;

tempstruct.reco = rmfield(tempstruct.reco, 'data');
BVf_Num_struct = tempstruct;
stO2_Num_struct = tempstruct;
R_Num_struct = tempstruct;
R2_Num_struct = tempstruct;

switch Dico_selected{:}
    
    case 'dico_full_se_60ms+ADC'
        Vessel_orient_struct = tempstruct;
        ADC_Num_struct = tempstruct;
    case {'dico_full_se_60ms', 'dico_full_se_60ms_extendeddiff'} 
        Vessel_orient_struct = tempstruct;
        ADC_Num_struct = [];
    otherwise
        Vessel_orient_struct = [];
        ADC_Num_struct = [];
end
% BVf-num
BVf_Num_struct.reco.data = DataOut(:,:,1,:);
BVf_Num_struct.reco.unit = {'%'};
BVf_Num_struct.reco.echo_label = {'BVf-num'};
BVf_Num_struct.reco.globalmax = 100;
BVf_Num_struct.reco.globalmin = 0;

% R-Num 
R_Num_struct.reco.data =DataOut(:,:,2,:);
R_Num_struct.reco.unit = {'um'};
R_Num_struct.reco.echo_label = {'R-Num '};
R_Num_struct.reco.globalmax = 1000;
R_Num_struct.reco.globalmin = 0;

% stO2-Num
stO2_Num_struct.reco.data =  DataOut(:,:,3,:);
stO2_Num_struct.reco.unit = {'%'};
stO2_Num_struct.reco.echo_label = {'stO2-Num'};
stO2_Num_struct.reco.globalmax = 100;
stO2_Num_struct.reco.globalmin = 0;

switch Dico_selected{:}
    case 'dico_full_se_60ms+ADC'
        % Vessel_orient_struct
        Vessel_orient_struct.reco.data =  DataOut(:,:,4,:);
        Vessel_orient_struct.reco.unit = {'degre'};
        Vessel_orient_struct.reco.echo_label = {'Vessel_orient-Num'};
        Vessel_orient_struct.reco.globalmax = 2;
        Vessel_orient_struct.reco.globalmin = 0;
        
        % ADC_Num_struct
        ADC_Num_struct.reco.data =  DataOut(:,:,5,:);
        ADC_Num_struct.reco.unit = {'%'};
        ADC_Num_struct.reco.echo_label = {'ADC-Num'};
        ADC_Num_struct.reco.globalmax = 3500;
        ADC_Num_struct.reco.globalmin = 0;
     case {'dico_full_se_60ms', 'dico_full_se_60ms_extendeddiff'} 
           % Vessel_orient_struct
        Vessel_orient_struct.reco.data =  DataOut(:,:,4,:);
        Vessel_orient_struct.reco.unit = {'degre'};
        Vessel_orient_struct.reco.echo_label = {'Vessel_orient-Num'};
        Vessel_orient_struct.reco.globalmax = 2;
        Vessel_orient_struct.reco.globalmin = 0;
end
% R2_Num_struct
switch Dico_selected{:}
    
    case {'dico_full_se_60ms', 'dico_full_se_60ms_extendeddiff'}  %'DenseLargeVessel'
        R2_Num_struct.reco.data =  DataOut(:,:,5,:);
     case 'dico_full_se_60ms+ADC'
        R2_Num_struct.reco.data =  DataOut(:,:,6,:);
    otherwise 
        R2_Num_struct.reco.data =  DataOut(:,:,4,:);
end

R2_Num_struct.reco.unit = {'a.u.'};
R2_Num_struct.reco.echo_label = {'R2-Num'};
R2_Num_struct.reco.globalmax = 1;
R2_Num_struct.reco.globalmin = 0;



