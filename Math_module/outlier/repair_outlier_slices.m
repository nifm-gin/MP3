function output = repair_outlier_slices(param1,VOI_filename,alpha,interpolation,plotDiagnostics,detrend,MeanMap)
%
% repair_outlier_slices: Repair individual corrupted slices in MR time series
%
%    Parameters:
%    Image Files   - Image files to be analyzed.
%    MRI Parameters - List of .brkhdr files for (f)MRI acquisitions
%    Mask File     - Optional mask for analysis. Default is to use mean and
%                    standard deviation of image intensities over the
%                    entire volume. If a mask is provided, only voxels with
%                    a mask value above the threshold are included.
%                    Default: ''.
%    Mask Threshold - Threshold for the mask image. Default: 0.5.
%    Outlier Threshold - Threshold to be applied for the outlier detection.
%                    This corresponds to the alpha used in the Student
%                    t-distribution at the basis of the calculation of the
%                    Thompson Tau. This correponds approximately to the
%                    fraction of samples that will be excluded on average
%                    from normally distributed data. Default: 0.0005.
%                    Any value other than [0.0005 0.001 0.005 0.01 0.05]
%                    requires changing the code as the statistics toolbox
%                    becomes necessary.
%    Output Prefix - Prefix to add to newly created image series. Default:
%                    'o'.
%    Data Is ASL   - 'Yes', 'No', 'true', 'false', '0' or '1'. If true,
%                    even and odd images will be analyzed separately.
%    Exclude Fname - Name of file to store list of excluded volumes.
%                    Default: excluded_frames.mat
%
%    EXAMPLE:
%
%
%    <pipeline>
%       ...
%                 <sequence select='pCASL_CO2' group='no'>
%                     <module name='GFB_get_dir' id='asl_scan_dir'>
%                         <param name='levelName' value='sequence'/>
%                         <param name='returntype' value='cell' />
%                     </module>
%
%                     <module name='GFB_repair_outlier_slices' id='m1'>
%                         <param name='Image Files' input='HIERARCHY' filter='.*\.(?i)imgFormat$'/>
%                         <param name='MRI Parameters' input='HIERARCHY' filter='.*\.(?i)brkhdr$'/>
%                         <param name='Mask File' value="{fullfile('TemplateDirectory','MaskName')}"/>
%                         <param name='Data Is ASL' value="Yes"/>
%                         <param name='Output Prefix' value="'o'"/>
%                         <param name='Exclude Fname' value='excluded_frames.mat' />
%                     </module>
%                 </sequence>
%       ...
%    </pipeline>
%
% This function is a module of the GIN-fMRI-Batch pipeline.

% Written by Jan Warnking, Inserm, 06/2013
%
% $LastChangedDate: 2012-09-30 22:47:41 +0200 (dim., 30 sept. 2012) $
% $LastChangedBy: warnking $
%
% $Id: GFB_remove_outliers.m 264 2012-09-30 20:47:41Z warnking $
%

DEBUG = false;

% load the parameter
map_filename = param1{1};
fid=fopen(map_filename ,'r');
if fid>0
    fclose(fid);
    load(map_filename);
else
    logbook{numel(logbook)+1} = sprintf('Somthing wrong with the parameter\n##$%s\n##$',...
        map_filename);
    output = [];
    return
end
if numel(param1) > 1
    uvascim.imageExit = uvascim.image;
    uvascim.image = param1{2};
end

% load the ROI

if isempty(VOI_filename{:})
    VOI_matrix = ones(size(uvascim.image.reco.data,1),size(uvascim.image.reco.data,2),size(uvascim.image.reco.data,4));
else
    fid=fopen(VOI_filename{:} ,'r');
    if fid>0
        fclose(fid);
        load(VOI_filename{:});
    else
        logbook{numel(logbook)+1} = sprintf('Somthing wrong with the VOI file\n##$%s\n##$',...
            VOI_filename);
        output = [];
        return
    end
    VOI_matrix = zeros(size(uvascim.image.reco.data,1),size(uvascim.image.reco.data,2),size(uvascim.image.reco.data,4));
    for iii = 1:uvascim.image.reco.no_slices
        for jjj = 1:numel(uvascroi)
            if abs(uvascim.image.reco.fov_offsets(3,1,iii) - uvascroi(jjj).fov_offsets(3)) < 0.0001
                VOI_matrix(:,:,iii)=  imresize(uvascroi(jjj).value(:,:),[size(VOI_matrix, 1) size(VOI_matrix,2)],'bilinear');
            end
        end
    end
%     VOI_matrix(VOI_matrix == 0) = NaN;
end
% VOI_matrix = permute(repmat(VOI_matrix,[1 1 1 2 30]),[1 2 4 3 5]);
% tmp = double(uvascim.image.reco.data).*VOI_matrix;

%% initialize parameters
acqTR = uvascim.image.acq.tr/1000; % convert ms to seconds

imgR = uvascim.image;
imgD = uvascim.image.reco.data;
% imgD = tmp;
nVoxels = size(imgD,1)*size(imgD,2);
nEchoes = size(imgD,3);
nSlices = size(imgD,4);
nExpts = size(imgD,5);
imgD = reshape(imgD,[nVoxels,nEchoes,nSlices,nExpts]); % flatten 2D slices to a single spatial dimension
origD = imgD; % keep a copy for interpolation
% offsetScan = imgR.acq.fov_offsets(:,3)';

maskD = cell(nSlices,1);
for itSlice=1:nSlices
    maskD{itSlice} = logical(VOI_matrix(:,:,itSlice));
end

refImage = zeros(nVoxels,nEchoes,nSlices);
for iEcho = 1:nEchoes
    for iSlice = 1:nSlices
        refImage(maskD{iSlice},iEcho,iSlice) = median(imgD(maskD{iSlice},iEcho,iSlice,:),4);
        % get signal difference wrt the reference image
        imgD(:,iEcho,iSlice,:) = imgD(:,iEcho,iSlice,:) - repmat(refImage(:,iEcho,iSlice),[1 1 1 nExpts]);
    end
end

% % plot for test
% testimgD = reshape(imgD,[size(imgR.reco.data,1),size(imgR.reco.data,2),nEchoes,nSlices,nExpts]);
% testrefImage = reshape(refImage,[size(imgR.reco.data,1),size(imgR.reco.data,2),nEchoes,nSlices]);
% figure,imagesc(testimgD(:,:,1,1,1));
% figure,imagesc(testrefImage(:,:,1,1));

% detect local outliers using a Hampel filter with a large window
DX          = 8;                               % Window Half size
T           = 6;                               % Threshold
% Threshold   = 0.1;                           % AdaptiveThreshold
idxOutlier = cell(nEchoes,nSlices);
if plotDiagnostics
    f1 = figure;
    % calculate full screen size
    fPos = get(0,'ScreenSize');
    tBarHeight = 40;
    mBarHeight = 80;
    fPos = fPos + [0 tBarHeight 0 -tBarHeight-mBarHeight];
    set(f1,'Position',fPos);
end
X = 1:nExpts;
K = struct('RT',acqTR,'HParam',nExpts*acqTR,'row',X);
K = spm_filter(K);

for iEcho = 1:nEchoes
    for iSlice = 1:nSlices
        thisD = squeeze(imgD(maskD{iSlice},iEcho,iSlice,:));
        
        % first cycle: detrend, detect local outliers, construct
        % reference image series
        if detrend
            filtD = thisD - (K.X0*(K.X0'*thisD'))'; % detrend
        else
            filtD = thisD;
        end
        refD = median(filtD,2);
        filtD = filtD - repmat(refD,[1 nExpts]); % subtract offset
        [~,outlM] = hampel(X,nanmean(filtD),DX,3);
        [~,outlS] = hampel(X,nanstd(filtD),DX,3);
        outlM = find(outlM);
        outlS = find(outlS);
        
        % second cycle: reconstruct more robust reference image and
        % detrend nanmean and nanstd timecourses based only on valid points
        thisKeep = setdiff(X,union(outlM,outlS));
        refD  = nanmean(filtD(:,thisKeep),2); % we excluded outliers, so nanmean should be fine
        filtD = filtD - repmat(refD,[1 nExpts]); % subtract offset
        thisM = nanmean(filtD);
        thisS = nanstd(filtD);
        if plotDiagnostics
            figure(f1);
            ax = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+iEcho);
            plot(ax,thisM,'c','LineWidth', 2)
            hold on
            ax = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+iEcho+nEchoes);
            plot(ax,thisS,'c','LineWidth', 2)
            hold on
        end
        
        if detrend
            thisKeep = setdiff(X,outlM);
            % orthogonalize X0 over the selected volumes
            X0orth = K.X0;
            for iN = 1:size(K.X0,2)
                X0orth(:,iN) = X0orth(:,iN) - (X0orth(:,1:iN-1)*(X0orth(thisKeep,1:iN-1)'*X0orth(thisKeep,iN))); % detrend;
                X0orth(:,iN) = X0orth(:,iN) / sqrt(sum(X0orth(thisKeep,iN).^2));
            end
            thisM = thisM - (X0orth*(X0orth(thisKeep,:)'*thisM(thisKeep)'))'; % detrend
            
            thisKeep = setdiff(X,outlS);
            % orthogonalize X0 over the selected volumes
            X0orth = K.X0;
            for iN = 1:size(K.X0,2)
                X0orth(:,iN) = X0orth(:,iN) - (X0orth(:,1:iN-1)*(X0orth(thisKeep,1:iN-1)'*X0orth(thisKeep,iN))); % detrend;
                X0orth(:,iN) = X0orth(:,iN) / sqrt(sum(X0orth(thisKeep,iN).^2));
            end
            % correct for non-zero average of the std
            avgS = nanmean(thisS(thisKeep));
            thisS = thisS - (X0orth*(X0orth(thisKeep,:)'*(thisS(thisKeep)'-avgS)))'; % detrend
            
            if plotDiagnostics
                figure(f1);
                ax = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+iEcho);
                plot(ax,thisM,'m','LineWidth', 2)
                ax = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+iEcho+nEchoes);
                plot(ax,thisS,'m','LineWidth', 2)
            end
        end
        
        % detect global outliers with Thompson tau method using robust biweight
        % statistics. Use very strict threshold here
        idxMg = find_outliers_Thompson_new(thisM,alpha,'biweight');
        idxSg = find_outliers_Thompson_new(thisS,alpha,'biweight');
        idxGlobalOutlier = union(idxMg,idxSg);
        % detect global outliers at a more moderate level, to intersect with
        % local outlier detection
        idxMg = find_outliers_Thompson_new(thisM,alpha,'biweight');
        idxSg = find_outliers_Thompson_new(thisS,alpha,'biweight');
        % detect local outliers
        try
            [~,idxMl] = hampel(X,thisM,DX,T,'Outliers',idxMg);
            % [~,idxMl] = hampel(X,imgM,DX,T,'Adaptive',Threshold);
            % same thing for standard deviations
            [~,idxSl] = hampel(X,thisS,DX,T,'Outliers',idxSg);
            % [~,idxSl] = hampel(X,imgS,DX,T,'Adaptive',Threshold);
            % exclude only frames that are both local outliers and global outliers
            % at a moderate significance level
            idxMl = intersect(find(idxMl),idxMg);
            idxSl = intersect(find(idxSl),idxSg);
        catch
            idxMl = [];
            idxSl = [];
        end
        idxLocalOutlier = union(idxMl,idxSl);
        % plot results if outliers have been detected
        if plotDiagnostics
            figure(f1);
            ax = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+iEcho);
            find_outliers_Thompson_new(thisM,alpha,'biweight',ax);
            plot(ax,idxMl,thisM(idxMl), 'og', 'MarkerSize', 6, 'LineWidth', 2)
            if iEcho==2 % make y scale of adjacent axes equal for better comparison
                tagAx = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+1);
                yl = get(tagAx,'YLim');
                set(ax,'YLim',yl);
            end
            ax = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+iEcho+nEchoes);
            find_outliers_Thompson_new(thisS,alpha,'biweight',ax);
            plot(ax,idxSl,thisS(idxSl), 'og', 'MarkerSize', 6, 'LineWidth', 2)
            if iEcho==2 % make y scale of adjacent axes equal for better comparison
                tagAx = subplot(nSlices,nEchoes*2,(iSlice-1)*2*nEchoes+1+nEchoes);
                yl = get(tagAx,'YLim');
                set(ax,'YLim',yl);
            end
        end
        
        idxOutlier{iEcho,iSlice} = union(idxOutlier{iEcho,iSlice},union(idxLocalOutlier,idxGlobalOutlier));
%         allOutlier = [allOutlier idxOutlier{iEcho,iSlice}'];
    end
end

if plotDiagnostics
    tHeight = 40/fPos(4);
    FilePath = map_filename;
    patternFile = 'Image_Analyses_data';
    IndxPath = strfind(FilePath,patternFile);
    SavePath = FilePath(1:IndxPath+numel(patternFile));
    FileNameMat = FilePath(IndxPath+numel(patternFile)+1:end-4);
    p1 = get(subplot(nSlices,nEchoes*2,1),'Position');
    p2 = get(subplot(nSlices,nEchoes*2,2*nEchoes),'Position');
    tPos = [(p1(1)+p2(1)+p2(3))/2 1-(1-p1(2)-p1(4)-tHeight)/2 0 0];
    ah = annotation(f1,'textbox',tPos,'String',FilePath,'Interpreter','none','FitBoxToText','on');
    tPos = get(ah,'Position');
    tPos(1) = tPos(1) - tPos(3)/2;
    tPos(2) = tPos(2) + tPos(4)/2;
    set(ah,'Position',tPos);
    outlFigName = fullfile(SavePath,sprintf('Outlier_file-%s_Detrend-%s.pdf',FileNameMat,num2str(detrend)));
    printpdf(f1, outlFigName)
%     set(f1,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','A1');
%     saveas(f1,outlFigName);
end

%% manage volumes to be excluded from the analysis

% % detect volumes to exclude
% NoutlSlicesPerVol = zeros(nExpts,1);
% for iEcho = 1:nEchoes
%     for iSlice = 1:nSlices
%         NoutlSlicesPerVol(idxOutlier{iEcho,iSlice}) = NoutlSlicesPerVol(idxOutlier{iEcho,iSlice}) + 1;
%     end
% end
% % exclude volumes if more than half of the slices are corrupt (ignore 2
% % added border slices)
% outlierVolumes = NoutlSlicesPerVol > ceil((nSlices-2)/2);
% exclude(outlierVolumes) = true;
% extrapVols = [(1:find(NoutlSlicesPerVol==0,1,'first')-1) (find(NoutlSlicesPerVol==0,1,'last')+1:nSlices)];
% exclude(extrapVols) = true;
% 
% % prune slices that would be extrapolated from list of outliers
% allOutlier = [];
% for iEcho = 1:nEchoes
%     for iSlice = 1:nSlices
%         idxOutlier{iEcho,iSlice} = setdiff(idxOutlier{iEcho,iSlice},extrapVols);
%         allOutlier = [allOutlier idxOutlier{iEcho,iSlice}'];
%     end
% end

% nanmean of data with exclusion of only outlier not all scan of experiment

if numel(param1) > 1
    imgR = uvascim.imageExit;
    imgD = uvascim.imageExit.reco.data;
    nVoxels = size(imgD,1)*size(imgD,2);
    nEchoes = size(imgD,3);
    nSlices = size(imgD,4);
    nExpts = size(imgD,5);
    imgD = reshape(imgD,[nVoxels,nEchoes,nSlices,nExpts]); % flatten 2D slices to a single spatial dimension
    origD = imgD; % keep a copy for interpolation
end

if MeanMap
    MeanData = NaN(size(imgR.reco.data,1),size(imgR.reco.data,2),size(imgR.reco.data,3),size(imgR.reco.data,4),1);
    imgMean = imgR;
    imgMean.reco.fov_offsets = imgR.reco.fov_offsets(:,:,:,1);
    imgMean.reco.fov_orientation = imgR.reco.fov_orientation(:,:,:,1);
    imgMean.reco.fov_phase_orientation = imgR.reco.fov_phase_orientation(:,:,1);
    imgMean.reco.no_expts = numel(1);
    imgMean.reco.label = imgR.reco.label(:,:,1);
    imgMean.reco.phaselabel = imgR.reco.phaselabel(:,:,1);
    imgMean.reco.scaling_factor = imgR.reco.scaling_factor(:,:,1);
    imgMean.reco.scaling_offset = imgR.reco.scaling_offset(:,:,1);
    
    for iEcho = 1:nEchoes
        if numel(param1) > 1
            iEchoOut = min(iEcho,size(idxOutlier,1));
        else
            iEchoOut = iEcho;
        end
            
        for iSlice = 1:nSlices
            thisOutl = idxOutlier{iEchoOut,iSlice};
            %             if (numel(thisOutl)/numel(X) > 0.6)
            validTimepoints = setdiff(X,thisOutl);
            % orient matrix the right way for interp1 to operate along time dimension
            thisD = nanmean(squeeze(origD(:,iEcho,iSlice,validTimepoints)),2);
            MeanData(:,:,iEcho,iSlice,1) = reshape(thisD,[size(imgR.reco.data,1),size(imgR.reco.data,2),1,1,1]);
        end
    end
    imgMean.reco.data = MeanData;
    output{2} = imgMean;
end

% interpolate data
if interpolation==1
    for iEcho = 1:nEchoes
        if numel(param1) > 1
            iEchoOut = min(iEcho,size(idxOutlier,1));
        else
            iEchoOut = iEcho;
        end
        %     origIdx = iEcho:nEchoes:nImages;
        for iSlice = 1:nSlices
            %         thisOutl = intersect(origIdx,idxOutlier{iSlice}); % select the tag or ctl volumes among the outliers
            thisOutl = idxOutlier{iEchoOut,iSlice};
            if numel(thisOutl)>0
                validTimepoints = setdiff(X,thisOutl);
                % orient matrix the right way for interp1 to operate along time dimension
                thisD = squeeze(origD(:,iEcho,iSlice,validTimepoints))';
                thisD = interp1(validTimepoints,thisD,thisOutl,'spline');
                origD(:,iEcho,iSlice,thisOutl) = reshape(thisD',[nVoxels,1,1,numel(thisOutl)]);
            end
        end
    end
    
    % reshape interpolate data for uvasc
    outputData = reshape(origD,[size(imgR.reco.data,1),size(imgR.reco.data,2),nEchoes,nSlices,nExpts]);
    imgR.reco.data = outputData;
    
    % % plot for debug
    % testorigD = reshape(origD,[size(imgR.reco.data,1),size(imgR.reco.data,2),nEchoes,nSlices,nExpts]);
    % figure,imagesc(testorigD(:,:,1,1,1));
    
elseif interpolation==0
    allOutlier = str2num(char(cellfun(@num2str,idxOutlier,'uni',0)))'; %#ok<ST2NM>
    OutlierExpt = unique(allOutlier);
    validTimepoints = setdiff(X,OutlierExpt);
    origD = origD(:,:,:,validTimepoints);
    outputData = reshape(origD,[size(imgR.reco.data,1),size(imgR.reco.data,2),nEchoes,nSlices,numel(validTimepoints)]);
    imgR.reco.data = outputData;
    imgR.reco.fov_offsets = imgR.reco.fov_offsets(:,:,:,validTimepoints);
    imgR.reco.fov_orientation = imgR.reco.fov_orientation(:,:,:,validTimepoints);
    imgR.reco.fov_phase_orientation = imgR.reco.fov_phase_orientation(:,:,validTimepoints);
    imgR.reco.no_expts = numel(validTimepoints);
    imgR.reco.label = imgR.reco.label(:,:,validTimepoints);
    imgR.reco.phaselabel = imgR.reco.phaselabel(:,:,validTimepoints);
    imgR.reco.scaling_factor = imgR.reco.scaling_factor(:,:,validTimepoints);
    imgR.reco.scaling_offset = imgR.reco.scaling_offset(:,:,validTimepoints);
end

output{1} = imgR;

end

function [filename] =  printpdf(fig, name)
% printpdf Prints image in PDF format without tons of white space

% The width and height of the figure are found
% The paper is set to be the same width and height as the figure
% The figure's bottom left corner is lined up with
% the paper's bottom left corner

% Set figure and paper to use the same unit
set(fig, 'Units', 'centimeters')
set(fig, 'PaperUnits','centimeters');

% Position of figure is of form [left bottom width height]
% We only care about width and height
pos = get(fig,'Position');

% Set paper size to be same as figure size
set(fig, 'PaperSize', [pos(3) pos(4)]);

% Set figure to start at bottom left of paper
% This ensures that figure and paper will match up in size
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);

% Print as pdf
print(fig, '-dpdf', name)

% Return full file name
filename = [name, '.pdf'];
end
