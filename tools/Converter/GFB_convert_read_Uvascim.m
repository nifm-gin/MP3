function [recdat, par] = GFB_convert_read_Uvascim(filename, par, applyRescale)
% function recdat = GFB_convert_read_rec(filename, par, applyRescale)
%
% Given the header data in par, load image data from filename (a REC file).
% If selectedScanes is not [], then only those scans specified will be
% loaded.
%
% (C)opyright 2005, Bennett Landman, bennett@bme.jhu.edu
% Revision History:
% Created: 2/11/2005
% Bug Fix: 5/04/2005 - Corrected function of "selectScans" parameter
% 26/8/10 remove option to load only selected dynamic scans (IT, JW)

if(~exist(filename,'file'))
    error(['loadREC: Unable to find file: ' filename]);
end

if ~isstruct(par) && isnan(par) % PAR file corrupt
    recdat = [];
    return
end

switch(par.scn.pix_bits)
    case 32
        type = 'float32';
        applyRescale=0; % not needed in this case
    case 16
        type = 'uint16'; % fixed?
    case 8
        type='uint8';
    otherwise
        error('Unsupported data type');
end

allInfo = [par.img.info];
% separate calculated from acquired images
convertCalculatedImages = false;

% The scan_seq field indicates the scanning sequence for each image,
% according to the following scheme (from rw_mr_datamodel__define.pro):
% ;
% ;   Image sequence is an ENUM which is translated into an integer.
% ;   Below is a list of the possible values
% ;
% ;   0   "IR"
% ;   1   "SE"
% ;   2   "FFE"
% ;   3   "DERIVED"
% ;   4   "PCA"
% ;   5   "UNSPECIFIED"
% ;   6   "SPECTRO"
% ;   7   "SI"
%
% Scanning sequence can be an actual sequence (values 0,1,2,6,7),
% or a calculation producing an image (values 3,4,5). The number of
% different sequences in an acquisition is the number of mixes.
% Unfortunately, there may be several sequence types with just one
% announced mix. Here we set num_mixes to the number of different scanning
% sequences for acquired data found.
acquiredSequences = [0,1,2,6,7];

% Look at announced scanning sequences
[imgMixIdcs,~,imgMlookup] = unique([allInfo.scan_seq]');
seqIsAcquired = ismember(imgMixIdcs,acquiredSequences);
par.max.num_mixes = sum(seqIsAcquired); % count acquired sequence types
imgIsAcquired = ismember([allInfo.scan_seq]',acquiredSequences);

% Look at types of acquired images (Real, Imag, Modulus, Phase)
[imgTypes,~,imgTlookup] = unique([allInfo(imgIsAcquired).image_type_mr]');
nImgTypes = length(imgTypes);
par.max.num_img_types = nImgTypes;

% Determine number of calculated images
if ~all(imgIsAcquired)
    notAcqInfo = allInfo(~imgIsAcquired);
    calcTypeSeq = [[notAcqInfo.image_type_mr]' [notAcqInfo.scan_seq]'];
    [calcTypes,~,calcTypeLookup] = unique(calcTypeSeq,'rows');
    nCalcImgTypes = size(calcTypes,1);
    calcTlookup = zeros(size(allInfo));
    calcTlookup(~imgIsAcquired) = calcTypeLookup;
    % check we understand how calculated images are organized
    if (nCalcImgTypes*par.max.num_slices==sum(~imgIsAcquired))
        convertCalculatedImages = true;
    else
        warning('Inconsistent number of calulated images. Skipping conversion of calculated images.')
    end
end

allSpecial = [par.img.special];
% figure out how many diffusion values/directions we have.
diffValues = [[allSpecial.number_of_diffusion_b_factor]' [allSpecial.gradient_orientation_number]'];
[imgDIdcs,~,imgDlookup] = unique(diffValues,'rows');
nDiffValues = size(imgDIdcs,1);

% figure out how many dynamic scans we really have (don't trust the header)
scnValues = [allInfo.dynamic_scan_num]';
[imgSIdcs,~,imgSlookup] = unique(scnValues);
nScanValues = size(imgSIdcs,1);
if nScanValues~=par.max.num_dynamics
    warning('Not all anounced dynamic scans found. Changing number and re-numerating. This is ugly.')
    par.max.num_dynamics = nScanValues;
end

% figure out how many label types we really have (don't trust the header)
if isfield(allInfo,'label_type')
    ltyValues = [allInfo.label_type]';
    [imgLIdcs,~,imgLlookup] = unique(ltyValues);
    nLtyValues = size(imgLIdcs,1);
    if nLtyValues~=par.max.num_label_types
        warning('Not all anounced label types found. Changing number and re-numerating. This is ugly.')
        par.max.num_label_types = nLtyValues;
    end
end

% % find number of slices in each output volume
% slices = zeros(size(scans));
% for j=1:length(scans)
%     slices(j) = max([allInfo(scanidx_backw==j).slice_num]);
% end

% See if we have all the info to rescale data, if requested
if(applyRescale),
    if(isfield(par.img,'vis'))
        try
            allVis = [par.img.vis];
            if(~isfield(allVis,'rescale_slope') || ~isfield(allVis,'rescale_intercept') || ~isfield(allVis,'scale_slope')),
                applyRescale=false;
            end
        catch ME
            warning(ME,'Availability of rescale info inconsistent across slices. Disabling rescaling.');
            applyRescale=false;
        end
    else
        applyRescale=false;
    end
end

% pre-fill data array in case slices are missing
recon_res = par.scn.recon_res(:)';
if(applyRescale)
    data_type = 'single';
else
    switch(par.scn.pix_bits)
        case 32
            data_type = 'single';
        case 16
            data_type = 'uint16';
        case 8
            data_type = 'uint8';
    end
end

% preallocate output variables
% recdat = cell(par.max.num_echo,par.max.card_phs,par.max.num_label_types,par.max.num_dynamics);
recdat.acquired = zeros([recon_res, par.max.num_slices, par.max.num_echo,...
    par.max.card_phs,par.max.num_label_types,par.max.num_dynamics,...
    nDiffValues,par.max.num_mixes,nImgTypes],data_type);
% create empty image structure of the same dimensions (without pixels)
% imgNew.acquired(par.max.num_slices, par.max.num_echo,...
%     par.max.card_phs,par.max.num_label_types,par.max.num_dynamics,...
%     nDiffValues,par.max.num_mixes,nImgTypes) = par.img(find(imgIsAcquired,1,'last'));
imgNew.acquired = reshape(par.img(imgIsAcquired),[par.max.num_slices, par.max.num_echo,...
    par.max.card_phs,par.max.num_label_types,par.max.num_dynamics,...
    nDiffValues,par.max.num_mixes,nImgTypes]);
if convertCalculatedImages
    recdat.calculated = zeros([recon_res, par.max.num_slices, nCalcImgTypes],data_type);
    imgNew.calculated = reshape(par.img(~imgIsAcquired),[par.max.num_slices, nCalcImgTypes]);
end
% Warning : missing dimension for diffusion data !

% fp = fopen(filename,'rb','ieee-le');
fp = fopen(filename);
if(fp==-1)
    error(['Unable to open: ' filename]);
end
img_len = prod(recon_res);
% img_size = img_len * par.scn.pix_bits/8; % length of one image in bytes
thisData = zeros(recon_res);
% load uvascim data
load(filename);

for j=1:length(par.img)
    if ~imgIsAcquired(j) && ~convertCalculatedImages
        continue % no need to waste time on data we won't use
    end
    

    % figure out where to put these data
    thisImg = par.img(j);
    thisSliceNum = thisImg.info.slice_num;
    if imgIsAcquired(j)
        thisEchoIdx = thisImg.info.echo_num;
        thisCdPhIdx = thisImg.info.cardiac_phase_num;
        if isfield(thisImg.info','label_type')
            thisLablIdx = imgLlookup(j);
        else
            thisLablIdx = 1;
        end
        thisScanIdx = imgSlookup(j); % thisImg.info.dynamic_scan_num;
        thisDValIdx = imgDlookup(j);
        thisMixIdx = imgMlookup(j);
        thisImgTypeIdx = imgTlookup(j);
        recdat.acquired(:,:,thisSliceNum,thisEchoIdx,thisCdPhIdx,thisLablIdx,thisScanIdx,thisDValIdx,thisMixIdx,thisImgTypeIdx) = thisData;
        imgNew.acquired(thisSliceNum,thisEchoIdx,thisCdPhIdx,thisLablIdx,thisScanIdx,thisDValIdx,thisMixIdx,thisImgTypeIdx) = thisImg;
    else
        thisTypeIdx = calcTlookup(j);
        recdat.calculated(:,:,thisSliceNum,thisTypeIdx) = thisData;
        imgNew.calculated(thisSliceNum,thisTypeIdx) = thisImg;
    end
    


    thisData(:) = squeeze(uvascim.image.reco.data(:,:,thisEchoIdx,thisSliceNum,thisLablIdx));
    thisData = thisData';
    thisData = flipdim(thisData, 1);
    if applyRescale %rescale?
        %             # === PIXEL VALUES =============================================================
        %             #  PV = pixel value in REC file, FP = floating point value, DV = displayed value on console
        %             #  RS = rescale slope,           RI = rescale intercept,    SS = scale slope
        %             #  DV = PV * RS + RI             FP = DV / (RS * SS)
        thisVis = thisImg.vis;
        SS = thisVis.scale_slope;
        RS = thisVis.rescale_slope;
        RI = thisVis.rescale_intercept;
        thisData = thisData/SS + (RI/(RS*SS));
    end
    
    if imgIsAcquired(j)
        recdat.acquired(:,:,thisSliceNum,thisEchoIdx,thisCdPhIdx,thisLablIdx,thisScanIdx,thisDValIdx,thisMixIdx,thisImgTypeIdx) = thisData;
        imgNew.acquired(thisSliceNum,thisEchoIdx,thisCdPhIdx,thisLablIdx,thisScanIdx,thisDValIdx,thisMixIdx,thisImgTypeIdx) = thisImg;
    else
        thisTypeIdx = calcTlookup(j);
        recdat.calculated(:,:,thisSliceNum,thisTypeIdx) = thisData;
        imgNew.calculated(thisSliceNum,thisTypeIdx) = thisImg;
    end
end
par.img = imgNew;
fclose(fp);

end