function outfiles=GFB_convert_ParUvascim2nii_BL(filelist, options)
% outfiles=GFB_convert_r2a(filelist, options)
%
% version 2.6.0, 17/07/09
%
% converts Philips Format (V3 and V4) to Analyze or Nifti format image
% volumes, used in spm2/spm99 and spm5, respectively. For more information,
% see http://www.fil.ion.ucl.ac.uk/spm. This function is part of the r2agui
% package with graphical interface (http://r2agui.sourceforge.net), but can
% be used standalone from the commandline or your own scripts (when all
% files of the r2agui package are in the matlab path). Files will
% be written to a seperate folder per PAR file (when indicated, see below),
% with optional prefix, to the same directiory as the PAR files are in or
% to an optional alternative directory. Filenames can be derived from PAR
% filename, or indicated by prefixes.
%
% requirements: none, (eg only matlab kernel).
%               previous versions of r2agui sometimes depended on
%               spm2/spm5, this one though is independent of spm_ functions
%               (but r2agui might try to use them for selection of files when
%               present).
%
% filelist: cell array (with {}) containing PAR file names (without path
% info) to convert into analyze/nifti
%
% options: structure containing options for conversion. Fields:
% options.prefix       : characters to prepend to all output folder and filenames,
%                         use '' when not needed. Cell array of prefixes
%                         can be used to use different prefix for each
%                         corresponding file in filelist.
% options.usefullprefix: when 1, do not append PARfilename to output files, use
%                        prefix only, plus filenumber
% options.pathpar      : complete path containing PAR files
% options.subaan       : when 1 checked, files will be written in a different
%                        subdirectory per PAR file, otherwise all to pathpar
% options.usealtfolder : when 1, files will be written to
%                        options.altfolder, including lowest level folder 
%                        containing parfile
% options.altfolder    : see above
% options.outputformat : 1 for Nifty output format (spm5), 2 for Analyze (spm2)
% options.angulation   : when 1: include affine transformation as defined in PAR
%                        file in hdr part of Nifti file (nifti only, EXPERIMENTAL!)
% options.rescale      : when 1: store intensity scale as found in PAR 
%                        file (assumed equall for all slices). Yields DV values. 
% options.dim          : when 3, single 3D nii files will be produced, when
%                        4, one 4D nii file will be produced, for example
%                        for time series or dti data
% options.dti_revertb0 : when 0 (default), philips ordering is used for DTI data 
%                       (eg b0 image last). When 1, b0 is saved as first image in 3D or 4D data
% options.label_as_dyn : when true (default), interleave label types and
%                        dynamic scans (sort by ascending L.ty)
% outfiles: list with all converted files

% Extensively modified by Jan Warnking, Inserm, 05/2012
%
% $LastChangedDate$
% $LastChangedBy$
%
% $Id$
%

SVNrev = '$Rev$';
SVNrev = SVNrev(7:end-2); %#ok<NASGU> % extract only string with revision number
thisFname = mfilename;
thisMID = [strrep(thisFname,'GFB_','GFB:'),':'];


%check some of the options
if ~isfield(options,'dim')
    options.dim=3;
end
if options.dim==4 && options.outputformat==2 % can't do 4D Analyze
    warning(thisMID,'Not possible to save to 4-D ANALYZE format. Falling back to 3D files.');
    options.dim = 3;
end
if ~isfield(options,'dti_revertb0')
    options.dti_revertb0=0;
end
% end of options checking

no_files=length(filelist);
prefixtemplate=options.prefix;
parPath=options.pathpar;

outfiles = cell(1,no_files(1));
for thisFileNo=1:no_files(1),
    if size(prefixtemplate,1)<=1
        if iscell(prefixtemplate)
            prefix=prefixtemplate{1};
        else
            prefix=prefixtemplate;
        end
    else
        prefix=prefixtemplate{thisFileNo};
    end
    
    % read all data
    parFile = fullfile(parPath,filelist{thisFileNo});
    parrec = GFB_convert_read_ParUvascim(parFile);
    if ~isstruct(parrec.hdr) && isnan(parrec.hdr) % PAR file corrupt
        outfiles{thisFileNo} = {};
        continue
    end
    
    % get some basic dimension info
    [nX, nY, nSlices, nEchoes, nCardiacPha, nLabelTypes, nDynamScans, nDiffValues, nMixes, nImgTypes] = size(parrec.scans.acquired);
    haveCalculatedImages = isfield(parrec.scans,'calculated');
    if haveCalculatedImages
        nCalcTypes = size(parrec.scans.calculated,4);
    end
    
    % save all parameters except image parameters separately
    globalParameters = rmfield(parrec.hdr,'img');
    imgParametersAcq = parrec.hdr.img.acquired;
    if haveCalculatedImages
        imgParametersCal = parrec.hdr.img.calculated;
    end
    % save some memory by removing original (now duplicated) header data
    parrec = rmfield(parrec,'hdr');
    
    % reorder to interleave label type with dynamic scans, if requested
    if nLabelTypes>1 && options.label_as_dyn,
        newDims = [nX,nY,nSlices,nEchoes,nCardiacPha,1,nLabelTypes*nDynamScans,nDiffValues,nMixes, nImgTypes];
        parrec.scans.acquired = reshape(parrec.scans.acquired,newDims);
        imgParametersAcq = reshape(imgParametersAcq,newDims(3:end));
        imgSpecial = [imgParametersAcq(1,1,1,1,:,1,1).special];
        dynTimes = [imgSpecial.dyn_scan_begin_time];
        oldTR = max(diff(dynTimes));
        newTR = oldTR/nLabelTypes;
        dynTimes = (0:nLabelTypes*nDynamScans-1)*newTR;
        nDynamScans = nDynamScans * nLabelTypes;
        nLabelTypes = 1;
        globalParameters.timing.dyn_scan_begin_times = dynTimes;
        % modify header to fake absence of label types
        for iImType = 1:nImgTypes
            for iMap = 1:nMixes
                for iDValue = 1:nDiffValues
                    for iDynScan = 1:nDynamScans
                        for iCardPha = 1:nCardiacPha
                            for iEcho = 1:nEchoes
                                for iSlice = 1:nSlices
                                    imgParametersAcq(iSlice,iEcho,iCardPha,1,iDynScan,iDValue,iMap,iImType).special.dyn_scan_begin_time = dynTimes(iDynScan);
                                    imgParametersAcq(iSlice,iEcho,iCardPha,1,iDynScan,iDValue,iMap,iImType).info.dynamic_scan_num = iDynScan;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % prepare matrix for list of files generated
    nOutfiles = numel(imgParametersAcq)/nSlices;
    if options.dim==4 && nDynamScans > 1 && nEchoes == 1 && nDiffValues == 1
        nOutfiles = nOutfiles/nDynamScans;
    elseif options.dim==4 && nDynamScans == 1 && nEchoes > 1 && nDiffValues == 1
         nOutfiles = nOutfiles/nEchoes;
    elseif options.dim==4 && nDynamScans == 1 && nEchoes == 1 && nDiffValues > 1
         nOutfiles = nOutfiles/nDiffValues;    
    end
    if haveCalculatedImages
        nOutfiles = nOutfiles + numel(imgParametersCal)/nSlices;
    end
    outfiles{thisFileNo} = cell(1,nOutfiles);
    outfcount=1;    
    
    % construct path and name of the data output
    if options.usealtfolder
        outfoldername = options.altfolder;
    else
        outfoldername = parPath;
    end
    if options.usefullprefix
        fNameRoot = prefix;
    else
        fNameRoot = [prefix,globalParameters.scn.protocol_name];
        % clean up file name somewhat
        fNameRoot = strrep(fNameRoot,'WIP','');
        fNameRoot = strrep(fNameRoot,'SENSE','');
        fNameRoot = strtrim(fNameRoot);
        fNameRoot = strrep(fNameRoot,' ','_');
    end
    if options.subaan,
        absDir = fullfile(outfoldername,fNameRoot);
        if ~exist(absDir,'dir')
            mkdir(absDir);
        end
    else
        absDir = outfoldername;
    end

    % construct description string of image
    Descrip=[globalParameters.scn.protocol_name,'; created by ',thisMID(1:end-1)];
    if length(Descrip)>80
        Descrip=Descrip(1:80); %truncate to 80
    else
        Descrip=[Descrip,char(ones(1,80-length(Descrip)))*32]; %#ok<AGROW> %add spaces to create string of 80
    end       
        
    % construct function to compute file names
    fnameFct1 = '@(s,i1,i2,i3,i4,i5,i6,i7)sprintf(''%s';
    fnameFct2 = ''',s';
    fnameFct3 = ')';
    if nEchoes > 1 && options.dim~=4
        fnameFct1 = [fnameFct1, '-e%02d']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i1']; %#ok<AGROW>
    end
    if nCardiacPha > 1
        fnameFct1 = [fnameFct1, '-c%02d']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i2']; %#ok<AGROW>
    end
    if nLabelTypes > 1
        fnameFct1 = [fnameFct1, '-l%01d']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i3']; %#ok<AGROW>
    end
    if nDiffValues > 1 && options.dim~=4
        fnameFct1 = [fnameFct1, '-b%03d']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i4']; %#ok<AGROW>
    end
    if nMixes > 1
        fnameFct1 = [fnameFct1, '-%s']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i5']; %#ok<AGROW>
    end
    if nImgTypes > 1
        fnameFct1 = [fnameFct1, '-%s']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i6']; %#ok<AGROW>
    end
    if nDynamScans > 1 && options.dim~=4
        fnameFct1 = [fnameFct1, '-d%04d']; %#ok<AGROW>
        fnameFct2 = [fnameFct2, ',i7']; %#ok<AGROW>
    end
    acqFnameFct = [fnameFct1 fnameFct2 fnameFct3];
    acqFnameFct = str2func(acqFnameFct);
    % convert acquired images
    % loop over all possible image dimensions and convert data, reading
    % header information individually per volume.
    for iType = 1:nImgTypes
        % assume that image type
        thisType = imgParametersAcq(1,1,1,1,1,1,1,iType).info.image_type_mr;
        for iMix = 1:nMixes
            thisMix = imgParametersAcq(1,1,1,1,1,1,iMix,iType).info.scan_seq;
            for iDValue = 1:nDiffValues
                % determine a common scale factor for all images of one image
                % type and one diffusion value.
                % First check if original scale factors are identical for all
                % images, with zero rescale intercept. In that case, use those
                % scale factors. There is no information to be gained by using
                % different ones.
                if globalParameters.scn.pix_bits == 32 % no rescaling required : floating point data originally
                    RS = 1;
                    RI = 0;
                    dataIsRescaled = false;
                else
                    thisVis = [imgParametersAcq(:,:,:,:,:,iDValue,iMix,iType).vis];
                    SSRSRI = [[thisVis.scale_slope]' [thisVis.rescale_slope]' [thisVis.rescale_intercept]'];
                    dSSRSRI = diff(SSRSRI);
                    % if we only have one image or all images have the same
                    % scaling parameters AND the intercept is zero (not
                    % sure nifti handles intercept correctly even though it
                    % is supported - analyze does not support it)
                    if (size(dSSRSRI,1)==0 || all(dSSRSRI(:)==0)) && all(SSRSRI(:,3)==0)
                        % go back to PV fom FV, to store data efficiently
                        % in the original integer format, without loss of
                        % information.
                        RS = 1/SSRSRI(1,1);
                        RI = 0;

                        parrec.scans.acquired(:,:,:,:,:,:,:,iDValue,iMix,iType) = round(parrec.scans.acquired(:,:,:,:,:,:,:,iDValue,iMix,iType)/RS);
                        dataIsRescaled = true;
                    else
                        % do not rescale differently, store as floating point
                        % to not loose any information
                        RS = 1;
                        RI = 0;
                        dataIsRescaled = false;
                    end
                end
                thisOptions = options;
                thisOptions.rescale = dataIsRescaled;
                for iLabType = 1:nLabelTypes
                    for iCardPha = 1:nCardiacPha
                        for iEcho = 1:nEchoes
                            if options.dim==4 && nDynamScans > 1
                                VolNameSinExt = fullfile(absDir,acqFnameFct(fNameRoot,iEcho,iCardPha,iLabType,iDValue,imgSstr(thisMix),imgTstr(thisType),1));
                                VolData = squeeze(parrec.scans.acquired(:,:,:,iEcho,iCardPha,iLabType,:,iDValue,iMix,iType));
                                VolName=[VolNameSinExt,'.nii'];
                                parameters = globalParameters;
                                parameters.img = imgParametersAcq(:,iEcho,iCardPha,iLabType,:,iDValue,iMix,iType);
                                % Adapt parameters for the scaling
                                if ~thisOptions.rescale % store any data that is not to be rescaled as floating point
                                    parameters.scn.pix_bits = 32;
                                end
                                % Apply rescaling determined above
                                parameters.img(1).vis.rescale_slope = RS;
                                parameters.img(1).vis.rescale_intercept = RI;
                                NHdr=CreateNiiHdr(parameters,thisOptions,Descrip); % for later use in Nii writing
                                WriteNii(VolName,NHdr,VolData);
                                outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
                            elseif options.dim==4 && nEchoes > 1
                                VolNameSinExt = fullfile(absDir,acqFnameFct(fNameRoot,iEcho,iCardPha,iLabType,iDValue,imgSstr(thisMix),imgTstr(thisType),1));
                                VolName=[VolNameSinExt,'.nii'];
                                if exist(VolName, 'file') == 2
                                    continue
                                end
                                VolData = squeeze(parrec.scans.acquired(:,:,:,:,iCardPha,iLabType,1,iDValue,iMix,iType));
                                parameters = globalParameters;
                                parameters.img = imgParametersAcq(:,:,iCardPha,iLabType,1,iDValue,iMix,iType);
                                % Adapt parameters for the scaling
                                if ~thisOptions.rescale % store any data that is not to be rescaled as floating point
                                    parameters.scn.pix_bits = 32;
                                end
                                % Apply rescaling determined above
                                parameters.img(1).vis.rescale_slope = RS;
                                parameters.img(1).vis.rescale_intercept = RI;
                                NHdr=CreateNiiHdr(parameters,thisOptions,Descrip); % for later use in Nii writing
                                WriteNii(VolName,NHdr,VolData);
                                outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
                                
                            elseif options.dim==4 && nDiffValues > 1
                                VolNameSinExt = fullfile(absDir,acqFnameFct(fNameRoot,iEcho,iCardPha,iLabType,iDValue,imgSstr(thisMix),imgTstr(thisType),1));
                                VolName=[VolNameSinExt,'.nii'];
                                if exist(VolName, 'file') == 2
                                    continue
                                end
                                VolData = squeeze(parrec.scans.acquired(:,:,:,iEcho,iCardPha,iLabType,1,:,iMix,iType));
                                parameters = globalParameters;
                                parameters.img = imgParametersAcq(:,iEcho,iCardPha,iLabType,1,:,iMix,iType);
                                % Adapt parameters for the scaling
                                if ~thisOptions.rescale % store any data that is not to be rescaled as floating point
                                    parameters.scn.pix_bits = 32;
                                end
                                % Apply rescaling determined above
                                parameters.img(1).vis.rescale_slope = RS;
                                parameters.img(1).vis.rescale_intercept = RI;
                                NHdr=CreateNiiHdr(parameters,thisOptions,Descrip); % for later use in Nii writing
                                WriteNii(VolName,NHdr,VolData);
                                outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
                                
                                
                            else
                                for iDynScan = 1:nDynamScans
                                    [~, parfilename, ~] = fileparts(parFile);
                                    %                                     VolNameSinExt = fullfile(absDir,acqFnameFct(fNameRoot,iEcho,iCardPha,iLabType,iDValue,imgSstr(thisMix),imgTstr(thisType),iDynScan));
                                    VolNameSinExt = fullfile(absDir,parfilename);
                                    VolData = parrec.scans.acquired(:,:,:,iEcho,iCardPha,iLabType,iDynScan,iDValue,iMix,iType);
                                    switch(options.outputformat)
                                        case 1 %nii
                                            VolName=[VolNameSinExt,'.nii'];
                                            parameters = globalParameters;
                                            parameters.img = imgParametersAcq(:,iEcho,iCardPha,iLabType,iDynScan,iDValue,iMix,iType);
                                            % Adapt parameters for the scaling
                                            if ~thisOptions.rescale % store any data that is not to be rescaled as floating point
                                                parameters.scn.pix_bits = 32;
                                            end
                                            % Apply rescaling determined above
                                            parameters.img(1).vis.rescale_slope = RS;
                                            parameters.img(1).vis.rescale_intercept = RI;
                                            NHdr=CreateNiiHdr(parameters,thisOptions,Descrip); % for later use in Nii writing
                                            WriteNii(VolName,NHdr,VolData);
                                            outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
                                        case 2 %analyze
                                            VolName=[VolNameSinExt,'.img'];
                                            
                                            %MatName=[VolNameSinExt];
                                            thisImgIdx = {1,iEcho,iCardPha,iLabType,iDynScan,iDValue,iMix,iType};
                                            if imgParametersAcq(thisImgIdx{:}).info.pix_bits <= 16
                                                Precision = sprintf('int%d', imgParametersAcq(thisImgIdx{:}).info.pix_bits);
                                            elseif imgParametersAcq(thisImgIdx{:}).info.pix_bits == 32
                                                Precision = 'float';
                                            elseif imgParametersAcq(thisImgIdx{:}).info.pix_bits == 64
                                                Precision = 'double';
                                            end
                                            Type = spm_type(Precision);
                                            Vox = [imgParametersAcq(thisImgIdx{:}).special.pix_spacing(:)' imgParametersAcq(thisImgIdx{:}).info.slicethk+imgParametersAcq(thisImgIdx{:}).info.slicegap];
                                            Dim = [imgParametersAcq(thisImgIdx{:}).info.recon_res(:)' parameters.max.num_slices];
                                            if length(Descrip)>80
                                                Descrip=Descrip(1:80); %truncate to 80
                                            else
                                                Descrip=[Descrip,char(ones(1,80-length(Descrip)))*32]; %#ok<AGROW> %add spaces to create string of 80
                                            end
                                            ID2 = fopen(VolName, 'w');
                                            fwrite(ID2, VolData,Precision);
                                            fclose(ID2);
                                            
                                            P = VolName;
                                            spm_hwrite(P,Dim,Vox,RS,Type,0,round(Orign),Descrip);
                                            %save(MatName,'M'); % add matfile holding angulation
                                            outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
                                    end % ... switch(options.outputformat)
                                end
                            end
                        end
                    end
                    disp(['Done writing file: ' VolName]);
                end
            end
        end
    end
    
    if haveCalculatedImages
        fnameFct1 = '@(s,i1)sprintf(''%s-map';
        fnameFct2 = ''',s';
        fnameFct3 = ')';
        if nCalcTypes > 1
            fnameFct1 = [fnameFct1, '%02d']; %#ok<AGROW>
            fnameFct2 = [fnameFct2, ',i1']; %#ok<AGROW>
        end
        calFnameFct = [fnameFct1 fnameFct2 fnameFct3];
        calFnameFct = str2func(calFnameFct);
        for iMap = 1:nCalcTypes
            % determine a common scale factor for all images of one image
            % type and one diffusion value.
            % First check if original scale factors are identical for all
            % images, with zero rescale intercept. In that case, use those
            % scale factors. There is no information to be gained by using
            % different ones.
            if globalParameters.scn.pix_bits == 32 % no rescaling required : floating point data originally
                RS = 1;
                RI = 0;
                dataIsRescaled = false;
            else
                thisVis = [imgParametersCal(:,iMap).vis];
                SSRSRI = [[thisVis.scale_slope]' [thisVis.rescale_slope]' [thisVis.rescale_intercept]'];
                dSSRSRI = diff(SSRSRI);
                if (size(dSSRSRI,1)==0 || all(dSSRSRI(:)==0)) && all(SSRSRI(:,3)==0)
                    RS = 1/SSRSRI(1,1);
                    RI = 0;
                    parrec.scans.calculated(:,:,:,iMap) = round(parrec.scans.calculated(:,:,:,iMap)/RS);
                    dataIsRescaled = true;
                else
                    % do not rescale differently, store as floating point
                    % to not loose any information
                    RS = 1;
                    RI = 0;
                    dataIsRescaled = false;
                end
            end
            thisOptions = options;
            thisOptions.rescale = dataIsRescaled;
            VolNameSinExt = fullfile(absDir,calFnameFct(fNameRoot,iMap));
            VolData = parrec.scans.calculated(:,:,:,iMap);
            if numel(size(VolData),4:5) == 2
                thisOptions.dim = 3;
            end
            switch(options.outputformat)
                case 1 %nii
                    VolName=[VolNameSinExt,'.nii'];
%                     disp(['Writing file: ' VolName]);
                    parameters = globalParameters;
                    parameters.img = imgParametersCal(:,iMap);
                    % Adapt parameters for the scaling
                    if ~thisOptions.rescale % store any data that is not to be rescaled as floating point
                        parameters.scn.pix_bits = 32;
                    end
                    % Apply rescaling determined above
                    parameters.img(1).vis.rescale_slope = RS;
                    parameters.img(1).vis.rescale_intercept = RI;
                    NHdr=CreateNiiHdr(parameters,thisOptions,Descrip); % for later use in Nii writing
                    WriteNii(VolName,NHdr,VolData);
                    outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
                case 2 %analyze
                    VolName=[VolNameSinExt,'.img'];
%                     disp(['Writing file: ' VolName]);
                    
                    %MatName=[VolNameSinExt];
                    if imgParametersCal(1,iMap).info.pix_bits <= 16
                        Precision = sprintf('int%d', imgParametersCal(1,iMap).info.pix_bits);
                    elseif imgParametersAcq(thisImgIdx{:}).info.pix_bits == 32
                        Precision = 'float';
                    elseif imgParametersAcq(thisImgIdx{:}).info.pix_bits == 64
                        Precision = 'double';
                    end
                    Type = spm_type(Precision);
                    Vox = [imgParametersCal(1,iMap).special.pix_spacing(:)' imgParametersCal(1,iMap).info.slicethk+imgParametersCal(1,iMap).info.slicegap];
                    Dim = [imgParametersCal(1,iMap).info.recon_res(:)' parameters.max.num_slices];
                    if length(Descrip)>80
                        Descrip=Descrip(1:80); %truncate to 80
                    else
                        Descrip=[Descrip,char(ones(1,80-length(Descrip)))*32]; %#ok<AGROW> %add spaces to create string of 80
                    end
                    ID2 = fopen(VolName, 'w');
                    fwrite(ID2, VolData,Precision);
                    fclose(ID2);
                    
                    P = VolName;
                    spm_hwrite(P,Dim,Vox,RS,Type,0,round(Orign),Descrip);
                    %save(MatName,'M'); % add matfile holding angulation
                    outfiles{thisFileNo}{outfcount}=VolName;outfcount=outfcount+1;
            end % ... switch(options.outputformat)
        end
        disp(['Done writing file: ' VolName]);
    end    
end % ...for thisFileNo=1:no_files(1),
end

function NHdr=CreateNiiHdr(pars,options,Descrip)
%create Nifti header from parameters as read from PAR file
GFB_convert_nifti_defines;
[M,realVoxSize,matrixSize]=calc_angulation(pars,options);
[qoffset_xyz quatern_bcd,qfac]=nifti_mat44_to_quatern(M); %calc quaternion (copied from dcm2nii)

NHdr.HdrSz=struct('val',348,'prec','int32');
NHdr.Data_Type=struct('val',char(ones(1,10)*32),'prec','char');%unused
NHdr.db_name=struct('val',char(ones(1,18)*32),'prec','char');%unused
NHdr.extents=struct('val',0,'prec','int32'); %unused
NHdr.session_error=struct('val',0,'prec','int16');%unused
NHdr.regular=struct('val',114,'prec','char'); %unused: in Analyze 7.5 this must be 114
NHdr.dim_info=struct('val',0,'prec','char'); %MRI slice order
if numel(pars.timing.dyn_scan_begin_times) > 1
    TR = diff(pars.timing.dyn_scan_begin_times(1:2));
else
    TR = 1; % for anatomicals etc, RT is set to NAN, probably doesnt make sense here, so set to 1
end
switch options.dim
    case 3
        NHdr.dim=struct('val',[options.dim matrixSize' 1 1 1 1],'prec','int16'); %Data array dimensions
        pixdim=struct('val',[qfac realVoxSize TR 1 1 1],'prec','float'); 
    case 4
        if pars.max.num_dynamics > 1
            NHdr.dim=struct('val',[options.dim matrixSize' pars.max.num_dynamics 1 1 1],'prec','int16'); %Data array dimensions
            pixdim=struct('val',[qfac realVoxSize TR 1 1 1],'prec','float');
        elseif pars.max.num_echo > 1
            NHdr.dim=struct('val',[options.dim matrixSize' pars.max.num_echo 1 1 1],'prec','int16'); %Data array dimensions
            pixdim=struct('val',[qfac realVoxSize TR 1 1 1],'prec','float');
        elseif pars.max.num_diffusion_values > 1
            NHdr.dim=struct('val',[options.dim matrixSize' pars.max.num_diffusion_values 1 1 1],'prec','int16'); %Data array dimensions
            pixdim=struct('val',[qfac realVoxSize TR 1 1 1],'prec','float');
        end
    otherwise
        error('N-D conversion not supported.');
end
NHdr.intent_p123=struct('val',[0 0 0],'prec','float'); %intent_p1, intent_p2, intent_p3: single;
NHdr.intent_code=struct('val',0,'prec','int16');
switch(pars.scn.pix_bits)
    case 8
        dt=2;
    case 16
        dt=4;
    case 32
        dt=16;
    case 64
        dt=64;
end
NHdr.datatype=struct('val',dt,'prec','int16');
NHdr.bitpix=struct('val',pars.scn.pix_bits,'prec','int16');
NHdr.slice_start=struct('val',0,'prec','int16');

% use real voxel dimensions as calculated from FOV/matrixsize in approp direction (CHECK!). 
% Because for older Philips releases, voxel dimensions in PAR file slice lines are rounded to 0.1!
NHdr.pixdim=pixdim;

NHdr.vox_offset=struct('val',352,'prec','float');%default value; r2agui does not plan to write out other than nii files, no magic
if options.rescale==1
    rs=pars.img(1).vis.rescale_slope;
    ri=pars.img(1).vis.rescale_intercept;
else
    rs=1;
    ri=0;
end
NHdr.scl_slope=struct('val',rs,'prec','float'); %scaling slope
NHdr.scl_inter=struct('val',ri,'prec','float');%scaling intercept
NHdr.slice_end=struct('val',0,'prec','int16');
NHdr.slice_code=struct('val',0,'prec','char'); %e.g. ascending
NHdr.xyzt_units=struct('val',10,'prec','char'); %e.g. mm and sec (=10)
NHdr.cal_maxmin=struct('val',[0 0],'prec','float');
NHdr.slice_duration=struct('val',0,'prec','float'); %time for one slice
NHdr.toffset=struct('val',0,'prec','float'); %time axis to shift
NHdr.glmaxmin=struct('val',[255 0],'prec','int32');
NHdr.descrip=struct('val',Descrip,'prec','char');
NHdr.aux_file=struct('val',char(ones(1,24)*32),'prec','char');
NHdr.qform_code=struct('val',kNIFTI_XFORM_SCANNER_ANAT,'prec','int16');
NHdr.sform_code=struct('val',kNIFTI_XFORM_SCANNER_ANAT,'prec','int16');
NHdr.quatern_bcd=struct('val',quatern_bcd,'prec','float');
NHdr.qoffset_xyz=struct('val',qoffset_xyz,'prec','float');
NHdr.srow_xyz=struct('val',M(1:3,:)','prec','float'); % 4D angulation matrix, ordered row-wise (without bottom row)
NHdr.intent_name=struct('val',char(ones(1,16)*32),'prec','char');
NHdr.magic=struct('val',kNIFTI_MAGIC_EMBEDDED_HDR,'prec','int32');
end


function WriteNii(fname,NHdr,Data3D)

fid=fopen(fname,'w');
%write header to nii binary
fn=fieldnames(NHdr);
for t=1:length(fn)
    %disp(fn{t});
    hdrfield=NHdr.(fn{t});
    fwrite(fid,hdrfield.val,hdrfield.prec);
end
switch NHdr.datatype.val
    case 2
        bitpixstr='int8';
    case 4
        bitpixstr='int16';
    case 16
        bitpixstr='float32';
    case 64
        bitpixstr='float64';
end
%now add 4 extra bytes in space between header and offset for data indicating
%that single .nii file ("n+1\0") rather than separate img/hdr files were written
%see http://nifti.nimh.nih.gov
niftitag=[hex2dec('6e') hex2dec('2b') hex2dec('31') hex2dec('00')];
fwrite(fid,niftitag,'char');

% add remaining 0s (probably not required)
if NHdr.vox_offset.val-NHdr.HdrSz.val-4>0
    fwrite(fid,zeros(NHdr.vox_offset.val-NHdr.HdrSz.val-4,1),'char');
end
%Data3D_rad=Data3D; %preallocate for speed


Data3D = flipdim(Data3D,2);

count=fwrite(fid,Data3D,bitpixstr); %actually write data
if count<length(Data3D(:))
    disp('Write failed: less bytes written than there voxels in image. Disk full?');
end
fclose(fid);
end


function [HdrMat,realVoxSize,matrixSize]=calc_angulation(pars,options)

imgOrientations = [pars.img.orient];
imgInfo = [pars.img.info];
imgSpecial = [pars.img.special];
if options.angulation && isstruct(imgOrientations)
    % trying to incorporate AP FH RL rotation angles: determined using some 
    % common sense, Chris Rordon's help + source code and trial and error, 
    % this is considered EXPERIMENTAL!
    angulations =  [imgOrientations.img_angulation]';
    if size(angulations,1)==1 || all(all(diff(angulations,1,1)==0,1),2)
        angAP = angulations(1,1) * (pi/180);
        angFH = angulations(1,2) * (pi/180);
        angRL = angulations(1,3) * (pi/180);
        r1 = [    1           0          0     ;     0      cos(angRL) -sin(angRL);      0      sin(angRL) cos(angRL)];
        r2 = [cos(angAP)      0      sin(angAP);     0          1           0     ; -sin(angAP)     0      cos(angAP)];
        r3 = [cos(angFH) -sin(angFH)     0     ; sin(angFH) cos(angFH)      0     ;      0          0          1     ];
        R_tot=[[r1*r2*r3,[0;0;0]];0 0 0 1];
    else
        warning('Angulations vary across slices. No angulation data converted');
        R_tot = eye(4);
    end
else
    R_tot=eye(4);
end
if numel(imgOrientations)>1 && ~all(diff([imgOrientations.slice_orientation]')==0)
    warning('Slice orientation varies across slices. Using orientation of first slice. This is likely WRONG!');
end
sliceOrientation = imgOrientations(1).slice_orientation;
if numel(imgInfo)>1 && ~all(diff([imgInfo.slicethk]')==0)
    warning('Slice thickness varies across slices. Using thickness of first slice.');
end
sliceThk = imgInfo(1).slicethk;
if numel(imgInfo)>1 && ~all(diff([imgInfo.slicegap]')==0)
    warning('Slice gap varies across slices. Using gap of first slice.');
end
sliceThk = sliceThk + imgInfo(1).slicegap;
if numel(imgSpecial)>1 && ~all(all(diff([imgSpecial.pix_spacing]')==0))
    warning('Pixel spacing varies across slices. Using voxel size of first slice.');
end
pixSpacing = imgSpecial(1).pix_spacing;
if ~all(all(diff([imgInfo.recon_res]')==0))
    warning('Recon resolution varies across slices. Using recon resolution of first slice.');
end
reconRes = imgInfo(1).recon_res;
offcMidslice = pars.orient.off_ctr_midslice;

switch sliceOrientation
    case 1	%transversal
        lmm= eye(4); %do not rotate
    case 2 	%sagittal
        lmm=[0  0 -1  0;
            1  0  0  0;
            0 -1  0  0;
            0  0  0  1];
    case 3 %coronal
        lmm=[1  0  0  0;  %rotate 90 degrees
            0  0  1  0;
            0 -1  0  0;
            0  0  0  1];
end
realVoxSize=[pixSpacing' sliceThk]; % return argument, used to fill in pixdim nifti header info
Zm = diag([realVoxSize 1]);
patient_to_tal   = diag([-1 -1 1 1]);
analyze_to_dicom = diag([1  -1 1 1]);
%A_tot=patient_to_tal*R_tot*Zm*lmm*analyze_to_dicom;
HdrMat=patient_to_tal*R_tot*lmm*Zm*analyze_to_dicom;
%A_tot=patient_to_tal*Zm*R_tot*lmm*analyze_to_dicom;
% get position of the center of the volume
p_orig = [(reconRes(1)-1)/2, (reconRes(2)-2)/2, (pars.max.num_slices-1)/2, 1];
matrixSize = [reconRes;pars.max.num_slices];

offsetA=HdrMat*p_orig';
if options.angulation
    % trying to incorporate AP FH RL translation: determined using some 
    % common sense, Chris Rordon's help + source code and trial and error, 
    % this is considered EXPERIMENTAL!
    HdrMat(1:3,4)=-offsetA(1:3) - [0 0 1;1 0 0;0 -1 0] * offcMidslice';
else
    HdrMat(1:3,4)=-offsetA(1:3);
end
end

function [qoffset_xyz,quatern_bcd,qfac]=nifti_mat44_to_quatern(A)
%procedure nifti_mat44_to_quatern( lR :TMatrix;
%                             var qb, qc, qd,
%                             qx, qy, qz,
%                             dx, dy, dz, qfac : single);


%var
%   r11,r12,r13 , r21,r22,r23 , r31,r32,r33, xd,yd,zd , a,b,c,d : double;
%   P,Q: TMatrix;  //3x3
%begin


% offset outputs are read write out of input matrix
qoffset_xyz = A(1:3,4);

% load 3x3 matrix into local variables
%FromMatrix(lR,r11,r12,r13,r21,r22,r23,r31,r32,r33);

%(* compute lengths of each column; these determine grid spacings  *)

d1=sqrt(sum( A(1:3,1:3).*A(1:3,1:3) ));
%xd := sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
%yd := sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
%zd := sqrt( r13*r13 + r23*r23 + r33*r33 ) ;

%(* if a column length is zero, patch the trouble *)
if (d1(1)==0 )
    A(:,1) = [1 0 0]'; d1(1) = 1;
end
if (d1(2)==0 )
    A(:,2) = [0 1 0]'; d1(2) = 1;
end
if (d1(3)==0 )
    A(:,3) = [0 0 1]'; d1(3) = 1;
end;


%(* normalize the columns *)

A(1:3,1)=A(1:3,1)/d1(1);
A(1:3,2)=A(1:3,2)/d1(2);
A(1:3,3)=A(1:3,3)/d1(3);

%(* At this point, the matrix has normal columns, but we have to allow
%   for the fact that the hideous user may not have given us a matrix
%  with orthogonal columns.

%   So, now find the orthogonal matrix closest to the current matrix.

%   One reason for using the polar decomposition to get this
%   orthogonal matrix, rather than just directly orthogonalizing
%   the columns, is so that inputting the inverse matrix to R
%   will result in the inverse orthogonal matrix at this point.
%   If we just orthogonalized the columns, this wouldn't necessarily hold. *)
Q =  A(1:3,1:3);

%//x  P = nifti_mat33_polar(Q) ;  (* P is orthog matrix closest to Q *)
%FromMatrix(P,r11,r12,r13,r21,r22,r23,r31,r32,r33);
[U,~,V]=svd(Q);
P=U*V;

%(*                            [ r11 r12 r13 ]               *)
%(* at this point, the matrix  [ r21 r22 r23 ] is orthogonal *)
%(*                            [ r31 r32 r33 ]               *)

%(* compute the determinant to determine if it is proper *)

zd = det(P); % should be -1 or 1 *)

if( zd > 0 )		%* proper *)
    qfac  = 1.0 ;
else 		%* improper ==> flip 3rd column *)
    qfac = -1.0 ;
    P(:,3)=-P(:,3);
end

%(* now, compute quaternion parameters *)

a = trace(P) + 1;

if( a > 0.5 )                 %(* simplest case *)
    a = 0.5 * sqrt(a) ;
    b = 0.25 * (P(3,2)-P(2,3)) / a ;
    c = 0.25 * (P(1,3)-P(3,1)) / a ;
    d = 0.25 * (P(2,1)-P(1,2)) / a ;
else                        %(* trickier case *)
    xd = 1.0 + P(1,1) - (P(2,2)+P(3,3)) ;  %(* 4*b*b *)
    yd = 1.0 + P(2,2) - (P(1,1)+P(3,3)) ;  %(* 4*c*c *)
    zd = 1.0 + P(3,3) - (P(1,1)+P(2,2)) ;  %(* 4*d*d *)
    if( xd > 1.0 )
        b = 0.5 * sqrt(xd) ;
        c = 0.25* (P(1,2)+P(2,1)) / b ;
        d = 0.25* (P(1,3)+P(3,1)) / b ;
        a = 0.25* (P(3,2)-P(2,3)) / b ;
        else if( yd > 1.0 )
            c = 0.5 * sqrt(yd) ;
            b = 0.25* (P(1,2)+P(2,1)) / c ;
            d = 0.25* (P(2,3)+P(3,2)) / c ;
            a = 0.25* (P(1,3)-P(3,1)) / c ;
        else
            d = 0.5 * sqrt(zd) ;
            b = 0.25* (P(1,3)+P(3,1)) / d ;
            c = 0.25* (P(2,3)+P(3,2)) / d ;
            a = 0.25* (P(2,1)-P(1,2)) / d ;
            end
    end
    if( a < 0.0 )
         b=-b ; c=-c ; d=-d;
    end
end
quatern_bcd = [b c d];
end


function [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP)
% writes a header
% (function copied from spm99, so spm99 does not have to be present)
% FORMAT [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
%
% P       - filename 	     (e.g 'spm' or 'spm.img')
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [sec])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - [i j k] of origin  (default = [0 0 0])
% DESCRIP - description string (default = 'spm compatible')
%
% s       - number of elements successfully written (should be 348)
%__________________________________________________________________________
%
% spm_hwrite writes variables from working memory into a SPM/ANALYZE
% compatible header file.  The 'originator' field of the ANALYZE format has
% been changed to ORIGIN in the SPM version of the header. funused1
% of the ANALYZE format is used for SCALE
%
% see also dbh.h (ANALYZE) spm_hread.m and spm_type.m
%
%__________________________________________________________________________
% @(#)spm_hwrite.m	2.2 99/10/29


% ensure correct suffix {.hdr} and open header file
%---------------------------------------------------------------------------
P = P(P ~= ' ');
q = length(P);
if q>=4 && (P(q - 3) == '.')
    P = P(1:(q - 4)); 
end;
P = [P '.hdr'];

% For byte swapped data-types, also swap the bytes around in the headers.
mach = 'native';
if spm_type(TYPE,'swapped'),
    if spm_platform('bigend'),
        mach = 'ieee-le';
    else
        mach = 'ieee-be';
    end;
    TYPE = spm_type(spm_type(TYPE));
end;
fid             = fopen(P,'w',mach);

if (fid == -1),
    error(['Error opening ' P '. Check that you have write permission.']);
end;
%---------------------------------------------------------------------------
data_type 	= ['dsr      ' 0];

P     		= [P '                  '];
db_name		= [P(1:17) 0];

% set header variables
%---------------------------------------------------------------------------
DIM		= DIM(:)'; if size(DIM,2) < 4; DIM = [DIM 1]; end
VOX		= VOX(:)'; if size(VOX,2) < 4; VOX = [VOX 0]; end
dim		= [4 DIM(1:4) 0 0 0];
pixdim		= [0 VOX(1:4) 0 0 0];
vox_offset      = OFFSET;
funused1	= SCALE;
glmax		= 1;
glmin		= 0;
bitpix 		= 0;
descrip         = zeros(1,80);
aux_file        = ['none                   ' 0];
origin          = [0 0 0 0 0];

%---------------------------------------------------------------------------
if TYPE == 1;   bitpix = 1;  glmax = 1;        glmin = 0;	end
if TYPE == 2;   bitpix = 8;  glmax = 255;      glmin = 0;	end
if TYPE == 4;   bitpix = 16; glmax = 32767;    glmin = 0;  	end
if TYPE == 8;   bitpix = 32; glmax = (2^31-1); glmin = 0;	end
if TYPE == 16;  bitpix = 32; glmax = 1;        glmin = 0;	end
if TYPE == 64;  bitpix = 64; glmax = 1;        glmin = 0;	end

%---------------------------------------------------------------------------
if nargin >= 7; origin = [ORIGIN(:)' 0 0];  end
if nargin <  8; DESCRIP = 'spm compatible'; end

d          	= 1:min([length(DESCRIP) 79]);
descrip(d) 	= DESCRIP(d);

fseek(fid,0,'bof');

% write (struct) header_key
%---------------------------------------------------------------------------
fwrite(fid,348,		'int32');
fwrite(fid,data_type,	'char' );
fwrite(fid,db_name,	'char' );
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int16');
fwrite(fid,'r',		'char' );
fwrite(fid,'0',		'char' );

% write (struct) image_dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');

fwrite(fid,dim,		'int16');
fwrite(fid,'mm',	'char' );
fwrite(fid,0,		'char' );
fwrite(fid,0,		'char' );

fwrite(fid,zeros(1,8),	'char' );
fwrite(fid,0,		'int16');
fwrite(fid,TYPE,	'int16');
fwrite(fid,bitpix,	'int16');
fwrite(fid,0,		'int16');
fwrite(fid,pixdim,	'float');
fwrite(fid,vox_offset,	'float');
fwrite(fid,funused1,	'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int32');
fwrite(fid,glmax,	'int32');
fwrite(fid,glmin,	'int32');

% write (struct) image_dimension
%---------------------------------------------------------------------------
fwrite(fid,descrip,	'char');
fwrite(fid,aux_file,    'char');
fwrite(fid,0,           'char');
fwrite(fid,origin,      'int16');
if fwrite(fid,zeros(1,85), 'char')~=85
    fclose(fid);
    spm_unlink(P);
    error(['Error writing ' P '. Check your disk space.']);
end

s   = ftell(fid);
fclose(fid);
end



function varargout=spm_platform(varargin)
% Platform specific configuration parameters for SPM
%
% FORMAT ans = spm_platform(arg)
% arg  - optional string argument, can be
%        - 'bigend'  - return whether this architecture is bigendian
%                      - Inf - is not IEEE floating point
%                      - 0   - is little end
%                      - 1   - big end
%        - 'filesys' - type of filesystem
%                      - 'unx' - UNIX
%                      - 'win' - DOS
%                      - 'mac' - Macintosh
%                      - 'vms' - VMS
%        - 'sepchar' - returns directory separator
%        - 'rootlen' - returns number of chars in root directory name
%        - 'user'    - returns username
%        - 'tempdir' - returns name of temp directory
%
% FORMAT PlatFontNames = spm_platform('fonts')
% Returns structure with fields named after the generic (UNIX) fonts, the
% field containing the name of the platform specific font.
%
% FORMAT PlatFontName = spm_platform('font',GenFontName)
% Maps generic (UNIX) FontNames to platform specific FontNames
%
% FORMAT SPM_PLATFORM = spm_platform('init',comp)
% Initialises platform specific parameters in global SPM_PLATFORM
% (External gateway to init_platform(comp) subfunction)
% comp         - computer to use [defaults to MatLab's `computer`]
% SPM_PLATFORM - copy of global SPM_PLATFORM
%
% FORMAT spm_platform
% Initialises platform specific parameters in global SPM_PLATFORM
% (External gateway to init_platform(computer) subfunction)
%
% FORMAT spm_platform('clear')
% Clears global SPM_PLATFORM containing platform specific parameters
%
%                           ----------------
% SUBFUNCTIONS:
%
% FORMAT init_platform(comp)
% Initialise platform specific parameters in global SPM_PLATFORM
% comp         - computer to use [defaults to MatLab's `computer`]
%
%-----------------------------------------------------------------------
%
% Since calls to spm_platform will be made frequently, most platform
% specific parameters are stored as a structure in the global variable
% SPM_PLATFORM. Subsequent calls use the information from this global
% variable, if it exists.
%
% Platform specific difinitions are contained in the data structures at
% the beginning of the init_platform subfunction at the end of this
% file.
%_______________________________________________________________________
% @(#)spm_platform.m	2.10 Matthew Brett 00/11/08


%-Initialise
%-----------------------------------------------------------------------
global SPM_PLATFORM
if isempty(SPM_PLATFORM), init_platform, end

if nargin==0, return, end


switch lower(varargin{1}), case 'init'                  %-(re)initialise
    %=======================================================================
    init_platform(varargin{2:end})
    varargout = {SPM_PLATFORM};

    case 'clear'                                       %-Clear SPM_PLATFORM
        %=======================================================================
        clear global SPM_PLATFORM

    case 'bigend'                      %-Return endian for this architecture
        %=======================================================================
        varargout = {SPM_PLATFORM.bigend};
        if ~isfinite(SPM_PLATFORM.bigend),
            if isnan(SPM_PLATFORM.bigend)
                error(['I don''t know if "',computer,'" is big-endian.'])
            else
                error(['I don''t think that "',computer,...
                    '" uses IEEE floating point ops.'])
            end
        end

    case 'filesys'                                      %-Return file system
        %=======================================================================
        varargout = {SPM_PLATFORM.filesys};

    case 'sepchar'                         %-Return file separator character
        %=======================================================================
        warning('use filesep instead (supported by MathWorks)')
        varargout = {SPM_PLATFORM.sepchar};

    case 'rootlen'           %-Return length in chars of root directory name
        %=======================================================================
        varargout = {SPM_PLATFORM.rootlen};

    case 'user'                                         %-Return user string
        %=======================================================================
        varargout = {SPM_PLATFORM.user};

    case 'tempdir'                              %-Return temporary directory
        %=======================================================================
        twd = getenv('SPMTMP');
        if isempty(twd)
            twd = tempdir;
        end
        varargout = {twd};


    case {'font','fonts'}    %-Map default font names to platform font names
        %=======================================================================
        if nargin<2, varargout={SPM_PLATFORM.font}; return, end
        switch lower(varargin{2})
            case 'times'
                varargout = {SPM_PLATFORM.font.times};
            case 'courier'
                varargout = {SPM_PLATFORM.font.courier};
            case 'helvetica'
                varargout = {SPM_PLATFORM.font.helvetica};
            case 'symbol'
                varargout = {SPM_PLATFORM.font.symbol};
            otherwise
                warning(['Unknown font ',varargin{2},', using default'])
                varargout = {SPM_PLATFORM.font.helvetica};
        end

    otherwise                                        %-Unknown Action string
        %=======================================================================
        error('Unknown Action string')

        %=======================================================================
end
end


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function init_platform(comp)             %-Initialise platform variables
%=======================================================================
if nargin<1, comp=computer; end
global SPM_PLATFORM

%-Platform definitions
%-----------------------------------------------------------------------
PDefs = {	'PCWIN',	'win',	0;...
    'MAC2',		'mac',	1;...
    'SUN4',		'unx',	1;...
    'SOL2',		'unx',	1;...
    'HP700',	'unx',	1;...
    'SGI',		'unx',	1;...
    'SGI64',	'unx',	1;...
    'IBM_RS',	'unx',	1;...
    'ALPHA',	'unx',	0;...
    'AXP_VMSG',	'vms',	Inf;...
    'AXP_VMSIEEE',	'vms',	0;...
    'LNX86',	'unx',	0;...
    'GLNX86',	'unx',  0;...
    'VAX_VMSG',	'vms',	Inf;...
    'VAX_VMSD',	'vms',	Inf	};

PDefs = cell2struct(PDefs,{'computer','filesys','endian'},2);


%-Which computer?
%-----------------------------------------------------------------------
ci = find(strcmp({PDefs.computer},comp));
if isempty(ci), error([comp,' not supported architecture for SPM']), end


%-Set bigend
%-----------------------------------------------------------------------
SPM_PLATFORM.bigend = PDefs(ci).endian;
% Commented out as ISIEEE is obsolete and will be removed in future
% versions of MATLAB:
%if ~isieee, SPM_PLATFORM.bigend = Inf; end	%-Last check for IEEE math


%-Set filesys
%-----------------------------------------------------------------------
SPM_PLATFORM.filesys = PDefs(ci).filesys;


%-Set filesystem dependent stuff
%-----------------------------------------------------------------------
%-File separators character
%-Length of root directory strings
%-User name finding
%-(mouse button labels?)
switch (SPM_PLATFORM.filesys)
    case 'unx'
        SPM_PLATFORM.sepchar = '/';
        SPM_PLATFORM.rootlen = 1;
        SPM_PLATFORM.user    = getenv('USER');
    case 'win'
        SPM_PLATFORM.sepchar = '\';
        SPM_PLATFORM.rootlen = 3;
        SPM_PLATFORM.user    = getenv('USERNAME');
        if isempty(SPM_PLATFORM.user)
            SPM_PLATFORM.user = spm_win32utils('username'); end
    case 'mac'
        SPM_PLATFORM.sepchar = ':';
        SPM_PLATFORM.rootlen = 1;			%-** Not sure!?
        SPM_PLATFORM.user    = '';			%-** Dunno!
    otherwise
        error(['Don''t know filesystem ',SPM_PLATFORM.filesys])
end

%-Fonts
%-----------------------------------------------------------------------
switch comp
    case {'SOL2'}	%-Some Sol2 platforms give segmentation violations with Helvetica
        SPM_PLATFORM.font.helvetica = 'Lucida';
        SPM_PLATFORM.font.times     = 'Times';
        SPM_PLATFORM.font.courier   = 'Courier';
        SPM_PLATFORM.font.symbol    = 'Symbol';
    case {'SUN4','HP700','SGI','SGI64','IBM_RS','ALPHA','LNX86','GLNX86'}
        SPM_PLATFORM.font.helvetica = 'Helvetica';
        SPM_PLATFORM.font.times     = 'Times';
        SPM_PLATFORM.font.courier   = 'Courier';
        SPM_PLATFORM.font.symbol    = 'Symbol';
    case {'PCWIN'}
        SPM_PLATFORM.font.helvetica = 'Arial Narrow';
        SPM_PLATFORM.font.times     = 'Times New Roman';
        SPM_PLATFORM.font.courier   = 'Courier New';
        SPM_PLATFORM.font.symbol    = 'Symbol';
end
end

function T = spm_type(x, arg)
% translates data type specifiers between SPM & Matlab representations
% FORMAT T = spm_type(x, arg)
% x    - specifier
% T    - type
% arg  - optional string argument, can be
%	 - 'swapped' - if type is byteswapped return 1.
%	 - 'maxval'  - return maximum allowed value.
%	 - 'minval'  - return minimum allowed value.
%	 - 'nanrep'  - return 1 if there is a NaN representation.
%	 - 'bits'    - return the number of bits per voxel.
%	 - 'intt'    - return 1 if values rounded to nearest integer.
%_______________________________________________________________________
%
% Original format specifiers are based on ANALYZE.  If the input is
% a number then the corresponding matlab string is returned by default.
% If the input is a string then the appropriate TYPE is returned.
% However, if the optional arg argument is supplied then other
% information will be returned instead.
%
% With no arguments, a list of data types is returned.
%
% Additional support was added for signed bytes, unsigned short and
% unsigned int (by adding 128 to the format specifiers for unsigned bytes
% signed short and signed int).  Byte swapped datatypes have the same
% identifiers as the non-byte-swapped versions, multiplied by a factor of
% 256.
%_______________________________________________________________________
% @(#)spm_type.m	2.3 John Ashburner, Andrew Holmes 99/04/27


prec = char('uint8','int16','int32','float','double','int8','uint16','uint32','uint8','int16','int32','float','double','int8','uint16','uint32');
types   = [    2      4      8   16   64   130    132    136,   512   1024   2048 4096 16384 33280  33792  34816];
swapped = [    0      0      0    0    0     0      0      0,     1      1      1    1     1     1      1      1];
maxval  = [2^8-1 2^15-1 2^31-1  Inf  Inf 2^7-1 2^16-1 2^32-1, 2^8-1 2^15-1 2^31-1  Inf   Inf 2^8-1 2^16-1 2^32-1];
minval  = [    0  -2^15  -2^31 -Inf -Inf  -2^7      0      0,     0  -2^15  -2^31 -Inf  -Inf  -2^7      0      0];
nanrep  = [    0      0      0    1    1     0      0      0,     0      0      0    1     1     0      0      0];
bits    = [    8     16     32   32   64     8     16     32,     8     16     32   32    64     8     16     32];
intt    = [    1      1      1    0    0     1      1      1,     1      1      1    0     0     1      1      1];

if nargin==0,
    T=types;
    return
end

if ischar(x),
    sel = [];
    msk = find(swapped==0);
    for i=msk,
        if strcmp(deblank(prec(i,:)),deblank(x)),
            sel = i;
            break
        end
    end
else
    sel = find(types == x);
end
if nargin == 1,
    if ischar(x),
        if isempty(sel), T = NaN;
        else T = types(sel); end;
    else
        if isempty(sel), T = 'unknown';
        else T = deblank(prec(sel,:)); end;
    end
elseif isempty(sel),
    T = NaN;
else
    switch lower(arg)
        case 'swapped', T = swapped(sel);
        case 'maxval',  T = maxval(sel);
        case 'minval',  T = minval(sel);
        case 'nanrep',  T = nanrep(sel);
        case 'bits',    T = bits(sel);
        case 'intt',    T = intt(sel);
        otherwise,      T = NaN;
    end;
end
end


% Return explict image type depending on image type number
function str=imgTstr(i) 
% From rw_mr_datamodel__define.pro:
% ;
% ; Image type values
% ;
% ;   Image type is an ENUM which is translated into an integer.
% ;   Below is a list of the possible values
% ;
% ;   0   "M"
% ;   1   "R"
% ;   2   "I"
% ;   3   "P"
% ;   4   "CR"
% ;   5   "T0"
% ;   6   "T1"
% ;   7   "T2"
% ;   8   "RHO"
% ;   9   "SPECTRO"
% ;  10   "DERIVED"
% ;  11   "ADC"
% ;  12   "RCBV"
% ;  13   "RCBF"
% ;  14   "MTT"
% ;  15   "TTP"
str = {'M','R','I','P','CR','T0','T1','T2','RHO','SPECTRO','DERIVED','ADC','RCBV','RCBF','MTT','TTP'};
if i>=0 && i<numel(str)-1
    str = str{i+1};
else
    warning('Illegal image type number.');
    str='UNKNOWN_TYPE';
end
end

% Return explict image sequence depending on image type number
function str=imgSstr(i) 
% From rw_mr_datamodel__define.pro:
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
str = {'IR','SE','FFE','DERIVED','PCA','UNSPECIFIED','SPECTRO','SI'};
if i>=0 && i<numel(str)-1
    str = str{i+1};
else
    warning('Illegal image sequence number.');
    str='UNKNOWN_SEQ';
end
end