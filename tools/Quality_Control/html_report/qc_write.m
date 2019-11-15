function qc_write(qcdir,inputdate,subject,img,overlay,contrast,command, z, cmap, subject_name)
if iscell(qcdir),   qcdir   = qcdir{1}; end
if iscell(overlay), overlay = overlay{1}; end
if isnumeric(command), command = num2str(command); end
if isnumeric(contrast), contrast = num2str(contrast); end
if exist('z','var'), z(z==0) = []; if isempty(z), clear('z'); Nslices=10; else, Nslices=length(z); end; end
if length(img)==1, div=min(Nslices,4);
elseif length(img)==2, div=min(Nslices,2);
else, div=1;
end

if ~exist(qcdir,'dir'), mkdir(qcdir); end
jsonfname = fullfile(qcdir, 'qc_results.json');
if exist(jsonfname,'file')
    qcjson = loadjson(jsonfname);    
    while iscell(qcjson)
        qcjson = cat(1,qcjson{:});
    end
else
    qcjson = [];
end


qcjson(end+1).moddate        = inputdate;
if ~exist(fullfile(qcdir,subject),'dir')
    mkdir(fullfile(qcdir,subject))
end
if ishandle(overlay)
    print(overlay,fullfile(qcdir,subject,'overlay_img.png'),'-dpng','-r0')
    qcjson(end).overlay_img    = strrep(fullfile(subject,'overlay_img.png'),'\','/');
elseif exist(overlay,'file') && (strcmp(overlay(max(1,end-3):end),'.nii') || strcmp(overlay(max(1,end-6):end),'.nii.gz'))
    [overlay_dat, h] = nii_load(overlay,1);
    overlay_dat = overlay_dat{end};
    overlay_dat = nanmean(overlay_dat,4);
    if ~exist('z','var')
        [~,~,z]=find3d(overlay_dat);
        z = round(linspace(min(z),max(z),Nslices));
    end
    U = uniquetol(overlay_dat(round(linspace(1,end,min(end,1000)))),.1/max(1,nanmax(overlay_dat(:)))); % Unique values (with a tolerance) in a vector made of 1000 points on the image. Can be a bit few ...
    U = U(~isnan(U));
    if length(U)<100 % this is a Mask
        ismask = true;
        segpng = makeimagestack(rot90(overlay_dat(:,:,z)),[],[],[floor(Nslices/div)+1 div],0);
    else % this is an image
        ismask = false;
        segpng = makeimagestack(rot90(overlay_dat(:,:,z)),-3,[],[floor(Nslices/div)+1 div],0);
    end
    % manage anisotropic pixels
    newdim = size(segpng).*(h.pixdim([3 2])/min(h.pixdim([3 2])));
    segpng = imresize(segpng,newdim,'nearest');
    
    % concat
    %segpng = repmat(segpng,[1 length(img)]);
    if length(img)==2
        segpng = cat(2,uint8(segpng==4)*4 + uint8(segpng==1) ,segpng);
    end
    segpng = imresize(segpng, 1000/size(segpng,2),'nearest');
    
    if ismask
        if ~exist('cmap','var') || isempty(cmap)
            cmap = jet(length(unique(segpng(:))));
        else
            cmap = circshift(cmap,1,1);
        end
        segpngrgb = ind2rgb(uint8(segpng),cmap);
        imwrite(segpngrgb,fullfile(qcdir,subject,'overlay_img.png'),'Alpha',double(segpng>0))
    else
        imwrite(segpng,fullfile(qcdir,subject,'overlay_img.png'))
    end
    qcjson(end).overlay_img    = strrep(fullfile(subject,'overlay_img.png'),'\','/');
elseif exist(overlay,'file')
    copyfile(overlay,fullfile(qcdir,subject,'overlay_img.png'))
    qcjson(end).overlay_img    = strrep(fullfile(subject,'overlay_img.png'),'\','/');
end
if ishandle(img)
    print(img,fullfile(qcdir,subject,'bkg_img.png'),'-dpng','-r0')
elseif exist(img{1},'file')  && (strcmp(img{1}(max(1,end-3):end),'.nii') || strcmp(img{1}(max(1,end-6):end),'.nii.gz'))
    [imgdat, h] = nii_load(img,1);
    for iim = 1:length(imgdat)
        tmpimg = imgdat{iim};
        tmpimg = nanmean(tmpimg,4);
        if ~exist('z','var'), z = round(linspace(1,size(tmpimg,3),min(6,size(tmpimg,3)))); end
        if iim == 1
            img2D = makeimagestack(rot90(tmpimg(:,:,z)),-3,[],[floor(Nslices/div)+1 div],0);
        else
            img2D = cat(2+double(length(imgdat)==3),makeimagestack(rot90(tmpimg(:,:,z)),-3,[],[floor(Nslices/div)+1 div],0),img2D);
        end
    end
    % resize image to square pixels
    newdim = size(img2D).*(h.pixdim([3 2])/min(h.pixdim([3 2])));
    img2D = imresize(img2D,newdim,'nearest');

    imwrite(imresize(img2D,1000/size(img2D,2)),fullfile(qcdir,subject,'bkg_img.png'))
else
    copyfile(img{1},fullfile(qcdir,subject,'bkg_img.png'))
end
qcjson(end).background_img = strrep(fullfile(subject,'bkg_img.png'),'\','/');
qcjson(end).contrast       = contrast;
qcjson(end).command        = command;
qcjson(end).subject        = subject_name;

% delete duplicated lines
[~,ia] = unique({qcjson.background_img},'last');
savejson([],qcjson(ia),jsonfname);

% Load json and add to html
qc_reload(qcdir)