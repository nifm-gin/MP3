
function cell_density = parametric_map_histology_vascular_staining_standalone(HE_filename, add_parameters)

clear all

% figure;imagesc(histo_image)
% final resolution
options.Resize = 'on';
answer = inputdlg({'res_x', 'res_y'}, 'Final resolution', 1, {'128', '128'},options);
if isempty(answer)
    return
end
res_x = str2double(answer{1});
res_y = str2double(answer{2});


[filename, directory]=uigetfile({'*.wfml;','Vascular staining Files (*.wfml)'}, 'Please select Vascular staining file(s)','MultiSelect','on');
if isequal(filename,0)
    return
end
if ischar(filename)
    filename = {filename};
end

file_name={};
legend_names={};
% list all the images to analyse
for x = 1:numel(filename)
    % Import the image information
    fid=fopen(fullfile(directory,filename{x}) ,'r');
    if fid>0
        fclose(fid);
        % Info of parent image information
        info_txt = fileread(fullfile(directory,filename{x}));
        % Info of each children images (ROI acquire usely usinf higher magnification)
        image_info_beg = strfind(info_txt, '<wide-field-image');
        image_info_end = strfind(info_txt, '</wide-field-image>');
        legend_beg = strfind(info_txt, '<wide-field-legend>');
        legend_end = strfind(info_txt, '</wide-field-legend>');
        % /4 because 4 folder do exist and interesting images are located
        % in the forth one ("/sw/"
        for  i=1:numel(image_info_beg)/4
            sw_image_nb = i*4;
            info_image = info_txt(image_info_beg(sw_image_nb)+length('<wide-field-image'):image_info_end(sw_image_nb)-1);
            image_height =  info_image(strfind(info_image, '<image-height>') + length('<image-height>') : strfind(info_image, '</image-height>')-1);
            % Check if this children has been acquired, if so save children
            % info (path and name)
            if str2double(image_height) > 0
                file_name = [file_name info_image(strfind(info_image, '<image-path>') + length('<image-path>') : strfind(info_image, '</image-path>')-1)];
                legend_names = [legend_names info_txt(legend_beg(i) + length('<wide-field-legend>') : legend_end(i)-1)];
            end
        end
    else
        warndlg('Somthing wrong with the data','Warning');
        return
    end
    
end
% the user selects the images he wants to analyze
if numel(file_name) > 1
    [image_selected, ok] = listdlg('PromptString','Select one ROI for to analyse: ',...
        'Name', 'Select a file',...
        'SelectionMode','multiple',...
        'ListSize', [400 300],...
        'ListString',legend_names');
    if ok == 0
        warndlg('Somthing wrong with the data','Warning');
        return
    end
else
    image_selected = 1;
end
legend_names_selected = legend_names(image_selected);
histo_to_analyse_pathnames = file_name(image_selected);

% PC/linux stuff
if ispc
    if numel(strfind(histo_to_analyse_pathnames, '/')) >0
        histo_to_analyse_pathnames = strrep(histo_to_analyse_pathnames, '/', '\');
    end
end
% color vector
He = [0.650 0.704 0.286]';
DAB = [0.268 0.570 0.776]';
Zeromatrix = [0.7110272 0.42318153 0.5615672]';

% Create Deconvolution matrix
M = [He/norm(He) DAB/norm(DAB) Zeromatrix/norm(Zeromatrix)];
D = inv(M);

% start the quantification
for z = 1:numel(histo_to_analyse_pathnames)
    info = imfinfo(fullfile(directory,histo_to_analyse_pathnames{z}));
    name = legend_names_selected{z};
    
    % defin x and y indices to order to
    final_size_x = ceil(info(1).Height/res_x)*res_x;
    final_size_y = ceil(info(1).Width/res_y)*res_y;
    
    ind_x = 1:final_size_x/res_x:final_size_x;
    ind_y = 1:final_size_y/res_y:final_size_y;
    
    % create empty matrix
    overlay = zeros(final_size_x, final_size_y, 3, 'uint8');
    Vascular_Surface = zeros(res_x, res_y);
    Vascular_Orientation = zeros(res_x, res_y);
    Vascular_minRadii = zeros(res_x, res_y);
    Vascular_maxRadii = zeros(res_x, res_y);
    
    
    for i = 1:res_x
        ['File : ', num2str(z), '/', num2str(numel(histo_to_analyse_pathnames)), ' and ', 'row ' num2str(i), '/', num2str(res_x)]
        for j = 1:res_y
              %% load full image in low res
            %             [histo_lowres,~] = imread(fullfile(directory,histo_to_analyse_pathnames{z}),2);
            % figure;imshow(histo_lowres)
            % [histo_Higres,~] = imread(fullfile(directory,histo_to_analyse_pathnames{z}));
            % figure;imshow(histo_Higres)             
            
            [sampleRGB,~] = imread(fullfile(directory,histo_to_analyse_pathnames{z}),1,'PixelRegion',{[ind_x(i),ind_x(i)+final_size_x/res_x-1], [ind_y(j),ind_y(j)+final_size_y/res_y-1]});
            [height, width, channel] = size(sampleRGB);
            %             figure;imshow(sampleRGB)
            % Apply Color Deconvolution
            sampleHEB_OD = zeros(height*width,1);
            sampleRGB_in_vector = reshape(single(sampleRGB), [height*width , 3,1]);
            for ii=1:size(sampleRGB_in_vector,1)
                    HEB = D * sampleRGB_in_vector(ii,:)';
                    sampleHEB_OD(ii) = HEB(2);
            end
             tmp = reshape(sampleHEB_OD, [height, width , 1]);
            % figure;image(tmp);
            
             mask_em= tmp < -40 == 1; 
            %   figure;imshow(mask_em);
            
            
            %%%% star the segmenting alog
%               mask_em = adapthisteq(reshape(sampleHEB_OD, [height, width , 1]));
%             mask_em = adapthisteq(single(sampleHEB_OD(:,:,2)));
            %             figure;imshow(mask_em)
%             mask_em = ~im2bw(mask_em);
            %             figure;imshow(mask_em)
            % 3.0518e-04 * the voxel resolution correspond to the smallest nuclei size
            
            se = strel('disk', 1);
            mask_erode = imerode(mask_em, se);
            %             figure;imshow(mask_erode)
             se = strel('disk', 3);
            mask_erode_dilate = imdilate(mask_erode, se);
              %             figure;imshow(mask_erode_dilate)
              
            mask_erode_dilate_fill = imfill(mask_erode_dilate, 'holes');
            %             figure;imshow(mask_erode_dilate_fill) 
            mask_erode_dilate_fill_despeckle  = bwareaopen(mask_erode_dilate_fill, 100);
            %             figure;imshow(mask_erode_dilate_fill_despeckle)
           

            mask_em_perim = bwperim(mask_erode_dilate_fill_despeckle);
            %             figure;imshow(mask_em_perim)
            sample_overlay = imoverlay(sampleRGB, mask_erode_dilate_fill_despeckle, [.3 1 .3]);
            %             figure;imagesc(sample_overlay)
            
            overlay(ind_x(i):ind_x(i)+size(sampleRGB,1)-1,ind_y(j):ind_y(j)+size(sampleRGB,2)-1, :) = sample_overlay;
%             vessels_info = regionprops(mask_erode_dilate_fill_despeckle, 'Area', 'Orientation', 'MinorAxisLength', 'MajorAxisLength');
%             
%             % output maps
%             Vascular_Surface(i,j)= sum([vessels_info(:).Area]) * 100 / (height * width);
%             Vascular_Orientation(i,j)= mean([vessels_info(:).Orientation]);
%             Vascular_minRadii(i,j)= mean([vessels_info(:).MinorAxisLength]);
%             Vascular_maxRadii(i,j)= mean([vessels_info(:).MajorAxisLength]);
        end
    end
%     figure; imagesc(overlay);
%     % display and save results
%     figure; imagesc(Vascular_Surface);
%     figure; imagesc(Vascular_Orientation);
%     figure; imagesc(Vascular_minRadii);
%     figure; imagesc(Vascular_maxRadii);
%     drawnow
%% save data
%     save([directory,  name, '-Vascular_Surface.mat'], 'Vascular_Surface');
%     save([directory,  name, '-Vascular_Orientation.mat'], 'Vascular_Orientation');
%     save([directory,  name, '-Vascular_minRadii.mat'], 'Vascular_minRadii');
%     save([directory,  name, '-Vascular_maxRadii.mat'], 'Vascular_maxRadii');
%     
    save([directory,  name, '-bw.mat'], 'mask_erode_dilate_fill_despeckle')
    save([directory,  name, '-Overlay.mat'], 'overlay')
    histo_image = importdata(fullfile(directory,histo_to_analyse_pathnames{z}));
    sampleRGB_resized_resized= imresize(histo_image, [res_x res_y]);
    
    imwrite(sampleRGB_resized_resized, [directory,  name, '-', answer{1}, 'x', answer{2}, '.tiff'],'tiff');
    clear sampleRGB_resized_resized
    
end


