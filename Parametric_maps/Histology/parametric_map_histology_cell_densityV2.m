
function cell_density = parametric_map_histology_cell_density(HE_filename, add_parameters)

% Import the image
% sampleRGB = imread('untitled.tif');
[filename, directory]=uigetfile({'*.wfml;','HE Files (*.wfml)'}, 'Please select HE file(s)','MultiSelect','on');
if isequal(filename,0)
    return
end
if ischar(filename)
    filename = {filename};
end

RGB_image =imread( 'C:\Users\lemassob\Desktop\Projet_Guerbet\HE\p\000001c6\sw\00063550.jpg');

% figure;imagesc(RGB_image)
% final resolution
res_y = 256;
res_x = 256;
% Standard values from literature
He = [0.550 0.758 0.351]';
Eo = [0.398 0.634 0.600]';
Bg = [0.754 0.077 0.652]';

% Create Deconvolution matrix
M = [He/norm(He) Eo/norm(Eo) Bg/norm(Bg)];
D = inv(M);


% Apply Color Deconvolution
% sampleRGB_OD = reshape(sampleRGB_OD, [x*y, channel]);
% sampleHEB_OD = zeros(x*y, 1);
% 
% pixel_nbr = x*y;
% 
% tic
% for i=1:pixel_nbr
%     RGB = sampleRGB_OD(i,:)';
% %      HEB = M\RGB;
%     HEB = D * RGB; %#ok<MINV>
%     sampleHEB_OD(i,1) = HEB(1);
% %     sampleHEB_OD(i,2) = HEB(2);
% %     sampleHEB_OD(i,3) = HEB(3);
% end

% resize the original file to be a multiple of the final resolution
RGB_image_resized = imresize(RGB_image, [ceil(size(RGB_image,1)/res_x)*res_x ceil(size(RGB_image,2)/res_y)*res_y]);
RGB_image_resized(1:size(RGB_image,1),1:size(RGB_image,2),:) = RGB_image;
overlay = RGB_image_resized;

% sampleRGB_resized = ones(ceil(size(sampleRGB,1)/res_x)*res_x, ceil(size(sampleRGB,2)/res_y)*res_y);
% sampleRGB_resized(1:size(sampleRGB,1),1:size(sampleRGB,2)) = sampleRGB;
% overlay = sampleRGB_resized;

ind_x = 1:size(RGB_image_resized,1)/res_x:size(RGB_image_resized,1);
ind_y = 1:size(RGB_image_resized,2)/res_y:size(RGB_image_resized,2);
%%
%% test function imclearborder 

for i = 1:res_x
    i
    for j = 1:res_y
        sampleRGB = RGB_image_resized(ind_x(i):ind_x(i)+size(RGB_image_resized,1)/res_x-1,ind_y(j):ind_y(j)+size(RGB_image_resized,2)/res_y-1,:);
%         figure;imshow(sampleRGB)
        % Convert RGB intensity to optical density (absorbance)
        sampleRGB_OD = -log((single(sampleRGB)+1)./256);

        
        %%%% star the segmenting alog
        I_eq = adapthisteq(sampleRGB_OD(:,:,1));
%           figure;imshow(I_eq)
        % clear sampleRGB_OD
        mask_em = imextendedmax(I_eq, 0.3); %0.21
%         figure;imshow(mask_em)
        mask_em = imfill(mask_em, 'holes');
        mask_em = bwareaopen(mask_em, round(3.0518e-04*res_y*res_x)); % 3.0518e-04 * the voxel resolution correspond to the smallest nuclei size
        % [~, count] = bwlabel(mask_em);
        mask_em_perim = bwperim(mask_em);
%          figure;imshow(mask_em_perim)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           D = bwdist(~mask_em);
%         figure, imshow(D,[],'InitialMagnification','fit')
%         title('Distance transform of ~bw')
 
    %3. Complement the distance transform, and force pixels that don't
    %   belong to the objects to be at -Inf.
 
%          D = -D;
%         D(~mask_em) = -Inf;
 
    %4. Compute the watershed transform, and display the resulting label
    %   matrix as an RGB image.
 
%         L = watershed(D); 
%         rgb = label2rgb(L,'jet',[.5 .5 .5]);
%         figure, imshow(rgb,'InitialMagnification','fit')
%         title('Watershed transform of D')
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
        
        sample_overlay = imoverlay(RGB_image_resized(ind_x(i):ind_x(i)+size(RGB_image_resized,1)/res_x-1,ind_y(j):ind_y(j)+size(RGB_image_resized,2)/res_y-1), mask_em_perim, [.3 1 .3]);
%         figure;imshow(sample_overlay);
        overlay(ind_x(i):ind_x(i)+size(RGB_image_resized,1)/res_x-1,ind_y(j):ind_y(j)+size(RGB_image_resized,2)/res_y-1, :) = sample_overlay;
%         figure;imshow(sampleRGB_OD);
        
        [~,  cell_density(i,j)]= bwlabel(mask_em);
    end
end

figure; image(cell_density);
figure; imshow(overlay);

imwrite(overlay,'54-Overlay.tiff','tiff');
imwrite(cell_density,'54-cell_density.tiff','tiff');

sampleRGB_resized256 = imresize(RGB_image_resized, [256 256]);
figure; imshow(sampleRGB_resized256);
imwrite(sampleRGB_resized256,'54-HE_256.tiff','tiff');


