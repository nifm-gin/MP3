
function cell_density = parametric_map_histology_cell_density(HE_filename, add_parameters)

% Import the image
% sampleRGB = imread('untitled.tif');
tic
HE_filename = 'C:\Users\lemassob\Desktop\Projet_Guerbet\HE\p\000001c6\00063550-52.svs';
% HE_filename = 'C:\Users\lemassob\Desktop\test_microscop\p\000001c6\sw\sampleRGB.mat';
sampleRGB =importdata(HE_filename);

[~, ~, channel] = size(sampleRGB);

% Convert RGB intensity to optical density (absorbance)
sampleRGB_OD = -log((single(sampleRGB)+1)./256);

% clear sampleRGB

% Construct color deconvolution matrix
% Take the average around region of interest
% H2 = ones(10,10) ./ 100;
% sampleRGB_OD_Blur = zeros(height,width,channel);
% 
% for k=1:channel
%     sampleRGB_OD_Blur(:,:,k) = filter2(H2,sampleRGB_OD(:,:,k),'same');
% end
% 
% % Manual sampling of three regions of interest from the micrograph to calculate 
% % relative optical intensity of each stains in three color channels
% He = reshape(sampleRGB_OD_Blur(44,348,:),3,1);
% Eo = reshape(sampleRGB_OD_Blur(361,290,:),3,1);
% Bg = reshape(sampleRGB_OD_Blur(425,77,:),3,1);
% 
% clear sampleRGB_OD_Blur
% Standard values from literature
He = [0.550 0.758 0.351]';
Eo = [0.398 0.634 0.600]';
Bg = [0.754 0.077 0.652]';

% Create Deconvolution matrix
M = [He/norm(He) Eo/norm(Eo) Bg/norm(Bg)];
D = inv(M);

% % Apply Color Deconvolution
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
% toc
% 
% sampleRGB_OD = reshape(sampleRGB_OD, [x, y, channel]);

%%%% star the segmenting alog
I_eq = adapthisteq(sampleRGB_OD(:,:,1));

% clear sampleRGB_OD
mask_em = imextendedmax(I_eq, 0.21);
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 10);
% [~, count] = bwlabel(mask_em);

mask_em_perim = bwperim(mask_em);
overlay = imoverlay(sampleRGB, mask_em_perim, [.3 1 .3]);
% figure;imshow(overlay);

res_y = 256;
res_x = 256;
mask_em_resized = false([ceil(size(mask_em,1)/res_x)*res_x ceil(size(mask_em,2)/res_y)*res_y]);
mask_em_resized(1:size(mask_em,1),1:size(mask_em,2)) = mask_em;

ind_x = 1:size(mask_em_resized,1)/res_x:size(mask_em_resized,1);
ind_y = 1:size(mask_em_resized,2)/res_y:size(mask_em_resized,2);

tic
for i = 1:res_x
    for j = 1:res_y
%         cell_density(i,j)= sum(sum(bwmorph(mask_em_resized(ind_x(i):ind_x(i)+size(mask_em_resized,1)/res_x-1,ind_y(j):ind_y(j)+size(mask_em_resized,2)/res_y-1), 'shrink', inf),1),2);
         [~,  cell_density(i,j)]= bwlabel(mask_em_resized(ind_x(i):ind_x(i)+size(mask_em_resized,1)/res_x-1,ind_y(j):ind_y(j)+size(mask_em_resized,2)/res_y-1));

    end

end
toc
figure;imagesc(cell_density);



