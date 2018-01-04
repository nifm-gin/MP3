
HE_filename = 'C:\Users\lemassob\Desktop\Projet_Guerbet\HE\p\000001c6\00063547.svs';
% HE_filename = 'C:\Users\lemassob\Desktop\test_microscop\p\000001c6\sw\sampleRGB.mat';
sampleRGB =importdata(HE_filename);
% figure;imshow(sampleRGB);
% Convert RGB intensity to optical density (absorbance)
% sampleRGB_OD = -log(single(sampleRGB+1)./256);
% figure;imshow(sampleRGB);


I = imread('http://blogs.mathworks.com/images/steve/60/nuclei.png');
I_cropped = I(400:900, 465:965);
imshow(I_cropped)
I_eq = adapthisteq(I_cropped);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x, y,channel] = size(sampleRGB_OD);
% 
% 
% % Standard values from literature
% He = [0.550 0.758 0.351]';
% Eo = [0.398 0.634 0.600]';
% Bg = [0.754 0.077 0.652]';
% 
% % Create Deconvolution matrix
% M = [He/norm(He) Eo/norm(Eo) Bg/norm(Bg)];
% D = inv(M);
% %Apply Color Deconvolution
% sampleRGB_OD = reshape(sampleRGB_OD, [x*y, channel]);
% sampleHEB_OD = zeros(x*y, 1);
% 
% pixel_nbr = x*y;
% 
% % Apply Color Deconvolution
% for i=1:pixel_nbr
%     RGB = sampleRGB_OD(i,:)';
% %      HEB = M\RGB;
%     HEB = D * RGB; %#ok<MINV>
%     sampleHEB_OD(i,1) = HEB(1);
% %      sampleHEB_OD(i,2) = HEB(2);
% %      sampleHEB_OD(i,3) = HEB(3);
% end
% I_eq = reshape(sampleHEB_OD(:,1), [x, y]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_eq = rgb2gray(sampleRGB);
I_eq = imcomplement(I_eq);
% figure;imshow(I_eq);

I_eq = adapthisteq(I_eq);
% figure;imshow(I_eq);
%%%% star the segmenting alog
bw =im2bw(I_eq, graythresh(I_eq));
% figure;imshow(bw)

bw2 = imfill(bw,'holes');
% bw3 = imopen(bw2, ones(5,5));
% bw4 = bwareaopen(bw3, 40);
 bw4 = bw;
% figure;imshow(bw3)
bw4_perim = bwperim(bw4);
overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
% imshow(overlay1)

mask_em = imextendedmax(I_eq, 5); 
% figure;imshow(mask_em)

% mask_em = imclose(mask_em, ones(5,5));
% mask_em = imfill(mask_em, 'holes');
  mask_em = bwareaopen(mask_em, 4);

% figure;imshow(mask_em)
overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
% figure;imshow(overlay2)

I_eq_c = imcomplement(I_eq);

% figure;imshow(mask_em)
I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);
% figure;imshow(I_mod)
L = watershed(I_mod);
% L2 = bwareaopen(L, 50);
% figure;imshow(label2rgb(L))

overlay3 = imoverlay(sampleRGB, ~bwlabel(L), [.3 1 .3]);
figure;imshow(overlay3)
