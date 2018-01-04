
function cell_density = parametric_map_histology_cell_density_V4(HE_filename)

HE_info = imfinfo(HE_filename);

% figure;imagesc(HE_data)
% final resolution
options.Resize = 'on';
answer = inputdlg({'res_x', 'res_y'}, 'Final resolution', 1, {'128', '128'},options);
if isempty(answer)
    return
end
res_x = str2double(answer{1});
res_y = str2double(answer{2});

% define x and y indices
final_size_x = ceil(HE_info(1).Height/res_x)*res_x;
final_size_y = ceil(HE_info(1).Width/res_y)*res_y;

ind_x = 1:final_size_x/res_x:final_size_x;
ind_y = 1:final_size_y/res_y:final_size_y;

% create empty matrix
overlay = zeros(final_size_x, final_size_y, 3, 'uint8');
cell_density = zeros(res_x, res_y);


for i = 1:res_x
       ['row ' num2str(i), '/', num2str(res_x)]
    for j = 1:res_y
%         tic
%          ['row ' num2str(i), '/', num2str(res_x) ' and column : ' num2str(j) '/', num2str(res_y)]
        [sampleRGB,~] = imread(HE_filename,1,'PixelRegion',{[ind_x(i),ind_x(i)+final_size_x/res_x-1], [ind_y(j),ind_y(j)+final_size_y/res_y-1]});
%         toc
        % Convert RGB intensity to optical density (absorbance)
        sampleRGB_OD = -log((single(sampleRGB)+1)./256);
        %%%% star the segmenting alog
        I_eq = adapthisteq(sampleRGB_OD(:,:,1));
        mask_em = imextendedmax(I_eq, 0.3); %0.21
        mask_em = imfill(mask_em, 'holes');
        % 3.0518e-04 * the voxel resolution correspond to the smallest nuclei size
        mask_em = bwareaopen(mask_em, round(3.0518e-04*res_y*res_x));
        % remove artifical edge when analysing empty picture
        if sum(~mask_em(:)) == 0
            mask_em_perim = false(size(mask_em));
        else
            mask_em_perim = bwperim(mask_em);
        end
        sample_overlay = imoverlay(sampleRGB, mask_em_perim, [.3 1 .3]);
        overlay(ind_x(i):ind_x(i)+size(sampleRGB,1)-1,ind_y(j):ind_y(j)+size(sampleRGB,2)-1, :) = sample_overlay;
        [~,  cell_density(i,j)]= bwlabel(mask_em);
    end
    
end


[pathstr,name,~] = fileparts(HE_filename);
% display and save results
figure; imagesc(cell_density);figure; imagesc(overlay); 
drawnow
save(fullfile(pathstr,  [name, '-cell_density.mat']), 'cell_density')

% save(fullfile(pathstr,  [name, '-Overlay.mat']), 'overlay')
histo_image = importdata(HE_filename);
sampleRGB_resized_resized= imresize(histo_image, [res_x res_y]);
 
imwrite(sampleRGB_resized_resized, fullfile(pathstr,  [name, '-resized.tiff']),'tiff');
imwrite(uint8(cell_density), fullfile(pathstr,  [name, '-cell_density.tiff']),'tiff');


