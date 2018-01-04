
function cell_density = parametric_map_histology_cell_density_standalone(HE_filename, add_parameters)

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


[filename, directory]=uigetfile({'*.wfml;','HE Files (*.wfml)'}, 'Please select HE file(s)','MultiSelect','on');
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

% start the quantification algo
for z = 1:numel(histo_to_analyse_pathnames)
    info = imfinfo(fullfile(directory,histo_to_analyse_pathnames{z}));
    name = legend_names_selected{z};

    % define x and y indices 
    final_size_x = ceil(info(1).Height/res_x)*res_x;
    final_size_y = ceil(info(1).Width/res_y)*res_y;

    ind_x = 1:final_size_x/res_x:final_size_x;
    ind_y = 1:final_size_y/res_y:final_size_y;
    
    % create empty matrix
    overlay = zeros(final_size_x, final_size_y, 3, 'uint8');
    cell_density = zeros(res_x, res_y);         
    for i = 1:res_x
        ['File : ', num2str(z), '/', num2str(numel(filename)), ' and ', 'row ' num2str(i), '/', num2str(res_x)]
        for j = 1:res_y
            tic
            [sampleRGB,~] = imread(fullfile(directory,histo_to_analyse_pathnames{z}),1,'PixelRegion',{[ind_x(i),ind_x(i)+final_size_x/res_x-1], [ind_y(j),ind_y(j)+final_size_y/res_y-1]});
            toc
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
    
    % display and save results
    figure; imagesc(cell_density);
    drawnow
    save([directory,  name, '-cell_density.mat'], 'cell_density')

    save([directory,  name, '-Overlay.mat'], 'overlay')
     histo_image = importdata(fullfile(directory,histo_to_analyse_pathnames{z}));
    sampleRGB_resized_resized= imresize(histo_image, [res_x res_y]);

    imwrite(sampleRGB_resized_resized, [directory,  name, '-', answer{1}, 'x', answer{2}, '.tiff'],'tiff');
    clear all
end


