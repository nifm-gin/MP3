function parametric_map_histology_vessel(Binary_vessel_staining_filename)

%  aa = imread(Binary_vessel_staining_filename);
% %   figure; imshow(aa)
%  BW = im2bw(aa);
% %   figure; imshow(BW)
% 
% test = bwlabeln(~BW);
% figure; imshow(test)
% L = labelmatrix(test);
% CC = bwconncomp(~BW);
% S = regionprops(CC);
% figure; imshow(S)
% 
%  BW2 = ~(bwareaopen(~BW, 20));
% CC2 = bwconncomp(~BW2);
% S2 = regionprops(CC2);
% S2 = struct2dataset(S2);
%  figure;hist(S2.Area(S2.Area < 200), 1000)
%  
%  figure; imshow(BW2)
%  
%  
 
 
% if jpg --> tiff
if ~isempty(strfind(Binary_vessel_staining_filename, '.jpg'))  
    if exist(strrep(Binary_vessel_staining_filename, '.jpg', '.tiff'), 'file') ~= 2
        imwrite(imread(Binary_vessel_staining_filename), strrep(Binary_vessel_staining_filename, '.jpg', '.tiff'));
    end
    Binary_vessel_staining_filename = strrep(Binary_vessel_staining_filename, '.jpg', '.tiff');
    file_info = imfinfo(Binary_vessel_staining_filename);
end



% figure;imagesc(HE_data)
% final resolution
options.Resize = 'on';
answer = inputdlg({'Acquisition: voxel_size (in µm)', 'Acquisition: thickness (in µm)', 'Final: voxel_size (in µm)', 'bwareaopen size'},...
    'Final resolution', 1, {'0.9084*0.9084', '20', '234*234', '20'},options);
if isempty(answer)
    return
end

acq_voxel_size_x = str2double(answer{1}(1:strfind(answer{1},'*')-1));
acq_voxel_size_y = str2double(answer{1}(strfind(answer{1},'*')+1:end));
t= str2double(answer{2});
final_res_x = str2double(answer{3}(1:strfind(answer{3},'*')-1));
final_res_y = str2double(answer{3}(strfind(answer{3},'*')+1:end));
bwareaopen_size = str2double(answer{4});
res_x = ceil(file_info(1).Height /ceil(final_res_x/acq_voxel_size_x));
res_y = ceil(file_info(1).Height /ceil(final_res_y/acq_voxel_size_y));

% define x and y indices


final_size_x = ceil(file_info(1).Height/res_x)*res_x;
final_size_y = ceil(file_info(1).Width/res_y)*res_y;

ind_x = 1:final_size_x/res_x:final_size_x;
ind_y = 1:final_size_y/res_y:final_size_y;

%% display grid on row image
% aa = imread(Binary_vessel_staining_filename);
% aa(ind_x,:,:) = 0;       %# Change every tenth row to black
% aa(:,ind_x,:) = 0;
% figure; imshow(aa)

%% create empty matrix
% overlay = zeros(final_size_x, final_size_y, 3, 'uint8');
vessel_density = zeros(res_x, res_y);
Vascular_surface = zeros(res_x, res_y);
BVhisto= zeros(res_x, res_y);
Vessel_diameter = zeros(res_x, res_y);
VSIhisto = zeros(res_x, res_y);


vessel_orientation = zeros(res_x, res_y);
%vessel_length = zeros(res_x, res_y);
vessel_max_minradii = zeros(res_x, res_y);

for i = 1:res_x
    disp(['row ' num2str(i), '/', num2str(res_x)]);
    parfor j = 1:res_y
        [sampleRGB,~] = imread(Binary_vessel_staining_filename,1,'PixelRegion',{[ind_x(i),ind_x(i)+final_size_x/res_x-1], [ind_y(j),ind_y(j)+final_size_y/res_y-1]});
        BW = im2bw(sampleRGB);
        %           figure;imshow(BW)
        BW2 = ~(bwareaopen(~BW, bwareaopen_size));
        %          figure;imshow(~BW2)      
        stats = regionprops('table',~BW2,'MajorAxisLength','MinorAxisLength','Area','Orientation');
        if ~isempty(stats)
            
            
            vessel_density(i,j) = size(stats,1);
            Surface_vessel = sum(stats.Area);
            Surface_tot = (final_size_x/res_x-1) * (final_size_y/res_y-1);
            VSurf = Surface_vessel/Surface_tot;
            % mean vessel radius (r in µm)
            r = mean(stats.MinorAxisLength)*((acq_voxel_size_x+acq_voxel_size_y)/2);  %(stats.MinorAxisLength+stats.MajorAxisLength)/2
            %             % mean vessel length (h)
            h = mean(stats.MajorAxisLength)*((acq_voxel_size_x+acq_voxel_size_y)/2);
            Vascular_surface(i,j) = (Surface_vessel/Surface_tot)*100;
            BVhisto(i,j) = (((r*h)/(r*h+(r+h)*t)) * (-2*log(1-VSurf)))*100;
            
            %vessel_length(i,j) = mean(stats.MajorAxisLength)*(acq_voxel_size_x*acq_voxel_size_y)/2;
            Vessel_diameter(i,j) = mean(stats.MinorAxisLength)*(acq_voxel_size_x*acq_voxel_size_y)/2;
            VSIhisto(i,j) =sum(stats.Area.*stats.MinorAxisLength)/sum(stats.Area)*(acq_voxel_size_x*acq_voxel_size_y)/2;
            vessel_orientation(i,j) = sum(stats.Area.*stats.Orientation)/sum(stats.Area);
            
            vessel_max_minradii(i,j) = max(stats.MinorAxisLength);
        end
        
        %         overlay(ind_x(i):ind_x(i)+size(sampleRGB,1)-1,ind_y(j):ind_y(j)+size(sampleRGB,2)-1, :) = sample_overlay;
    end
end
% figure('Name','vessel_density') ;imagesc(fliplr(rot90(vessel_density,-1)))
% figure('Name','Vascular_surface') ;imagesc(fliplr(rot90(Vascular_surface,-1)));
% figure('Name','BVhisto') ;imagesc(fliplr(rot90(BVhisto,-1)));
% 
% figure('Name','VD') ;imagesc(fliplr(rot90(VD,-1)));
% figure('Name','VSIhisto') ;imagesc(fliplr(rot90(VSIhisto,-1)));
% 
% figure('Name','vessel_orientation') ;imagesc(fliplr(rot90(vessel_orientation,-1))); 
% figure('Name','vessel_length') ;imagesc(fliplr(rot90(vessel_length,-1))); 
% figure('Name','vessel_max_minradii') ;imagesc(fliplr(rot90(vessel_max_minradii,-1))); 

%% display Binary data used for quantification
[sampleRGB,~] = imread(Binary_vessel_staining_filename);
% figure; imshow(sampleRGB)
BW = im2bw(sampleRGB);
BW2 = ~(bwareaopen(~BW, bwareaopen_size));
% BW2(ind_x,:,:) = 0;       %# Change every tenth row to black
% BW2(:,ind_y,:) = 0;
% figure;imshow(BW2)
imwrite(BW2, strrep(Binary_vessel_staining_filename, '.tiff', ['_bwareaopen_', num2str(bwareaopen_size), '.jpg']))

tagstruct.ImageLength = res_x;
tagstruct.ImageWidth = res_y;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
% tagstruct.RowsPerStrip    = 16; %%%%
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

%% write vessel_density
t = Tiff( strrep(Binary_vessel_staining_filename, '_Binary', '_Vessel_density'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(vessel_density,-1))));
t.close();

%% write BVhisto
t = Tiff(strrep(Binary_vessel_staining_filename, '_Binary', '_BVhisto'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(BVhisto,-1))));
t.close();


%% write VSIhisto
t = Tiff(strrep(Binary_vessel_staining_filename, '_Binary', '_VSIhisto'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(VSIhisto,-1))));
t.close();

%% write Vessel diameter
t = Tiff(strrep(Binary_vessel_staining_filename, '_Binary', '_Vessel_diameter'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(Vessel_diameter,-1))));
t.close();

%% write Vessel Orientation
t = Tiff( strrep(Binary_vessel_staining_filename, '_Binary', '_Vessel_orientation'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(vessel_orientation,-1))));
t.close();

%% write vessel_max_minradii
t = Tiff(strrep(Binary_vessel_staining_filename, '_Binary', '_Vessel_max_minradii'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(vessel_max_minradii,-1))));
t.close();

%% write Vascular_surface
t = Tiff( strrep(Binary_vessel_staining_filename, '_Binary', '_Vascular_surface'),'w');
t.setTag(tagstruct);
t.write(single(fliplr(rot90(Vascular_surface,-1))));
t.close();

delete(Binary_vessel_staining_filename)

