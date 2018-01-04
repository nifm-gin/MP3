function new_hdr = update_nifti_hdr(original_hdr, New_data, new_filename)
%% Following a Brick of analysis, the hdr needs to be updated before 
%% beeing used to save the result the the Brick

new_size = size(New_data);
new_hdr = original_hdr; 
new_hdr.ImageSize = new_size;
new_hdr.PixelDimensions = new_hdr.PixelDimensions(1:length(new_size));
new_hdr.Filename = new_filename;
new_hdr.DisplayIntensityRange = [min(New_data(:)) max(New_data(:))];

% update the raw structure
new_hdr.raw.dim = ones([1,8]);
new_hdr.raw.dim(1)=length(new_size);
new_hdr.raw.dim(2:length(new_size)+1) = new_size;
new_hdr.raw.pixdim =  ones([1,length(new_hdr.raw.pixdim)]);
new_hdr.raw.pixdim = original_hdr.raw.pixdim(1:length(new_size));

