% Register the moving image to the fixed image using the tform struct
function moving_reg = transformMovingImage_perso(moving, fixed, tform)

if(length(size(fixed))==2)       
    [M, N] = size(fixed);    
    moving_reg = imtransform(moving, tform,...
        'XData', [1 N], 'YData', [1 M], ...
        'Size', [M N]);
else    
    nDims      = ndims(fixed);
    resampler  = makeresampler('linear','fill');
    fillValue  = 0;
    tdims_a    = 1:nDims;
    tdims_b    = 1:nDims;    
    outSize    = size(fixed);
    moving_reg = tformarray(moving,tform,resampler,...
        tdims_a, tdims_b,...
        outSize, [], fillValue);
end


end