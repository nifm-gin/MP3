function [M_out, C_out]=sortMemC(M, C)

[C_out IDX]=sort(C);
if size(M,4) == 3   % for 3D image
    for k = 1 : length(C)
        M_out(:,:,:,k) = M(:,:,:,IDX(k));
    end
elseif size(M,4) ==1   % for 2D image
    for k = 1 : length(C)
        M_out(:,:,k) = M(:,:,IDX(k));
    end
else
    error('sortMemC: wrong dimension of the membership function');
end


