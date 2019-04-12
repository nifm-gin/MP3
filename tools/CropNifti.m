function [OutputVol,OutputMat] = CropNifti(InputVol,InputMat)
%CropNifti Crop the empty slices before and after the nifti slices and adapt
%the transform matrix.
%   Detailed explanation goes here

%The input and output matrix are supposed to be as:
% [a, b, c, d]
% [e, f, g, h]
% [i, j, k, l]
% [0, 0, 0, m]

NiftiSlice = [];
for i=1:size(InputVol, 3)
    A = InputVol(:,:,i);
    if any(A(:))
        NiftiSlice = [NiftiSlice, i];
    end
end
assert(min(NiftiSlice) == NiftiSlice(1));
assert(max(NiftiSlice) == NiftiSlice(end));
FirstSlice = NiftiSlice(1);
FinalSlice = NiftiSlice(end);

OutputVol = InputVol(:,:,FirstSlice:FinalSlice,:,:);
OutputMat = InputMat;
Movement = [0;0;FirstSlice-1;1];
OutputMat(:,4) = InputMat * Movement;
%OutputMat(3,4) = InputMat(3,4) + (FirstSlice-1)*InputMat(3,3);

end

