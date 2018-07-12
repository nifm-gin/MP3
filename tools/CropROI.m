function [OutputVol,OutputMat] = CropROI(InputVol,InputMat)
%CropROI Crop the empty slices before and after the ROI slices and adapt
%the transform matrix.
%   Detailed explanation goes here

%The input and output matrix are supposed to be as:
% [a, b, c, d]
% [e, f, g, h]
% [i, j, k, l]
% [0, 0, 0, m]

ROISlices = [];
for i=1:size(InputVol, 3)
    A = InputVol(:,:,i);
    if any(A(:))
        ROISlices = [ROISlices, i];
    end
end
assert(min(ROISlices) == ROISlices(1));
assert(max(ROISlices) == ROISlices(end));
FirstSlice = ROISlices(1);
FinalSlice = ROISlices(end);

OutputVol = InputVol(:,:,FirstSlice:FinalSlice);
OutputMat = InputMat;
Movement = [0;0;FirstSlice-1;1];
OutputMat(:,4) = InputMat * Movement;
%OutputMat(3,4) = InputMat(3,4) + (FirstSlice-1)*InputMat(3,3);

end

