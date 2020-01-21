function [snr_map] = computeSNRmap(Img, backgroundROI)
% Img is a 3D or 4D matrix X x Y x time x Z
% the third dimension is time

tmp = Img .* backgroundROI;
tmp(tmp == 0) = nan;
snr_map = Img(:,:,1,:) ./ nanstd(tmp(:));