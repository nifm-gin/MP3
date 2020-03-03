function [snr_map] = computeSNRmap(Img, backgroundROI)
% Img is a 3D or 4D matrix X x Y x time x Z
% the third dimension is time

tmp = Img .* backgroundROI;
tmp(tmp == 0) = nan;
snr_map = max(Img,[],3) ./ nanstd(tmp(:));