function clims = imgadjustlim(img,lowlimuc)
%function clims = imgadjustlim(img)
%function clims = imgadjustlim(img,lowlimuc)
%   give low and high clims for displaying images with imagesc
%   with 1% of pixels saturated at high and low values
%
% INPUT:
% img : 2D matrix
% lowlimuc : low limit user choice (overwrite the clim calculated by the
% algorithm). Useful if the algorithm return a negative value as low limit
% and a positive (or null) value is needed for displaying.
%
% OUTPUT:
% clims : [lowlim highlim]

[I,J] = size(img);
img(isinf(img) | isnan(img)) = 0;
trie = sort(img(:));
if nargin == 2
    lowlim = lowlimuc;
else
    lowlim = trie(round(I*J/100));
end
highlim = trie(round(99*I*J/100));

clims = [lowlim highlim]
