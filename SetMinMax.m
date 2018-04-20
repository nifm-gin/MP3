function [hmin,hmax] = SetMinMax(image)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

AUTO_THRESHOLD = 5000;
nbpix = size(image, 1)*size(image, 2);
limit = nbpix/10;
threshold = nbpix/AUTO_THRESHOLD;
nbbins = 256;
figure('Name', 'Hist');
h = histogram(image(:), nbbins);

i=0;
found = 0;
while ~found && i<nbbins
    counts = h.Values(i+1);
    if counts>limit
        counts=0;
    end
    found = counts>threshold;
    i=i+1;
end
hmin = h.BinEdges(i);

i=nbbins;
while ~found && i>0
    counts = h.Values(i-1);
    if counts>limit
        counts=0;
    end
    found = counts>threshold;
    i=i-1;
end
hmax = h.BinEdges(i);
close('Hist')
end

