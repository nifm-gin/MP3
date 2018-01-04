function [I1_TF,I2_TF]=MutualTransform(I1,I2,SizeKernel,SampleDistance,Bx,By,Bz)
% This function MutualTransform transforms two pictures with different modalities 
% into the grey levels of the opposite modality. 
% The function divides the images in small overlapping regions, and calculates 
% mutual histograms from the regions. A (gaussian) kernel is used so that 
% pixels far away from the region counts less in the histogram than nearby 
% pixels. The mutual histograms are then used to find the maximum
% correlation between the intensity of the pixels in image 1 and image 2.
% This maximum correlation is then used to paint a new picture of the opposite modality.
%
% [I1_TF,I2_TF]=MutualTransform(I1,I2,SizeKernel,SampleDistance,Bx,By,Bz)
%
% inputs,
%  I1 : Input Image 1
%  I2 : Input Image 2, with a different (MRI) modality than image I1.
%  SizeKernel: The size of the gaussian kernel which is used to make a
%           local mutual histogram.
%  SampleDistance: The distance between the spacing of the 
%           kernels  / image regions
%  Bx,By,Bz : Translation images from demon registration in x, y (and z) 
%           which deform the images to get a better modality transformation.
%
% outputs,
%  I1_TF : Modality transformed painting of image 1
%  I2_TF : Modality transformed painting of image 2
%
% Example,
%
% Imoving=im2double(imread('D:\Matlab\demons\version 7\images\modtest3.png'));
% Istatic=im2double(imread('D:\Matlab\demons\version 7\images\modtest2.png'));
% 
% Bx=zeros(size(Imoving))-38;
% By=zeros(size(Imoving));
% 
% [Imoving_TF,Istatic_TF]=MutualTransform(Imoving,Istatic,14,7,Bx,By);
% 
% figure,
% subplot(2,2,1), imshow(Imoving)
% subplot(2,2,2), imshow(Istatic)
% subplot(2,2,3), imshow(Imoving_TF)
% subplot(2,2,4), imshow(Istatic_TF)

SizeKernel=floor(SizeKernel/2)*2+1;
SampleDistance=floor(SampleDistance/2)*2+1;

% Histogram number of bins
HistogramBins=255;

% Range of input image
range=getrangefromclass(I1);

% Limit images to the range specified by the input class
I1(I1>range(2))=range(2); I1(I1<range(1))=range(1);
I2(I2>range(2))=range(2); I2(I2<range(1))=range(1);

if(ndims(I1)==2)
    % Make a gaussian filter kernel
    Hkernel=kernel_gaussian([SizeKernel SizeKernel], SizeKernel/6);

    % Make transformed images
    if(exist('Bx','var'))
        % Estimate the forward transformation field from the backward
        % transformation field
        [Fx,Fy]=backwards2forwards(Bx,By); 

        % Make the transformed I1 image
        I1_moved=movepixels(I1,Bx,By);

        % Make the transformed I2 image
        I2_moved=movepixels(I2,Fx,Fy);
    else
        I1_moved=I1; I2_moved=I2;
    end
    [I1_TF,I2_TF]=mutual_transform_2d_double(double(I1),double(I2),double(I1_moved),double(I2_moved),double(range),double(HistogramBins),double(Hkernel),double(SampleDistance));
else
    % Make a gaussian filter kernel
    Hkernel=kernel_gaussian([SizeKernel SizeKernel SizeKernel], SizeKernel/6);

    % Make transformed images
    if(exist('Bx','var'))
        % Estimate the forward transformation field from the backward
        % transformation field
        [Fx,Fy,Fz]=backwards2forwards(Bx,By,Bz); 

        % Make the transformed I1 image
        I1_moved=movepixels(I1,Bx,By,Bz);

        % Make the transformed I2 image
        I2_moved=movepixels(I2,Fx,Fy,Fz);
    else
        I1_moved=I1; I2_moved=I2;
    end
    if(isa(I1,'double'))
        [I1_TF,I2_TF]=mutual_transform_3d_double(double(I1),double(I2),double(I1_moved),double(I2_moved),double(range),double(HistogramBins),double(Hkernel),double(SampleDistance));
    else
        keyboard;
        [I1_TF,I2_TF]=mutual_transform_3d_single(single(I1),single(I2),single(I1_moved),single(I2_moved),single(range),single(HistogramBins),single(Hkernel),single(SampleDistance));
    end
end


function h=kernel_gaussian(Hsize,sigma)
siz   = (Hsize-1)/2; std = sigma;
% Make gaussian kernel
if(length(Hsize)==2) % Detect if 2D or 3D kernel
    [x,y] = ndgrid(-siz(1):siz(1),-siz(2):siz(2));
    arg   = -(x.*x + y.*y)/(2*std*std);
else
    [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
    arg   = -(x.*x + y.*y + z.*z)/(2*std*std);
end
% No negative values
h = exp(arg);
h(h<eps*max(h(:))) = 0;
% Total sum of kernel must be one
sumh = sum(h(:));
if (sumh ~= 0), h  = h/sumh; end
     






