function I=imgaussian(I,sigma,siz)
% IMGAUSSIAN filters an 1D, 2D or 3D image with an gaussian filter.
% This function uses IMFILTER, for the filtering but instead of using
% a multidimensional gaussian kernel, it uses the fact that a gaussian
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D, or 3D input image
%   SIGMA: The sigma used for the gaussian
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filterd image
%
% example,
%   I = im2double(rgb2gray(imread('peppers.png')));
%   figure, imshow(imgaussian(I,3));
% 
% Function is written by D.Kroon University of Twente (October 2008)

if(~exist('siz','var')), siz=sigma*6; end

% Make 1D gaussian kernel
x=-ceil(siz/2):ceil(siz/2);
H = exp(-(x.^2/(2*sigma^2))); 
H = H/sum(H(:));

% Filter each dimension with the 1D gaussian kernels
if(ndims(I)==1)
    I=imfilter(I,H,'replicate');
elseif(ndims(I)==2)
    Hx=reshape(H,[length(H) 1]); 
    Hy=reshape(H,[1 length(H)]); 
    I=imfilter(imfilter(I,Hx,'replicate'),Hy,'replicate');
elseif(ndims(I)==3)
    Hx=reshape(H,[length(H) 1 1]); 
    Hy=reshape(H,[1 length(H) 1]); 
    Hz=reshape(H,[1 1 length(H)]);
    I=imfilter(imfilter(imfilter(I,Hx,'replicate'),Hy,'replicate'),Hz,'replicate');
else
    error('imgaussian:input','unsupported input dimension');
end

        

