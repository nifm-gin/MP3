function [e,egrad]=affine_registration_error(par,scale,I1,I2,type,mode)
% This function affine_registration_error, uses affine transfomation of the
% 3D input volume and calculates the registration error after transformation.
%
% [e,egrad]=affine_registration_error(parameters,scale,I1,I2,type,Grid,Spacing,MaskI1,MaskI2,Points1,Points2,PStrength,mode);
%
% input,
%   parameters (in 2D) : Rigid vector of length 3 -> [translateX translateY rotate]
%                        or Affine vector of length 7 -> [translateX translateY
%                                           rotate resizeX resizeY shearXY
%                                           shearYX]
%
%   parameters (in 3D) : Rigid vector of length 6 : [translateX translateY translateZ
%                                           rotateX rotateY rotateZ]
%                       or Affine vector of length 15 : [translateX translateY translateZ,
%                             rotateX rotateY rotateZ resizeX resizeY
%                             resizeZ,
%                             shearXY, shearXZ, shearYX, shearYZ, shearZX, shearZY]
%
%   scale: Vector with Scaling of the input parameters with the same lenght
%               as the parameter vector.
%   I1: The 2D/3D image which is rigid or affine transformed
%   I2: The second 2D/3D image which is used to calculate the
%       registration error
%   type: The type of registration error used see registration_error.m
% (optional)
%   mode: If 0: linear interpolation and outside pixels set to nearest pixel
%            1: linear interpolation and outside pixels set to zero
%            2: cubic interpolation and outsite pixels set to nearest pixel
%            3: cubic interpolation and outside pixels set to zero
%             
% outputs,
%   e: registration error between I1 and I2
%   egrad: error gradient of input parameters
% example,
%   see example_3d_affine.m
%
% Function is written by D.Kroon University of Twente (April 2009)

if(~exist('mode','var')), mode=0; end


% Scale the inputs
par=par.*scale;

% Delta
delta=1e-5;

% Special case for simple squared difference (speed optimized code)
if((size(I1,3)>3)&&(strcmp(type,'sd')))
    M=getransformation_matrix(par);
    if(isa(I1,'double'))
        if(nargout>1)
            [e,mgrad]=affine_error_3d_double(double(I1),double(I2),double(M),double(mode));
            Me=[mgrad(1) mgrad(2) mgrad(3) mgrad(4);
                mgrad(5) mgrad(6) mgrad(7) mgrad(8);
                mgrad(9) mgrad(10) mgrad(11) mgrad(12);
                0        0        0         0];
            egrad=zeros(1,length(par));
            for i=1:length(par)
                par2=par;
                par2(i)=par(i)+delta*scale(i);
                Mg=getransformation_matrix(par2);
                diffM=(Mg-M).*Me;
                egrad(i)=sum(diffM(:))/delta;
            end
        else
            e=affine_error_3d_double(double(I1),double(I2),double(M),double(mode));
        end
    else
        if(nargout>1)
            [e,mgrad]=affine_error_3d_single(single(I1),single(I2),single(M),single(mode));
            Me=[mgrad(1) mgrad(2) mgrad(3) mgrad(4);
                mgrad(5) mgrad(6) mgrad(7) mgrad(8);
                mgrad(9) mgrad(10) mgrad(11) mgrad(12);
                0        0        0         0];
            egrad=zeros(1,length(par));
            for i=1:length(par)
                par2=par;
                par2(i)=par(i)+delta*scale(i);
                Mg=getransformation_matrix(par2);
                diffM=(Mg-M).*Me;
                egrad(i)=sum(diffM(:))/delta;
            end
        else
            e=affine_error_3d_single(single(I1),single(I2),single(M),single(mode));
        end
    end
    return;
end

% Normal error calculation between the two images, and error gradient if needed
% by final differences
if(size(I1,3)<4)
    e=affine_registration_error_2d(par,I1,I2,type,mode);
    if(nargout>1)
        egrad=zeros(1,length(par));
        for i=1:length(par)
            par2=par; par2(i)=par(i)+delta*scale(i);
            egrad(i)=(affine_registration_error_2d(par2,I1,I2,type,mode)-e)/delta;
        end
    end
else
    e=affine_registration_error_3d(par,I1,I2,type,mode);
    if(nargout>1)
        egrad=zeros(1,length(par));
        for i=1:length(par)
            par2=par; par2(i)=par(i)+delta*scale(i);
            egrad(i)=(affine_registration_error_3d(par2,I1,I2,type,mode)-e)/delta;
        end
    end
end
 

function e=affine_registration_error_2d(par,I1,I2,type,mode)
M=getransformation_matrix(par);
I3=affine_transform(I1,M,mode);

% registration error calculation.
e = image_difference(I3,I2,type);

function e=affine_registration_error_3d(par,I1,I2,type,mode)
M=getransformation_matrix(par);
I3=affine_transform(I1,M,mode);
% registration error calculation.
e = image_difference(I3,I2,type);


function M=getransformation_matrix(par)
switch(length(par))
    case 6  %3d
        M=make_transformation_matrix(par(1:3),par(4:6));
    case 9  %3d
        M=make_transformation_matrix(par(1:3),par(4:6),par(7:9));
    case 15 %3d
        M=make_transformation_matrix(par(1:3),par(4:6),par(7:9),par(10:15));
    case 3 % 2d
        M=make_transformation_matrix(par(1:2),par(3));
    case 5 % 2d
        M=make_transformation_matrix(par(1:2),par(3),par(4:5));
    case 7 % 2d
        M=make_transformation_matrix(par(1:2),par(3),par(4:5),par(6:7));
end

