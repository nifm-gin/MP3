function I=affine_registration_image(par,scale,I1,I2,type)
% This function affine_registration_image, uses affine transfomation of the
% 3D input volume and calculates the registration error after transformation.
%
% I=affine_registration_image(parameters,scale,I1,I2,type);
%
% input,
%   parameters (in 2D) : Rigid vector of length 3 -> [translateX translateY rotate]
%                        or Affine vector of length 7 -> [translateX translateY  
%                                           rotate resizeX resizeY shearXY shearYX]
%
%   parameters (in 3D) : Rigid vector of length 6 : [translateX translateY translateZ
%                                           rotateX rotateY rotateZ]
%                       or Affine vector of length 15 : [translateX translateY translateZ,
%                             rotateX rotateY rotateZ resizeX resizeY resizeZ, 
%                             shearXY, shearXZ, shearYX, shearYZ, shearZX, shearZY]
%   
%   scale: Vector with Scaling of the input parameters with the same lenght
%               as the parameter vector.
%   I1: The 2D/3D image which is affine transformed
%   I2: The second 2D/3D image which is used to calculate the
%       registration error
%   type: The type of registration error used see image_difference.m
%
% outputs,
%   I: An volume image with the registration error between I1 and I2
%
% example,
%   see example_3d_rigid.m
%
% Function is written by D.Kroon University of Twente (July 2008)

if(size(I1,3)<4)
    par=par.*scale;
    if(length(par)==3)
        M=make_transformation_matrix(par(1:2),par(3));
    elseif(length(par)==5)
        M=make_transformation_matrix(par(1:2),par(3),par(4:5));
    else
        M=make_transformation_matrix(par(1:2),par(3),par(4:5),par(6:7));
    end
    I3=affine_transform(I1,M);
    [t,I] = image_difference(I3,I2,type);
else
    par=par.*scale;
    if(length(par)==6)
        M=make_transformation_matrix(par(1:3),par(4:6));
    elseif(length(par)==9)
        M=make_transformation_matrix(par(1:3),par(4:6),par(7:9));
    else
        M=make_transformation_matrix(par(1:3),par(4:6),par(7:9),par(10:15));
    end
    I3=affine_transform(I1,M);
    [t,I] = image_difference(I3,I2,type);
end

