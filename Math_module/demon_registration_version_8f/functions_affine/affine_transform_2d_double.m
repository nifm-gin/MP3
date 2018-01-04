function Iout=affine_transform_2d_double(Iin,M,black)
% Affine transformation function (Rotation, Translation, Resize)
% This function transforms a volume with a 3x3 transformation matrix 
%
% Iout=affine_transform_2d_double(Iin,Minv,black)
%
% inputs,
%   Iin: The greyscale input image
%   Minv: The (inverse) 3x3 transformation matrix
%   black: If true pixels from outside the image are set to zero 
%           if false to the nearest old pixel.
% output,
%   Iout: The transformed image
%
% example,
%   % Read image
%   I=im2double(imread('lenag2.png'))
%   % Make a transformation matrix
%   M=make_transformation_matrix([2 3],[1.0 1.1],2);
%   % Transform the image
%   Iout=affine_transform_2d_double(I,M)
%   % Show the image
%   figure, imshow(Iout);
%
% Function is written by D.Kroon University of Twente (February 2009)
  
% Make all x,y indices
[x,y]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);

% Calculate center of the image
mean=size(Iin)/2;

% Make center of the image coordinates 0,0
xd=x-mean(1); 
yd=y-mean(2);

% Calculate the Transformed coordinates
Tlocalx = mean(1) + M(1,1) * xd + M(1,2) *yd + M(1,3) * 1;
Tlocaly = mean(2) + M(2,1) * xd + M(2,2) *yd + M(2,3) * 1;

% All the neighborh pixels involved in linear interpolation.
xBas0=floor(Tlocalx); 
yBas0=floor(Tlocaly);
xBas1=xBas0+1;           
yBas1=yBas0+1;

% Linear interpolation constants (percentages)
xCom=Tlocalx-xBas0; 
yCom=Tlocaly-yBas0;
perc0=(1-xCom).*(1-yCom);
perc1=(1-xCom).*yCom;
perc2=xCom.*(1-yCom);
perc3=xCom.*yCom;

% limit indexes to boundaries
check_xBas0=(xBas0<0)|(xBas0>(size(Iin,1)-1));
check_yBas0=(yBas0<0)|(yBas0>(size(Iin,2)-1));
xBas0(check_xBas0)=0; 
yBas0(check_yBas0)=0; 
check_xBas1=(xBas1<0)|(xBas1>(size(Iin,1)-1));
check_yBas1=(yBas1<0)|(yBas1>(size(Iin,2)-1));
xBas1(check_xBas1)=0; 
yBas1(check_yBas1)=0; 

% Get all neigborh intensities
intensity_xyz0=Iin(1+xBas0+yBas0*size(Iin,1));
intensity_xyz1=Iin(1+xBas0+yBas1*size(Iin,1)); 
intensity_xyz2=Iin(1+xBas1+yBas0*size(Iin,1));
intensity_xyz3=Iin(1+xBas1+yBas1*size(Iin,1));

% Make pixels before outside Ibuffer black
if(black>0)
    intensity_xyz0(check_xBas0|check_yBas0)=0;
    intensity_xyz1(check_xBas0|check_yBas1)=0;
    intensity_xyz2(check_xBas1|check_yBas0)=0;
    intensity_xyz3(check_xBas1|check_yBas1)=0;
end
Iout=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3;

% From linear array to square matrix
Iout=reshape(Iout,size(Iin));
 
    