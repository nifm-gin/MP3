function M=make_transformation_matrix(t,r,s,h)
% This function make_transformation_matrix.m creates an affine 
% 2D or 3D transformation matrix from translation, rotation, resize and shear parameters
%
% M=make_transformation_matrix.m(t,r,s,h)
%
% inputs (2D),
%   t: vector [translateX translateY]
%   r: vector [rotate]   (rotation clockwise about the origin)
%   s: vector [resizeX resizeY]
%   h: vector [ShearXY, ShearYX]
%
% inputs (3D),
%   t: vector [translateX translateY translateZ]
%   r: vector [rotateX rotateY rotateZ]
%   s: vector [resizeX resizeY resizeZ]
%   h: vector [ShearXY, ShearXZ, ShearYX, ShearYZ, ShearZX, ShearZY]
%
% outputs,
%   M: 2D or 3D affine transformation matrix
%
% examples,
%   % 2D
%   M=make_transformation_matrix([2 3],[1.0 1.1],2);
%   % 3D
%   M=make_transformation_matrix([0.5 0 0],[1 1 1.2],[0 0 0])
% 
% Function is written by D.Kroon University of Twente (October 2008)

% Process inputs
if(length(t)==2)
	if(~exist('r','var')||isempty(r)), r=0; end
	if(~exist('s','var')||isempty(s)), s=[1 1]; end
	if(~exist('h','var')||isempty(h)), h=[0 0]; end
else
	if(~exist('r','var')||isempty(r)), r=[0 0 0]; end
	if(~exist('s','var')||isempty(s)), s=[1 1 1]; end
	if(~exist('h','var')||isempty(h)), h=[0 0 0 0 0 0]; end
end

% Calculate affine transformation matrix
if(length(t)==2)
    % Make the transformation matrix
    M=mat_tra_2d(t)*mat_siz_2d(s)*mat_rot_2d(r)*mat_shear_2d(h); 
else
    % Make the transformation matrix    
    M=mat_tra_3d(t)*mat_siz_3d(s)*mat_rot_3d(r)*mat_shear_3d(h); 
end


function M=mat_rot_2d(r)
	M=[ cos(r) sin(r) 0;
	   -sin(r) cos(r) 0;
	   0 0 1];
   
function M=mat_siz_2d(s)
	M=[s(1) 0    0;
	   0    s(2) 0;
	   0    0    1];

function M=mat_shear_2d(h)
	M=[1    h(1) 0;
	   h(2) 1    0;
	   0    0    1];
	      
function M=mat_tra_2d(t)
	M=[1 0 t(1);
	   0 1 t(2);
	   0 0 1];


function M=mat_rot_3d(r)
    r=r*(pi/180);
    Rx=[1 0 0 0;
        0 cos(r(1)) -sin(r(1)) 0;
        0 sin(r(1)) cos(r(1)) 0;
        0 0 0 1];

    Ry=[cos(r(2)) 0 sin(r(2)) 0;
        0 1 0 0;
        -sin(r(2)) 0 cos(r(2)) 0;
        0 0 0 1];

    Rz=[cos(r(3)) -sin(r(3)) 0 0;
        sin(r(3)) cos(r(3)) 0 0;
        0 0 1 0;
        0 0 0 1];
    M=Rx*Ry*Rz;

function M=mat_siz_3d(s)
	M=[s(1) 0    0    0;
	   0    s(2) 0    0;
	   0    0    s(3) 0;
	   0    0    0    1];

function M=mat_shear_3d(h)
	M=[1    h(1) h(2) 0;
	   h(3) 1    h(4) 0;
	   h(5) h(6) 1    0;
	   0 0 0 1];
	   
function M=mat_tra_3d(t)
	M=[1 0 0 t(1);
	   0 1 0 t(2);
	   0 0 1 t(3);
	   0 0 0 1];
