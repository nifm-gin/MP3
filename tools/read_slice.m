function Y = read_slice(Vi,Vref, echo_nbr, expt_nbr, view_mode)
% Fonction used to rotate Vi to the voxel space of Vo
% Vi = image to move
% Vref = reference image (voxel space)
%-Loop over planes reading to Y
interpolation = 0;


%tic

Y   = zeros(Vref(1).dim(1:3));       % initialize output volume

xy = Vref(1).dim(1:2);
index_3D_vol = echo_nbr*expt_nbr;
% compute the transformation to apply between the Vi and the Vref
mat_tmp = Vref(1).mat\Vi(index_3D_vol).mat;

Vi_index_3D_vol=Vi(index_3D_vol);

for p = 1:Vref(1).dim(3)
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
    M = inv(B*mat_tmp);
%     d = spm_slice_vol(Vi_index_3D_vol,M,xy,3);
%     Y(:,:,p) = reshape(d,xy);
    Y(:,:,p) = spm_slice_vol(Vi_index_3D_vol,M,xy,interpolation)/Vi_index_3D_vol.pinfo(1);
end

%t1=toc

% ALTERNATIVE sans spm
% tic
% %N1 = read_volume(Vi, Vi, 0, 'Axial');
% Ii = niftiinfo(Vi.fname);
% Ni = niftiread(Ii);
% Iref = niftiinfo(Vref.fname);
% Mati = Ii.Transform.T;
% Matref = Iref.Transform.T;
% D = Mati/Matref;
% %D = Mat2.'/Mat1.';
% %D(:,4) = [0;0;0;1];
% R = imref3d(size(Ni));
% tform1 = affine3d(D);
% [N1_new,ref1] = imwarp(Ni,tform1, 'OutputView', R, 'interp', 'nearest');
% t2 = toc
% %figure;imshow3D(Y);figure;imshow3D(permute(N1_new, [2,1,3]))
% Y = N1_new;


%% reorient the matrix (sag, cor, trans)
tolerance = 1;
R = Vref(1).mat(1:3,1:3);

if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
    R_sort = sort(abs(R(:)));
    R(find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
end
%% weird code to still display non orthogonal images 
% need to check it !!
% There is too much distortion in the loaded image for any non-orthogonal rotation or shearing
if sum(sum(R==0) == [2 2 2]) ~= 3
    Rnew = zeros([3,3]);
    R = Vref(1).mat(1:3,1:3);
    Rnew(abs(R(1:3,1)) == max(abs(R(1:3,1))),1)  =  R(abs(R(1:3,1)) == max(abs(R(1:3,1))),1);%max(abs(R(1:3,1)));
    Rnew(abs(R(1:3,2)) == max(abs(R(1:3,2))),2)  =  R(abs(R(1:3,2)) == max(abs(R(1:3,2))),2);%max(abs(R(1:3,2)));
    Rnew(abs(R(1:3,3)) == max(abs(R(1:3,3))),3)  =  R(abs(R(1:3,3)) == max(abs(R(1:3,3))),3);%max(abs(R(1:3,3)));
    R =   Rnew;
end


inv_R = inv(R);
orient = get_orient(inv_R);

if ~isequal(orient, [1 2 3])
    rot_orient = mod(orient + 2, 3) + 1;
    flip_orient = orient - rot_orient;
    
    for i = 1:3
        if flip_orient(i)
            Y = flipdim(Y, i);
        end
    end
    %  get index of orient (rotate inversely)
    permutation(1:ndims(Y)) = 1:ndims(Y);
    [~, permutation(1:3)] = sort(rot_orient);
    Y = permute(Y, permutation);
end
switch view_mode
    case 'Axial'
        
    case 'Coronal'
        Y = permute(Y, [3 2 1]);
        Y = flip(Y,1);
        Y = flip(Y,2);

    case 'Saggital'
        Y = permute(Y, [3 1 2]);
        Y = flip(Y,1);
        Y = flip(Y,2);
        
end

% info = niftiinfo(Vi(1).fname);
% Y = cast(Y, info.Datatype);




function orient = get_orient(R)

orient = zeros(1,3);
for i = 1:3
    switch find(R(i,:)) * sign(sum(R(i,:))) 
        case 1
            orient(i) = 5;		% Left to Right
        case 2
            orient(i) = 4;		% Posterior to Anterior
        case 3
            orient(i) = 3;		% Inferior to Superior
        case -1
            orient(i) = 2;		% Right to Left
        case -2
            orient(i) = 1;		% Anterior to Posterior
        case -3
            orient(i) = 6;		% Superior to Inferior
    end
end

return;					
