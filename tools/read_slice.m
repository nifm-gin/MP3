function Y = read_slice(Vi,Vref, echo_nbr, expt_nbr, view_mode, num_slice, datatype)
% Fonction used to rotate Vi to the voxel space of Vo
% Vi = image to move
% Vref = reference image (voxel space)
%-Loop over planes reading to Y
tic
interpolation = 0;
			



%Y   = zeros(Vref(1).dim(1:3));       % initialize output volume

xy = Vref(1).dim(1:2);
index_3D_vol = echo_nbr*expt_nbr;
% compute the transformation to apply between the Vi and the Vref
mat_tmp = Vref(1).mat\Vi(index_3D_vol).mat;

Vi_index_3D_vol=Vi(index_3D_vol);
t = toc;
p = num_slice;
B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
M = inv(B*mat_tmp);
tbis = toc;
tmmp = spm_slice_vol(Vi_index_3D_vol,M,xy,interpolation)/Vi_index_3D_vol.pinfo(1);


%Y = nan(size(tmmp,1), size(tmmp,2), Vref(1).dim(3));
%Y(:,:,num_slice) = tmmp;
Y = tmmp;
% 
% for p = 1:Vref(1).dim(3)
%     disp('Boucle part1')
%     clem = toc
%     B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
%     disp('Boucle part2')
%     clem = toc
%     M = inv(B*mat_tmp);
%     disp('Boucle part3')
%     clem = toc
% %     d = spm_slice_vol(Vi_index_3D_vol,M,xy,3);
% %     Y(:,:,p) = reshape(d,xy);
%     Y(:,:,p) = spm_slice_vol(Vi_index_3D_vol,M,xy,interpolation)/Vi_index_3D_vol.pinfo(1);
%     disp('Boucle part4')
%     clem = toc
% end
t1 = toc;

%% reorient the matrix (sag, cor, trans)
tolerance = 1;
R = Vref(1).mat(1:3,1:3);

if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
    R_sort = sort(abs(R(:)));
    R(find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
end
t2 = toc;
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

t3 = toc;
inv_R = inv(R);
orient = get_orient(inv_R);
t4 = toc;
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
t5 = toc;
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

t6 = toc;
%info = niftiinfo(Vi(1).fname);
%Y = cast(Y, info.Datatype);
Y = cast(Y, datatype);
t7 = toc;
disp(['Read_slice : ', num2str(t), '  ', num2str(tbis), '  ', num2str(t1),'  ', num2str(t2), '  ', num2str(t3),'  ', num2str(t4), '  ', num2str(t5),'  ', num2str(t6), '  ', num2str(t7)]);


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
