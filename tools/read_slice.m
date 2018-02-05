function Y = read_slice(Vi,Vref, echo_nbr, expt_nbr)
% Fonction used to rotate Vi to the voxel space of Vo
% Vi = image to move
% Vo = reference image (voxel space)
%-Loop over planes reading to Y

			

Y   = zeros(Vref(1).dim(1:3));       % initialize output volume

xy = Vref(1).dim(1:2);
index_3D_vol = echo_nbr*expt_nbr ;

mat_tmp = Vref(1).mat\Vi(index_3D_vol).mat;
Vi_index_3D_vol=Vi(index_3D_vol);

parfor p = 1:Vref(1).dim(3)
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
%     M = inv(B*(Vref(1).mat\Vi(index_3D_vol).mat));
    M = inv(B*mat_tmp);
%     d = spm_slice_vol(Vi(index_3D_vol),M,xy,3);
    d = spm_slice_vol(Vi_index_3D_vol,M,xy,3);
    Y(:,:,p) = reshape(d,xy);
end


%% reorient the matrix (sag, cor, trans)
tolerance = 1;
R = Vref(1).mat(1:3,1:3);


if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
    R_sort = sort(abs(R(:)));
    R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
end
%% weird code to still display non orthogonal images 
% need to check it !!
%  There is too much distortion in the loaded image for any non-orthogonal rotation or shearing
if sum(sum(R==0) ~= [2 2 2]) ~= 3
    R = zeros([3,3]);
    R(1,1) = Vref(1).mat(1,1);
    R(2,2) = Vref(1).mat(2,2);
    R(3,3) = Vref(1).mat(3,3);
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




function orient = get_orient(R)

% % orient = [];
% for i = 1:3
%     switch find(R(i,:)) * sign(sum(R(i,:))) 
%         case 1
%             orient = [orient 5];		% Left to Right
%         case 2
%             orient = [orient 4];		% Posterior to Anterior
%         case 3
%             orient = [orient 3];		% Inferior to Superior
%         case -1
%             orient = [orient 2];		% Right to Left
%         case -2
%             orient = [orient 1];		% Anterior to Posterior
%         case -3
%             orient = [orient 6];		% Superior to Inferior
%     end
% end


orient = zeros(1,3);
parfor i = 1:3
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
