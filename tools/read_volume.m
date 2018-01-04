function Y = read_volume(Vi,Vo)
% Fonction used to rotate Vi to the voxel space of Vo
% Vi = image to move
% Vo = reference image (voxel space)
%-Loop over planes reading to Y

Y   = zeros(Vo(1).dim(1:3));       % initialize output volume
Vi_size = Vi(1).private.dat.dim;

switch  length(Vi_size)
    case 3  %% 3D data
        for p = 1:Vo(1).dim(3)
            B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
            M = inv(B*(Vo(1).mat\Vi.mat));
            d = spm_slice_vol(Vi,M,Vo(1).dim(1:2),3);
            Y(:,:,p) = reshape(d,Vo(1).dim(1:2));
        end
    case 4  %% 4D data
        for n = 1:Vi_size(4)
            for p = 1:Vo(1).dim(3)
                B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
                M = inv(B*(Vo(1).mat\Vi(n).mat));
                d = spm_slice_vol(Vi(n),M,Vo(1).dim(1:2),3);
                Y(:,:,p,n) = reshape(d,Vo(1).dim(1:2));
            end
        end
    case 5 %% 5D data
        for n = 1:Vi_size(4)*Vi_size(5)
            for p = 1:Vo(1).dim(3)
                B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
                M = inv(B*(Vo(1).mat\Vi(n).mat));
                d = spm_slice_vol(Vi(n),M,Vo(1).dim(1:2),3);
                Y(:,:,p,n) = reshape(d,Vo(1).dim(1:2));
            end
        end
        Y = reshape(Y, [size(Y, 1), size(Y, 2), size(Y, 3), Vi_size(4:5)]);
        
end




%% reorient the matrix (sag, cor, trans)
tolerance = 1;
R = Vo(1).mat(1:3,1:3);


if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
    R_sort = sort(abs(R(:)));
    R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
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

orient = [];

for i = 1:3
    switch find(R(i,:)) * sign(sum(R(i,:))) 
        case 1
            orient = [orient 5];		% Left to Right
        case 2
            orient = [orient 4];		% Posterior to Anterior
        case 3
            orient = [orient 3];		% Inferior to Superior
        case -1
            orient = [orient 2];		% Right to Left
        case -2
            orient = [orient 1];		% Anterior to Posterior
        case -3
            orient = [orient 6];		% Superior to Inferior
    end
end

return;					
