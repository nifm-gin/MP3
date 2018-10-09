function Y = read_volume(Vi,Vo,interpolation, view_mode)
% Fonction used to rotate Vi to the voxel space of Vo
% Vi = image to move
% Vo = reference image (voxel space)
%-Loop over planes reading to Y
% interpolation method for the resampling:
%              0         : Zero-order hold (nearest neighbour)
%              1         : First-order hold (trilinear interpolation)
%              2->127    : Higher order Lagrange (polynomial) interpolation
%                          using different holds (second-order upwards)
%              -127 - -1 : Different orders of sinc interpolation
% View_mode --> 'Axial', 'Saggital', or 'Coronal' (default : 'Axial')
if ~exist('view_mode', 'var')
    view_mode = 'Axial';
end
    
if ~exist('interpolation', 'var')
    interpolation = 3;
end

Y   = zeros(Vo(1).dim(1:3));       % initialize output volume
Vi_size = Vi(1).private.dat.dim;

switch  length(Vi_size)
    case {2,3}  %% 3D data
        for p = 1:Vo(1).dim(3)
            B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
            M = inv(B*(Vo(1).mat\Vi.mat));
            d = spm_slice_vol(Vi,M,Vo(1).dim(1:2),interpolation);
            Y(:,:,p) = reshape(d,Vo(1).dim(1:2));
        end
    case 4  %% 4D data
        for n = 1:Vi_size(4)
            for p = 1:Vo(1).dim(3)
                B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
                M = inv(B*(Vo(1).mat\Vi(n).mat));
                d = spm_slice_vol(Vi(n),M,Vo(1).dim(1:2),interpolation);
                Y(:,:,p,n) = reshape(d,Vo(1).dim(1:2));
            end
        end
    case 5 %% 5D data
        for n = 1:Vi_size(4)*Vi_size(5)
            for p = 1:Vo(1).dim(3)
                B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
                M = inv(B*(Vo(1).mat\Vi(n).mat));
                d = spm_slice_vol(Vi(n),M,Vo(1).dim(1:2),interpolation);
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

%% weird code to still display non orthogonal images
% need to check it !!
% There is too much distortion in the loaded image for any non-orthogonal rotation or shearing
if sum(sum(R==0) == [2 2 2]) ~= 3
    %     R = zeros([3,3]);
    %     R(1,find(abs(Vo(1).mat(:,1)) ==max(abs(Vo(1).mat(:,1))))) = Vo(1).mat(1,find(abs(Vo(1).mat(:,1)) == max(abs(Vo(1).mat(:,1)))));
    %     R(2,find(abs(Vo(1).mat(:,2)) ==max(abs(Vo(1).mat(:,2))))) = Vo(1).mat(2,find(abs(Vo(1).mat(:,2)) == max(abs(Vo(1).mat(:,2)))));
    %     R(3,find(abs(Vo(1).mat(:,3)) ==max(abs(Vo(1).mat(:,3))))) = Vo(1).mat(3,find(abs(Vo(1).mat(:,3)) == max(abs(Vo(1).mat(:,3)))));
    Rnew = zeros([3,3]);
    R = Vo(1).mat(1:3,1:3);
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
