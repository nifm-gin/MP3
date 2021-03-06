function Data = write_volume(Data,Vo, view_mode)
%% apply the inverse transformation (rotation/translation) in order 
%% to go back the the Vo referencial
%% Coded by BL 10102017
% View_mode --> 'Axial', 'Saggital', or 'Coronal' (default : 'Axial')
if ~exist('view_mode', 'var')
    view_mode = 'Axial';
end
    

tolerance = 1;
R = Vo(1).mat(1:3,1:3);


if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
    R_sort = sort(abs(R(:)));
    R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
end

%% weird code to still display non orthogonal images 
% need to check it !!
% There is too much distortion in the loaded image for any non-orthogonal rotation or shearing
% if sum(sum(R==0) == [2 2 2]) ~= 3
%     R = zeros([3,3]);
%     R(1,find(abs(Vo(1).mat(:,1)) ==max(abs(Vo(1).mat(:,1))))) = Vo(1).mat(1,find(abs(Vo(1).mat(:,1)) == max(abs(Vo(1).mat(:,1)))));
%     R(2,find(abs(Vo(1).mat(:,2)) ==max(abs(Vo(1).mat(:,2))))) = Vo(1).mat(2,find(abs(Vo(1).mat(:,2)) == max(abs(Vo(1).mat(:,2)))));
%     R(3,find(abs(Vo(1).mat(:,3)) ==max(abs(Vo(1).mat(:,3))))) = Vo(1).mat(3,find(abs(Vo(1).mat(:,3)) == max(abs(Vo(1).mat(:,3)))));
% end
if sum(sum(R==0) == [2 2 2]) ~= 3
    Rnew = zeros([3,3]);
    R = Vo(1).mat(1:3,1:3);
    Rnew(abs(R(1:3,1)) == max(abs(R(1:3,1))),1)  =  R(abs(R(1:3,1)) == max(abs(R(1:3,1))),1);%max(abs(R(1:3,1)));
    Rnew(abs(R(1:3,2)) == max(abs(R(1:3,2))),2)  =  R(abs(R(1:3,2)) == max(abs(R(1:3,2))),2);%max(abs(R(1:3,2)));
    Rnew(abs(R(1:3,3)) == max(abs(R(1:3,3))),3)  =  R(abs(R(1:3,3)) == max(abs(R(1:3,3))),3);%max(abs(R(1:3,3)));
    R =   Rnew;
end


switch view_mode
    case 'Axial'
        
    case 'Coronal'
        Data = flip(Data,2);
        Data = flip(Data,1);
        Data = permute(Data, [3 2 1]);

    case 'Saggital'
         Data = flip(Data,2);
         Data = flip(Data,1);
        Data = permute(Data, [2 3 1]);
        
end

inv_R = inv(R);
orient = get_orient(inv_R);

if ~isequal(orient, [1 2 3])
    rot_orient = mod(orient + 2, 3) + 1;
    flip_orient = orient - rot_orient;
    
    V0_size = Vo(1).private.dat.dim;
    %V0_dim = 1:length(V0_size);
    V0_dim = 1:length(size(Data));
    V0_dim(1:3) = rot_orient;
    Data = permute(Data, V0_dim);
    for i = [3 2 1]
        if flip_orient(i)
            Data = flipdim(Data, i);
        end
    end
    
end





function orient = get_orient(R)

orient = [];

for i = 1:3
    switch find(R(i,:)) * sign(sum(R(i,:))) 
        case 1
            orient = [orient 5];		% 5 Left to Right
        case 2
            orient = [orient 4];		% 4 Posterior to Anterior
        case 3
            orient = [orient 3];		% 3 Inferior to Superior
        case -1
            orient = [orient 2];		% 2 Right to Left
        case -2
            orient = [orient 1];		% 1 Anterior to Posterior
        case -3
            orient = [orient 6];		% 6 Superior to Inferior
    end
end