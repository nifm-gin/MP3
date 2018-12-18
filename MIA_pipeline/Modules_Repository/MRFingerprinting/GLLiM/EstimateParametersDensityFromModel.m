function [density] = EstimateParametersDensityFromModel(X, f_estim, y_samples, verb)

narginchk(2, 4);
if nargin == 2, verb = 0; end

[~,dens] = gllim_inverse_dens(X', f_estim, ones(size(X')), y_samples, verb);

% for c = 1:size(y_samples)
%     clear g
%     for k = 1:20
%         x       = y_samples(c,:);
%         g(k,:)  = dens.alpha(k) * normpdf(x, dens.mu(c,k), dens.S(1,1,k));
%     end
%     density(c,:) = sum(g,1);
% end

[X1,X2] = meshgrid(y_samples(1,:), y_samples(2,:));
for k = 1:20
    tmp = mvnpdf([X1(:) X2(:)], dens.mu(:,k)', dens.S(:,:,k));
    F(:,:,k) = reshape(tmp, length(y_samples(2,:)), length(y_samples(1,:)));    
end
density = sum(F, 3);

