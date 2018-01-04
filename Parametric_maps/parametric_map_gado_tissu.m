function Ct = parametric_map_gado_tissu(xdata,BETA,Cp)
% developpee par Marine Beaumont au cours de sa these

% modele de l'evolution de la concentration de gado dans un voxel
% global Cp
xdata = xdata.';
dt(2:length(xdata)) = diff(xdata);%(2)-xdata(1); % determination de l'intervalle de temps pour l'integrale
dt(1) = dt(2);
dt = dt.';
Ct = zeros(numel(xdata),1,'single');

if length(BETA)==3
    BETA(1) = abs(BETA(1)); % on interdit les valeurs negatives pour le volume sanguin
    % determination de Ct par convolution
    % code Jan
    Exp_x_b3 = exp(xdata*BETA(3));
    Ct = BETA(1)*Cp + BETA(2)*(cumsum(Cp.*Exp_x_b3.*dt)./Exp_x_b3);
%     x_b3 = xdata*BETA(3);
%     Ct = BETA(1)*Cp + BETA(2)*exp(-x_b3).*cumsum(Cp.*exp(x_b3).*dt);
%     for i=1:numel(xdata)
%         Ct(i) = BETA(1).*Cp(i).' +  BETA(2).*sum(Cp(1:i).'.*exp(-(xdata(i)-xdata(1:i)).*BETA(3)).*dt(1:i));
% %         Ct(i,1)=  BETA(1).*Cp(i).' +  BETA(2).*sum(Cp(1:i).'.*exp(-BETA(2) * (xdata(i)-xdata(1:i))./BETA(3)).*dt(1:i));
%     end
end

return


