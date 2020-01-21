function [Rmse, Nrmse, Mae, Nmae] = EvaluateEstimation(Ytrue, Yestim)

N       = size(Ytrue,2);
Rmse    = zeros(1, N);
Nrmse   = Rmse;

% TODO: can be optimized = remove this loop
for i = 1:N
    Rmse(i)  = nanmean( ( Ytrue(:,i) - Yestim(:,i) ).^2 )^.5;
    Nrmse(i) = Rmse(i) / nanmean( ( Ytrue(:,i) - nanmean(Ytrue(:,i)) ).^2 ).^.5;
    
    Mae(i)   = nanmean( abs(Ytrue(:,i) - Yestim(:,i)) );
    Nmae(i)  = Mae(i) ./ nanmean( Ytrue(:,i) - nanmean(Ytrue(:,i)) );
end
