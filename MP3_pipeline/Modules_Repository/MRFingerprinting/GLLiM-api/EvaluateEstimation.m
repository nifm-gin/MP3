function [Rmse, Nrmse, Mae, Nmae] = EvaluateEstimation(Ytrue, Yestim)

Rmse    = nanmean( (Ytrue - Yestim).^2 ).^.5;
Nrmse   = Rmse ./ nanmean( (Ytrue - nanmean(Ytrue)).^2 ).^.5;

Mae     = nanmean( abs(Ytrue - Yestim) );
Nmae    = Mae ./ nanmean( Ytrue - nanmean(Ytrue) );

end