function [Theta] = updateSigma(Theta, var_noise)

Theta.Sigma = Theta.Sigma + var_noise * eye(size(Theta.Sigma(:,:,1)));
end

