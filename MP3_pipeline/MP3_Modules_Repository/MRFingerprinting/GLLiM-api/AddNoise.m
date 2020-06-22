function [X_noisy, real_snr] = AddNoise(X, snr, verb)
% Add gaussian white noise to X. SNR is the signal-to-noise ratio
%
% inputs: 
%   - X: N*T matrix signals to noise
%   - snr: number ~= 0 
%   - verb: ~= 0 to see comments
%
% output:
%   - X_noisy: N*T matrix signals = X + noise

narginchk(2, 3);
if nargin == 2, verb = 0; end

if ~any(imag(X) ~= 0)
    X_noisy = abs(X + randn(size(X)) .* repmat(max(abs(X),[],2)./snr, 1,size(X,2)));
    
    real_snr = max(X,[],2) ./ std(X - X_noisy, [],2);
    
else
    X_noisy = complex(real(X) + randn(size(X)) .* max(abs(X),[],2)./snr, ...
                      imag(X) + randn(size(X)) .* max(abs(X),[],2)./snr);
                  
    real_snr = max(abs(X),[],2) ./ std(abs(X) - abs(X_noisy), [],2);
end
        