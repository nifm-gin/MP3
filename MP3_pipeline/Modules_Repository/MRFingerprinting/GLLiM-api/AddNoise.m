function [X_noisy, real_snr] = AddNoise(X, snr, verb)
% Add gaussian white noise to X. snr is the signal-to-noise ratio
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

X_noisy = abs(X + randn(size(X)) .* repmat(max(X, [], 2) ./ snr, 1,size(X,2)));

real_snr = max(X,[],2) ./ std(X - X_noisy,[],2);

        