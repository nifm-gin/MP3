function [Kopti] = FindOptimalK(Xtrain, Ytrain, K, Lw, maxiter, cstr, verb)

if ~exist('Kmax','var'),    K       = 2:2:50; end
if ~exist('Lw','var'),      Lw      = 0; end
if ~exist('maxiter','var'), maxiter = 200; end
if ~exist('cstr','var'),    cstr.Sigma  = 'i';
                            cstr.Gammat = '';
                            cstr.Gammaw = ''; end
if ~exist('verb','var'),    verb    = 1; end


%% Compute the Kmax regression

Nrmse   = zeros(length(K), size(Ytrain,2));
bic     = zeros(1,length(K));
[N,Lt]  = size(Ytrain);
D       = size(Xtrain,2);

parfor_progress(length(K));
parfor k = 1:length(K)
    
    try
        [~,~,ll]	= EstimateInverseFunction(Ytrain, Xtrain, K(k), Lw, maxiter, cstr, 0);
        
        switch cstr.Sigma
            case 'i*'
                M = k* (D*(Lw+Lt+1)     + Lt*(Lt+3)/2 + 1);
            case 'i'
                M = k* (D*(Lw+Lt+1)     + Lt*(Lt+3)/2 + 1)  + k-1;
            case 'd*'
                M = k* (D*(Lw+Lt+2)     + Lt*(Lt+3)/2) + D  + k-1;
            case 'd'
                M = k* (D*(Lw+Lt+2)     + Lt*(Lt+3)/2)      + k-1;
            otherwise
                M = k* (D*(Lw+Lt+D+1)   + Lt*(Lt+3)/2)      + k-1;
        end
        
        bic(k) = -2*ll + M*log(N)
    end
    
    parfor_progress;
end
parfor_progress(0);


%% Find the best K

[~,l]   = min(bic);
Kopti   = K(l);


if verb == 1
    figure
    
    plot(K,bic, 'o-', 'LineWidth', 1.5)
    hold on
    line('YData', ylim, 'XData', [Kopti Kopti], 'Color','r')
end


