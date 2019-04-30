function [Lwopti] = FindOptimalLw(Xtrain, Ytrain, K, Lw, maxiter, cstr, verb)

if ~exist('Kmax','var'),    K       = 30; end
if ~exist('Lw','var'),      Lw      = 0:1:10; end
if ~exist('maxiter','var'), maxiter = 200; end
if ~exist('cstr','var'),    cstr.Sigma  = 'i';
                            cstr.Gammat = '';
                            cstr.Gammaw = ''; end
if ~exist('verb','var'),    verb    = 1; end


%% Compute the Kmax regression

Nrmse   = zeros(length(Lw), size(Ytrain,2));
bic     = zeros(1,length(Lw));
[N,Lt]  = size(Ytrain);
D       = size(Xtrain,2);

parfor_progress(length(Lw));
parfor k = 1:length(Lw)
    
    try
        [~,~,ll]	= EstimateInverseFunction(Ytrain, Xtrain, K, Lw(k), maxiter, cstr, 0);
        
        switch cstr.Sigma
            case 'i*'
                M = K* (D*(Lw(k)+Lt+1)     + Lt*(Lt+3)/2 + 1);
            case 'i'
                M = K* (D*(Lw(k)+Lt+1)     + Lt*(Lt+3)/2 + 1)  + K-1;
            case 'd*'
                M = K* (D*(Lw(k)+Lt+2)     + Lt*(Lt+3)/2) + D  + K-1;
            case 'd'
                M = K* (D*(Lw(k)+Lt+2)     + Lt*(Lt+3)/2)      + K-1;
            otherwise
                M = K* (D*(Lw(k)+Lt+D+1)   + Lt*(Lt+3)/2)      + K-1;
        end
        
        bic(k) = -2*ll + M*log(N)
    end
    
    parfor_progress;
end
parfor_progress(0);


%% Find the best K

[~,l]   = min(bic);
Lwopti  = Lw(l);


if verb == 1
    figure
    
    plot(Lw,bic, 'o-', 'LineWidth', 1.5)
    hold on
    line('YData', ylim, 'XData', [Lwopti Lwopti], 'Color','r')
end


