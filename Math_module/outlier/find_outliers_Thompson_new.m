function [Index] = find_outliers_Thompson_new(X, alpha, method, pri)
% find_outliers_Thompson
% Find outliers in a data series using the "modified Thompson's Tau method"
% 
% find_outliers_Thompson(X)
%   For vectors or matrix with at least 3 valid values (the routine manages 
%   NaNs), find_outliers_Thompson(X) gives the indexes of elements in X  
%   that are considered outliers as defined by the modified Thompson Tau 
%   method.
%   If no outliers are found an empty matrix is returned.
%   The indexes of the outliers are ordered with respect to their severity.
%   When X is matrix it is taken as a whole sample. To operate on single
%   columns or rows, a series of calls to the routine are needed.
% 
% find_outliers_Thompson(X, alpha)
%   As additional input alpha, the critical value for the test, can be
%   specified (a value between 0 and 0.5 is expected). Otherwise a default
%   value of 0.01(1%) will be used.
%   To skip this use the empty value [].
% 
% find_outliers_Thompson(X, alpha, method)
%   As third input the methods for the calculation of the sample statistics
%   on wich the test is performed can be specified. 
%   The three methods gives usually the same result even if some data can 
%   be an outlier or not depending on the selected method. In case of 
%   outliers the 'biweight'(default) is expected to be the most reilable 
%   method. Available values for 'method' are:
%    'median'  the statistics used are median and pseudo-deviation, very 
%              robust to outliers. The pseudo-deviation is computed 
%              with respect to the normal distribution.
%    'mean'    in this case the statistic is not robust to outliers but it 
%              can be useful to compare the results of robust statistics 
%              with the more used versions of the Tau test.
%              In particular, with the mean statistic if a 3 sample data is
%              given and it includes a very big outlier it will be accepted
%              anyway.
%    'biweight' (default statistics) is only a bit less robust than the 
%               median but it is more efficient for the definition of the 
%               center of the distribution. It is more time consuming then 
%               the other two methods. In this case biweight mean and 
%               biweight standard deviation are used.
%   To skip this use the empty value [].
% 
% find_outliers_Thompson(X, ..., ..., pri)
%   With pri = 1 the routine plots a graph showing the original X values, 
%   the outliers identified with the test, if any, the final deviation 
%   around the center of the sample and the final critical bounders. 
%   If pri = 0 (as default) the routine doesn't plot any figure.
% 
% 
% WARNING: if the statistical toolbox is not installed the routine works
%          aniway but alpha will be fixed to 0.01. A linear interpolation 
%          of table values is used.
% 
%   Example 1: 
%   If 
%     X = [44 -5 45 48 41 45 42 43 44 46 46 150]'
% 
%   then the routine will return the vector: 
% 
%     Index = [12; 2]
% 
%
%   Example 2: 
%   If 
%     X = [34 -5 35 35 33 34; 37 38 35 35 36 150]'
% 
%   then the routine will anyway return the vector: 
% 
%     Index = [12; 2]
%   
%   To obtain a separate analysis of the two columns write:
%   
%     Index1 = find_outliers_Thompson(X(:,1), alpha)
%     Index2 = find_outliers_Thompson(X(:,2), alpha)
%
%   To obtain a separate analysis on many columns:
% 
%     Nc = size(X, 2); % number of columns of X
%     for n = 1:Nc; a = find_outliers_Thompson(X(:,n), alpha);
%     Index1(1:numel(a),n) = a; end
% 
% 
% References:
%
% on the test
%    Measurement Uncertainty, Part I, ASME PTC 19.1 1998
% recommended for the individuation of outliers in a set of repeted
% measurements.
%    Thompson R. 1985. A note on restricted maximum likelihood estimation
% with an alternative outlier model. Journal of the Royal Statistical 
% Society. Series B vol 47, 53-55
%    Christensen R., Pearson LM., Johnson W. 1992. Case deletion
% diagnostics for mixed models. Technometyrics. Vol 34, 38-45.
% 
% on the robust methods and the biweight mean:
%    Lanzante JR. 1996. Resistant, robust and non-parametric techniques for
% the analysis of climate data: theory and examples, including applications
% to historical radiosonde station data. International Journal of
% Climatology, vol.16, 1197-1226.
%
% 
% Author:  Michele Rienzner
% e-mail:  michele.rienzner@unimi.it
% Release: 3.1
% Release date: 5/17/2010

%% checkup script: use it to verify how the test works
% The four lines below generates a vector (X) containing at the most 3
% outliers in the positions 2, 5 and 10 in a sample of 20 data values. and 
% 3 values are missing (index 13:18). Than the routine is called and an 
% estimation of the outliers is donem their indexes are displayed and a
% figure is plotted. This is repeated 5 times giving different results.
% 
% alpha = 0.01; method = 'biweight';
% for n= 1:5; X = randn(20,1); X([2 5 10]) = X([1 5 10])*10; X(13:18) = NaN;
% [Index] = find_outliers_Thompson(X, alpha, method, 1), end
% 

%% Notes on the algorithm:
% 
% tau = t*(N-1) / (sqrt(N)*sqrt(N-2+t^2));
% where t is the value of the Student's-t with cumulate non-passing 
% probability of 1-alpha/2 and N-2 degrees of freedom
% 
% Suspicious = max(abs(x-mean(x)));
% S = std(x);
% if Suspicious > tau*S
% Suspicious is an outlier, its index is recorded and the value is erased 
% in x
% 
% the former steps are executed again until no new outlier is found.


%% Default statements
Def_alpha  = 0.01;          % usual value for alpha
Tab_alpha  = [1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 0.001 0.005 0.01 0.05];  % tabulated tau values available for some alpha values
Def_method = 'biweight';    % method for the computation of mean and deviation
Def_pri    = NaN;             % default is "do not plot anything"
X = X(:);                   % making X a column array

%% print figure
if nargin>3
    if ~isempty(pri) && (pri==1 || ishandle(pri))
        X0 = X; N0 = numel(X0); maxX = max(abs(X0));
        if ~ishandle(pri) || pri==1
            figure;
            pri = gca;
        end
        plot(pri,X0, '.b'); 
        hold(pri,'on');
        xlabel(pri,'Index of the values'); 
        ylabel(pri,'Values');
    end
end

%% initializations
I = find(~isnan(X)); 
X = X(I);                      % only the valid data are analyzed
N = numel(X);                  % determine the number of valid values in X
PrAllInd = 10;                 % preallocation size for Index
Index    = zeros(PrAllInd,1);  % preallocation for Index
h        = 0;                  % last writing position on Index;
statTool = license('test', 'statistics_toolbox'); % check the presence of the statistic toolbox
% statTool = false; % never use statistics license for this routine

%% Input control
if N<3
    error('Input data must have at least 3 valid values.');
end
if nargin<2
    alpha = [];
end
if isempty(alpha)    
    alpha = Def_alpha;
    if any(abs(alpha-Tab_alpha) < eps)
        alpha = find(abs(alpha-Tab_alpha) < eps,1); % make alpha an index into tau table values
        statTool = false; % do not use stat tool if tabulated alpha is available: save license usage
    end
elseif alpha<0 || alpha>0.5
    error('The second input "alpha" must be within 0 and 0.5 .')
elseif ~any(abs(alpha-Tab_alpha) < eps) % see if alpha is tabulated
    statTool = true;
else
    alpha = find(abs(alpha-Tab_alpha) < eps,1); % make alpha an index into tau table values
    statTool = false; % do not use stat tool if tabulated alpha is available: save license usage
end
if statTool
    [statTool,~]=license('checkout','statistics_toolbox'); % actually try to check the license out
    if ~statTool
        warning(['No statistics toolbox license available. To work with non-tabulated values of alpha the statistics toolbox is needed.\n',...
            'Reverting to default alpha of %f'],Def_alpha);           
        alpha = Def_alpha;
        if any(abs(alpha-Tab_alpha) < eps)
            alpha = find(abs(alpha-Tab_alpha) < eps,1); % make alpha an index into tau table values
            statTool = false; % do not use stat tool if tabulated alpha is available: save license usage
        end
    end
end
if nargin<3
    method = [];
end
if isempty(method)
    method = Def_method;
end
if nargin<4
    pri = [];
end
if isempty(pri)
    pri = Def_pri;
elseif pri~=0 && pri~=1 && ~ishandle(pri)
    error('Allowed values for ''pri'' are only 0 (do not plot), 1 (plot), or a valid graphics handle (axes to plot to)')
end

[X, indx] = sort(X);        % sorting of valid data values (this makes the routine simpler)
indI      = I(indx);        % indexing of the original data values in the sorted vector (needed to point the values on the original vector)
[Mr, Sr]  = meanStd_sub(X, method);      % robust values for mean and standard deviation
tau       = tau_sub(N, alpha, statTool); % value of the tau statistic

%% Find out the outliers
dev = abs(X([1,N])-Mr); % deviation from the mean of fist and last samples
while max(dev) > tau*Sr && N>=3
    if dev(1)-dev(2) > 0
        % the lower value is an outlier, its index is recorded while the
        % value is erased
        h = h+1;        Index(h) = indI(1);     
        X(1) = [];      indI(1) = [];   N = N-1;
    else
        % the higher value is an outlier, its index is recorded while the
        % value is erased
        h = h+1;        Index(h) = indI(N);     
        X(N) = [];      indI(N) = [];   N = N-1;
    end
    if N>=3
        [Mr, Sr]  = meanStd_sub(X, method);      % update the values for mean and standard deviation
        tau       = tau_sub(N, alpha, statTool); % update the value of the tau statistic
        dev       = abs(X([1,N])-Mr);            % update peak deviation from the mean
    end
end

%% Output cleanup
if h<PrAllInd % null elements of Index must be removed
    Index(h+1:PrAllInd) = [];
end

%% Plot of the outliers 
if ishandle(pri)
    n_max = fix(abs(maxX-Mr)/Sr);
    if numel(I)~=numel(X0)
        I1 = find(isnan(X0));
        plot(pri,I1, ones(numel(I1))*Mr, 's', 'MarkerSize', 3, 'MarkerFaceColor', ...
            [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5])
    end
    plot(pri,[1; N0], Mr+Sr*([1; 1]*(-n_max:n_max)), ':', 'Color', [0.75 0.75 0.75])            % plot the final deviations borders
    plot(pri,[1; N0], [Mr; Mr], 'k-')
    plot(pri,[1; N0], [Mr+tau*Sr, Mr-tau*Sr; Mr+tau*Sr, Mr-tau*Sr], 'r-')
    plot(pri,Index, X0(Index), 'or', 'MarkerSize', 8, 'LineWidth', 2) % highlight the outliers
    if statTool
        title(pri,['Find Outliers with Thompson Tau (\alpha=', num2str(alpha), ', Method=', method, ')']);
    else
        title(pri,['Find Outliers with Thompson Tau (\alpha=', num2str(Tab_alpha(alpha)), ', Method=', method, ')']);
    end
end

end

function tau = tau_sub(N, alpha, statTool)
% sub function computing the right value of Tau

if statTool                                 % the statistics toolbox is present 
    t = tinv(alpha/2, N-2);                 % finds the exact critical t-value
    tau = -t*(N-1)/(sqrt(N)*sqrt(N-2+t^2)); % computates tau
    
else % (alpha fixed = 0.01)
    % table of values 
    N_v = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 ...
        26 27 28 29 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 1000];
    
    tau_v = [ ...
        1.1547 1.5000 1.7887 2.0396 2.2608 2.4566 2.6299 2.7834 2.9197 ...
        3.0409 3.1493 3.2464 3.3339 3.4129 3.4846 3.5500 3.6097 3.6645 ...
        3.7149 3.7614 3.8045 3.8444 3.8816 3.9163 3.9487 3.9790 4.0075 ...
        4.0342 4.1471 4.2339 4.3027 4.3585 4.4046 4.4434 4.4766 4.5051 ...
        4.5300 4.5519 4.5713 4.5886 4.6041 4.6181 4.8637 ;...
        1.1547 1.5000 1.7884 2.0375 2.2545 2.4436 2.6083 2.7521 2.8780 ...
        2.9888 3.0868 3.1739 3.2518 3.3217 3.3848 3.4419 3.4939 3.5414 ...
        3.5850 3.6251 3.6620 3.6962 3.7279 3.7575 3.7850 3.8107 3.8349 ...
        3.8575 3.9525 4.0251 4.0824 4.1287 4.1670 4.1990 4.2263 4.2498 ...
        4.2703 4.2883 4.3042 4.3183 4.3311 4.3425 4.5422 ;...
        1.1547 1.5000 1.7881 2.0360 2.2503 2.4354 2.5954 2.7341 2.8548 ...
        2.9604 3.0534 3.1358 3.2091 3.2748 3.3339 3.3873 3.4358 3.4801 ...
        3.5206 3.5577 3.5920 3.6237 3.6530 3.6803 3.7057 3.7295 3.7517 ...
        3.7725 3.8599 3.9265 3.9789 4.0213 4.0561 4.0854 4.1103 4.1317 ...
        4.1503 4.1666 4.1811 4.1940 4.2055 4.2159 4.3968 ;...
        1.1547 1.4999 1.7867 2.0294 2.2344 2.4072 2.5533 2.6776 2.7841 ...
        2.8760 2.9561 3.0264 3.0884 3.1436 3.1929 3.2373 3.2774 3.3138 ...
        3.3470 3.3773 3.4053 3.4310 3.4548 3.4768 3.4973 3.5164 3.5343 ...
        3.5510 3.6207 3.6736 3.7150 3.7484 3.7758 3.7987 3.8182 3.8349 ...
        3.8494 3.8621 3.8734 3.8834 3.8924 3.9005 4.0400 ;...
        1.1547 1.4999 1.7854 2.0246 2.2237 2.3895 2.5281 2.6450 2.7444 ...
        2.8298 2.9037 2.9683 3.0252 3.0755 3.1205 3.1608 3.1972 3.2301 ...
        3.2601 3.2875 3.3127 3.3358 3.3572 3.3770 3.3954 3.4125 3.4285 ...
        3.4435 3.5058 3.5529 3.5898 3.6194 3.6437 3.6640 3.6813 3.6961 ...
        3.7089 3.7202 3.7301 3.7389 3.7469 3.7540 3.8769 ;...
        1.1547 1.4993 1.7789 2.0039 2.1836 2.3280 2.4453 2.5420 2.6228 ...
        2.6912 2.7497 2.8003 2.8444 2.8832 2.9177 2.9484 2.9759 3.0008 ...
        3.0234 3.0439 3.0627 3.0799 3.0958 3.1105 3.1241 3.1367 3.1485 ...
        3.1595 3.2052 3.2395 3.2662 3.2876 3.3051 3.3197 3.3321 3.3427 ...
        3.3519 3.3599 3.3670 3.3734 3.3790 3.3841 3.4711 ;...
        1.1547 1.4985 1.7730 1.9883 2.1564 2.2890 2.3954 2.4821 2.5539 ...
        2.6143 2.6658 2.7100 2.7485 2.7823 2.8121 2.8386 2.8624 2.8838 ...
        2.9032 2.9208 2.9369 2.9517 2.9653 2.9778 2.9894 3.0002 3.0103 ...
        3.0196 3.0584 3.0875 3.1102 3.1282 3.1430 3.1554 3.1658 3.1747 ...
        3.1825 3.1892 3.1952 3.2005 3.2052 3.2095 3.2824 ;...
        1.1547 1.4925 1.7424 1.9222 2.0536 2.1525 2.2291 2.2900 2.3394 ...
        2.3803 2.4147 2.4440 2.4693 2.4913 2.5106 2.5277 2.5429 2.5566 ...
        2.5689 2.5801 2.5902 2.5995 2.6081 2.6160 2.6232 2.6300 2.6362 ...
        2.6420 2.6661 2.6840 2.6979 2.7090 2.7180 2.7255 2.7319 2.7373 ...
        2.7420 2.7461 2.7497 2.7529 2.7558 2.7584 2.8022 ;...
        1.1546 1.4850 1.7150 1.8722 1.9832 2.0649 2.1271 2.1761 2.2155 ...
        2.2478 2.2749 2.2979 2.3176 2.3347 2.3497 2.3629 2.3747 2.3853 ...
        2.3948 2.4034 2.4112 2.4183 2.4249 2.4309 2.4365 2.4416 2.4464 ...
        2.4509 2.4692 2.4829 2.4934 2.5018 2.5087 2.5144 2.5192 2.5233 ...
        2.5268 2.5299 2.5326 2.5351 2.5372 2.5392 2.5722 ;...
        1.1511 1.4250 1.5712 1.6563 1.7110 1.7491 1.7770 1.7984 1.8153 ...
        1.8290 1.8403 1.8498 1.8579 1.8649 1.8710 1.8764 1.8811 1.8853 ...
        1.8891 1.8926 1.8957 1.8985 1.9011 1.9035 1.9057 1.9078 1.9096 ...
        1.9114 1.9186 1.9240 1.9281 1.9314 1.9340 1.9362 1.9381 1.9397 ...
        1.9411 1.9423 1.9433 1.9443 1.9451 1.9459 1.9586];    
    
    tau_v = tau_v(alpha,:); % choose appropriate row for alpha value
    
    % pick the better value from tau_v
    if N<=30 % the exact (rounded) value is present
        tau = tau_v(N-2);
    elseif N<=N_v(end) % an estimation is possible interpolating the table values
        I1 = find(N_v <  N, 1, 'last');
        I2 = find(N_v >= N, 1, 'first');
        t1 = tau_v(I1);
        t2 = tau_v(I2);
        N1 = N_v(I1);
        N2 = N_v(I2);
        tau = t1 + (t2-t1)/(N2-N1)*(N-N1);
%     elseif N<1000 % an estimation is possible interpolating between tau(100) and tau(1000)
%         tau = tau_v(end) + (tau1000-tau_v(end))/900*(N-100);
    else % after N=1000 the value of tau is quite flat
        tau = tau_v(end);
    end
end
end

function [Mr, Sr] = meanStd_sub(X, method)
% sub function computing Mean and standard deviation of the sample
% in a context of outliers a robust estimation of mean and standard
% deviation is usually suggested.

switch method 
    case 'mean'
        Mr = mean(X);
        Sr = std(X);
    case 'median'
        N   = numel(X);            % number of values in the sample
        % Mr  = median(X);           % robust estimation of the mean
        half = floor(N/2);         % do not use median function : X is already sorted, this is faster
        Mr = X(half+1);
        if 2*half == N             % Average if even number of elements
            Mr = (X(half)+Mr)/2;
        end
        
        % a   = quantile(X, [0.25, 0.75]); 
        % Sr  = (a(2)-a(1))/1.349;   % robust estimation of the standard deviation

        % robust estimation of the standard deviation using 25% and 75%
        % quantiles
        q  = [0 100*(0.5:(N-0.5))./N 100]';
        X  = [X(1); X(1:N); X(N)];
        X  = interp1q(q,X,[25 75]');
        Sr =  1.4826*diff(X)/2;   
    
    case 'biweight'
        N   = numel(X);            % number of values in the sample
	    c   = 7.5;                 % censor value, it corresponds to a certain number of standard 
        %                            deviations for a normal distribution:
                                   % c=6 is 4 std; c=7.5 is 5 std; c=9 is 6 std.
        % M   = median(X);         % median of the sample
        half = floor(N/2);         % do not use median function : X is already sorted, this is faster
        M = X(half+1);
        if 2*half == N             % Average if even number of elements
            M = (X(half)+M)/2;
        end
        
        XmM = X-M;                 % pre-calculate and store for speed
        % MAD = median(abs(XmM));    % median absolute deviation that is the median of the sample 
        %                            of the absolute values of the differences from the median
        x = sort(abs(XmM));        % do not use median function, only perform necessary operations here
        MAD = x(half+1);
        if 2*half == N             % Average if even number of elements
            MAD = (x(half)+MAD)/2;
        end
        w   = XmM/(c*MAD);         % weights for the computation of the biweight mean
        w(abs(w)>1) = 1;           % censoring of the wheights
        wsq = w.^2;                % pre-calculate these once and store for speed
        w2  = (1-wsq).^2;
        XmMw = XmM.*w2;
        Mr  = M+ sum(XmMw)/sum(w2); % computation of biwheight mean
        Sr  = sqrt(N*sum(XmMw.^2))/abs(sum((1-wsq).*(1-5*wsq))); % computation of biwheight std
    otherwise
        error('the method parameter has a wrong value; permitted values are: ''mean'', ''median'', ''biweight'' ')
end
end
