function [B, D, STATUS]=levenbergmarquardt(funct,xdat,ydat,B_init)
% function [B, D, STATUS]=levenbergmarquardt(funct,xdat,ydat,B_init)
% function [B, D, STATUS]=levenbergmarquardt(pars,xdat,ydat)
% 
% Non-linear least-squares fit based on the Levenberg-Marquardt algorithm
% Input parameters:
%    funct    : function representing the model to fit the data:
%                  Y = funct(xdat,B)
%               must return an estimate Y of data given X and the column
%               vector B containing all model parameters. Y must be of the
%               same size as ydat.
%               The parameter 'funct' can be given either as a string
%               containing the name of the function to use, or as an inline
%               function (warning: SLOW!), or as a 'function_handle' to a
%               Matlab function (fastest).
%   pars      : alternative first input parameter for advanced use. If this
%               is a structure, it is assumed to be a completely
%               initialized structure of internal optimization parameters,
%               for complete control over the optimization process. 
%               The input parameter 'B_init' is superfluous in this case.
%               See the code for documentation on the fields of the
%               structure and an example of how to initialize it.
%   xdat      : X data corresponding to the observed data 'ydat'. This can
%               be any type of Matlab variable. It will be passed unchanged
%               as a parameter to 'funct'.
%   ydat      : observed data. Numeric array of the same size as the output
%               of 'funct'.
%   B_init    : initial estimate of the model parameters. Numeric array of
%               arbitrary shape, but with a number of elements identical to
%               that of the column vector 'B' passed to 'funct'.
%
% Output parameters:
%   B         : optimized model parameters. Numerical array of the same
%               size as as B_init.
%   D         : error on the parameter estimates. Numerical array of the
%               same size as B_init.
%   STATUS    : exit status. Numerical scalar with values
%                  -1 : the algorithm converged
%                  >0 : abort due to either failure to converge within the
%                       allowed number of iterations or a badly
%                       conditioned matrix of partial derivatives.
%
% Code by J. Warnking, INSERM, Grenoble, 2007
% based on code provided by E. Barbier, INSERM, Grenoble, 2007,
% which seems to be based on code by A. Fahlman, Naval Medical Research
% Center, Bethesda, MD and Carleton University Ottawa 
% (see http://handle.dtic.mil/100.2/ADA407508)

if isstruct(funct),
    pars = funct;
    B_init = pars.B;
else
    b_num      = numel(B_init);
    if ischar(funct),
        F      = eval(['@',funct]);        % create a 'function_handle' from the function name
    else                                   % this might already be an inline function or function_handle
        F      = funct;
    end;
    Y          = F(xdat,B_init(:));     % get model data w/ the initial guess of parameters

    % initialize the structure of optimization parameters. We'll carry this
    % around to exchange information across the recursive function calls
    pars       = struct(...
        'N',     numel(ydat),...       % number of Y data points
        'b_num', b_num,...             % number of parameters of the model
        'A',     zeros(b_num),...      % variance covariance matrix of partial derivatives
        'G',     zeros(b_num,1),...    % gradient vector to update paramters
        'D',     zeros(b_num,1),...    % vector of parameter std's
        'q',     [],...                % ??? (some sort of normalized covariance matrix)
        'B',     B_init(:),...         % current parameter guess
        'INC',   1+eye(b_num)*0.001,...% perturbation of B to estimate partial derivatives
        'L',     0.1*eye(b_num),...    % lambda
        'dL',    10,...                % factorial increment or decrement of lambda
        'CONV',  0.001,...             % convergence criterion
        'T1',    30,...                % limit of number of iterations (and output status) 
        'T2',    0,...                 % iteration counter
        'F',     F,...                 % model function
        'E',     [],...                % error vector for current parameter estimate
        'SO',    [],...                % sum of squares error of current parameter guess
        'verbose',0);                  % flag indicating whether progress messages should be printed
end;

if pars.verbose,
    format short g;
end;

% initialize errors for the initial parameter guess
Y          = pars.F(xdat,pars.B);      % get model data w/ the initial guess of parameters
pars.E     = ydat(:)-Y(:);             % error vector for current parameter estimate
pars.SO    = sum(pars.E.^2);           % sum of squares error of current parameter guess

pars       = mod1(pars,ydat,xdat);     % do the parameter search

B          = reshape(pars.B,size(B_init)); % prepare output parameters
D          = reshape(pars.D,size(B_init));
STATUS     = pars.T1;
return
%%%%%%%%%%%%%%END LEVENBERGMARQUARDT%%%%%%%%%%%%%%

%%%%%%%%%%%%%%BEGIN MOD1%%%%%%%%%%%%%%
function pars = mod1(pars,ydat,xdat)

pars.T2    = pars.T2+1;                % increase iteration counter
p          = zeros(pars.N,pars.b_num); % partial derivative is size= [N,b_num]
for i=1:pars.b_num
    B      = pars.B.*pars.INC(:,i);    % B gets changed temporarily
    e2     = ydat-pars.F(xdat,B);      % e2 is the error term w/ new parameters B
    %figure(100),plot(pars.F(xdat,pars.B)), hold on,plot(ydat,'+'), hold off
    p(:,i) = (pars.E-e2(:))/(B(i)-pars.B(i)); % creating the partial derivative for each y and beta
end

    
pars.A     = p'*p;                     % variance-covariance matrice %curvature matrix
pars.D     = sqrt(diag(pars.A));       % vector of std's
warning('off','MATLAB:divideByZero');  % We're (almost) always getting div' by zero here
pars.q     = pars.A./(pars.D*pars.D');
pars.G     = (p'*pars.E)./pars.D;      % gradient vector determining the next parameter guess
warning('on','MATLAB:divideByZero');
pars       = mod2(pars,ydat,xdat);
return
%%%%%%%%%%%%%%END MOD1%%%%%%%%%%%%%%

%%%%%%%%%%%%%%BEGIN MOD2%%%%%%%%%%%%%%
function pars = mod2(pars,ydat,xdat)

pars.q     = pars.q .* (1+pars.L);
if rcond(pars.q) > 5*eps
    p      = (inv(pars.q)*pars.G)./pars.D;
else
    % TODO: consider setting T1 to a value indicating malconditioned matrix
    pars       = mod3(pars,ydat,xdat); % still try to estimate parameter uncertainties
    if pars.verbose,
        disp('BADLY CONDITIONED MATRIX --- ABORT');
    end;
    return;                            % we can't go on without p. Abort.
end
% now decide how to continue
if abs(p./pars.B) < pars.CONV          % we have converged! Prepare output parameters
    pars.T1    = -1;                   % output status: converged
    pars       = mod3(pars,ydat,xdat);
    if pars.verbose,
        disp('CONVERGENCE --- parameter uncertainties :');
        disp(pars.D');
    end;
elseif pars.T2 < pars.T1               % some parameters have not yet converged and we still have iterations left.
    B          = pars.B + p;           % try new parameters by adding their partial derivative
    E          = ydat-pars.F(xdat,B);  % error of the model with new parameter estimates
    SO         = sum(E(:).^2);         % summing up the squared errors
    if pars.verbose,
        disp(['New SSE = ',num2str(SO),' for B = ']);
        disp(B');
    end;
    if SO > pars.SO                    % the new estimate is worse (ssq error increased)
        pars.L = pars.L*pars.dL;       % increment lambda by a factor dL
        if pars.verbose,
            disp(['Increasing lambda. New value: ',num2str(pars.L(1))]);
        end;
        pars   = mod2(pars,ydat,xdat); % call module2 again. This is not counted as an iteration.
    else                               % error got smaller
        pars.E = E(:);                 % save these values (no need to recompute them)
        pars.SO= SO;
        pars.B = B;                    % save improved model parameters
        pars.L = pars.L/pars.dL;
        if pars.verbose,
            disp(['Decreasing lambda. New value: ',num2str(pars.L(1))]);
        end;
        pars   = mod1(pars,ydat,xdat); % call mod1 again and start over from scratch, i.e. new iteration
    end
else                                   % abort due to exceeded number of iterations
    pars.T1    = pars.T1-1;            % keep for backwards compatibility (useless IMHO)
    pars       = mod3(pars,ydat,xdat); % still try to estimate parameter uncertainties
    if pars.verbose,
        disp('NUMBER OF ITERATIONS EXCEEDED');
    end;
end
return
%%%%%%%%%%%%%%END MOD2%%%%%%%%%%%%%%

%%%%%%%%%%%%%%BEGIN MOD3%%%%%%%%%%%%%%
function pars = mod3(pars,ydat,xdat)
% estimate uncertainty on final parameters

df         = pars.N-pars.b_num;    % the degrees of freedom
V          = pars.SO./df;          % summing up the squared errors, i.e. variance
if rcond(pars.A) > 5*eps
    pars.A = V*inv(pars.A);        % Var-covar matrix
    pars.D = sqrt(diag(pars.A));
else
    pars.D = zeros(pars.b_num,1);
end
return
%%%%%%%%%%%%%%END MOD3%%%%%%%%%%%%%%
