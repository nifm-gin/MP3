function y=AB_expt1_v1(xdat,BETA)
y=abs( BETA(2) * ( 1 - 2 * exp(-xdat/BETA(1))) );
return

% y         Mz
% BETA(1)   T1
% BETA(2)   Mz a l'equilibre