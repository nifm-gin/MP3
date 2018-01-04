function y=fit_T1_3param(xdat,BETA)
y=abs( BETA(2) * ( 1 - 2 * BETA(3) * exp(-xdat/BETA(1))) );
return

% y         Mz
% BETA(1)   T1
% BETA(2)   M0z : Mz a l'equilibre
% BETA(3)   alpha : efficacite d'inversion