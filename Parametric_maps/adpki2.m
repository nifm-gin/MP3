function X=adpki2(xfix,xvar)

Mb0=xvar(1);
f=xvar(2);
tau=xvar(3);
T1t=xvar(4);

l=xfix(1);
a0=xfix(2);
T=xfix(3);
step=xfix(4);
size=xfix(5);
T1b = xfix(6);

%T1m=xfix(4);

Ts2=T/2;

fl=f/l;
Z0=Mb0/T1t+fl*Mb0;
A= -1/T1t-fl    ;
eATs2=expm(A*Ts2);
B=-2*fl*Mb0 *a0*exp(-tau/T1b);
F=eATs2*(eATs2-1)*B/A; %A negatif donc signe inverse et A=-1/T1sat
AZ0=Z0/A;
Xf=F/(1-expm(A*T)); %T=2*delta (publi)
Xm=eATs2*Xf + (eATs2-1)*B/A;

X=[];

% for t=T+step:step:T+tau
%     X=[X expm(A*(t-Ts2-tau))*Xm];
%     trest = t;
% end
% 
% for t=tau + step:step:Ts2+tau
%     X=[X expm(A*(t-tau))*Xf + (expm(A*(t-tau))-1)*B/A];
%     trest = t;
% end
% 
% for t=Ts2+tau+step:step:T
%     X=[X expm(A*(t-Ts2-tau))*Xm];
%     trest = t;
% end

for t=T:step:T+tau -step
    X=[X expm(A*(t-Ts2-tau))*Xm];
    trest = t;
end

for t=tau:step:Ts2+tau-step
    X=[X expm(A*(t-tau))*Xf + (expm(A*(t-tau))-1)*B/A];
    trest = t;
end

for t=Ts2+tau:step:T-step
    X=[X expm(A*(t-Ts2-tau))*Xm];
    trest = t;
end


% for t=T:step:T+tau-step
%     X=[X expm(A*(t-Ts2-tau))*Xm];
%     trest = t;
% end
% 
% for t=tau :step:Ts2+tau-step
%     X=[X expm(A*(t-tau))*Xf + (expm(A*(t-tau))-1)*B/A];
%     trest = t;
% end
% 
% for t=Ts2+tau:step:T-step
%     X=[X expm(A*(t-Ts2-tau))*Xm];
%     trest = t;
% end

%% Fin
if numel(X)< size
    nb = size-numel(X);
   for t=trest+step:step:trest+nb*step
       X=[X expm(A*(t-Ts2-tau))*Xm];
   end
end

if numel(X)>size
   X=X(1:size); 
end

X=X-AZ0;
%X=X+Mb0;

%figure(200), plot(0:0.05:4-0.05,X);