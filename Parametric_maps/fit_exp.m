function [t2et,m01,convergence] = fit_exp(xdata,ydata,n)

xdata=xdata(n:end);
ydata=ydata(n:end);
% [var1, bbb1]=polyfit(xdata,log(abs(ydata)),1);
% m01=exp(var1(2));
% t2s1=-1/var1(1);

t2s=(xdata(1)-xdata(end-1))/log(ydata(end-1)/ydata(1));
if t2s<=0 | isnan(t2s),
    ts2=30;
end
[aaa, bbb,  convergence]=levenbergmarquardt('AB_t2s',xdata, abs(ydata),[t2s max(abs(ydata))*1.5]);
% [aaa, bbb,  convergence]=levenbergmarquardt('AB_t2s',xdata, (ydata),[40 max(ydata)*1.5]);
t2et=(aaa(1));
m01=aaa(2);