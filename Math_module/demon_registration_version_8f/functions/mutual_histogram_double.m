function [hist12, hist1, hist2]=mutual_histogram_double(I1,I2,Imin,Imax,nbins)
% This function makes a 2D joint histogram of 1D,2D...ND images
% and also calculates the seperate histograms of both images.
%
% [hist12, hist1, hist2]=mutual_histogram_double(I1,I2,Imin,Imax,nbins)
%
% Function is written by D.Kroon University of Twente (February 2009)
%

% Number of bins must be integer
nbins=round(nbins);

% scaling value 
scav=nbins/(Imax-Imin);

% Indices (1D) of all pixels
index=1:numel(I1);

% Calculate histogram positions
xd=scav*(I1(index)-Imin);
yd=scav*(I2(index)-Imin);

% Calculate both neighbors and interpolation percentages
xm=floor(xd); xp=xm+1;
xmd=xp-xd; xpd=xd-xm;
ym=floor(yd); yp=ym+1;
ymd=yp-yd; ypd=yd-ym;

% Fit to range ...
xm(xm<0)=0; xp(xp<0)=0;
xm(xm>(nbins-1))=(nbins-1); xp(xp>(nbins-1))=(nbins-1);
ym(ym<0)=0; yp(yp<0)=0;
ym(ym>(nbins-1))=(nbins-1); yp(yp>(nbins-1))=(nbins-1);

xm=xm+1; ym=ym+1;
xp=xp+1; yp=yp+1;


hist12=zeros(nbins,nbins);
hist1=zeros(nbins,1);
hist2=zeros(nbins,1);

for i=1:numel(I1),
    hist1(xm(i))=hist1(xm(i))+xmd(i); 
    hist1(xp(i))=hist1(xp(i))+xpd(i);
    hist2(ym(i))=hist2(ym(i))+ymd(i); 
    hist2(yp(i))=hist2(yp(i))+ypd(i);
    hist12(xm(i),ym(i))=hist12(xm(i),ym(i))+xmd(i).*ymd(i);
    hist12(xp(i),ym(i))=hist12(xp(i),ym(i))+xpd(i).*ymd(i);
    hist12(xm(i),yp(i))=hist12(xm(i),yp(i))+xmd(i).*ypd(i);
    hist12(xp(i),yp(i))=hist12(xp(i),yp(i))+xpd(i).*ypd(i);
end

hist1=hist1./numel(I1); 
hist2=hist2./numel(I1); 
hist12=hist12./numel(I1);
  
 
