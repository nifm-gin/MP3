function [Fx,Fy]=backwards2forwards_2d_double(Bx,By,H)
% This function will turn a backwards transformation field into 
% a forwards transformation field.
%
% [Fx,Fy]=backwards2forwards_2d_double(Bx,By,H); 
%
% inputs,
%   Bx,By : The backward transformation fields
%   H : The Splatting kernel 
% outputs,
%   Fx,Fy : The forward transformation fields
%
% Function is written by D.Kroon University of Twente (Februari 2009)

% Gaussian kernel
hx_center=-floor(size(H,1)/2)-1; 
hy_center=-floor(size(H,2)/2)-1; 

Fx=zeros(size(Bx)); Fy=zeros(size(Bx));
Num=zeros(size(Bx));
[x,y]=ndgrid(1:size(Bx,1),1:size(Bx,2));
x=x(:); y=y(:);
nx=x+Bx(:); ny=y+By(:);
valx=-Bx(:); valy=-By(:);
for i=1:length(nx);
        Tlocalx=nx(i); Tlocaly=ny(i);
        % All the neighborh pixels involved in linear interpolation.
        xBas0=floor(Tlocalx);  yBas0=floor(Tlocaly);
        xBas1=xBas0+1; yBas1=yBas0+1;
        % Linear interpolation constants (percentages)
        xCom=Tlocalx-xBas0; yCom=Tlocaly-yBas0;
        perc0=(1-xCom).*(1-yCom); perc1=(1-xCom).*yCom; perc2=xCom.*(1-yCom); perc3=xCom.*yCom;
        for iHx=1:size(H,1)
            for iHy=1:size(H,2)
                for t=1:4
                    switch(t),
                        case(1), tx=xBas0+iHx+hx_center; ty=yBas0+iHy+hy_center; perc=perc0*H(iHx,iHy);
                        case(2), tx=xBas0+iHx+hx_center; ty=yBas1+iHy+hy_center; perc=perc1*H(iHx,iHy);
                        case(3), tx=xBas1+iHx+hx_center; ty=yBas0+iHy+hy_center; perc=perc2*H(iHx,iHy);
                        case(4), tx=xBas1+iHx+hx_center; ty=yBas1+iHy+hy_center; perc=perc3*H(iHx,iHy);
                    end
                    if(~((tx<1)||(tx>size(Bx,1))||(ty<1)||(ty>size(Bx,2)))), 
                        Fx(tx,ty)=Fx(tx,ty)+valx(i)*perc; Fy(tx,ty)=Fy(tx,ty)+valy(i)*perc;
                        Num(tx,ty)=Num(tx,ty)+perc;
                    end
                end
            end
        end
end
Fx=Fx./(Num+eps);
Fy=Fy./(Num+eps);
