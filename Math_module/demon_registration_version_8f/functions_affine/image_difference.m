function [t,I]=image_difference(V,U,type,Mask,MaskNum)
% This function gives a registration error and error image I between the two
% images or volumes.
%
% [t,I]=image_difference(I1,I2,type,Mask)
%
% inputs,
%   I1: Input image 1
%   I2: Input image 2
%   type: Type of similarity / error measure
% (optional)
%   Mask: Image/volume which is multiplied with the individual pixel errors
%         before calculation of the te total (mean) similarity error.
%
% if type,
% 'd'  : differences between I1 and I2
% 'sd' : squared differences
% 'mi' : normalized mutual information
% 'mip': normalized mutual information with image split in multiple small regions
% 'gd' : gradient differences
% 'gc' : gradient correlation
% 'cc' : normalized cros correlation
% 'pi' : pattern intensity
% 'ld' : log absolute difference
%
%  Example,
%    I1=im2double(imread('lenag1.png'));
%    I2=im2double(imread('lenag2.png'));
%    [t,I] = image_difference(I1,I2,'sd');
%    disp(t);
%    imshow(I,[])
%
% This function is written by D.Kroon University of Twente (April 2009)

if(exist('type','var')==0), type='p'; end
if(exist('Mask','var')==0), Mask=[]; end

if(exist('MaskNum','var')==0), MaskNum=numel(V); end

if(MaskNum==0), t=2; I=2; return; end
if(isempty(V)), t=2; I=2; return; end

% If color image make mask for RGB
if(~isempty(Mask)&&(size(V,3)==3)&&(size(Mask,3)==1)) 
	Mask_rgb(:,:,1)=Mask; Mask_rgb(:,:,2)=Mask; Mask_rgb(:,:,3)=Mask;
	Mask=Mask_rgb;
end

switch(type)
    case 'd'
        if ( nargout > 1 )
            [t,I]=registration_error_differences(V,U,Mask,MaskNum);
        else
            t=registration_error_differences(V,U,Mask,MaskNum);
        end
    case 'sd'
        if ( nargout > 1 )
            [t,I]=registration_error_squared_differences(V,U,Mask,MaskNum);
        else
            t=registration_error_squared_differences(V,U,Mask,MaskNum);
        end
        
    case 'mi'
        if ( nargout > 1 )
            [t,I]=registration_error_mutual_info(V,U,Mask);
        else
            t=registration_error_mutual_info(V,U,Mask);
        end
    case 'mip'
        if ( nargout > 1 )
            [t,I]=registration_error_local_mutual_info(V,U,Mask);
        else
            t=registration_error_local_mutual_info(V,U,Mask);
        end
        
    case 'gd'
        if ( nargout > 1 )
            [t,I]=registration_error_gradient_difference(V,U,Mask,MaskNum);
            I=2-I;
        else
            t=registration_error_gradient_difference(V,U,Mask,MaskNum);
        end
        t=2-t;
    case 'gc'
        if ( nargout > 1 )
            [t,I]=registration_error_gradient_correlation(V,U,Mask,MaskNum);
            I=1-I;
        else
            t=registration_error_gradient_correlation(V,U,Mask,MaskNum);
        end
        t=1-t;
    case 'cc'
        if ( nargout > 1 )
            [t,I]=registration_error_normalized_cross_correlation(V,U,Mask,MaskNum);
            I=1-I;
        else
            t=registration_error_normalized_cross_correlation(V,U,Mask,MaskNum);
        end
        t=1-t;
    case 'pi'
        if ( nargout > 1 )
            [t,I]=registration_error_pattern_intensity(V,U,Mask,MaskNum);
            I=1-I;
        else
            t=registration_error_pattern_intensity(V,U,Mask,MaskNum);
        end
        t=1-t;
    case 'ld'
        if ( nargout > 1 )
            [t,I]=registration_error_log_absolute_distance(V,U,Mask,MaskNum);
        else
            t=registration_error_log_absolute_distance(V,U,Mask,MaskNum);
        end
    otherwise
        error('Unknown error type')
end
if(isnan(t)), warning('imagedifference:NaN','NaN in error image'); t=2; I=2; end


function [t,I]=registration_error_log_absolute_distance(V,U,Mask,MaskNum)
I=log(abs(V-U)+1);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [t,I]=registration_error_normalized_cross_correlation(V,U,Mask,MaskNum)
Vvar=V-mean(V(:)); Uvar=U-mean(U(:));
I=(Vvar.*Uvar)/((sqrt(sum(Vvar(:).^2))*sqrt(sum(Uvar(:).^2)))+eps);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [t,I]=registration_error_gradient_correlation(V,U,Mask,MaskNum)
if(size(U,3)<4)
    [Gx,Gy]=sobel2();
    I=zeros(size(U));
    for i=1:size(U,3)
        if(i==1)
            I(:,:,i)=(1/2)*(registration_error_normalized_cross_correlation(conv2(V(:,:,i),Gx,'same'),conv2(U(:,:,i),Gx,'same'))...
                +registration_error_normalized_cross_correlation(conv2(V(:,:,i),Gy,'same'),conv2(U(:,:,i),Gy,'same')));
        end
    end
else
    [Gx,Gy,Gz]=sobel3();
    I=(1/3)*(registration_error_normalized_cross_correlation(convn(V,Gx,'same'),convn(U,Gx,'same'))...
        +registration_error_normalized_cross_correlation(convn(V,Gy,'same'),convn(U,Gy,'same'))...
        +registration_error_normalized_cross_correlation(convn(V,Gz,'same'),convn(U,Gz,'same')));
end
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [Gx,Gy]=sobel2()
Gx=[1 0 -1;2 0 -2;1 0 -1];  Gy=[1 2 1;0 0 0;-1 -2 -1];

function [Gx,Gy,Gz]=sobel3()
Gx=zeros(3,3,3);Gy=zeros(3,3,3); Gz=zeros(3,3,3);
Gx(:,:,1)=[-1 -3 -1;-3 -6 -3;-1 -3 -1]; Gx(:,:,2)=[ 0  0  0; 0  0  0; 0  0  0]; Gx(:,:,3)=[ 1  3  1; 3  6  3; 1  3  1];
Gy(:,:,1)=[ 1  3  1; 0  0  0;-1 -3 -1]; Gy(:,:,2)=[ 3  6  3; 0  0  0;-3 -6 -3]; Gy(:,:,3)=[ 1  3  1; 0  0  0;-1 -3 -1];
Gz(:,:,1)=[-1  0  1;-3  0  3;-1  0  1]; Gz(:,:,2)=[-3  0  3;-6  0  6;-3  0  3]; Gz(:,:,3)=[-1  0  1;-3  0  3;-1  0  1];

function [t,I]=registration_error_squared_differences(V,U,Mask,MaskNum)
if(isempty(Mask)&&(nargout==1))
    if(isa(V,'double'))
        t=squared_difference_double(double(V),double(U));
    else
        t=squared_difference_single(single(V),single(U));
    end
else
    I=(V-U).^2;
    if(~isempty(Mask)), I=I.*Mask;  end
    t=sum(I(:))/(MaskNum);
end

function [t,I]=registration_error_differences(V,U,Mask,MaskNum)
I=(V-U);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [t,I]=registration_error_pattern_intensity(V,U,Mask,MaskNum)
Idiff=V./(mean(V(:))+1e-5)-U./(mean(U(:))+1e-5);
o=0.3; %determines if grey-value varion is a structure (must be laster than noise)
r=5; numr=0; listr=[];
if(size(U,3)<4)
    if(size(U,3)>1), U=rgb2gray(U); V=rgb2gray(V); Idiff=V./(mean(V(:))+1e-5)-U./(mean(U(:))+1e-5); end
    
    for u=-r:r
        for v=-r:r
            if((u^2+v^2)<=r^2), numr=numr+1; listr(numr,:)=[u v]; end
        end
    end
    for u=-r:r
        for v=-r:r
            if((u^2+v^2)<=r^2), numr=numr+1; listr(numr,:)=[u v]; end
        end
    end
    SP1=zeros(size(U)-2*r,class(Idiff));
    for i=1:size(listr,1),
        u=listr(i,1); v=listr(i,2);
        SP1=SP1+o^2./(o^2+(Idiff(1+r:end-r,1+r:end-r)-Idiff(1+r+u:end-r+u,1+r+v:end-r+v)).^2);
    end
else
    for u=-r:r,
        for v=-r:r
            for w=-r:r
                if((u^2+v^2+w^2)<=r^2), numr=numr+1; listr(numr,:)=[u v w]; end
            end
        end
    end
    SP1=zeros(size(U)-2*r,class(Idiff));
    for i=1:size(listr,1),
        u=listr(i,1); v=listr(i,2); w=listr(i,3);
        SP1=SP1+o^2./(o^2+(Idiff(1+r:end-r,1+r:end-r,1+r:end-r)-Idiff(1+r+u:end-r+u,1+r+v:end-r+v,1+r+w:end-r+w)).^2);
    end
end
I=SP1./(size(listr,1)+eps);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;


function [t,I]=registration_error_gradient_difference(V,U,Mask,MaskNum)
if(size(U,3)<4)
    Sgdiff=zeros(size(U));
    for i=1:size(U,3)
        a=mean(V(:))/mean(U(:));
        [Gx,Gy]=sobel2();
        Idiffv=conv2(V(:,:,i),Gx,'same')-a*conv2(U(:,:,i),Gx,'same');
        Idiffh=conv2(V(:,:,i),Gy,'same')-a*conv2(U(:,:,i),Gy,'same');
        Av=var(Idiffv(:));
        Ah=var(Idiffh(:));
        Sgdiff(:,:,i)=(Av./(Av+Idiffv.^2+eps))+(Ah./(Ah+Idiffh.^2+eps));
    end
else
    a=mean(V(:))/mean(U(:));
    [Gx,Gy,Gz]=sobel3();
    Idiffv=convn(V,Gx,'same')-a*convn(U,Gx,'same');
    Idiffh=convn(V,Gy,'same')-a*convn(U,Gy,'same');
    Idiffz=convn(V,Gz,'same')-a*convn(U,Gz,'same');
    Av=var(Idiffv(:));
    Ah=var(Idiffh(:));
    Az=var(Idiffz(:));
    Sgdiff=(Av./(Av+Idiffv.^2+eps))+(Ah./(Ah+Idiffh.^2+eps))+(Az./(Az+Idiffz.^2+eps));
end
I=Sgdiff;
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;


function [t,I]=registration_error_local_mutual_info(V,U,Mask)
% Split the image in multiple regions to allow a more local mutual
% information measure.
numreg=8;
t=0;
x1=round(linspace(1,size(V,1),numreg+1)); x2=x1(2:end);
y1=round(linspace(1,size(V,2),numreg+1)); y2=y1(2:end);
z1=round(linspace(1,size(V,3),numreg+1)); z2=z1(2:end);
if(size(V,3)<4)
    for i=1:numreg
        for j=1:numreg
            Vpart=V(x1(i):x2(i),y1(j):y2(j),:); Upart=U(x1(i):x2(i),y1(j):y2(j),:);
            if(~isempty(Mask)), Maskpart=Mask(x1(i):x2(i),y1(j):y2(j),:); else Maskpart=[]; end
            [tpart,temp]=registration_error_mutual_info(Vpart,Upart,Maskpart);
            t=t+tpart;
        end
    end
    t=t/(numreg^2);
else
    for i=1:numreg
        for j=1:numreg
            for k=1:numreg
                Vpart=V(x1(i):x2(i),y1(j):y2(j),z1(k):z2(k)); Upart=U(x1(i):x2(i),y1(j):y2(j),z1(k):z2(k));
                if(~isempty(Mask)), Maskpart=Mask(x1(i):x2(i),y1(j):y2(j),z1(k):z2(k)); else Maskpart=[]; end
                [tpart,temp]=registration_error_mutual_info(Vpart,Upart,Maskpart);
                t=t+tpart;
            end
        end
    end
    t=t/(numreg^3);
end
I=[];

function [t,I]=registration_error_mutual_info(V,U,Mask)
% This function t=registration_error_mutual_info(V,U) gives a registration error
% value based on mutual information (H(A) + H(B)) / H(A,B)

% Make a joint image histogram and single image histograms
%bins=round(numel(V)^(1/ndims(V)));
bins=128;

% Remove unmasked pixels
if(~isempty(Mask)),
    sizes=size(V);
    V=V(:); V(~Mask)=[]; U=U(:); U(~Mask)=[]; 
    if(isempty(V)), t=0; I=[]; return, end
    V=reshape(V,sizes); U=reshape(U,sizes);
end

if(size(U,3)==3),
    % Colored image split into Hue en Intensity
    V_Hue=atan(sqrt(3)*(V(:,:,2)-V(:,:,3))./(2*(V(:,:,1)-V(:,:,2)-V(:,:,3)+eps)));
    V_Hue(~isfinite(V_Hue))=0;
    V_Int=(V(:,:,1)+V(:,:,2)+V(:,:,3))/3;
    
    U_Hue=atan(sqrt(3)*(U(:,:,2)-U(:,:,3))./(2*(U(:,:,1)-U(:,:,2)-U(:,:,3)+eps)));
    U_Hue(~isfinite(U_Hue))=0;
    U_Int=(U(:,:,1)+U(:,:,2)+U(:,:,3))/3;
    
    range=getrangefromclass(V);
    if(isa(V,'double'))
        [hist12, hist1, hist2]=mutual_histogram_double(double(V_Hue),double(U_Hue),double(range(1)),double(range(2)),double(bins));
    else
        [hist12, hist1, hist2]=mutual_histogram_single(single(V_Hue),single(U_Hue),single(range(1)),single(range(2)),single(bins));
    end
    % Calculate probabilities
    p1=hist1./numel(V); p2=hist2./numel(V); p12=hist12./numel(V);
    p1log=p1 .* log(p1+eps); p2log=p2 .* log(p2+eps); p12log=p12.* log(p12+eps);
    % Calculate amount of Information
    HA = -sum(p1log); HB = -sum(p2log); HAB = -sum(p12log(:));
    % Studholme, Normalized mutual information
    t_Hue=2-(HA+HB)/HAB;
    
    range=getrangefromclass(V);
    if(isa(V,'double'))
        [hist12, hist1, hist2]=mutual_histogram_double(double(V_Int),double(U_Int),double(range(1)),double(range(2)),double(bins));
    else
        [hist12, hist1, hist2]=mutual_histogram_single(single(V_Int),single(U_Int),single(range(1)),single(range(2)),single(bins));
    end
    % Calculate probabilities
    %p1=hist1./numel(V); p2=hist2./numel(V); p12=hist12./numel(V);
    p1=hist1; p2=hist2; p12=hist12;
    p1log=p1 .* log(p1+eps); p2log=p2 .* log(p2+eps); p12log=p12.* log(p12+eps);
    % Calculate amount of Information
    HA = -sum(p1log); HB = -sum(p2log); HAB = -sum(p12log(:));
    % Studholme, Normalized mutual information
    t_Int=2-(HA+HB)/HAB;
    
    t=t_Hue+t_Int;
    I=[];
else
    range=getrangefromclass(V);
    if(isa(V,'double'))
        [hist12, hist1, hist2]=mutual_histogram_double(double(V),double(U),double(range(1)),double(range(2)),double(bins));
    else
        [hist12, hist1, hist2]=mutual_histogram_single(single(V),single(U),single(range(1)),single(range(2)),single(bins));
    end
    % Calculate probabilities
    p1=hist1; p2=hist2; p12=hist12;
    %p1=hist1./numel(V); p2=hist2./numel(V); p12=hist12./numel(V);
    p1log=p1 .* log(p1+eps); 
    p2log=p2 .* log(p2+eps);
    p12log=p12.* log(p12+eps);
    % Calculate amount of Information
    HA = -sum(p1log);
    HB = -sum(p2log); 
    HAB = -sum(p12log(:));
    % Studholme, Normalized mutual information
    t=2-(HA+HB)/HAB; I=[];
end
if(isnan(t)), t=0; end

