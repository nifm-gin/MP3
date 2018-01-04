function [E,Egrad]=demons_energy(M,F,M_TF,F_TF,T,sizes,options)
% This function DEMONS_ENERGY calculates the registration energy and
% gradient of the registration energy at each pixel. Which can be used by an 
% (steepest gradient) optimizer to register two images.
%
% The energy gradient on each pixel equals the the demon velocity field 
% which is  described by the paper of Thirion 1998 and extended 
% by Cachier 1999 and He Wang 2005. 
%  
% GradientPixels =-(ImageStatic-ImageMoving) * 
%    ((GradientImageStatic./(||GradientImageStatic||^2+alpha^2*(ImageStatic-ImageMoving)^2))+
%    (GradientImageMoving./(||GradientImageMoving||^2+alpha^2*(ImageStatic-ImageMoving)^2)));
%
% This gradient is gaussian smoothed by this function to get fluid regularization.
% 
% Also the energy value itself at a certain transformation is calculated
%    E =(ImageStatic-ImageMoving)^2 + 
%         alpha^2*(ImageStatic-ImageMoving)^2.*||UpdateTransformationField||^2;
% 
% This energy value equation is not described by Thirion but needed by the
% line search of the optimizer. The paper by Tom Vercauteren et Al.
% "Non-parametric Diffeomorphic Image ...", gives more insight on how the
% Thirions equation can be derived from it.
%
% The gradient of the transformation field it self can also be smoothed by
% this function to get diffusion regularization (not part of original 
% demon registration).
%
% [E,Egrad]=DEMONS_ENERGY(M,F,M_TF,F_TF,T,SIZES,OPTIONS)
%
% inputs,
%  M: 2D or 3D moving input image
%  F: 2D or 3D static image
%  M_TF: image M but then modality transformed into grey values of F. Only
%         needed in case of multiple modalities, otherwise set to []
%  F_TF: image F but then modality transformed into grey values of M. Only
%         needed in case of multiple modalities, otherwise set to []
%  T: Current transformation field as a long vector
%  SIZES: The dimensions of the transformation field variable
%  OPTIONS: The parameters.
%     OPTIONS.sigma_fluid : Sigma used to smooth gradient(/update) field
%                           (Fluid regularization)
%     OPTIONS.sigma_diff : Sigma used to smooth current transformation
%                           field (Diffusion regularization)
%     OPTIONS.alpha : Kind of noise suppression constant, see equations
%                      above.
%     OPTIONS.interpolation : linear (default) or cubic.
%
% Outputs,
%  E: Current registration energy (registration error)
%  Egrad: Energy gradient on all pixels
%
% note: This function works better with steepest gradient optimizers,
%       than other Matlab optimizers.
%
% Example 1,
%   % Read both images 
%    M=im2double(imread('lenag1.png'));
%    F=im2double(imread('lenag2.png'));
%   % Parameters
%    options.sigma_fluid=8; options.sigma_diff=0; options.alpha=1.5;
%   % Transformation field
%    T=zeros([size(M) 2]);
%   % Store the dimensions of transformation field, and make, a long vector from T
%    sizes=size(T); T=T(:);
%   % Start the demon energy registration optimizer
%    optim=struct('Display','iter','GradObj','on','MaxIter',500,'OutputFcn', @store_transf);
%    T=fminlbfgs(@(x)demons_energy(M,F,[],[],x,sizes,options),T,optim);
%   % Reshape O_trans from a vector to a matrix
%    T=reshape(T,sizes);
%   % If diffusion regularization is used: smooth transformation field
%    if(options.sigma_diff>0)
%      Hsmooth=fspecial('gaussian',[options.sigma_diff*6 options.sigma_diff*6],options.sigma_diff);
%      T(:,:,1)=imfilter(T(:,:,1),Hsmooth); T(:,:,2)=imfilter(T(:,:,2),Hsmooth);
%    end
%   % Transform the input image
%    Ms=movepixels(M,T(:,:,1),T(:,:,2));
%   % Show the registration results
%    figure,
%    subplot(1,3,1), imshow(M); title('input image 1');
%    subplot(1,3,2), imshow(F); title('input image 2');
%    subplot(1,3,3), imshow(Ms); title('transformed image 1');
%
% See also FMINSD, DEMONREGISTRATION
%
% Function is written by D.Kroon University of Twente (October 2008)

% Get transformation field of previous itteration
global last_transformation_field updatemovie;
if(isempty(last_transformation_field)), last_transformation_field=T; 
else
    if(nnz(size(last_transformation_field)-size(T))>0)
        last_transformation_field=T; 
    end
end

% Check/set input options
defaultoptions=struct('sigma_fluid',8, 'sigma_diff',0,'alpha',1.5,'interpolation','linear');
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('BsplineRegistrationGradient:unknownoption','unknown options found');
    end
end
% Set all parameters
alpha=options.alpha;
sigma_diff=options.sigma_diff;
sigma_fluid=options.sigma_fluid;
if(strcmpi(options.interpolation(1),'l')), interpolation_mode=0; else interpolation_mode=2; end 
  
if(ndims(M)==2)
    C=reshape(last_transformation_field,sizes);
    cx=C(:,:,1); cy=C(:,:,2);

    % Get the current transformation update
    T=reshape(T,sizes);
    ux=T(:,:,1)-C(:,:,1); uy=T(:,:,2)-C(:,:,2);

    if(sigma_diff>0)
        % Smooth transformation field
        sx=imgaussian(cx,sigma_diff);
        sy=imgaussian(cy,sigma_diff);
    else
        sx=cx; sy=cy;
    end

     % Calculate current moving image
    Ms=movepixels(M,sx+ux,sy+uy,[],interpolation_mode);
    if(~isempty(M_TF)), 
        Ms_TF=movepixels(M_TF,sx+ux,sy+uy,[],interpolation_mode); 
    end
    
    % If no multiple modalities, set F_TF to F and Ms_TF to Ms
    if(isempty(F_TF)), F_TF=F; Ms_TF=Ms; end
    
    % Calculate difference image
    Idiff_M=Ms-F_TF;
    Idiff_F=Ms_TF-F;
            

    % Calculate error after current transformation
    E_M=Idiff_M.^2+alpha^2*Idiff_M.^2.*(ux.^2+uy.^2);
    E_F=Idiff_F.^2+alpha^2*Idiff_F.^2.*(ux.^2+uy.^2);
    
    E=0.5*(sum(E_M(:))+sum(E_F(:)));
        
    % Make movie
    if(updatemovie),
        makemovie(movepixels(M_TF,sx,sy,[],interpolation_mode));
        updatemovie=false;
    end
  
    
    % If gradient needed, also determine the gradient.
    if ( nargout > 1 )
        % Calculate gradients in both images
        [My,Mx] = gradient(Ms);
        [Fy,Fx] = gradient(F);

        % Calculate gradient (update field) in x and y direction
        Ma=1./(eps+(Mx.^2+My.^2)+alpha^2*Idiff_M.^2);
        Fa=1./(eps+(Fx.^2+Fy.^2)+alpha^2*Idiff_F.^2);
        ux=Idiff_F.*(Fx.*Fa)+Idiff_M.*(Mx.*Ma); 
        uy=Idiff_F.*(Fy.*Fa)+Idiff_M.*(My.*Ma);

        % Smooth gradient incase of fluid-like regularization
        if(sigma_fluid>0)
            ux=imgaussian(ux,sigma_fluid);
            uy=imgaussian(uy,sigma_fluid);
        end

        % Both x and y gradient to one output variable
        Egrad=zeros(size(T));
        Egrad(:,:,1)=ux;
        Egrad(:,:,2)=uy;
        Egrad=Egrad(:);
    end
else
    C=reshape(last_transformation_field,sizes);
    cx=C(:,:,:,1); cy=C(:,:,:,2); cz=C(:,:,:,3);

    % Get the current transformation update
    T=reshape(T,sizes);
    ux=T(:,:,:,1)-C(:,:,:,1); 
    uy=T(:,:,:,2)-C(:,:,:,2);
    uz=T(:,:,:,3)-C(:,:,:,3);

    if(sigma_diff>0)
        sx=imgaussian(cx,sigma_diff);
        sy=imgaussian(cy,sigma_diff);
        sz=imgaussian(cz,sigma_diff);
    else
        sx=cx; sy=cy; sz=cz;
    end 

    % Calculate current moving image
    Ms=movepixels(M,sx+ux,sy+uy,sz+uz,interpolation_mode);
    if(~isempty(M_TF)), 
        Ms_TF=movepixels(M_TF,sx+ux,sy+uy,sz+uz,interpolation_mode); 
    end
    
    % If no multiple modalities, set F_TF to F and Ms_TF to Ms
    if(isempty(F_TF)), F_TF=F; Ms_TF=Ms; end
    
    % Calculate difference image
    Idiff_M=Ms-F_TF;
    Idiff_F=Ms_TF-F;
    
    % Calculate error after current transformation
    E_M=Idiff_M.^2+alpha^2*Idiff_M.^2.*(ux.^2+uy.^2+uz.^2);
    E_F=Idiff_F.^2+alpha^2*Idiff_F.^2.*(ux.^2+uy.^2+uz.^2);
    
    E=0.5*(sum(E_M(:))+sum(E_F(:)));
        

    % If gradient needed, also determine the gradient.
    if ( nargout > 1 )
        % Calculate gradients in both images
        [My,Mx,Mz] = gradient(Ms);
        [Fy,Fx,Fz] = gradient(F);

        % Calculate gradient (update field) in x and y direction
        Ma=1./(eps+(Mx.^2+My.^2+Mz.^2)+alpha^2*Idiff_M.^2);
        Fa=1./(eps+(Fx.^2+Fy.^2+Fz.^2)+alpha^2*Idiff_F.^2);
        ux=Idiff_F.*(Fx.*Fa)+Idiff_M.*(Mx.*Ma); 
        uy=Idiff_F.*(Fy.*Fa)+Idiff_M.*(My.*Ma);
        uz=Idiff_F.*(Fz.*Fa)+Idiff_M.*(Mz.*Ma);
        
        % Smooth gradient incase of fluid-like regularization
        if(sigma_fluid>0)
            ux=imgaussian(ux,sigma_fluid);
            uy=imgaussian(uy,sigma_fluid);
            uz=imgaussian(uz,sigma_fluid);
        end

        % Both x and y and z gradient to one output variable
        Egrad=zeros(size(T));
        Egrad(:,:,:,1)=ux;
        Egrad(:,:,:,2)=uy;
        Egrad(:,:,:,3)=uz;
        Egrad=Egrad(:);
    end
end
    
function makemovie(I)
global Mmovie Mindex;
if(isempty(Mindex)), Mindex=1; end
I=imresize(I,[256 256]);
if(isempty(Mmovie)), 
    Mmovie = im2frame(uint8(I*256),gray(256));
else
    Mmovie(Mindex) = im2frame(uint8(I*256),gray(256));
end
Mindex=Mindex+1;



