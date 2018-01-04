function [Ireg,Bx,By,Fx,Fy] = register_images(Imoving,Istatic,Options)
% This function register_images is the most easy way to register two
% images both affine and nonrigidly.
%
% Features:
% - It can be used with images from different type of scans or modalities.
% - It uses both a rigid transform and a nonrigid registration.
% - It uses multilevel refinement
% - It can be used with images of different sizes.
% - The function will automaticaly detect single modality or multiple
%   modalities, and choose the right registration method.
%
% [Ireg,Bx,By,Fx,Fy] = register_images(Imoving,Istatic,Options);
%
% Inputs,
%   Imoving : The image which will be registerd
%   Istatic : The image on which Imoving will be registered
%   Options : Registration options, see help below
%
% Outputs,
%   Ireg : The registered moving image
%   Bx, By : The backwards transformation fields of the pixels in 
%       x and y direction seen from the  static image to the moving image.
%   Fx, Fy : The (approximated) forward transformation fields of the pixels in 
%       x and y direction seen from the moving image to the static image.
%       (See the function backwards2forwards)
% Options,
%   Options.SigmaFluid : The sigma of the gaussian smoothing kernel of the pixel
%                   velocity field / update field, this is a form of fluid
%                   regularization, (default 4)
%   Options.SigmaDiff : The sigma for smoothing the transformation field
%                   is not part of the orignal demon registration, this is 
%                   a form of diffusion regularization, (default 1)
%   Options.Interpolation : Linear (default) or Cubic.
%   Options.Alpha : Constant which reduces the influence of edges (and noise)
%                   and limits the update speed (default 4). 
%   Options.Similarity : Choose 'p' for single modality and 'm' for
%                   images of different modalities. (default autodetect)
%   Options.Registration: Rigid, Affine, NonRigid  
%   Options.MaxRef : Maximum number of grid refinements steps.
%   Options.Verbose: Display Debug information 0,1 or 2 
%                   
% Notes,
%   In case of Multiple Modalities affine registration is done with mutual
%   information. The non-rigid registration is done by first doing a
%   modality transformation (paints regions in image 1 with the intensity
%   pallette of those regions in image2 and visa versa), and than 
%   using "normal" pixel based demon registration. See MutualTransform.m
% 
%
% Example,
%   % Read two greyscale images of Lena
%   Imoving=imread('images/lenag1.png'); 
%   Istatic=imread('images/lenag3.png');
% 
%   % Register the images
%   [Ireg,Bx,By,Fx,Fy] = register_images(Imoving,Istatic,struct('Similarity','p'));
%
%   % Show the registration result
%   figure,
%   subplot(2,2,1), imshow(Imoving); title('moving image');
%   subplot(2,2,2), imshow(Istatic); title('static image');
%   subplot(2,2,3), imshow(Ireg); title('registerd moving image');
%   % Show also the static image transformed to the moving image
%   Ireg2=movepixels(Istatic,Fx,Fy);
%   subplot(2,2,4), imshow(Ireg2); title('registerd static image');
%
%  % Show the transformation fields
%   figure,
%   subplot(2,2,1), imshow(Bx,[]); title('Backward Transf. in x direction');
%   subplot(2,2,2), imshow(Fx,[]); title('Forward Transf. in x direction');
%   subplot(2,2,3), imshow(By,[]); title('Backward Transf. in y direction');
%   subplot(2,2,4), imshow(Fy,[]); title('Forward Transf. in y direction');
%
% % Calculate strain tensors
%   E = strain(Fx,Fy);
% % Show the strain tensors
%   figure,
%   subplot(2,2,1), imshow(E(:,:,1,1),[]); title('Strain Tensors Exx');
%   subplot(2,2,2), imshow(E(:,:,1,2),[]); title('Strain Tensors Exy');
%   subplot(2,2,3), imshow(E(:,:,2,1),[]); title('Strain Tensors Eyx');
%   subplot(2,2,4), imshow(E(:,:,2,2),[]); title('Strain Tensors Eyy');
%
% Example Multi-Modalities
%   % Read two brain images 
%   Imoving=im2double(imread('images/brain_T1_wave.png')); 
%   Istatic=im2double(imread('images/brain_T2.png'));
%
%   % Register the images
%   [Ireg,Bx,By] = register_images(Imoving,Istatic,struct('SigmaFluid',4));
%
%   figure,
%   subplot(1,3,1), imshow(Imoving); title('moving image');
%   subplot(1,3,2), imshow(Istatic); title('static image');
%   subplot(1,3,3), imshow(Ireg); title('registerd moving image');
%
%   % Read normal T1 image and transformation field
%   Inormal=im2double(imread('images/brain_T1.png'));
%   load('images/wave_field.mat');
%
%   % Show the difference with ideal image
%   figure, imshow(Imoving-Inormal,[-0.5 0.5]); title('unregistered')
%   figure, imshow(Ireg-Inormal,[-0.5 0.5]); title('registered');
%   disp(['pixel abs difference : ' num2str(sum(abs(Imoving(:)-Inormal(:))))])
%   disp(['pixel abs difference : ' num2str(sum(abs(Imoving(:)-Ireg(:))))])
%
%   % Show Warp field
%   figure,
%   subplot(2,2,1), imshow(BxNormal,[-20 20]); title('Bx Normal');
%   subplot(2,2,2), imshow(Bx,[-20 20]); title('Bx');
%   subplot(2,2,3), imshow(ByNormal,[-20 20]); title('By Normal');
%   subplot(2,2,4), imshow(By,[-20 20]); title('By');


% Function is written by D.Kroon University of Twente (March 2009)


% add all needed function paths
try
    functionname='register_images.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions'])
    addpath([functiondir '/functions_affine'])
	addpath([functiondir '/functions_nonrigid'])
catch me
    disp(me.message);
end

% Disable warning
warning('off', 'MATLAB:maxNumCompThreads:Deprecated')

% Process inputs
defaultoptions=struct('Similarity',[],'Registration','NonRigid','MaxRef',[],'Verbose',2,'SigmaFluid',4,'Alpha',4,'SigmaDiff',1,'Interpolation','Linear');
if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})),  Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('register_images:unknownoption','unknown options found');
    end
end

% Set parameters
MaxRef=Options.MaxRef;

% Start time measurement
if(Options.Verbose>0), tic; end

% Store the class of the inputs
Iclass=class(Imoving);

% Convert the inputs to double
Imoving=im2double(Imoving);
Istatic=im2double(Istatic);

% Resize the moving image to fit the static image
if(sum(size(Istatic)-size(Imoving))~=0)
    Imoving = imresize(Imoving,size(Istatic),'bicubic');
end

% Make smooth images for histogram and fast affine registration
ISmoving=imgaussian(Imoving,2.5,[10 10]);
ISstatic=imgaussian(Istatic,2.5,[10 10]);

% Detect if the mutual information or pixel distance can be used as 
% similarity measure. By comparing the histograms.
if(isempty(Options.Similarity))
    Hmoving= hist(ISmoving(:),60)./numel(Imoving);
    Hstatic = hist(ISstatic(:),60)./numel(Istatic);
    Hmoving(1)=0; Hstatic(1)=0;
    if(sum(log(abs(Hmoving-Hstatic)+1))>0.3), 
        Options.Similarity='m'; 
        if(Options.Verbose>0), disp('Multi Modalities, Mutual information is used'); drawnow; end
    else
        Options.Similarity='p';
        if(Options.Verbose>0), disp('Same Modalities, Pixel Distance is used'); drawnow; end
    end
end
if(Options.Similarity(1)=='p'), type_affine='sd'; else type_affine='mi'; end
% Register the moving image affine to the static image

% Affine register the smoothed images to get the registration parameters
if(strcmpi(Options.Registration(1),'R'))
    if(Options.Verbose>0), disp('Start Rigid registration'); drawnow; end
    % Parameter scaling of the Translation and Rotation
    scale=[1 1 1];
    % Set initial affine parameters
    x=[0 0 0];
elseif(strcmpi(Options.Registration(1),'A'))
    if(Options.Verbose>0), disp('Start Affine registration'); drawnow; end
    % Parameter scaling of the Translation, Rotation, Resize and Shear
    scale=[1 1 1 0.01 0.01 1e-4 1e-4];
    % Set initial affine parameters
    x=[0 0 0 100 100 0 0];
elseif(strcmpi(Options.Registration(1),'N'))
    if(Options.Verbose>0), disp('Start Affine part of Non-Rigid registration'); drawnow; end
    % Parameter scaling of the Translation, Rotation, Resize and Shear
    scale=[1 1 1 0.01 0.01 1e-4 1e-4];
    % Set initial affine parameters
    x=[0 0 0 100 100 0 0];
else
     warning('register_images:unknownoption','unknown registration method');    
end

for refine_itt=1:2
    if(refine_itt==2)
        ISmoving=Imoving; ISstatic=Istatic;
    end
	% Use struct because expanded optimset is part of the Optimization Toolbox.
    optim=struct('GradObj','off','GoalsExactAchieve',1,'Display','off','MaxIter',100,'MaxFunEvals',1000,'TolFun',1e-14,'DiffMinChange',1e-6);
    if(Options.Verbose>0), optim.Display='iter'; end
	x=fminlbfgs(@(x)affine_registration_error(x,scale,ISmoving,ISstatic,type_affine),x,optim);            
end

% Scale the translation, resize and rotation parameters to the real values
x=x.*scale;

if(strcmpi(Options.Registration(1),'R'))
    % Make the rigid transformation matrix
    M=make_transformation_matrix(x(1:2),x(3));
else
    % Make the affine transformation matrix
    M=make_transformation_matrix(x(1:2),x(3),x(4:5));    
end


% Make center of the image transformation coordinates 0,0
[x,y]=ndgrid(0:(size(Imoving,1)-1),0:(size(Imoving,2)-1));
xd=x-(size(Imoving,1)/2); yd=y-(size(Imoving,2)/2);

% Calculate the backwards transformation fields
Bx = ((size(Imoving,1)/2) + M(1,1) * xd + M(1,2) *yd + M(1,3) * 1)-x;
By = ((size(Imoving,2)/2) + M(2,1) * xd + M(2,2) *yd + M(2,3) * 1)-y;

% Initialize the modality transformed image variables
M_TF=[]; F_TF=[];
    
% The nonrigid part of the registration
if(strcmpi(Options.Registration(1),'N'))
    
    % Demon registration parameters
    refinements=floor(log2(min(size(Imoving))/16));
    if(refinements>MaxRef), refinements=MaxRef; end
    parameters.sigma_diff=Options.SigmaDiff;
     
    % Non-rigid registration
    if(Options.Verbose>0), disp('Start non-rigid demon registration'); drawnow; end

    % Do every refinements step twice if modality transformation enabled
    if(Options.Similarity(1)=='m'),  loop=2; else loop=1; end
    
    % Loop trough all refinements steps.
    for j=0:refinements
        for l=1:loop
            % Set scaling parameters.resizepercentageentage
            resizepercentage=1/2^(refinements-j);
            if(resizepercentage>1), resizepercentage=1; end

            parameters.alpha=Options.Alpha*sqrt(resizepercentage);
            parameters.sigma_fluid=Options.SigmaFluid;

            if(Options.Verbose>0), disp(['Scaling resizepercentageentage : ' num2str(resizepercentage)]), end

            % Incase of multiple modalities, transform both images to their
            % opposite modalities.
            if(Options.Similarity(1)=='m')
                if(Options.Verbose>0), disp('Start modality transformation'); drawnow; end
                Bx_large=imresize(Bx,size(Imoving),'bicubic')*(size(Imoving,1)/size(Bx,1));
                By_large=imresize(By,size(Imoving),'bicubic')*(size(Imoving,2)/size(By,2));
                [Imoving_TF,Istatic_TF]=MutualTransform(Imoving,Istatic,15*sqrt(1/resizepercentage),4,Bx_large,By_large);
                if(Options.Verbose>0), disp('Finished modality transformation'); drawnow; end
            end

            sigma = 0.3/resizepercentage;
            % Set and resize the moving image and static image
            M=imresize(imgaussian(Imoving,sigma,[sigma*6 sigma*6]),resizepercentage,'bicubic'); 
            F=imresize(imgaussian(Istatic,sigma,[sigma*6 sigma*6]),resizepercentage,'bicubic');
            
            % Resize the modality transformed images
            if(Options.Similarity(1)=='m')
                M_TF=imresize(imgaussian(Imoving_TF,sigma,[sigma*6 sigma*6]),resizepercentage,'bicubic'); 
                F_TF=imresize(imgaussian(Istatic_TF,sigma,[sigma*6 sigma*6]),resizepercentage,'bicubic');
            end

            % Resize the transformation field to current image size
            Bx=imresize(Bx,size(M),'bicubic')*(size(M,1)/size(Bx,1));
            By=imresize(By,size(M),'bicubic')*(size(M,2)/size(By,2));

            % Put transformation fields in x and y direction in one variable
            B=zeros([size(M) 2]); B(:,:,1)=Bx; B(:,:,2)=By;
            % Store the dimensions of transformation field, and make a long vector from T
            sizes=size(B); B=B(:);

            % Parameters
            options.sigma_fluid=parameters.sigma_fluid;
            options.sigma_diff=parameters.sigma_diff;
            options.alpha=parameters.alpha;
            options.interpolation=Options.Interpolation;
            
            % Optimizer parameters
            optim=struct('Display','off','StoreN',10,'GoalsExactAchieve',0,'HessUpdate','lbfgs','GradObj','on','OutputFcn', @store_transf,'MaxIter',200,'TolFun',1e-14,'DiffMinChange',1e-5);
            if(l==loop), 
                optim.TolX = 0.02; 
            else
                optim.TolX = 0.1; 
            end
            if(Options.Verbose>1), optim.Display='iter'; end

            % Start the demon energy registration optimizer
            B=fminlbfgs(@(x)demons_energy(M,F,M_TF,F_TF,x,sizes,options),B,optim);

            % Reshape B from a vector to an x and y transformation field
            B=reshape(B,sizes);
            Bx=B(:,:,1); By=B(:,:,2);
        end
    end

    % Scale everything back if not already
    if(resizepercentage~=1)
        Bx=imresize(Bx,size(Imoving),'bicubic')*(size(Imoving,1)/size(Bx,1));
        By=imresize(By,size(Imoving),'bicubic')*(size(Imoving,2)/size(By,2));
    end
end
    
% Transform the input image
Ireg=movepixels(Imoving,Bx,By,[], 3);

if ( nargout>3 )
    % Make the forward transformation fields from the backwards
    [Fx,Fy]=backwards2forwards(Bx,By);
end

% Set the class of output to input class
if(strcmpi(Iclass,'uint8')), Ireg=im2uint8(Ireg); end
if(strcmpi(Iclass,'uint16')), Ireg=im2uint16(Ireg); end
if(strcmpi(Iclass,'uint32')), Ireg=im2uint32(Ireg); end
if(strcmpi(Iclass,'int8')), Ireg=im2int8(Ireg); end
if(strcmpi(Iclass,'int16')), Ireg=im2int16(Ireg); end
if(strcmpi(Iclass,'int32')), Ireg=im2int32(Ireg); end
if(strcmpi(Iclass,'single')), Ireg=im2single(Ireg); end

% End time measurement
if(Options.Verbose>0), toc, end

