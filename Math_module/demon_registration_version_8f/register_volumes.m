function [Ireg,Bx,By,Bz,Fx,Fy,Fz] = register_volumes(Imoving,Istatic,Options)
% This function register_volumes is the most easy way to register two
% 3D images both affine and nonrigidly.
%
% Features:
% - It can be used with images from different type of scans or modalities.
% - It uses both a rigid transform and a nonrigid registration.
% - It uses multilevel refinement
% - It can be used with images of different sizes.
% - The function will automaticaly detect single modality or multiple
%   modalities, and choose the right registration method.
%
% [Ireg,Bx,By,Bz,Fx,Fy,Fz] = register_volumes(Imoving,Istatic,Options);
%
% Inputs,
%   Imoving : The 3D image which will be registerd
%   Istatic : The 3D  image on which Imoving will be registered
%   Options : Registration options, see help below
%
% Outputs,
%   Ireg : The registered moving image
%   Bx, By, Bz : The backwards transformation fields of the pixels in 
%       x,y and z direction seen from the  static image to the moving image.
%   Fx, Fy, Fz : The (approximated) forward transformation fields of the pixels in 
%       x, y and z direction seen from the moving image to the static image.
%       (See the function backwards2forwards)
%
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
%   % add needed function paths
%   functiondir=which('register_volumes.m');
%   addpath([functiondir(1:end-length('register_volumes.m')) '/images'])
%
%   % Get the volume data
%   [Imoving,Istatic]=get_example_data;
%
%   % Register the images
%   Ireg = register_volumes(Imoving,Istatic);
%
%   % Show the results
%   showcs3(Imoving);
%   showcs3(Istatic);
%   showcs3(Ireg);
%
% Function is written by D.Kroon University of Twente (March 2009)


% add all needed function paths
% try
%     functionname='register_images.m';
%     functiondir=which(functionname);
%     functiondir=functiondir(1:end-length(functionname));
%     addpath([functiondir '/functions'])
%     addpath([functiondir '/functions_affine'])
% 	addpath([functiondir '/functions_nonrigid'])
% catch me
%     disp(me.message);
% end

% Disable warning
warning('off', 'MATLAB:maxNumCompThreads:Deprecated')

% Process inputs
defaultoptions=struct('Similarity',[],'Registration','NonRigid','MaxRef',[],'Verbose',2,'SigmaFluid',2,'Alpha',2,'SigmaDiff',1,'Interpolation','Linear');
% defaultoptions=struct('Similarity',[],'Registration','Rigid','MaxRef',[],'Verbose',2,'SigmaFluid',4,'Alpha',4,'SigmaDiff',1,'Interpolation','Linear');

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

% Convert input uint8, uint32 etc. to single.
if(~strcmpi(Iclass,'single')&&~strcmpi(Iclass,'double'))
    range=getrangefromclass(Imoving);
    Imoving=single(Imoving)./range(2);
    range=getrangefromclass(Istatic);
    Istatic=single(Istatic)./range(2);
end
% Do the inverse registration if the resolution of the moving image is
% greater thant the static one! --> add by B. Lemasson
% if size(Imoving,1) > size(Istatic,1) ||...
%     size(Imoving,2) > size(Istatic,2)
%     old_Imoving = Imoving;    old_Istatic = Istatic;
%     Imoving = old_Istatic;
%     Istatic = old_Imoving;
%     inverse_registration = 1;
% else
     inverse_registration = 0;
% end
% Resize the moving image to fit the static image
if(sum(size(Istatic)-size(Imoving))~=0)
    Imoving=imresize3d(Imoving,[],size(Istatic),'cubic');
end

% Smooth both images for faster registration
ISmoving=imgaussian(Imoving,2.5,[10 10 10]);
ISstatic=imgaussian(Istatic,2.5,[10 10 10]);

% Detect if the mutual information or pixel distance can be used as 
% similarity measure. By comparing the histograms.
if(isempty(Options.Similarity))
    Hmoving= hist(ISmoving(:),60)./numel(Imoving);
    Hstatic = hist(ISstatic(:),60)./numel(Istatic);
    if(sum(log(abs(Hmoving-Hstatic)+1))>0.5), 
        Options.Similarity='m'; 
        if(Options.Verbose>0), disp('Multi Modalities, Mutual information is used'); drawnow; end
    else
        Options.Similarity='p';
        if(Options.Verbose>0), disp('Same Modalities, Pixel Distance is used'); drawnow; end
    end
end
%%%%test forcing Same Modalities!!
 Options.Similarity='p';
        if(Options.Verbose>0), disp('Same Modalities, Pixel Distance is used'); drawnow; end
%%%%

if(Options.Similarity(1)=='p'), type_affine='sd'; else type_affine='mi'; end

% Register the moving image affine to the static image

% Affine register the smoothed images to get the registration parameters
if(strcmpi(Options.Registration(1),'R'))
    if(Options.Verbose>0), disp('Start Rigid registration'); drawnow; end
    % Parameter scaling of the Translation and Rotation
    scale=[1 1 1 1 1 1];
    % Set initial rigid parameters
    x=[0 0 0 0 0 0];
elseif(strcmpi(Options.Registration(1),'A'))
    if(Options.Verbose>0), disp('Start Affine registration'); drawnow; end
    % Parameter scaling of the Translation, Rotation, Resize and Shear
    scale=[1 1 1 1 1 1 0.01 0.01 0.01 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4];
    % Set initial rigid parameters
    x=[0 0 0 0 0 0 100 100 100 0 0 0 0 0 0];
elseif(strcmpi(Options.Registration(1),'N'))
    if(Options.Verbose>0), disp('Start Affine part of Non-Rigid registration'); drawnow; end
    % Parameter scaling of the Translation, Rotation, Resize and Shear
    scale=[1 1 1 1 1 1 0.01 0.01 0.01 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4];
    % Set initial rigid parameters
    x=[0 0 0 0 0 0 100 100 100 0 0 0 0 0 0];
else
     warning('register_volumes:unknownoption','unknown registration method');    
end

% try to register the unsmothed data (directly)!!
for refine_itt=2:2
    if(refine_itt==2)
        ISmoving=Imoving; ISstatic=Istatic;
    end
	% Use struct because expanded optimset is part of the Optimization Toolbox.
    optim=struct('GradObj','off','GoalsExactAchieve',1,'Display','off','MaxIter',500,'MaxFunEvals',5000,'TolFun',1e-14,'DiffMinChange',1e-6, 'TolX',1e-6);
    if(Options.Verbose>0), optim.Display='iter'; end
	x=fminlbfgs(@(x)affine_registration_error(x,scale,ISmoving,ISstatic,type_affine),x,optim);            
end

%%% Modified by B. Lemasson
% 
% % first do a rigid registation with ISmoving and ISstatic data
% % Use struct because expanded optimset is part of the Optimization Toolbox.
% optim=struct('GradObj','off','GoalsExactAchieve',1,'Display','off','MaxIter',100,'MaxFunEvals',1000,'TolFun',1e-14,'DiffMinChange',1e-7, 'TolX',1e-7);
% if(Options.Verbose>0), optim.Display='iter'; end
% x=fminlbfgs(@(x)affine_registration_error(x,scale,ISmoving,ISstatic,type_affine),x,optim);
% % the same with the real data
% ISmoving=Imoving; ISstatic=Istatic;
% x=fminlbfgs(@(x)affine_registration_error(x,scale,ISmoving,ISstatic,type_affine),x,optim);
% 
% 
% Options.Registration = 'NonRigid';
% x(7:15) = [100 100 100 0 0 0 0 0 0];
% scale(7:15)=[0.01 0.01 0.01 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4];
% optim=struct('GradObj','off','GoalsExactAchieve',1,'Display','off','MaxIter',100,'MaxFunEvals',1000,'TolFun',1e-14,'DiffMinChange',1e-6, 'TolX',1e-6);
% if(Options.Verbose>0), optim.Display='iter'; end
% x=fminlbfgs(@(x)affine_registration_error(x,scale,ISmoving,ISstatic,type_affine),x,optim);



% Scale the translation, resize and rotation parameters to the real values
x=x.*scale;

if(strcmpi(Options.Registration(1),'R'))
    % Make the rigid transformation matrix
    M=make_transformation_matrix(x(1:3),x(4:6));
else
    % Make the affine transformation matrix
    M=make_transformation_matrix(x(1:3),x(4:6),x(7:9),x(10:15));    
end


% Make center of the image transformation coordinates 0,0
[x,y,z]=ndgrid(0:(size(Imoving,1)-1),0:(size(Imoving,2)-1),0:(size(Imoving,3)-1));
xd=x-(size(Imoving,1)/2); yd=y-(size(Imoving,2)/2);  zd=z-(size(Imoving,3)/2);

% Calculate the backwards transformation fields
Bx = ((size(Imoving,1)/2) + M(1,1) * xd + M(1,2) *yd + M(1,3) *zd + M(1,4)* 1)-x;
By = ((size(Imoving,2)/2) + M(2,1) * xd + M(2,2) *yd + M(2,3) *zd + M(2,4)* 1)-y;
Bz = ((size(Imoving,3)/2) + M(3,1) * xd + M(3,2) *yd + M(3,3) *zd + M(3,4)* 1)-z;

% display result
% Ireg=movepixels(Imoving,Bx,By,Bz,3);
% figure;imshow3D(Iregtest);figure;imshow3D(Istatic);

% Initialize the modality transformed image variables
% M_TF=[]; F_TF=[];
%     
% % The nonrigid part of the registration
% if(strcmpi(Options.Registration(1),'N'))
%     
%     % Demon registration parameters
%     refinements=floor(log2(min(size(Imoving))/16));
%     if(refinements>MaxRef), refinements=MaxRef; end
%     parameters.sigma_diff=Options.SigmaDiff;
%      
%     % Non-rigid registration
%     if(Options.Verbose>0), disp('Start non-rigid demon registration'); drawnow; end
% 
%     % Do every refinements step twice if modality transformation enabled
%     if(Options.Similarity(1)=='m'),  loop=2; else loop=1; end
%     
%     % add by BL
%     if refinements < 0, refinements = 0; end
% %     refinements = 0; loop = 1;
%     % Loop trough all refinements steps.
%     for j=0:refinements
%         for l=1:loop
%             % Set scaling parameters.resizepercentageentage
%             resizepercentage=1/2^(refinements-j);
%             if(resizepercentage>1), resizepercentage=1; end
% 
%             parameters.alpha=Options.Alpha*sqrt(resizepercentage);
%             parameters.sigma_fluid=Options.SigmaFluid;
% 
%             if(Options.Verbose>0), disp(['Scaling resizepercentageentage : ' num2str(resizepercentage)]), end
% 
%             % Incase of multiple modalities, transform both images to their
%             % opposite modalities.
%             if(Options.Similarity(1)=='m')
%                 if(Options.Verbose>0), disp('Start modality transformation'); drawnow; end
% 				Bx_large=imresize3d(Bx,[],size(Imoving),'cubic')*(size(Imoving,1)/size(Bx,1));
% 				By_large=imresize3d(By,[],size(Imoving),'cubic')*(size(Imoving,2)/size(By,2));
% 				Bz_large=imresize3d(Bz,[],size(Imoving),'cubic')*(size(Imoving,3)/size(Bz,3));
%                 [Imoving_TF,Istatic_TF]=MutualTransform(Imoving,Istatic,15*sqrt(1/resizepercentage),4,Bx_large,By_large,Bz_large);
%                 if(Options.Verbose>0), disp('Finished modality transformation'); drawnow; end
%             end
% 
%             sigma = 0.3/resizepercentage;
%             
% 			% Set and resize the moving image and static image
%             M=imresize3d(imgaussian(Imoving,sigma,round([sigma*6 sigma*6 sigma*6])),resizepercentage,[],'linear'); 
%             F=imresize3d(imgaussian(Istatic,sigma,round([sigma*6 sigma*6 sigma*6])),resizepercentage,[],'linear');
% 
%             % Resize the modality transformed images
%             if(Options.Similarity(1)=='m')
% 				M_TF=imresize3d(imgaussian(Imoving_TF,sigma,round([sigma*6 sigma*6 sigma*6])),resizepercentage,[],'linear'); 
% 				F_TF=imresize3d(imgaussian(Istatic_TF,sigma,round([sigma*6 sigma*6 sigma*6])),resizepercentage,[],'linear');
%             end
% 
%             % Resize the transformation field to current image size
%             Bx=imresize3d(Bx,[],size(M),'cubic')*(size(M,1)/size(Bx,1));
% 			By=imresize3d(By,[],size(M),'cubic')*(size(M,2)/size(By,2));
% 			Bz=imresize3d(Bz,[],size(M),'cubic')*(size(M,3)/size(Bz,3));
% 		
% 		    % Put transformation fields in x and y direction in one variable
% 			B=zeros([size(M) 3]); B(:,:,:,1)=Bx; B(:,:,:,2)=By; B(:,:,:,3)=Bz;
% 			% Store the dimensions of transformation field, and make a long vector from T
% 			sizes=size(B); B=B(:);
% 			
%             % Parameters
%             options.sigma_fluid=parameters.sigma_fluid;
%             options.sigma_diff=parameters.sigma_diff;
%             options.alpha=parameters.alpha;
%             options.interpolation=Options.Interpolation;
% 
%             % Optimizer parameters
%             optim=struct('Display','off','StoreN',10,'GoalsExactAchieve',0,'HessUpdate','lbfgs','GradObj','on','OutputFcn', @store_transf,'MaxIter',200,'TolFun',1e-14,'DiffMinChange',1e-5);
%             if(l==loop), 
%                 optim.TolX = 0.02; 
%             else
%                 optim.TolX = 0.02;%0.1 
%             end
%             if(Options.Verbose>1), optim.Display='iter'; end
% 
%            % Start the demon energy registration optimizer
%             B=fminlbfgs(@(x)demons_energy(M,F,M_TF,F_TF,x,sizes,options),B,optim);
% 
%             % Reshape B from a vector to an x and y transformation field
%             B=reshape(B,sizes);
%  			Bx=B(:,:,:,1); By=B(:,:,:,2); Bz=B(:,:,:,3); 
% 
%         end
%     end
% 
%     % Scale everything back if not already
%     if(resizepercentage~=1)
% 	    Bx=imresize3d(Bx,[],size(Imoving),'cubic')*(size(Imoving,1)/size(Bx,1));
%         By=imresize3d(By,[],size(Imoving),'cubic')*(size(Imoving,2)/size(By,2));
%         Bz=imresize3d(Bz,[],size(Imoving),'cubic')*(size(Imoving,3)/size(Bz,3));
%     end
% end
    
% Transform the input image

% Do the inverse registration if the resolution of the moving image is
% greater thant the static one! cf above --> add by B. Lemasson
if inverse_registration == 1
	% Make the forward transformation fields from the backwards with the
	% same geometry of the static image
    Bx=imresize3d(Bx,[],size(old_Istatic),'cubic')*(size(old_Istatic,1)/size(Bx,1));
    By=imresize3d(By,[],size(old_Istatic),'cubic')*(size(old_Istatic,2)/size(By,2));
    Bz=imresize3d(Bz,[],size(old_Istatic),'cubic')*(size(old_Istatic,3)/size(Bz,3));
	[Fx,Fy,Fz]=backwards2forwards(Bx,By,Bz);
    
    old_Imoving_resized=imresize3d(old_Imoving,[],size(old_Istatic),'cubic');
    Ireg=movepixels(old_Imoving_resized,Fx,Fy,Fz,3);
else
    Ireg=movepixels(Imoving,Bx,By,Bz,3); 
    % no needed and too long
%     [Fx,Fy,Fz]=backwards2forwards(Bx,By,Bz);
    Fx = [];
    Fy = [];
    Fz = [];
end


% Set the class of output to input class
if(strcmpi(Iclass,'uint8')), Ireg=uint8(Ireg*((2^8)-1)); end
if(strcmpi(Iclass,'uint16')), Ireg=uint16(Ireg*((2^16)-1)); end
if(strcmpi(Iclass,'uint32')), Ireg=uint32(Ireg*((2^32)-1)); end
if(strcmpi(Iclass,'int8')), Ireg=int8(Ireg*((2^7)-1)); end
if(strcmpi(Iclass,'int16')), Ireg=int16(Ireg*((2^15)-1)); end
if(strcmpi(Iclass,'int32')), Ireg=int32(Ireg*((2^31)-1)); end

% End time measurement
if(Options.Verbose>0), toc, end

