function [files_in,files_out,opt] = Module_ASL_InvEff(files_in,files_out,opt)
% This is a template file for "brick" functions in NIAK.
%
% SYNTAX:
% [IN,OUT,OPT] = PSOM_TEMPLATE_BRICK(IN,OUT,OPT)
%
% _________________________________________________________________________
% INPUTS:
%
% IN        
%   (string) a file name of a 3D+t fMRI dataset .
%
% OUT
%   (structure) with the following fields:
%       flag_test
%   CORRECTED_DATA
%       (string, default <BASE NAME FMRI>_c.<EXT>) File name for processed 
%       data.
%       If OUT is an empty string, the name of the outputs will be 
%       the same as the inputs, with a '_c' suffix added at the end.
%
%   MASK
%       (string, default <BASE NAME FMRI>_mask.<EXT>) File name for a mask 
%       of the data. If OUT is an empty string, the name of the 
%       outputs will be the same as the inputs, with a '_mask' suffix added 
%       at the end.
%
% OPT           
%   (structure) with the following fields.  
%
%   TYPE_CORRECTION       
%      (string, default 'mean_var') possible values :
%      'none' : no correction at all                       
%      'mean' : correction to zero mean.
%      'mean_var' : correction to zero mean and unit variance
%      'mean_var2' : same as 'mean_var' but slower, yet does not use as 
%      much memory).
%
%   FOLDER_OUT 
%      (string, default: path of IN) If present, all default outputs 
%      will be created in the folder FOLDER_OUT. The folder needs to be 
%      created beforehand.
%
%   FLAG_VERBOSE 
%      (boolean, default 1) if the flag is 1, then the function prints 
%      some infos during the processing.
%
%   FLAG_TEST 
%      (boolean, default 0) if FLAG_TEST equals 1, the brick does not do 
%      anything but update the default values in IN, OUT and OPT.
%           
% _________________________________________________________________________
% OUTPUTS:
%
% IN, OUT, OPT: same as inputs but updated with default values.
%              
% _________________________________________________________________________
% SEE ALSO:
% NIAK_CORRECT_MEAN_VAR
%
% _________________________________________________________________________
% COMMENTS:
%
% _________________________________________________________________________
% Copyright (c) <NAME>, <INSTITUTION>, <START DATE>-<END DATE>.
% Maintainer : <EMAIL ADDRESS>
% See licensing information in the code.
% Keywords : PSOM, documentation, template, brick

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)
  
     % define every option needed to run this module
      %%   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'output_filename_ext','ASL_InvEff'};
    module_option(:,4)   = {'OutputSequenceName','AllName'};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
% list of everything displayed to the user associated to their 'type'
     % --> user_parameter(1,:) = user_parameter_list
     % --> user_parameter(2,:) = user_parameter_type
     % --> user_parameter(3,:) = parameter_default
     % --> user_parameter(4,:) = psom_parameter_list
     % --> user_parameter(5,:) = Scans_Input_DOF (degree-of-freedom)
     % --> user_parameter(6,:) = Help : text data which describe the parameter (it
     % will be display to help the user)
     user_parameter(:,1)   = {'Description','Text','','','','',...
        {'ASL_InvEff Description'}'};
user_parameter(:,2)   = {'Select a Pcasl scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
user_parameter(:,3)   = {'Parameters','','','','', '',''};
user_parameter(:,4)   = {'   .Output filename extension','char','ASL_InvEff','output_filename_ext','','',...
    {'Specify the string to be added to the filename input.'
    'Default filename extension is ''ASL_InvEff''.'}'};


VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'ASL_InvEff'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%
opt.NameOutFiles = {opt.output_filename_ext};




if isempty(files_out)
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('ASL_InvEff:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
if ~ischar(files_in.In1{1}) 
    error('files in should be a char');
end

[path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
if ~strcmp(ext_nii, '.nii')
     error('First file need to be a .nii');  
end


if isfield(opt,'threshold') && (~isnumeric(opt.threshold))
    opt.threshold = str2double(opt.threshold);
    if isnan(opt.threshold)
        disp('The threshold used was not a number')
        return
    end
end


%% Building default output names
if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
    opt.folder_out = path_nii;    
end

if isempty(files_out)
   files_out.In1 = {cat(2,opt.folder_out,filesep,name_nii,'_',opt.output_filename_ext,ext_nii)};
end

%% If the test flag is true, stop here !
if opt.flag_test == 1
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load input Nii file
N = niftiread(files_in.In1{1});
% load nifti info
info = niftiinfo(files_in.In1{1});

%% load input JSON file
J = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));



%ASL=input_image
ASL = N;


% DataControl = ASL(:,:,:,2,:);
% DataLabel = ASL(:,:,:,1,:);
% alpha = abs( (DataControl - DataLabel) ./ (2 * DataControl) )*100;
% 
% ASL_InvEff = [];
% ASL_InvEff(:,:,:,1,:) = alpha;


DataControl = ASL(:,:,2,:,:);
DataLabel = ASL(:,:,1,:,:);
alpha = abs( (DataControl - DataLabel) ./ (2 * DataControl) )*100;

% I need to squeeze the data to be able to display it in MIA, but I believe
% it might be a fraud.
ASL_InvEff = squeeze(alpha);


%% Initial function
% if(isempty(regexpi(ASL.acq.ppl_name,'\w*casl\w*')) ||...
%         (isempty(regexpi(ASL.acq.ppl_name,'\w*gefc\w*'))&&isempty(regexpi(ASL.acq.ppl_name,'\w*FcFlash\w*'))) || ...
%         isreal(ASL.reco.data))
%     msgbox(sprintf('%s \n \t - aquisition:\t %s \n \t - labeling:\t %s \n \t - reco:\t %s','Scan must combine:', 'GEFC','CASL or pCASL','complex'),'Wrong scan!', 'ASL_InvEff map warning');
%     ASL_InvEff = [];
%     return
% else
%     ASL_InvEff=ASL;
%    
%     % calcul de l'efficacite d'inversion PV6.0 (parties reelles et
%     % imaginaires sur images separees)
%     
%     
%     DataControl = ASL.reco.data(:,:,2,:,:);
%     DataLabel = ASL.reco.data(:,:,1,:,:);
%     alpha = abs( (DataControl - DataLabel) ./ (2 * DataControl) )*100;
%     
%     
%     % Mise en forme donnees
%     ASL_InvEff.reco.data = [];
%     ASL_InvEff.reco.data(:,:,1,:,:)=alpha;
%     
%     ParamConfig=sprintf('##$QuantifMethod=ASL_InvEff\n##$ASL=%s\n##$ASL scan info\n%s\n##END=',...
%         ASL_filename,[ASL.reco.iminfos{:}]);
%     
%     ASL_InvEff.reco.paramQuantif = ParamConfig;
%     ASL_InvEff.reco.texte = 'ASL_InvEff';
%     ASL_InvEff.reco.date = date;
%     ASL_InvEff.reco.no_echoes = 1;
%     ASL_InvEff.reco.globalmin=min(ASL_InvEff.reco.data(:));
%     ASL_InvEff.reco.globalmax=max(ASL_InvEff.reco.data(:));
%     ASL_InvEff.reco.fov_offsets = zeros([3,ASL_InvEff.reco.no_echoes,ASL_InvEff.reco.no_slices,ASL_InvEff.reco.no_expts]);
%     ASL_InvEff.reco.fov_orientation = zeros([9,ASL_InvEff.reco.no_echoes,ASL_InvEff.reco.no_slices,ASL_InvEff.reco.no_expts]);
%     ASL_InvEff.reco.label = {''};
%     ASL_InvEff.reco.phaselabel  = {''};
%     ASL_InvEff.reco.fov_phase_orientation  = zeros([ASL_InvEff.reco.no_echoes,ASL_InvEff.reco.no_slices,ASL_InvEff.reco.no_expts]);
%     ASL_InvEff.reco.scaling_factor = [];
%     ASL_InvEff.reco.scaling_offset = [];
%     for m_expt=1:ASL_InvEff.reco.no_expts,
%         for m_slice=1:ASL_InvEff.reco.no_slices,
%             for m_echo=1:ASL_InvEff.reco.no_echoes
%                 ASL_InvEff.reco.fov_offsets(:,m_echo,m_slice,m_expt) = ASL.reco.fov_offsets(:,1,m_slice,m_expt);
%                 ASL_InvEff.reco.fov_orientation(:,m_echo,m_slice,m_expt) = ASL.reco.fov_orientation(:,1,m_slice,m_expt);
%                 ASL_InvEff.reco.label(m_echo,m_slice,m_expt) = ASL.reco.label(1,m_slice,m_expt);
%                 ASL_InvEff.reco.phaselabel(m_echo,m_slice,m_expt) = ASL.reco.phaselabel(1,m_slice,m_expt);
%                 ASL_InvEff.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = ASL.reco.fov_phase_orientation(1,m_slice,m_expt);
%                 ASL_InvEff.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
%                 ASL_InvEff.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
%             end
%         end
%     end
%     if isfield(ASL_InvEff.reco, 'fov')
%         ASL_InvEff.reco.displayedecho=ASL_InvEff.reco.fov;
%     end
%     if isfield(ASL_InvEff.reco, 'reco_number')
%         ASL_InvEff.reco.reco_number=ASL_InvEff.reco.reco_number;
%     else
%         ASL_InvEff.reco.reco_number=1;
%     end
%     if isfield(ASL_InvEff.reco, 'scan_number')
%         ASL_InvEff.reco.scan_number=ASL_InvEff.reco.scan_number;
%     end
%     
%     ASL_InvEff.reco.displayedecho=ASL_InvEff.reco.no_echoes;
%     ASL_InvEff.reco.displayedslice=ASL_InvEff.reco.no_slices;
%     ASL_InvEff.reco.displayedexpt=1;
%     ASL_InvEff.reco=orderfields(ASL_InvEff.reco);
%     
%     ASL_InvEff.clip=[0 100 1];
% end

OutputImages = ASL_InvEff;


% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
%info2.PixelDimensions = [info.PixelDimensions, 0];
info2.ImageSize = size(OutputImages);
info2.Description = [info.Description, 'Modified by ASL_InvEff Module'];

% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 
