function [files_in,files_out,opt] = Module_Fit_T2_T2star(files_in,files_out,opt)
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
    module_option(:,3)   = {'trash_below',0};
    module_option(:,4)   = {'trash_after',Inf};
    module_option(:,5)   = {'threshold',5};
    module_option(:,6)   = {'output_filename_ext','_Fit_T2_T2star'};
    module_option(:,7)   = {'OutputSequenceName','Extension'};
    module_option(:,8)   = {'RefInput',1};
    module_option(:,9)   = {'InputToReshape',1};
    module_option(:,10)   = {'Table_in', table()};
    module_option(:,11)   = {'Table_out', table()};
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
        {'This module performs a voxel-by-voxel fitting using the following formula.'
    'y=abs(x(2) * exp(-xdata/x(1)))'
    'If a value equal to -1: the voxel signal is below the threshold or the fitting algorithm did not converge'}'};
user_parameter(:,2)   = {'Select a Multi Gradient Echo or a Multi Spin Echo scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
user_parameter(:,3)   = {'Parameters','','','','', '',''};
user_parameter(:,4)   = {'   .Output filename extension','char','_Fit_T2_T2star','output_filename_ext','','',...
    {'Specify the string to be added to the filename input.'
    'Default filename extension is ''_Fit_T2_T2star''.'}'};
user_parameter(:,5)   = {'   .Threshold','numeric', 5,'threshold','','',...
    {'Specify the minimal signal necessary to perform the fitting (voxel-wise)'
    'Default value is 5 (which corresponds to a threshold of 5 percent of the maximal signal)'}'};
user_parameter(:,6)   = {'   .Delete echoes before (ms)','numeric', 0,'trash_below','','',...
    {'User can use this setting in order to apply a fitting only on a certain part of the data'
    'Every data acquired before the user defined value (in millisecond) will be trash'
    'Default valuie is 0 --> the fitting is performed on every data'}'};
user_parameter(:,7)   = {'   .Delete echoes after (ms)','numeric', Inf,'trash_after','','',...
    {'User can use this setting in order to apply a fitting only on a certain part of the data'
    'Every data acquired after the user defined value (in millisecond) will be trash'
    'Default valuie is Inf --> the fitting is performed on every data'}'};

VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);

% So far no input file is selected and therefore no output
%
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'T2star_map'};
    
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
    error('T2star_map:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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



% define the threshold and variables
echotime = J.EchoTime;
firstecho = sum(echotime < opt.trash_below) + 1;
first_echotime=echotime(firstecho);

% Get information from the JSON data
if isfield(J, 'SpinEchoTime')
    lastecho = sum(echotime<J.SpinEchoTime/2);
else
    lastecho = length(echotime);
end

last_echotime= echotime(lastecho);
% update the last_echotime if the user add a specific value
if last_echotime > opt.trash_after
    lastecho = sum(echotime < opt.trash_after);
end
echotimes_used = echotime(firstecho:lastecho);

% reshape the data to a vector matric (speed the fitting process)
%data_to_fit = reshape(double(data.img), [size(data.img,1)*size(data.img, 2)*size(data.img,3) numel(EchoTime)]);
data_to_fit = reshape(double(N), [size(N,1)*size(N, 2)*size(N,3) numel(echotime)]);
maxim=max(data_to_fit(:)) * opt.threshold/100;

%% create empty structures
T2map_tmp = NaN(size(data_to_fit,1),1);
%fit_err = NaN(size(data_to_fit,1),1);
%wrongpix = NaN(size(data_to_fit,1),1);

parfor voxel_nbr=1:size(data_to_fit,1)
    tempydata=data_to_fit(voxel_nbr,:);
    if tempydata(1)>maxim
        % fit initializing
        t2s=(first_echotime-last_echotime)/log(tempydata(lastecho)/tempydata(firstecho));
        if t2s<=0 || isnan(t2s)
            t2s=30;
        end
        % apply the fit
        %[aaa, bbb,  ~]=levenbergmarquardt('AB_t2s',echotime_used, tempydata(firstecho:lastecho)',[t2s max(tempydata(firstecho:lastecho))*1.5]);
        [aaa, ~,  ~]=levenbergmarquardt('AB_t2s',echotimes_used, tempydata(firstecho:lastecho)',[t2s max(tempydata(firstecho:lastecho))*1.5]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 %#ok<AND2>
            T2map_tmp(voxel_nbr)=aaa(1);
           % fit_err(voxel_nbr)=bbb(1);
        %else
            %wrongpix(voxel_nbr)=2; % the fit does not work
        end
   % else % else below the fit threshold
    %    wrongpix(voxel_nbr)=1;
    end
end

% reshape matrix
OutputImages=reshape(T2map_tmp,[size(N,1) size(N, 2) size(N,3)]);
OutputImages(OutputImages < 0) = -1;
OutputImages(OutputImages > 5000) = -1;
OutputImages(isnan(OutputImages)) = -1;

% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);
info2.Description = [info.Description, 'Modified by T2star_map Module'];

% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);

% so far copy the .json file of the first input
copyfile(strrep(files_in.In1{1}, '.nii', '.json'), strrep(files_out.In1{1}, '.nii', '.json'))
% 
