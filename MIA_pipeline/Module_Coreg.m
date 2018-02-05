function [files_in,files_out,opt] = Module_Coreg(files_in,files_out,opt)
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
    %fields   = {'Type', 'HSize', 'Sigma', 'flag_test' , 'folder_out', 'output_filename_ext'};
    fields   = {'folder_out', 'flag_test', 'output_filename_ext', 'FinalResolution', 'Function', 'Separation', 'Tolerence', 'Hist_Smooth', 'Interpolation', 'Warpping', 'Masking'};
    defaults = {'', true, '_Coreg', 'Unchanged', 'mi', 'Auto= [slice thickness voxel_size voxel_size/2]', '0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001', '7 7', '4th Degree B-Spline', 'No wrap', false};
    opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
    
    % list of everything displayed to the user associated to their 'type'
    user_parameter_list = {'Select one scan or more as input to compute the coreg on'; 'Select one scan or more as input reference image'; 'Select one scan or more as input image to apply the coreg on'; 'Parameters'; '   .Output filename extension';  '   .FinalResolution';  '   .Function';  '   .Separation'; '   .Tolerence'; '   .Hist_Smooth'; '   .Interpolation'; '   .Warpping'; '   .Masking'; ''; ''};
    user_parameter_type = {'Scan'; 'Scan'; 'Scan'; ''; 'char'; 'cell'; 'cell'; 'char';'numeric';'numeric';'cell';'cell';'logical'; 'logical'; 'char'};
    parameter_default = {'';'';''; ''; '_Coreg'; {'Unchanged', 'Same as Ref', '64', '112', '128', '192', '256', '384', '512'}; {'mi','ncc', 'nmi', 'ecc'}; 'Auto= [slice thickness voxel_size voxel_size/2]'; '0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001'; '7 7'; {'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'}; {'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'}; false; 1;''};
    psom_parameter_list = {'';'';''; ''; 'output_filename_ext'; 'FinalResolution'; 'Function'; 'Separation';'Tolerence' ; 'Hist_Smooth';'Interpolation'; 'Warpping'; 'Masking';'flag_test'; 'folder_out' };
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields'};
    %opt.table = table(categorical(user_parameter_list), categorical(user_parameter_type), categorical(parameter_default), categorical(psom_parameter_list), 'VariableNames', VariableNames);
    opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, 'VariableNames', VariableNames);
    %opt.parameter_link_psom = {'output_filename_ext', '   .Output filename extension'; 'Type', '   .Type'; 'HSize','   .HSize'; 'Sigma', '   .Sigma'};
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
  
end
%%%%%%%%

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Smoothing:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
if ~ischar(files_in) 
    error('files in should be a char');
end

[path_nii,name_nii,ext_nii] = fileparts(char(files_in));
if ~strcmp(ext_nii, '.nii')
     error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, ext_nii, files_in);  
end

if isfield(opt,'threshold') && (~isnumeric(opt.threshold))
    opt.threshold = str2double(opt.threshold);
    if isnan(opt.threshold)
        disp('The threshold used was not a number')
        return
    end
end

%% Options
% %fields   = {'Type', 'HSize', 'Sigma', 'flag_test' , 'folder_out', 'output_filename_ext'};
% %defaults = {5, 3, 1, true, '', '_Smooth'};
% fields   = {'folder_out', 'flag_test', 'output_filename_ext', 'Type', 'HSize', 'Sigma'};
% defaults = {'', false, '_Smooth', {'gaussian'}, 3, 1};
% if nargin < 5
%     opt = psom_struct_defaults(struct(),fields,defaults);
% else
%     opt = psom_struct_defaults(opt,fields,defaults);
% end

%% Check the output files structure
%fields    = {'filename'};
%defaults  = {'output defaults' };
%files_out = psom_struct_defaults(files_out,fields,defaults);

%% Building default output names
if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
    opt.folder_out = path_nii;    
end

%if isempty(files_out.filename)
%    files_out.filename = cat(2,opt.folder_out,filesep,name_nii,opt.output_filename_ext,ext_nii);
%end

if strcmp(files_out, '')
    files_out = cat(2,opt.folder_out,filesep,name_nii,opt.output_filename_ext,ext_nii);
end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N = niftiread(files_in);
info = niftiinfo(files_in);
[path, name, ext] = fileparts(files_in);
jsonfile = [path, '/', name, '.json'];
fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
J = jsondecode(raw);

Informations = whos('N');
FilteredImages = zeros(size(N), Informations.class);
NbDim = length(size(N));
if NbDim>4
    warning('Too much dimensions. This module deals with at most 4 dimensions.')
end
h = fspecial(opt.Type,str2double(opt.HSize), str2double(opt.Sigma));
for i=1:size(N,3)
    for j=1:size(N,4)
        FilteredImages(:,:,i,j) = imfilter(N(:,:,i,j), h, 'replicate');
    end
end

info2 = info;
info2.Filename = files_out;
info2.Filemoddate = char(datetime('now'));
info2.Description = [info.Description, 'Modified by Smoothing Module'];

niftiwrite(FilteredImages, files_out, info2)
JMod = jsonencode(J);
[path, name, ext] = fileparts(files_out);
jsonfile = [path, '/', name, '.json'];
fidmod = fopen(jsonfile, 'w');
fwrite(fidmod, JMod, 'uint8');
fclose(fidmod);