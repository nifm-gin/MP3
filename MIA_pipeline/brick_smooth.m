function [files_in,files_out,opt] = brick_smooth(files_in,files_out,opt)
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
%       
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('T2map:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
if ischar(files_in) 
    error('files in should be a char');
end
if size(files_in,1) == 2
    error('2 files in are required');
end
[path_nii,name_nii,ext_nii] = fileparts(char(files_in(1)));
if ~strcmp(ext_nii, '.nii')
     error('First file need to be a .nii');  
end
[path_json,name_json,ext_json] = fileparts(char(files_in(2)));
if ~strcmp(ext_json, '.json')
     error('First file need to be a .json');  
end
% 
% if ~isnumeric(opt.threshold)
%     opt.threshold = str2double(opt.threshold);
%     if isnan(opt.threshold)
%         disp('The threshold used was not a number')
%         return
%     end
% end

%% Options
fields   = {'Type', 'Sigma', 'Hsize', 'flag_test' , 'folder_out', 'output_filename_ext'};
defaults = {'gaussian', 1, 3, false, '', '_smoothed'};
if nargin < 3
    opt = psom_struct_defaults(struct(),fields,defaults);
else
    opt = psom_struct_defaults(opt,fields,defaults);
end

%% Check the output files structure
fields    = {'filename'};
defaults  = {'output defaults' };
files_out = psom_struct_defaults(files_out,fields,defaults);

%% Building default output names
if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
    opt.folder_out = path_nii;    
end

if isempty(files_out.filename)
    files_out.filename{1} = cat(2,opt.folder_out,filesep,name_nii,opt.output_filename_ext,ext_nii);
    files_out.filename{2} = cat(2,opt.folder_out,filesep,name_nii,opt.output_filename_ext,'.json');
end

%% If the test flag is true, stop here !
if opt.flag_test == 1
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load input Nii file
data.hdr = spm_vol(files_in{1});
data.img = spm_read_vols(data.hdr);

%% load input JSON fil
data.json = spm_jsonread(files_in{2});
opt.smooth.Type = 'gaussian';
opt.smooth.Sigma = 1;
opt.smooth.Hsize = 3;

h = fspecial(opt.smooth.Type,opt.smooth.Hsize,opt.smooth.Sigma);
[~, ~, E, Z] = size(data.img);
for s=1:Z
    for e=1:E
        data_smooth.img(:,:,e,s) = imfilter(data.img(:,:,e,s),h,'replicate');
    end
end
%% if spm_function is used
% data_smooth.hdr = spm_vol([files_in{1}, ', 1']);
% data_smooth.hdr.fname = files_out.filename;
% spm_write_vol(data_smooth.hdr, squeeze(data_smooth.img(:,:,10,:)));

% if matlab function is used
data_smooth.hdr = niftiinfo(files_in{1});
data_smooth.hdr.Filename = files_out.filename{1};
data_smooth.hdr = update_nifti_hdr(data_smooth.hdr, data_smooth.img, files_out.filename{1});
niftiwrite(single(data_smooth.img), files_out.filename{1}, data_smooth.hdr)

%% need to update the json structure here before saving it with the T2map
spm_jsonwrite(files_out.filename{2}, data.json);

