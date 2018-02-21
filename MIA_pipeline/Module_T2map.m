function [files_in,files_out,opt] = Module_T2map(files_in,files_out,opt)
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
    fields   = {'threshold'  , 'flag_test' , 'folder_out', 'output_filename_ext', 'OutputSequenceName'};
    defaults = {5, true, '', 'T2map', 'AllName'};
    opt.Module_settings = psom_struct_defaults(struct(),fields,defaults);
    
    % list of everything displayed to the user associated to their 'type'
    user_parameter_list = {'Select a Multi Spin Echo scan as input'; 'Parameters'; '   .Output filename extension' ; '   .Threshold'};
    user_parameter_type = {'1Scan'; ''; 'char'; 'numeric'};
    parameter_default = {''; ''; 'T2map'; 5};
    psom_parameter_list = {''; ''; 'output_filename_ext'; 'threshold'};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields'};
    opt.table = table(user_parameter_list, user_parameter_type, parameter_default, psom_parameter_list, 'VariableNames', VariableNames);
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('T2map:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%% Options
% fields   = {'threshold'  , 'flag_test' , 'folder_out', 'output_filename_ext'};
% defaults = {5, false, '', '_T2_map'};
% if nargin < 3
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
%data.hdr = spm_vol(files_in{1});
%data.img = spm_read_vols(data.hdr);

%% load input JSON fil
%data.json = spm_jsonread(files_in{2});

% Get information from the JSON data
%EchoTime = data.json.EchoTime;

N = niftiread(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
[path, name, ext] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
%raw = reshape(raw, 1,length(raw));
J = jsondecode(raw);

EchoTime = J.EchoTime;

% reshape the data to a vector matric (speed the fitting process)
%data_to_fit = reshape(double(data.img), [size(data.img,1)*size(data.img, 2)*size(data.img,3) numel(EchoTime)]);
data_to_fit = reshape(double(N), [size(N,1)*size(N, 2)*size(N,3) numel(EchoTime)]);

%% create empty structures
T2map_tmp = NaN(size(data_to_fit,1),1);
% M0map_tmp = NaN(size(data_to_fit,1),1);
% T2_Error_map_tmp = NaN(size(data_to_fit,1),1);
% M0_Error_map_tmp = NaN(size(data_to_fit,1),1);

% define the threshold and variables
maxim=max(data_to_fit(:)) * opt.threshold/100;
t2init_Cte = EchoTime(1) - EchoTime(end-1);


%init matlabpool
% schd = parcluster();
% poolobj = parpool('local', schd.NumWorkers);

parfor voxel_nbr = 1:size(data_to_fit,1)
    tmp_voxel_data=data_to_fit(voxel_nbr,:);
    if max(tmp_voxel_data(:))>= maxim
        %% fit data
        t2init=(t2init_Cte)/log(tmp_voxel_data(end-1)/tmp_voxel_data(1));
        if t2init<=0 || isnan(t2init)
            t2init=30;
        end
        [aaa, ~,  convergence]=levenbergmarquardt('AB_t2s',EchoTime', abs(tmp_voxel_data),[t2init max(abs(tmp_voxel_data))*1.5]);
        % the algorithm converged
        if convergence == -1
            % to remove when good input data
            if isreal(aaa(1))
                T2map_tmp(voxel_nbr)=aaa(1);
%                 M0map_tmp(voxel_nbr)=aaa(2);
%                 T2_Error_map_tmp(voxel_nbr)=bbb(1);
%                 M0_Error_map_tmp(voxel_nbr)=bbb(2);
            end
        end
    end
end
% delete(poolobj);

% [~,filename,~] = fileparts(MSE_map_filename);

% reshape matrix
OutputImages=reshape(T2map_tmp,[size(N,1) size(N, 2) size(N,3)]);
OutputImages(OutputImages < 0) = -1;
OutputImages(OutputImages > 5000) = -1;
OutputImages(isnan(OutputImages)) = -1;
% M0map.img=reshape(M0map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% T2_Error_map.img=reshape(T2_Error_map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% M0_Error_map.img=reshape(M0_Error_map_tmp,[size(data.img,1) size(data.img, 2) size(data.img,3)]);
% save the T2 map
%% if spm_function is used
% T2map.hdr = spm_vol([files_in{1}, ', 1']);
% T2map.hdr.fname = files_out.filename;
% spm_write_vol(T2map.hdr, T2map.img);

%% if matlab function is used
%T2map.hdr = niftiinfo(files_in{1});
%T2map.hdr = update_nifti_hdr(T2map.hdr, T2map.img, files_out.filename);
%niftiwrite(single(T2map.img), files_out.filename, T2map.hdr)

%% need to update the json structure here before saving it with the T2map
%spm_jsonwrite(strrep(files_out.filename, '.nii', '.json'), data.json);

info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
%info2.raw.datatype = 16;
%info2.BitsPerPixel = 32;
info2.ImageSize = size(OutputImages);
%info2.raw.dim(1) = 3;
%info2.raw.dim(5) = 1;
%info2.raw.bitpix = 32;
%info2.raw = struct();
info2.Description = [info.Description, 'Modified by T2map Module'];



niftiwrite(OutputImages, files_out.In1{1}, info2);
JMod = jsonencode(J);
[path, name, ext] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
fidmod = fopen(jsonfile, 'w');
fwrite(fidmod, JMod, 'uint8');
fclose(fidmod);
% 
% % save the M0map map
% M0map.hdr = spm_vol([MSE_map_filename, ', 1']);
% M0map.hdr.fname = char(strcat(filename, '-M0map.nii'));
% M0map.img(isnan(M0map.img)) = -1;
% spm_write_vol(M0map.hdr, M0map.img);
% 
% % save the T2_Error_map 
% T2_Error_map.hdr = spm_vol([MSE_map_filename, ', 1']);
% T2_Error_map.hdr.fname = char(strcat(filename, '-T2_Error.nii'));
% T2_Error_map.img(T2_Error_map.img < 0) = -1;
% T2_Error_map.img(T2_Error_map.img > 50) = -1;
% T2_Error_map.img(isnan(T2_Error_map.img)) = -1;
% spm_write_vol(T2_Error_map.hdr, T2_Error_map.img);
% 
% % save the M0_Error_map map
% M0_Error_map.hdr = spm_vol([MSE_map_filename, ', 1']);
% M0_Error_map.hdr.fname = char(strcat(filename, '-M0_Error.nii'));
% M0_Error_map.img(isnan(M0_Error_map.img)) = -1;
% spm_write_vol(M0_Error_map.hdr, M0_Error_map.img);
