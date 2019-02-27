function opt_up = psom_defaults(def,opt,flag_warning)
% Update a structure based on a template of default values.
% SYNTAX: OPT_U = PSOM_DEFAULTS(DEF,OPT,FLAG_WARNING)
% 
% DEF (structure) define the default values for a series of fields.
% OPT (structure, default empty structure) user-specified values for a subset
%   of fields in DEF.
% FLAG_WARNING (boolean, default true) issue warnings for fields present in 
%   OPT and not in DEF.
% 
% NOTES
%
%   If a default value NaN is assigned to any field, the function will return an 
%   error if that field is absent in OPT.
%
%   The behaviour of PSOM_DEFAULTS is recursive, in the sense that it applies 
%   to all fields of DEF which are themselves structures.
%
% EXAMPLE
%
%   def.warp.fwhm = 2;
%   def.warp.nb_iter = 4;
%   def.verbose = true;
%   def.write = NaN;
% 
%   opt.warp.fwhm = 5;
%   opt.write = true;
% 
%   opt_u = psom_defaults(def,opt);
% 
% Returns:
%   opt_u.warp.fwhm = 5;
%   opt_u.warp.nb_iter = 4;
%   opt_u.verbose = true;
%   opt_u.write = true;
% 
% Copyright (c) Pierre Bellec, Centre de recherche de l'institut de 
% Geriatrie de Montreal, Departement d'informatique et de recherche 
% operationnelle, Universite de Montreal, 2017
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.

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

if ~isstruct(def)
  error('DEF needs to be a structure')
end

if nargin<2
  opt = struct();
end

if ~isstruct(opt)
  error('OPT needs to be a structure')
end

if nargin < 3
  flag_warning = true;
end

%% Check for unexpected fields
list_def = fieldnames(def);
list_opt = fieldnames(opt);
if flag_warning&&any(~ismember(list_opt,list_def))
  warning('The following fields are not expected in the structure %s \n    %s',inputname(1),list_opt{~ismember(list_opt,list_def)})
end

%% Check for missing fields
list_val = struct2cell(def);
is_missing = false(size(list_val));
for num_f = 1:length(list_def)
    is_missing(num_f) =   isreal(list_val{num_f}) ...
                       &&(length(list_val{num_f}(:))==1) ... 
                       &&isnan(list_val{num_f}) ...
                       &&~isfield(opt,list_def{num_f});
end
if any(is_missing)
    error(sprintf('A value must be specified in %s for the following fields (%s)',inputname(1),list_def{is_missing}));
end

%% Set defaults
opt_up = def;
for num_f = 1:length(list_opt)
    opt_up.(list_opt{num_f}) = opt.(list_opt{num_f});    
end

%% Recursively set defaults inside substructures
for ff = 1:length(list_def)
  if isstruct(def.(list_def{ff}))&&isfield(opt,list_def{ff})
    opt_up.(list_def{ff}) = psom_defaults(def.(list_def{ff}),opt.(list_def{ff}),flag_warning);
  end
end
