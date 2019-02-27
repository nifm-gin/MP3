function cell_files = psom_files2cell(files,flag_rec)
% Convert a string, cell of strings or a structure where each field is a
% string into one single cell of strings.
%
% SYNTAX :
% CELL_FILES = PSOM_FILES2CELL(FILES)
%
% _________________________________________________________________________
% INPUTS :
%
% FILES         
%       (string, cell of strings or a structure where each terminal field 
%       is a string or a cell of strings) All those strings are file names.
%
% _________________________________________________________________________
% OUTPUTS :
%
% CELL_FILES    
%       (cell of strings) all file names in FILES stored in a cell of 
%       strings (wheter it was initially string, cell of strings or 
%       structure does not matter).
%
% _________________________________________________________________________
% COMMENTS : 
%
% Empty file names, or file names equal to 'gb_niak_omitted' are ignored.
%
% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008-2010.
% DIRO, CRIUGM, University of Montreal, 2010-2017.
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.
% Keywords : string

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

if isstruct(files)     % That's a structure
    
    cell_files = struct2cell(files);
    cell_files = cell_files(:)';
    if ~iscellstr(cell_files)
        for ii = 1:length(cell_files)
            if ~ischar(cell_files{ii})
                cell_files{ii} = psom_files2cell(cell_files{ii},true);
            end
        end
        cell_files = [cell_files{:}];
    end

elseif iscellstr(files) %% That's a cell
    
    cell_files = files(:)';
    
    
elseif ischar(files) % That's a string

    cell_files = {files};
    
else    
    
    if ~isempty(files)
        error('FILES should be a string or a cell of strings, or a structure with arbitrary depths whos terminal fields are strings or cell of strings');
    end
    
end

if (nargin == 1)&&~isempty(cell_files)
    mask = strcmp(cell_files,'gb_niak_omitted');
    mask = mask | strcmp(cell_files,'gb_psom_omitted');
    mask = mask | strcmp(cell_files,'');
    cell_files = cell_files(~mask);
    cell_files = regexprep(cell_files,[filesep '+'],filesep);
end


%!test
%! a_str = '/path/to/file.ext';
%! a_funky_str = '/path/to//file.ext';
%! an_omited_str = 'gb_psom_omitted';
%! an_niak_omited_str = 'gb_niak_omitted';
%! an_empty_str = '';
%! a_cell_o_str = {a_str,a_funky_str;an_omited_str,an_empty_str} ;
%! assert(psom_files2cell(a_str), {a_str}) ;
%! assert(psom_files2cell(a_funky_str), {a_str}) ;
%! assert(psom_files2cell(an_omited_str), {}) ;
%! assert(psom_files2cell(an_niak_omited_str), {}) ;
%! assert(psom_files2cell(an_empty_str), {}) ;
%! assert(psom_files2cell(a_cell_o_str), {a_str,a_str}) ;
%! a_struct_o_str.a_str = a_str ; 
%! a_struct_o_str.a_funky_str = a_funky_str ;
%! a_struct_o_str.an_omited_str = an_omited_str ;
%! a_struct_o_str.an_empty_str = an_empty_str ;
%! a_struct_o_str.a_cell_o_str = a_cell_o_str ;
%! assert(psom_files2cell(a_struct_o_str), {a_str, a_str, a_str, a_str}) ;