function tf = startsWith(s, pattern, varargin)
%STARTSWITH True if text starts with pattern.
%   TF = startsWith(S,PATTERN) returns true if any element of string array S
%   starts with PATTERN. TF is the same size as S.
%
%   S can be a string array, a character vector, or a cell array of
%   character vectors. So can PATTERN. PATTERN and S need not be the same
%   size. If PATTERN is nonscalar, startsWith returns true if S starts with
%   any element of PATTERN.
%
%   TF = startsWith(S,PATTERN,'IgnoreCase',IGNORE) ignores case when searching 
%   for PATTERN at the start of S if IGNORE is true. The default value of IGNORE 
%   is false.
%
%   Examples
%       S = string('data.tar.gz');
%       P = string('data');
%       startsWith(S,P)                     returns  1
%
%       S = string({'abstracts.docx','data.tar.gz'});
%       P = 'data';         
%       startsWith(S,P)                     returns  [0 1]
%
%       S = string('data.tar.gz');
%       P = {'ab','data'};
%       startsWith(S,P)                     returns  1
%
%       S = string({'DATA.TAR.GZ','SUMMARY.PPT'});
%       P = string('data');
%       startsWith(S,P,'IgnoreCase',true)   returns  [1 0]
%
%   See also endsWith, contains.

%   Copyright 2015-2016 The MathWorks, Inc.

    narginchk(2, inf);

    if ~ischar(s) && ~iscellstr(s) && ~isstring(s)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringS = char(s);
        tf = stringS.startsWith(pattern, varargin{:});
    catch E
        throw(E);
    end
end
