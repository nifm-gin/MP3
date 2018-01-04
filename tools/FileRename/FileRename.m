function [Status, Msg] = FileRename(Source, Dest, Mode)
% Rename file or folder
% This function renames the existing file or folder specified by the string
% Source to the name given by the string Dest. You can use FileRename to move
% a file from one folder to another folder or drive, but folders can be renamed
% only, not moved.
%
% Files and folders can be renamed by Matlab's MOVEFILE also, but this C-Mex is
% faster (timings vary with the size and number of the files due to the
% caching of write operations by the hard disk and the OS):
%    Matlab 2009a: 4 to 50 times faster,
%    Matlab 6.5:   1600 times faster (!).
%
% [Status, Msg] = FileRename(Source, Dest, [Mode])
% INPUT:
%   Source: String, name of the source file or folder.
%           Unicode and UNC paths are considered.
%   Dest:   String, name of the destination file or folder.
%   Mode:   String, if 'forced' an existing Dest file is overwritten,
%           if it is not write protected. Folders are *not* overwritten.
%           Optional, default: 'DoNotOverwrite'.
%
% OUTPUT:
%   Status: Scalar DOUBLE. Optional.
%            0: Success
%           -1: Source is not existing
%           -2: Dest is existing already
%           -3: Dest is write protected, in forced [Mode] only
%           -4: Unknown problems:
%               Source or Dest is accessed from another program,
%               Source is a folder and Dest is on another drive.
%   Msg: String, empty on success, some information in case of problems.
%
% COMPILE: The fast C-Mex file must be compiled before using.
%   See FileRename.c for details.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
%         Compiler: LCC2.4, OWC1.8, BCC5.5, MSVC2008
% Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
% Author: Jan Simon, Heidelberg, (C) 2006-2010 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R0c V:002 Sum:k2h6PfSIX+16 Date:29-Nov-2010 01:15:58 $
% $License: BSD $
% $UnitTest: uTest_FileRename $
% $File: Tools\GLFile\FileRename.m $
% History:

% Initialize: ==================================================================
persistent Warned
if isempty(Warned)
   warning(['JSimon:', mfilename, ':NoMex'], ...
      'Cannot find compiled Mex. MOVEFILE is used instead.');
end

% Do the work: =================================================================
% Fast alternative, but slower than the C-Mex:
%   java.io.File(Source).renameTo(java.io.File(Dest));

if nargin == 2
   [Status, Msg] = movefile(Source, Dest);
elseif nargin == 3
   [Status, Msg] = movefile(Source, Dest, Mode);
else
   error(['JSimon:', mfilename, ':BadNInput'], ...
      [mfilename, ': 2 or 3 inputs required.']);
end

% Handle problems:
if Status ~= 1
   if ~exist(Source, 'file')
      Status = -1;
   elseif exist(Dest, 'file')
      Status = -2;
   elseif exist(Dest, 'dir')  % Or write protected
      Status = -3;
   else
      Status = -4;
   end
end

return;
