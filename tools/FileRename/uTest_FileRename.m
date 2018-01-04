function uTest_FileRename(doSpeed)
% Automatic test: FileRename
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_FileRename(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed tested are defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R0q V:011 Sum:umvYd+x3Mk+2 Date:28-Nov-2010 03:32:07 $
% $License: BSD $
% $File: Tools\UnitTests_\uTest_FileRename.m $
% History:
% 001: 10-Nov-2009 13:18, Need a test for new features.
% 005: 20-Aug-2010 10:36, BUGFIX: DIR replies the date in the local language.
%      In Germany the date of DIR can be "01-Dez-2010" and the comparison with
%      DATESTR(GetFileTime) failed, because it replies English month names.
%      Now DATENUM is used to compare the times.
% 008. 01-Oct-2010 14:24, BUGFIX: The same problem appeared again.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
ErrID = ['JSimon:', mfilename, ':'];

% Initial values: --------------------------------------------------------------
nFile = 400;  % Number of test files for speed measurements

% Program Interface: -----------------------------------------------------------
if nargin == 0
   doSpeed = true;
end

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
disp(['== Test FileRename  ', datestr(now, 0)]);

bakCD      = cd;
Temp       = tempdir;
TestFolder = fullfile(Temp, 'Test_FileRename_');
if strcmpi(cd, TestFolder);
   cd(Temp);
   bakCD = Temp;
end

% Cleanup test folder, if it was not removed due to a former crash:
if DirExist(TestFolder);
   Cleanup(TestFolder);
end

% Create the test folder:
[Success, Msg] = mkdir(Temp, 'Test_FileRename_');  %#ok<*NASGU>
if Success ~= 1
   error([ErrID, 'NoFolder'], 'Cannot create folder?!');
end

cd(TestFolder);

Folder1 = fullfile(TestFolder, 'Folder1');
Folder2 = fullfile(TestFolder, 'Folder2');

% Create some test files: ------------------------------------------------------
fprintf('\nCreate test files...\n');
TestFile = cell(1, nFile);
x = rand(100, 100);
for i = 1:nFile
   Name = sprintf('File%.3d.test', i);
   FID = fopen(Name, 'w');
   if FID < 0
      error('Cannot create file?!');
   end
   fwrite(FID, x);
   fclose(FID);
   TestFile{i} = Name;
end
TestFile = strcat(TestFolder, filesep, TestFile);
FileDir  = dir(TestFile{1});
FileSize = FileDir.bytes;
fprintf('\n');

% At first without using output arguments: -------------------------------------
% FileRename stops with an error on problems!
fprintf('== Without catching outputs:\n');

% Rename one file:
File_1 = TestFile{1};
File_2 = [File_1, '2'];
try
   FileRename(File_1, File_2);
catch
   error('Crashed: FileRename(file, not existing)\n%s', lasterr);
end
if FileExist(File_1) || ~FileExist(File_2)
   error('Failed: FileRename(file, not existing)');
end
fprintf('  ok: FileRename(file, not existing)\n');

tooLazy = false;
try
   FileRename(File_1, File_2);
   tooLazy = true;
catch
   fprintf('  ok: FileRename(not existing file, existing file) rejected\n');
end
if tooLazy
   error('Failed: FileRename(not existing file, existing file) not rejected');
end

try
   FileRename(File_2, File_1);
catch
   error('Crashed: FileRename(file, not existing)\n%s', lasterr);
end
if FileExist(File_2) || ~FileExist(File_1)
   error('Failed: FileRename(file, not existing)');
end
fprintf('  ok: FileRename(file, not existing)\n');

% Try to overwrite existing file:
[Status, Msg] = copyfile(File_1, File_2);
if Status ~= 1
   error('Cannot copy file?!\n%s', Msg);
end

try
   FileRename(File_1, File_2);
   tooLazy = true;
catch
   fprintf('  ok: FileRename(file, existing file) rejected\n');
end
if tooLazy
   error('Failed: FileRename(file, existing file) not rejected');
end

try
   FileRename(File_1, File_2, 'forced');
   fprintf('  ok: FileRename(file, existing file, forced)\n');
catch
   error('Crashed: FileRename(file, existing file, forced):\n%s', lasterr);
end

% Recreate File_1:
[Status, Msg] = copyfile(File_2, File_1);
if Status ~= 1
   error('Cannot copy file?!\n%s', Msg);
end

[Status, Msg] = mkdir(TestFolder, 'Folder1');
if Status ~= 1
   error('Cannot create folder?!\n%s', Msg);
end

try
   FileRename(File_1, Folder1);
   tooLazy = true;
catch
   fprintf('  ok: FileRename(file, existing folder) rejected\n');
end
if tooLazy
   error('Failed: FileRename(file, existing folder) not rejected\n');
end

try
   FileRename(File_1, Folder1, 'forced');
   tooLazy = true;
catch
   fprintf(['  ok: ', ...
      'FileRename(file, existing folder, forced) rejected\n']);
end
if tooLazy
   error(['Failed: FileRename(file, existing folder, ', ...
      'forced) not rejected\n']);
end

% Write protected destination:
fileattrib(File_2, '-w');
try
   FileRename(File_1, File_2);
   tooLazy = true;
catch
   fprintf('  ok: FileRename(file, protected file) rejected\n');
end
if tooLazy
   error('FileRename(file, protected file) not rejected?!');
end

try
   FileRename(File_1, File_2, 'forced');
   tooLazy = true;
catch
   fprintf('  ok: FileRename(file, protected file, force) rejected\n');
end
if tooLazy
   error('Failed: FileRename(file, protected file, force) not rejetced?!');
end

fileattrib(File_2, '+w');

% Write protected source:
fileattrib(File_1, '-w');
FileRename(File_1, 'Dummy');
if FileExist('Dummy') && ~FileExist(File_1)
   fprintf('  ok: FileRename(protected file, not existing)\n');
else
   error('Failed: FileRename(protected file, not existing) not successful');
end

FileRename('Dummy', File_1);
fileattrib(File_1, '+w');

% Rename a folder:
FileRename(Folder1, Folder2);
if DirExist(Folder2) && ~DirExist(Folder1)
   fprintf('  ok: FileRename(folder, not existing)\n');
else
   error('Failed: FileRename(folder, not existing) did not work');
end

try
   FileRename(Folder1, Folder2, 'forced');
   tooLazy = true;
catch
   fprintf('  ok: FileRename(not existing, folder, forced) rejected\n');
end
if tooLazy
   error('Failed: FileRename(not existing, folder, forced) not rejected');
end
   
if DirExist(Folder2) || ~DirExist(Folder1)
   fprintf('  ok: Existence of folders as expected\n');
else
   error('Failed: Strange existence of folders?!');
end

[Status, Msg] = mkdir(TestFolder, 'Folder1');
try
   FileRename(Folder1, Folder2, 'forced');
   tooLazy = true;
catch
   fprintf('  ok: FileRename(folder, existing folder, forced) rejected\n');
end
if tooLazy
   error('Failed: FileRename(folder, existing folder, forced) not rejected');
end

% Now with using an output: ----------------------------------------------------
fprintf('== With catching outputs:\n');
delete(File_2);
try
   [Status, Msg] = FileRename(File_1, File_2);
catch
   error('Crashed: [S,M]=FileRename(file, not existing)\n%s', lasterr);
end
if FileExist(File_1) || ~FileExist(File_2) || Status ~= 1
   error('Failed: [S,M]=FileRename(file, not existing)\n%s\%s', File_1, File_2);
end
fprintf('  ok: [S,M]=FileRename(file, not existing)\n');

try
   [Status, Msg] = FileRename(File_1, File_2);
catch
   error('Crashed: [S,M]=FileRename(not existing, file)\n%s', lasterr);
end
if Status == 1
   error('Failed: [S,M]=FileRename(not existing, file)');
else
   fprintf('  ok: [S,M]=FileRename(not existing, file) => Status=%d\n', Status);
end

try
   [Status, Msg] = FileRename(File_2, File_1);
catch
   error('Crashed: [S,M]=FileRename(file, not existing)\n%s', lasterr);
end
if FileExist(File_2) || ~FileExist(File_1) || Status ~= 1
   error('Failed: [S,M]=FileRename(file, not existing)');
end
fprintf('  ok: [S,M]=FileRename(file, not existing)\n');

% Try to overwrite existing file:
[Status, Msg] = copyfile(File_1, File_2);
if Status ~= 1
   error('Cannot copy file?!\n%s', Msg);
end

try
   [Status, Msg] = FileRename(File_1, File_2);
catch
   error('Crashed: [S,M]=FileRename(file, existing file)\n%s', lasterr);
end
if Status < 0
   fprintf('  ok: [S,M]=FileRename(file, existing file) => Status=%d\n', ...
      Status);
else
   error('Failed: [S,M]=FileRename(file, existing file) not rejected');
end

try
   [Status, Msg] = FileRename(File_1, File_2, 'forced');
catch
   error('Crash: [S,M]=FileRename(file, existing file, forced): \n%s', lasterr);
end
if Status == 1
   fprintf('  ok: [S,M]=FileRename(file, existing file, forced)\n');
else
   error(['Failed: ', ...
         '[S,M]=FileRename(file, existing file, forced) => Status=%d'], Status);
end

% Recreate File_1:
[Status, Msg] = copyfile(File_2, File_1);
if Status ~= 1
   error('Cannot copy file?!\n%s', Msg);
end

[Status, Msg] = mkdir(TestFolder, 'Folder1');
try
   [Status, Msg] = FileRename(File_1, Folder1);
catch
   error('Crash: [S,M]=FileRename(file, existing folder)\n%s', lasterr);
end
if Status == 1 || ~FileExist(File_1) || ~DirExist(Folder1)
   error('Failed: [S,M]=FileRename(file, existing folder) worked');
else
   fprintf(['  ok: [S,M]=FileRename(file, existing folder) => ', ...
      'Status=%d\n'], Status);
end

try
   [Status, Msg] = FileRename(File_1, Folder1, 'forced');
catch
   error('Crash: [S,M]=FileRename(file, existing folder, forced)\n%s', lasterr);
end
if Status == 1 || ~FileExist(File_1) || FileExist(Folder1)
   error('Failed: [S,M]=FileRename(file, existing folder, forced) worked');
else
   fprintf(['  ok: [S,M]=FileRename(file, existing folder) => ', ...
      'Status=%d\n'], Status);
end

% Write protected destination:
fileattrib(File_2, '-w');
try
   [Status, Msg] = FileRename(File_1, File_2);
catch
   error('Crash: [S,M]=FileRename(file, protected file)\n%s', lasterr);
end
if Status == 1 || ~FileExist(File_1) || ~FileExist(File_2)
   error('Failed: [S,M]=FileRename(file, protected file) worked');
else
   fprintf(['  ok: [S,M]=FileRename(file, protected file) => ', ...
      'Status=%d\n'], Status);
end

try
   [Status, Msg] = FileRename(File_1, File_2, 'forced');
catch
   error('Crash: [S,M]=FileRename(file, protected file, force)\n%s', lasterr);
end
if Status == 1 || ~FileExist(File_1) || ~FileExist(File_2)
   error('Failed: [S,M]=FileRename(file, protected file, force) worked');
else
   fprintf(['  ok: [S,M]=FileRename(file, protected file, force) => ', ...
      'Status=%d\n'], Status);
end

fileattrib(File_2, '+w');

% Write protected source:
fileattrib(File_1, '-w');
try
   [Status, Msg] = FileRename(File_1, 'Dummy');
catch
   error('Crash: [S,M]=FileRename(protected file, not existing)\n%s', lasterr);
end
if FileExist('Dummy') && ~FileExist(File_1) && Status == 1
   fprintf('  ok: [S,M]=FileRename(protected file, not existing)\n');
else
   error('Failed: [S,M]=FileRename(protected file, not existing) did not work');
end

[Status, Msg] = FileRename('Dummy', File_1);
fileattrib(File_1, '+w');

% Rename a folder:
rmdir(Folder2);
[Status, Msg] = FileRename(Folder1, Folder2);
if DirExist(Folder2) && ~DirExist(Folder1) && Status == 1
   fprintf('  ok: [S,M]=FileRename(folder, not existing)\n');
else
   error('Failed: [S,M]=FileRename(folder, not existing)');
end

try
   [Status, Msg] = FileRename(Folder1, Folder2, 'forced');
catch
   fprintf('Crash: [S,M]=FileRename(not existing, folder)\n%s', lasterr);
end
if Status < 0
   fprintf('  ok: [S,M]=FileRename(folder, not existing) rejected\n');
else
   error('Failed: [S,M]=FileRename(folder, not existing) did work?!');
end
   
if DirExist(Folder2) || ~DirExist(Folder1)
   fprintf('  ok: Existence of folders as expected\n');
else
   error('Failed: Strange existence of folders?!');
end

[Status, Msg] = mkdir(TestFolder, 'Folder1');
try
   [Status, Msg] = FileRename(Folder1, Folder2, 'forced');
catch
   fprintf('Crash: [S,M]=FileRename(not existing, folder, force)\n%s', lasterr);
end
if Status < 0
   fprintf('  ok: [S,M]=FileRename(folder, existing, force) rejected\n');
else
   error('Failed: [S,M]=FileRename(folder, existing, force) did work?!');
end

% Speed: -----------------------------------------------------------------------
if doSpeed
   fprintf('\n== Speed tests:\n  Rename %d files of %dkB\n', ...
      nFile * 2, round(FileSize / 1000));

   % Matlab 6.5 opens a DOS box for each MOVEFILE!!! Very slow.
   if sscanf(version, '%d', 1) < 7.0
      reduce = 100;
      nFileR = ceil(nFile / reduce);
      xMsg   = ' (extrapolated)';
   else
      nFileR = nFile;
      reduce = 1;
      xMsg   = '';
   end
      
   A = TestFile;
   B = strcat(A, '_');
   
   % Silent warm up for a fair comparison:
   for i = 1:nFile
      [Status, Msg] = FileRename(A{i}, B{i});
   end
   
   fprintf('  ');
   fprintf('%-16s', 'MOVEFILE', 'FileRename', 'java.io.File.renameTo');
   fprintf('\n');
   for i = 1:4
      fprintf('  ');
      
      tic;
      for i = 1:nFileR
         [Status, Msg] = movefile(B{i}, A{i});
      end
      for i = 1:nFile
         [Status, Msg] = movefile(A{i}, B{i});
      end
      mTime = toc * reduce;
      fprintf('%-16s', sprintf('%.4f sec', mTime));
      
      tic;
      for i = 1:nFile
         [Status, Msg] = FileRename(B{i}, A{i});
      end
      for i = 1:nFile
         [Status, Msg] = FileRename(A{i}, B{i});
      end
      mexTime = toc;
      fprintf('%-16s', sprintf('%.4f sec', mexTime));

      tic;
      for i = 1:nFile
         S = java.io.File(B{i});
         S.renameTo(java.io.File(A{i}));
      end
      for i = 1:nFile
         S = java.io.File(A{i});
         S.renameTo(java.io.File(B{i}));
      end
      jTime = toc;
      fprintf('%-16s', sprintf('%.4f sec', jTime));
      fprintf('\n');
   end
end

% Success! Goodbye: ------------------------------------------------------------
cd(bakCD);
if strcmpi(cd, TestFolder)
   cd(Temp);
end

Cleanup(TestFolder);

fprintf('\nFileRename seems to work fine.\n');
   
return;

% ******************************************************************************
function Cleanup(TestFolder)

% Delete the test folder:
fprintf('\nCleanup test files...\n');
if ispc
   [internStatus, Msg] = dos(['rmdir /s /q "', TestFolder, '"']);
else
   [internStatus, Msg] = unix(['rm -rf "', TestFolder, '"']);
end
if DirExist(TestFolder)
   warning(['JSimon:', mfilename, ':FolderNotDeleted'], ...
      'Cannot delete folder: %s', TestFolder);
end

return;

% ******************************************************************************
function  E = FileExist(File)
E = any(exist(File, 'file')) && not(exist(File, 'dir'));
return;

% ******************************************************************************
function  E = DirExist(File)
E = any(exist(File, 'dir'));
return;
