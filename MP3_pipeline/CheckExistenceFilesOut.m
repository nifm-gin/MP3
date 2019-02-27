function [flag_exist] = CheckExistenceFilesOut(Job)
%UNTITLED Summary of this function goes here
%   flag_exist = 1 means all the files exist.
%   flag_exist = 0 means at least one file is missing.
if isempty(Job.files_out)
    flag_exist = 0;
    return
end
Files_out = fieldnames(Job.files_out);
Flags = zeros(1, length(Files_out));
for i=1:length(Files_out)
    file = Job.files_out.(Files_out{i});
    Fl = zeros(1,length(file));
    for j=1:length(file)
        [~, name, ext] = fileparts(file{j});
        filemod = strrep(file{j}, ['/Tmp/', name, ext], ['/Derived_data/', name, ext]);
        Fl(j) = exist(filemod, 'file');
    end
    Flags(i) = all(Fl);
end
flag_exist = all(Flags);
end
