function [Status, Message, Wrong_File] = Check_files(files_in)
%% Check validity of each files of the PSOM 'files_in' structure and return :
%       Status : Booleean - 1 if all the files are valid, 0 if at least one file is incorrect.
%       Message : String - Explaining the potential errors.
%       Wrong_File : String - Name of the first invalid file.
% This function checks if the files contained in files_in are strings, if
% they are nifti files, and if matlab succeed to open them.

if ~isstruct(files_in)
    error('files_in is not a structure.')
end


Fields = fieldnames(files_in);
Wrong_File = '';
Message = 'All files are valid !';
for i=1:length(Fields)
    Input = files_in.(Fields{i});
    for j=1:length(Input)
        Status = 1;
        file = Input{j};
        if ~ischar(file)
            Status = 0;
            Message = 'Files_in should be a char';
            Wrong_File = file;
            return
        end
        [~,~,ext_nii] = fileparts(file);
        if ~strcmp(ext_nii, '.nii')
            Status = 0;
            Message = 'Files need to be .nii';
            Wrong_File = file;
            return
        end
        [fid, message]=fopen(file ,'r');
        if fid == -1
            if exist(strrep(file, '.nii', '.nii.gz'),'file')
                gunzip(strrep(file, '.nii', '.nii.gz'))
                assert(exist(file, 'file')==2)
                delete(strrep(file, '.nii', '.nii.gz'))
                continue
            end
            Status = 0;
            Message = message;
            Wrong_File = file;
            return
        else
            fclose(fid);
        end
    end
end

end