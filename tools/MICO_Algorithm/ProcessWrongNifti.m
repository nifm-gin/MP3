function ProcessWrongNifti(filewrong)
    fid = fopen(filewrong, 'r+', 'l');
    %fseek(fid, 0, 'bof');
    % fread(fid, 1, 'int32');
    % fread(fid, 10, 'uchar'); % data_type (unused)
    % fread(fid, 18, 'uchar'); % db_name (unused)
    % fread(fid, 1, 'int32'); % extents (unused)
    % fread(fid, 1, 'short'); % session_error (unused)
    % fread(fid, 1, 'uchar'); % regular (unused)
    % dim_info = fread(fid, 1, 'uchar=>char')';
    fseek(fid, 40, 'bof');
    dim = fread(fid, 8, 'short')';
    A = dim;
    A(~dim)=1;
    fseek(fid, 40, 'bof');
    fwrite(fid, A, 'short')
    fseek(fid, 40, 'bof');
    dim = fread(fid, 8, 'short')'
    fclose(fid);
end