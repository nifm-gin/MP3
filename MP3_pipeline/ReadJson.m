function J = ReadJson(jsonfile)
%ReadJson Open Json file and return the data in a Matlab ariable

fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
if verLessThan('matlab','9.5')
    % -- Code to run in MATLAB R2018a and earlier here --
  J = jsondecode(raw);
else
    % -- Code to run in MATLAB R2018a and later here --
   J = jsondecode(raw');
end
