function J = ReadJson(jsonfile)
%ReadJson Open Json file and return the data in a Matlab ariable

fid = fopen(jsonfile, 'r');
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
J = jsondecode(raw);
end

