function [status] = CreateJsonFromNifti(Nifti_filename, SequenceName, Tp, template_json)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('template_json')
    template_json = fullfile([fileparts(mfilename('full')),'template.json']);
end
if ~exist('SequenceName')
    SequenceName = 'Undefined';
end
if ~exist('Tp')
    Tp = 'Undefined';
end

fid = fopen(template_json);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
JSON = jsondecode(str);
new_JSON = JSON;

info = niftiinfo(Nifti_filename);

json_filename = strrep(Nifti_filename, '.nii', '.json');
new_JSON.DatasetDataFile = {Nifti_filename};
new_JSON.DatasetHeaderFile = {Nifti_filename};

new_JSON.SequenceName.value = char(SequenceName);
new_JSON.ProtocolName.value = char(SequenceName);
new_JSON.AcquisitionDate.value = char(Tp);

new_JSON.CreationDate.value = {date};
new_JSON.x1stRowAffineTransform.value = info.Transform.T(:,1); % [4*1 double]
new_JSON.x2ndRowAffineTransform.value = info.Transform.T(:,2); % [4*1 double]
new_JSON.x3rdRowAffineTransform.value = info.Transform.T(:,3); % [4*1 double]
new_JSON.Datatype = num2str(info.raw.datatype);
text = jsonencode(new_JSON);
fid2 = fopen(json_filename, 'w');
fprintf(fid2, text);
fclose(fid2);
end

