function data_selected = findpatient_selected(handles)

% If the VOI list is displayed and the selected patient has no VOI, the
% "scan_seleted" variable is set to 1 and the returned value is the first
% scan of the same patient and the same TP.


patient_seleted = get(handles.MP3_name_list, 'String');
patient_seleted = patient_seleted(get(handles.MP3_name_list, 'Value'),:);
time_point_seleted = get(handles.MP3_time_points_list, 'String');
time_point_seleted = time_point_seleted(get(handles.MP3_time_points_list, 'Value'),:);
scan_seleted = get(handles.MP3_scans_list, 'String');
if isempty(scan_seleted) && ~isempty(time_point_seleted)
    scan_seleted = 1;
else
    scan_seleted= scan_seleted(get(handles.MP3_scans_list, 'Value'),:);
end

data_selected = nan(size(patient_seleted,1) + size(time_point_seleted,1) +size(scan_seleted,1) - 2,1);
for i = 1:size(patient_seleted,1)
    tmp = find(handles.database.Patient == patient_seleted(i,:));
    data_selected(i) = tmp(1);
end