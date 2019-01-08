function data_selected = get_data_selected(handles)


patient_seleted = get(handles.MP3_name_list, 'Value');
time_point_seleted = get(handles.MP3_time_points_list, 'Value');
scan_seleted = get(handles.MP3_scans_list, 'Value');

id_listing = unique(handles.database.Patient);

Patient_filter = handles.database.Patient== id_listing(patient_seleted);

tp_listing = unique(handles.database.Tp(Patient_filter));
% set(handles.MIA_time_points_list, 'String', string(tp_listing));
tp_filter = handles.database.Tp== tp_listing(time_point_seleted);

if get(handles.MP3_scan_VOIs_button, 'Value') == 0  % search inside the scan listing
    is_scan =  handles.database.Type == 'Scan';
    nii_listing = handles.database.SequenceName(Patient_filter & tp_filter & is_scan);
    
else % search inside the ROI listing
    is_ROI =  handles.database.Type == 'ROI';
    nii_listing = handles.database.SequenceName(Patient_filter & tp_filter & is_ROI);
    if isempty(nii_listing)
        data_selected = {};
        return
    end
end

for i=1:numel(nii_listing(scan_seleted))
    sequence_filter =  handles.database.SequenceName== nii_listing(scan_seleted(i));
    data_selected(i) = find(Patient_filter & tp_filter & sequence_filter == 1);
    
end
