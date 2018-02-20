function data_selected = finddata_selected(handles)


patient_seleted = get(handles.MIA_name_list, 'String');
patient_seleted = patient_seleted(get(handles.MIA_name_list, 'Value'),:);
time_point_seleted = get(handles.MIA_time_points_list, 'String');
time_point_seleted = time_point_seleted(get(handles.MIA_time_points_list, 'Value'),:);
scan_seleted = get(handles.MIA_scans_list, 'String');
scan_seleted= scan_seleted(get(handles.MIA_scans_list, 'Value'),:);

data_selected = nan(size(patient_seleted,1) + size(time_point_seleted,1) +size(scan_seleted,1) - 2,1);
ii=1;
for i = 1:size(patient_seleted,1)
    for j= 1:size(time_point_seleted,1)
        for z = 1:size(scan_seleted,1)
            data_selected(ii) = find(handles.database.Patient == patient_seleted(i,:) &...
                handles.database.Tp == time_point_seleted(j,:) &...
                handles.database.SequenceName == scan_seleted(z,:));
            ii=ii+1;
        end
    end    
end
