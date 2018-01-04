function [Status,ValidSlices] = check_compatibility(FirstMap,SecondMap)
%% Function for check compatibilty of scans and VOI
% return a logical: true if scan are compatible, false if this is not the case

% FirstMap_filename = sprintf('%s_filename',FirstMap);
% SecondMap_filename = sprintf('%s_filename',SecondMap);

SliceOffsetFirstMap = FirstMap.acq.fov_offsets(:,3);
NumSlicesF = FirstMap.reco.no_slices;
if isempty(SecondMap)
    Status = true;
    ValidSlices = {true(NumSlicesF,1);[]};
    return
end
SliceOffsetSecondMap = SecondMap.acq.fov_offsets(:,3);
NumSlicesS = SecondMap.reco.no_slices;

ValidSlicesF = false(NumSlicesF,1);
ValidSlicesS = false(NumSlicesS,1);
for itSlices = 1 : NumSlicesF
    ValidSlicesTemp = abs(SliceOffsetFirstMap(itSlices) - SliceOffsetSecondMap) < 10e-10;
    if max(ValidSlicesTemp)~=0
        ValidSlicesF(itSlices) = true;
        ValidSlicesS = ValidSlicesS + logical(ValidSlicesTemp);
    end
end
ValidSlices{1} = ValidSlicesF;
ValidSlices{2} = ValidSlicesS;

if FirstMap.reco.thickness ~= SecondMap.reco.thickness
    warning_text = sprintf('##$ Can not calculate the map because there is\n##$ slice thickness missmatch between\n##$ %s\n##$ and \n##$ %s\n',...
        FirstMap.acq.ppl_name,SecondMap.acq.ppl_name);
    msgbox(warning_text, sprintf('%s warning',SecondMap.acq.ppl_name));
    Status = false;
    ValidSlices = [];
    return
elseif max(ValidSlicesF) == 0
    warning_text = sprintf('##$ Can not calculate the map because there is\n##$ slice offset missmatch between\n##$ %s\n##$ and \n##$ %s\n',...
        FirstMap.acq.ppl_name,SecondMap.acq.ppl_name);
    msgbox(warning_text, sprintf('%s warning',SecondMap.acq.ppl_name));
    Status = false;
    ValidSlices = [];
    return
elseif FirstMap.reco.no_samples ~= SecondMap.reco.no_samples
    warning_text = sprintf('##$ Can not calculate the map because there is\n##$ number of sample missmatch between\n##$ %s\n##$ and \n##$ %s\n',...
        FirstMap.acq.ppl_name,SecondMap.acq.ppl_name);
    msgbox(warning_text, sprintf('%s warning',SecondMap.acq.ppl_name));
    Status = false;
    ValidSlices = [];
    return
else
    Status = true;
    return
end
 