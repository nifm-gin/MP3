function [Map,MapStatus] = extract_map(Map_filename,finalMap,Type)
%% Function for retreive map
% input: filename of map or roi, the map to return in script, and type (map or roi)
% ouput: Map -> return map image or roi, Status -> logical which is true if map is loaded or if filename is empty

if(~isempty(Map_filename))
    fid=fopen(Map_filename ,'r');
    if fid>0
        fclose(fid);
        Map = load(Map_filename);
        if strcmp(Type,'VOI')
            Map = Map.uvascroi;
        elseif strcmp(Type,'map')
            Map = Map.uvascim.image;
        end
        MapStatus = true;
    else
        warning_text = sprintf('##$ Can not calculate the %s map because there is\n##$ Something wrong with the data\n##$%s=%s\n##$',...
            finalMap,Type,Map_filename);
        msgbox(warning_text, sprintf('%s map warning',finalMap));
        Map = [];
        MapStatus = false;
    end
else
    Map = [];
    MapStatus = true;
end
