function output = char2logical(input)
%% Function for convert char to logical
% wrtitten by C.Debacker
if strcmpi(input,'yes') || strcmpi(input,'ok') || strcmpi(input,'1') || strcmpi(input,'true') || strcmpi(input,'on')
    output = true;
elseif strcmpi(input,'no') || strcmpi(input,'cancel') || strcmpi(input,'0') || strcmpi(input,'false') || strcmpi(input,'off')
    output = false;
else
    output = [];
end