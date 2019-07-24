function Mat = ReplaceNaNCell(nestCell)

% Function to treat nested cells returned by ReadJson in some cases
% Replaces char elements by NaN, applies cell2mat on inner cells, then
% cell2mat on the main cell
% ********************************************
% Usage: Mat = ReplaceNaNCell(nestCell)
% ********************************************
% Input: nestCell ={M*1}{N*1}
% Output: Mat = M*N matrix

Mat = cellfun(@fun, nestCell, 'UniformOutput', 0);
Mat = cell2mat(Mat);
end

function [mat1D] = fun(Cell1D)
    charPos = cellfun(@ischar, Cell1D);
    Cell1D{charPos} = str2double(Cell1D{charPos});
    mat1D = cell2mat(Cell1D');
end
