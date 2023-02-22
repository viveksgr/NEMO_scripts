function joined_array = arr2vec(vec1,vec2, sort_switch)
% Creates a 1D cell array of arrays from an array of vectors
% sort_switch = 1, sort each cell of the cell array in ascending order.

if nargin<3
    sort_switch = 1;
end

if iscell(vec1)
    vec1 = cell2mat(vec1);
end

if iscell(vec2)
    vec2 = cell2mat(vec2);
end

if sort_switch==1
    joined_array = mat2cell(sort([vec1, vec2],2), ones(length(vec1),1), [2]);
else
    joined_array = mat2cell([vec1, vec2], ones(length(vec1),1), [2]);
end

