function [kidx] = k_index_sorter(idx,target)
% Resort indices of idx to target. Idx and kidx are nx1 and target is mx1
% (m classes)

[~,target] = sort(target); % Inverse of permutation target
kidx = zeros(length(idx),1);
for ii = 1:length(target)
    ids = idx==ii;
    kidx(ids) = target(ii);
end
    

