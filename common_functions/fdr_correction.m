% FDR correction on correlation maps. Run 2-step 2. Outdated. New version
% in adaptive permutation.

p_vals  = zeros(size(corr_voxel_final));
for ii = 1:length(p_vals)
    p_vals(ii)=r2p(corr_voxel_final(ii),160);
end

[sort_p, argsort] = sort(p_vals);
p_list = 1:1:length(p_vals);
p_list = p_list*0.01/length(p_list);

p_res = p_list<sort_p';
p_thresh_ind = find(p_res,1);
p_thresh = sort_p(p_thresh_ind);