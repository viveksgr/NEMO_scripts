function [corrs_l,dist_l] = corrs_and_dist(corrs_mat_l,i_red,j_red,k_red)
%Correlation matrix and 3d locations of voxels provided. Returns list of
%correlation and distances for each possible pair.
n_iter = sum(sum(triu(ones(length(i_red)),1)));
dist_l = zeros(1,n_iter);
corrs_l = zeros(1,n_iter);
kk = 1;
for ii = 1:length(i_red) % First voxel
    for jj = ii:length(i_red)
        dist = sqrt((i_red(ii)-i_red(jj))^2+(j_red(ii)-j_red(jj))^2+(k_red(ii)-k_red(jj))^2);
        if dist>0
            corrs_l(kk) = corrs_mat_l(ii,jj);
            dist_l(kk) = dist;
            kk = kk+1;
        end
    end
end
dist_l = 2*sqrt(3)*dist_l;
end

