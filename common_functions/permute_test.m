function [corr_voxel_final,p_val, p_correct] = permute_test(pred_vox,test_vox,statpath,iter,smooth_)
% Do a permutation test and find the p_val for the correlation between
% pred_series and actual_series under the null hypothesis that there is 0
% correlation. Resolution = 1/iter.
% p_val = Nx1, pred_series = NxM, actual_series = NxM. N = number of
% independent vectors. M = size of each vector.

% Varargin inputs = number of iterations, smoothing and work directory with
% mask file
% VivekSagar2016@u.northwestern.edu, Feb 20, 2019.

try
    assert(and(size(pred_vox,1)==size(test_vox,1),size(pred_vox,2)==size(test_vox,2)));
catch
    error('Sizes of input matrices should be same');
end


if nargin<5
    smooth_ = false;
end

if nargin<4
    iter = 2000; % Min p value = 1/iter
end

if smooth_
    corr_voxel_ = iter_corr(pred_vox,test_vox);
    write_nifty(corr_voxel_,statpath);
    mask = spm_read_vols(spm_vol(fullfile(pwd, 'mask.nii')));
    mask_unpacked = logical(reshape(mask,1,[]));
    corr_voxel_packed = spm_read_vols(spm_vol(fullfile(statpath, 'scorr_voxel_.nii')));
    corr_voxel_unpacked = reshape(corr_voxel_packed,1,[]);
    corr_voxel_final = corr_voxel_unpacked(mask_unpacked)';
else
    corr_voxel_final = iter_corr(pred_vox,test_vox);
end

compare_matrix = zeros(size(test_vox,1),iter)'; 

for tt = 1:iter
    if rem(tt,iter/100)==0
        fprintf('%d percent complete\n',tt*100/iter)
    end
    pred_vox_shuff = pred_vox(:,randperm(size(pred_vox,2)));
    compare_matrix(tt,:) = iter_corr(pred_vox_shuff,test_vox);
end

h5write('smooth_corr.h5', pwd, compare_matrix)

p_val = zeros(size(test_vox,1),1);
for ii = 1:length(p_val)
    p_val(ii)=sum(compare_matrix(:,ii)>corr_voxel_final(ii))/iter;
end

p_correct = mafdr(p_val);
end

