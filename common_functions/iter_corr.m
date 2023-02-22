function [corr_voxel_final] = iter_corr(pred_vox,test_vox)
%Two matrices are NxM where N are independent vectors and M is the size of
%each vector. corr_voxel gives the correlation between pairs of N vectors
%of size M.
% if nargin<3
%     speed = true; % Must have .mex corr file
% end
% 
% corr_voxel_final = zeros(size(test_vox,1),1);
% if speed
%     for ii = 1:size(test_vox,1)
%         corr_voxel_final(ii) = fastcorr(test_vox(ii,:)',pred_vox(ii,:)');
%     end
% else
%     for ii = 1:size(test_vox,1)
%         temp = corrcoef(test_vox(ii,:)',pred_vox(ii,:)');
%         corr_voxel_final(ii) = temp(1,2);
%     end
% end

% corr_voxel_final = zeros(size(test_vox,1),1);
% for ii = 1:size(test_vox,1)
%     corr_voxel_final(ii) = fastcorr(test_vox(ii,:)',pred_vox(ii,:)');
% end

% Make each column zero-mean
pred_vox = bsxfun(@minus,pred_vox,mean(pred_vox,2));
test_vox = bsxfun(@minus,test_vox,mean(test_vox,2));

% L2 normalize each column
pred_vox = bsxfun(@times,pred_vox,1./sqrt(sum(pred_vox.^2,2)));
test_vox = bsxfun(@times,test_vox,1./sqrt(sum(test_vox.^2,2)));

% Take the dot product of the columns and then sum
corr_voxel_final=sum(pred_vox.*test_vox,2);
end

