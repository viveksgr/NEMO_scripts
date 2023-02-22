function smoothed = nansmooth3(img,bins)
% Compute moving average for each point in img over bins ignoring nans.

if nargin<2
    bins = 1; 
end

smoothed = zeros(size(img));
for ii = 1+bins:size(img,1)-bins
    for jj = 1+bins:size(img,2)-bins
        for kk = 1+bins:size(img,3)-bins
            tiny_vox = img(ii-bins:ii+bins,jj-bins:jj+bins,kk-bins:kk+bins);        
%             kernel_ = gaussian_kernel;
%             tiny_vox = tiny_vox.*kernel_;
            smoothed(ii,jj,kk) = nanmedian(tiny_vox(:));
        end
    end
end

