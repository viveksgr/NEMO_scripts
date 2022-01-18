function write_nifty(var_,statpath,name_,maskpath)
% Change the dimensions of var_ and save it as a nifty file in statpath
% Statpath must have mask.nii such that sum(mask) == length(var_). 
if nargin<4
    maskpath = fullfile(statpath,'mask.nii');
end

if nargin<3
    name_ = sprintf('%s.nii',inputname(1));
end

% Load Mask
mask_header = spm_vol(maskpath);
[mask,XYZmm] = spm_read_vols(mask_header);
mask(isnan(mask))=0;
mask_pred = logical(reshape(mask,1,[]));

% Make sure that the sizes match
assert(length(var_)==sum(mask_pred));

% Create image_mat
XYZvx = round(mask_header.mat\[XYZmm; ones(1,size(XYZmm,2))]);
XYZvx(end,:)=[];
voxel_id = XYZvx(:,mask_pred);
image_mat = zeros(size(mask));
for ii = 1:length(var_)
    image_mat(voxel_id(1,ii),voxel_id(2,ii),voxel_id(3,ii))=var_(ii);
end

% Nan mask
% image_mat(~logical(mask))=nan;

corr_v = mask_header(1);
corr_v.fname = fullfile(statpath,name_);
corr_v.private.dat.fname = corr_v.fname;
corr_v.dt = [64,0];
spm_write_vol(corr_v,image_mat)

% if smooth_
%     matlabbatch = [];
%     matlabbatch{1}.spm.spatial.smooth.data = {corr_v.fname};
%     matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
%     matlabbatch{1}.spm.spatial.smooth.dtype = 0;
%     matlabbatch{1}.spm.spatial.smooth.im = 1;
%     matlabbatch{1}.spm.spatial.smooth.prefix = 's';    
%     spm_jobman('run', matlabbatch);
% end

end

