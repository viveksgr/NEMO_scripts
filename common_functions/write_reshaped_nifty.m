function write_reshaped_nifty(image_mat,statpath,smooth_,f_name,f_name_save)
% Change the dimensions of var_ and save it as a nifty file in statpath
% Statpath must have mask.nii such that sum(mask) == length(var_). 

if nargin<5
    f_name_save = inputname(1);
end

if nargin<4
    f_name = 'mask.nii';
end

if nargin<3
    smooth_ = false;
end


% Load Mask
mask_header = spm_vol(fullfile(statpath,f_name));
mask = spm_read_vols(mask_header);

% Nan mask
% image_mat(~logical(mask))=nan;

corr_v = mask_header(1);
corr_v.fname = fullfile(statpath,sprintf('%s.nii',f_name_save));
corr_v.private.dat.fname = corr_v.fname;
corr_v.dt = [64,0];
spm_write_vol(corr_v,image_mat)

if smooth_
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.smooth.data = {corr_v.fname};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';    
    spm_jobman('run', matlabbatch);
end

end

