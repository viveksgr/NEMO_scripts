function mask_extractor(img,mask,impath,maskpath)

% Maskpath
if nargin<4
    maskpath = pwd;
end

% Image path
if nargin<3
    impath = pwd;
end

fullim = fullfile(impath,img);
img_mat = spm_read_vols(spm_vol(fullim));
mask_mat = spm_read_vols(spm_vol(fullfile(maskpath,mask)));
mask_mat(isnan(mask_mat))=0;
mask_mat = mask_mat>0.01;
img_mat(~mask_mat)=0;

filename = sprintf('m_%s',img(1:end-4));
write_reshaped_nifty(img_mat,impath,false,img,filename)
