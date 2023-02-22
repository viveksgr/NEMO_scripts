function binarizer_mask(maskname,maskthresh,mpath)
% Binarizes nii mask in the mpath (default pwd) with values <maskthresh = 0 and
% >=maskthresh=1. Saves output with a prefix_m in the same directory.
if nargin<3
    mpath = pwd;
end

mask = spm_read_vols(spm_vol(fullfile(mpath,maskname)));
mask(mask<maskthresh)=0;
mask(mask>0)=1;

filename = sprintf('m_%s',maskname(1:end-4));
write_reshaped_nifty(mask,pwd,false,maskname,filename)