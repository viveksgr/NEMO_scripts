function apply_all_mask(fname,s)

dir_ = pwd;
file_n = fullfile(pwd,fname);
file_ = spm_read_vols(spm_vol(file_n));

if s==1
   maskpath = 'C:\Data\NEMO\NEMO_01\imaging\nii\masks\ARC2_anatgw.nii';
elseif s==2
    maskpath = 'C:\Data\NEMO\NEMO_02\imaging\nii\masks\ARC2_anatgw.nii';
elseif s==4
    maskpath = 'C:\Data\NEMO\NEMO_04\imaging\nii\masks\ARC2_anatgw.nii';
end

all_mask = spm_read_vols(spm_vol(maskpath));
all_mask(isnan(all_mask)) = 0;
mask_ = logical(all_mask);
file_(~mask_) = 0;

eval(sprintf('%s_masked=%s;',fname(1:end-4),'file_'))
eval(sprintf('write_reshaped_nifty(%s_masked,pwd,false,fname);',fname(1:end-4)))
end