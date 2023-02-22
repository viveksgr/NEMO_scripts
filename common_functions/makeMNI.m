function makeMNI(nums,name_,wdir)
% Make masks in the MNI space from AAL atlas. ROIs are specified by name
% and the resultant mask is named <name_.nii>
if nargin<3
    wdir = 'C:\Data\MNI+masks';
end
aal_mat = spm_read_vols(spm_vol(fullfile(wdir,'aal.nii')));
aals = zeros([size(aal_mat),length(nums)]);
for ii = 1:length(nums)
    aals(:,:,:,ii) = aal_mat == nums(ii);
end
aals(isnan(aals)) = 0;
aals = sum(aals,4);

% Binarization:
aals(aals>0) = 1;
write_reshaped_nifty(aals,wdir,false,'aal.nii',name_)