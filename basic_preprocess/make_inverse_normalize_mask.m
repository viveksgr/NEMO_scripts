% Construction of anatomical masks in subject's native space based on MNI template.
% Masks need to be further optimized manually.

% Path and files etc.
root = 'C:\Data\NEMO\NEMO_02'; % For functional files
% addpath('C:\spm12\')

% Mask image in normalized space to be inverse nomralized into native space
maskpath = 'C:\Data\NEMO\NEMO_02\imaging\nii\masks';
maskfile = 'FFA.nii'; % MNI mask
anatpath = 'C:\Data\NEMO\NEMO_02\imaging\nii\anat'; % Path of anatomical T1 image

% anatomical in native space
n = dir(fullfile(anatpath, 'sNEMO02.nii'));
anatfile = fullfile(anatpath, n(1).name);
n = dir(fullfile(anatpath, 'y_sNEMO02.nii'));
deffile = fullfile(anatpath, n(1).name);
% deformation field for normalization 

% Create inverse normalization parameters
cd(anatpath); 
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {deffile};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {anatfile};
matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'inverse.nii';
matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.savepwd = 1;
spm_jobman('run', matlabbatch)

% write inverse normalized masks
matlabbatch = [];
invdeffile = fullfile(anatpath, 'y_inverse.nii');
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {invdeffile};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {fullfile(maskpath, maskfile)};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -85
    %        78 76 85];
    78 95 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
spm_jobman('run', matlabbatch)

% % copy written file to anat folder
source = fullfile(maskpath, sprintf('w%s', maskfile));
% target = fullfile(anatpath, sprintf('w%s', maskfile));
% movefile(source,target)

% Reslicing the estimated files in the T1 space
matlabbatch = [];
datapath = fullfile(root, 'imaging', 'nii', 'set_01', 'sess_02', 'run_01');
n = dir(fullfile(datapath, 'rf*.nii'));
fimages{1,1} = sprintf('%s,1', fullfile(datapath, n(1).name));
% fimages{1,1} = sprintf('%s,1', anatfile);
matlabbatch{1}.spm.spatial.coreg.write.ref = fimages(1,1);
matlabbatch{1}.spm.spatial.coreg.write.source = {source};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run', matlabbatch)
