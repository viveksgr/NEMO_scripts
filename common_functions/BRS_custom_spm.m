function BRS_custom_spm(substatpath,filer,hrfs,box_sw)
% Run encoding model with settings in filer using hrf specified by hrfs(1) = hrf_max
% (time to peak in sec) and hrfs(2) = hrf_min (time to valley in sec) in folder substatpath.
% Box_sw = use a window of the hrf to ensure the start and end of
% convolution kernel is 0. Default = true.

if nargin<4
    box_sw = true;
end

% Hrf bases function
hrf_kernel = spm_hrf(filer.RT/filer.Mtime_resol, [hrfs(1) hrfs(2) 1 1 6 0 32]);
if box_sw
    boxwind = tukeywin(366,0.2);
    hrf_temp = hrf_kernel.*boxwind;
    hrf_kernel = hrf_temp/sum(hrf_temp);
end

% Create design matrix
[X,~,~] = spm_Volterra(filer.U,hrf_kernel,1); % Design matrix at microtime resolution
X2 = X((0:(filer.nscan(1) - 1))*filer.fMRI_T + filer.fMRI_T0 + 32,:); % Downsample DM          
X3= [X2 filer.X_noise]; % Design matrix 

% Treat design matrix ad 
dlmwrite(fullfile(substatpath,sprintf('regfile%02d%02d.txt',hrfs(1),hrfs(2))),X3)
n = dir(fullfile(substatpath,'regfile*'));
regfile = fullfile(substatpath,n(1).name);

% Model specification
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_spec.dir = {substatpath};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.4;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = filer.m_cell;
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = filer.file_list;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {regfile};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
spm_jobman('run', matlabbatch);

% SPM concatenate
n_volumes = round(filer.nscan/(filer.nruns*filer.nsess));
scans = n_volumes*ones(1,filer.nruns*filer.nsess);
model_loc = fullfile(substatpath,'SPM.mat');
spm_fmri_concatenate(model_loc,scans)

% SPM model estimation
matlabbatch = [];
fname = fullfile(substatpath, 'SPM.mat');
matlabbatch{1}.spm.stats.fmri_est.spmmat =  {fname};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);

% Collate all beta files
collate_betafiles(filer,substatpath,hrfs)
end

% SPM collate beta-files
function collate_betafiles(filer,substatpath,hrfs)
% Collate all the betafiles in substatpath
mask = spm_read_vols(spm_vol(filer.maskfile));

file_list = dir(fullfile(substatpath,sprintf('beta*.nii')));
filename = cell(size(file_list));
for i= 1:length(file_list)
    fname = fullfile(substatpath, file_list(i).name);
    filename{i,1}= sprintf('%s,1',fname); % Datafiles
end
file_V = spm_vol(filename); % File headers for betafiles
file_V = cell2mat(file_V);
file_V = file_V(1:filer.n_odors*filer.nruns,:); % Only choose the ones for odors.

[Y, ~] = spm_read_vols(file_V(1));
dims = size(Y);

odor_responses = zeros([dims filer.n_odors*filer.nruns]); % Zscored
odor_responses_unn = zeros([dims filer.n_odors*filer.nruns]); % Unnormalized

dims1 = dims(1);
dims2 = dims(2);
dims3 = dims(3);

for ii = 1:dims1
    fprintf('.')
    odor_temp = zeros(dims2,dims3,filer.n_odors*filer.nruns);
    odor_temp_unn = zeros(dims2,dims3,filer.n_odors*filer.nruns);
    for jj = 1:dims2
        for kk = 1:dims3
            if mask(ii,jj,kk)~=0
                fprintf('.')
                odor_resp = spm_get_data(file_V,[ii;jj;kk]);
                odor_temp(jj,kk,:) = zscore(odor_resp); % Odor amplitude per voxel 
                odor_temp_unn(jj,kk,:) = odor_resp;
            end
        end
    end
    odor_responses(ii,:,:,:)=odor_temp;
    odor_responses_unn(ii,:,:,:)=odor_temp_unn;
end

% Change 4D to 2D array
mask_1d = logical(reshape(mask,1,[]));
voxel_ids = find(mask_1d);
two_D_ = zeros(size(voxel_ids,1),filer.n_odors*filer.nruns); % Two dimensional representation of odor amplitude
for ii = 1:length(voxel_ids)
    [xx,yy,zz]  = ind2sub(size(mask),voxel_ids(ii)); % Index of the mask
    two_D_(ii,:) = odor_responses(xx,yy,zz,:);
end

% Add third dimension (of bases function used)
odor_responses = zeros(size(two_D_,1),1,size(two_D_,2)); 
odor_responses(:,1,:) = two_D_;
name_chunk = sprintf('full_zscored%02d%02d.mat',hrfs(1),hrfs(2));
save(fullfile(substatpath,name_chunk),'odor_responses','odor_responses_unn')
end


