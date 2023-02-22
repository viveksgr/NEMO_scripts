%% The Encoding Model
% Construct encoding model to predict neural responses as a linear function
% of behavioral descriptors

% Inputs:
% Fullmet_FIR.mat: Matrix of voxels x HRF bases x Odors
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% argmaxes2.mat: If not fully crossvalidating, the index of which hrf bases to choose. 
% ARC3_anatgw.nii: A binary mask of gray matter voxels such that number of
% voxels in Fullmet_FIR.mat = sum(ARC3_anatgw,'all')

% Outputs:
% corr_voxel_3d: A 3d voxel map in subject's native space specifying the prediction accuracy of the
% model.
% train_cell: Matrix specifying training weights for a voxel, given perceptual feature, and fold of crossvalidation

% Vivek Sagar (VivekSagar2016@u.northwestern.edu). August 31, 2022

%% General Settings
s = 3; % Subject
modelname = 'S2_EM_LOOCV';
root = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\NEMO_scripts\Data'; % Path to the Data folder
respath = fullfile(root,sprintf('EM_%02d',s),modelname);
mkdir(respath)
statpath = fullfile(root,'supporting_files',sprintf('NEMO_s%02d',s));

mask_anat_name = 'EEM'; % 'EEM' for no mask. Or use 'APC','PPC','Amygdala','OFC' for restricting analysis to a single ROI
linux_config = false;
warning('off','all') 

% Behavior analysis settings.
fname_behav = sprintf('behav_ratings_NEMO%02d',s); % C for chemical and 'P' for perceptual basis
behav_options = [];
% behav_options.pca_anal = 0.4;
% behav_options.intercept = true; % Turn off the intercept with ridge regression
% behav_options.orthog = true;
% behav_options.normalization = true;% behav_options.polybase = 2;
% behav_options.fourier = 3;
% behav_options.regress_orth = 2;
% behav_options.binarize = [0 1]; % binarize(1) = threshold of binarizing. binarize(2) = first descriptor to binarize.                                      % = 3 for s_1 and 5 for s_2               
% behav_options.rem_low = true; 

normalize_vox = true;
lesions = false;
timing_map = false;

% Visualization settings
read_mask = false; % Read functional mask from a previous model instead of writing a new one.
accu = true;

% Regularization settings
regul = 2; % 0 for no regression, 1 for lasso, 1.5 for elastic net, 2 for ridge, 2.5 for tikhonov
lambda = 0.05; % Specify lambda or alpha for elastic net if no regul
cross_vald_regul = false; % If cross validating regularization
krange = [0.001 0.01, 0.05, 0.2];%linspace(0.01,0.4,10); % Lambda values. Alpha values in elastic net are always 0.1:0.1:1. Change in function.

% General
% File path
if linux_config
    addpath('/home/vsh3681/spm12')
    statpath = fullfile('/home/vsh3681/results', modelname);
    addpath('/home/vsh3681/Scripts/imaging_analysis/common_functions')
    parpool(20)
else
    sn = sprintf('NEMO_%02d', s);
    root = 'C:\Data\NEMO';
    addpath('C:\spm12')
end 

% Voxel_responses
Odor_D = fullfile(statpath,'Fullmet_FIR.mat');
load(Odor_D,'odor_responses')

behav_file = fullfile(statpath,fname_behav);
load(behav_file,'behav')
behav.ratings(isnan(behav.ratings))=0;
behav = analyse_behav(behav,behav_options,lesions);

% Voxel and anat mask
mask = spm_read_vols(spm_vol(fullfile(statpath, 'ARC3_anatgw.nii')));
mask(isnan(mask))=0;
mask_pred = logical(reshape(mask,1,[]));
if ~strcmp(mask_anat_name,'EEM') % Use an anat mask
    mask_anat_full = sprintf('%s.nii',mask_anat_name);
    mask_anat = spm_read_vols(spm_vol(fullfile(statpath,mask_anat_full)));
    mask_anat_stack = logical(reshape(mask_anat,1,[]));
    
 
    mask_pred_anat = mask_anat_stack(mask_pred);
    odor_responses = odor_responses(mask_pred_anat,:,:);

else
    mask_pred_anat = mask_pred;
end

% Initialize all matrices to save results
nfolds = size(behav.ratings,1);
nvox = size(odor_responses,1);
ncomp = size(behav.ratings,2);

test_cell = zeros(nvox,nfolds);
pred_cell = zeros(nvox,nfolds);
train_cell = zeros(nvox,ncomp,nfolds); % These are training weights
load(fullfile(statpath,'argmaxes2.mat'))
argmaxes = argmaxes{1}; % Choose these time-bins from odor responses
argmaxes(argmaxes==0) = 6;

%% Folds for LOOCV
lambda_mat = zeros(sum(mask_pred_anat),nfolds);
fold_ind_inner = crossvalind('Kfold',nfolds-1,4);

odor_resp_fold = zeros(size(odor_responses,1),size(odor_responses,3));
for ii = 1:size(odor_responses,1)
    odor_resp_fold(ii,:)=squeeze(odor_responses(ii,argmaxes(ii),:));
end

fold_ind = 1:nfolds;
for fold_ = fold_ind
    %% Training and weight estimation
    fprintf('Running fold %d \n', fold_);
    test_mask = fold_ind==fold_; % fold_ind comes from cv_fold.mat
    behav_tt = behav.ratings(test_mask,:);
     
    % Validation step - estimate the argmax across second dim for maximum
    % prediction accuracy.
    % --------------------------------------------------------------------
    tr_sets = setdiff(unique(fold_ind),fold_);
  
    odor_resp_tt = odor_resp_fold(:,test_mask);
    odor_resp_tr = odor_resp_fold(:,~test_mask);
    behav_tr = behav.ratings(~test_mask,:);
    test_cell(:,fold_) =  odor_resp_tt;
    
    if normalize_vox
        mean_tr = nanmean(odor_resp_tr,2);
        std_tr = nanstd(odor_resp_tr,[],2);
        odor_resp_tr = (odor_resp_tr-repmat(mean_tr,1,size(odor_resp_tr,2)))./repmat(std_tr,1,size(odor_resp_tr,2)); % Broadcasting not available for old version of MATLAB on quest
        odor_resp_tt = (odor_resp_tt-repmat(mean_tr,1,size(odor_resp_tt,2)))./repmat(std_tr,1,size(odor_resp_tt,2));
    end
    
    switch regul
        case 2 % Ridge regularization
            for ii = 1:size(odor_resp_tr,1)
                if range(odor_resp_tr(ii,:)')>0
                if cross_vald_regul
                    lambda = max_lambda(odor_resp_tr(ii,:)',behav_tr,fold_ind_inner,krange,2);
                    lambda_mat(ii,fold_) = lambda;
                end
                train_cell(ii,:,fold_) = ridge(odor_resp_tr(ii,:)',behav_tr,lambda);
                pred_cell(ii,fold_) = squeeze(train_cell(ii,:,fold_))*behav_tt';
                end
            end
    end
end

%% Visualization
corr_voxel_final = iter_corr(pred_cell,test_cell);
mask_header = spm_vol(fullfile(statpath, 'ARC3_anatgw.nii'));
[mask,XYZmm] = spm_read_vols(mask_header);
mask(isnan(mask))=0;
XYZvx = round(mask_header.mat\[XYZmm; ones(1,size(XYZmm,2))]);
XYZvx(end,:)=[];
voxel_id = XYZvx(:,mask_pred);
sizer_ = size(mask);
corr_voxel_3d = nii_reshaper(corr_voxel_final,voxel_id,sizer_);
corr_voxel_3d(~logical(mask))=nan;
scorr_voxel_3d = nansmooth3(corr_voxel_3d);
scorr_voxel_3d = atanh(scorr_voxel_3d); % Smoothed map for visualization

% Correct for multiple comparisons and store p_value
corr_pval = corr_voxel_3d(logical(mask));
for ii = 1:length(corr_pval)
    corr_pval(ii) = r2p(corr_pval(ii),size(behav.ratings,1));
end
% p_thresh = fdr_benjhoc(corr_pval);
% r_thresh = atanh(r2t(p_thresh,size(behav.ratings,1)));
r_thresh = [0.1771, 0.1829, 0.2079]; % Manually put r_thresh computed from a previous run of this script in which
% FDR correction is implemented over the whole brain and not just olfactory areas.
r_thresh = r_thresh(s);
fprintf('Threshold for 0.05 fdr: %f', r_thresh);

save(fullfile(respath,mask_anat_name))
write_reshaped_nifty(corr_voxel_3d,statpath,false,'ARC3_anatgw.nii'); % Actual correlation
movefile(fullfile(statpath,'corr_voxel_3d.nii'),respath)
write_reshaped_nifty(scorr_voxel_3d,statpath,false,'ARC3_anatgw.nii'); % Smoothened
movefile(fullfile(statpath,'scorr_voxel_3d.nii'),respath)
