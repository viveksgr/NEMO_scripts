%% Representational similarity analysis
% RSA analysis to study representation of 
% 1. perceptual vs chemical properties
% 2. odor category vs odor identity.

% Major inputs:
% fullmat.mat: Matrix of voxels x HRF bases x Odors from FIR_model.m
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% chem_corr.mat: correlation matrix of odors based on chemical properties
% Gray matter masks of ROIs
% anat_gw.nii: A binary mask of gray matter voxels such that number of
% voxels in fullmat.mat = sum(anat_gw,'all')
% maskfile: mask of gray matter voxels in the whole brain
% anat_masks: masks of ROIs in native space of subjects.
% fmaskfile: functional mask of voxels with significant odor evoked
% activity (obtained from FIR_model.m, sniff_contrast = true)

% Options:
% dister = To toggle perceptual_vs_chemical (false) or category_vs_identity (true)
% chem_dister_regress = regress chemical similarity out from category and
% identity (true)
% sz_cntrl = correct for size differences in ROIs (true)

% Major outputs:
% p_mat_tt: representational similarity (bootstrapped, and for different
% ROIs) for perceptual properties or odor category
% c_mat_tt: representational similarity (bootstrapped, and for different
% ROIs) for chemical properties or odor identity
% Vivek Sagar (VivekSagar2016@u.northwestern.edu). April 9, 2022

%% Model path
s = 1; % Subject: 1,2,3
modelname = 'S1_test';
root = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\NEMO_scripts\Data'; % Path to the Data folder
respath = fullfile(root,sprintf('RSA_%02d',s),modelname);
mkdir(respath)
statpath = fullfile(root,'supporting_files',sprintf('NEMO_s%02d',s));
maskpath = statpath;

valint = false; % Probe intensity and valence separately
dister = true; % Change correlation to max correlation. For checking coarse vs fine identity representation.
chem_dister_regress = false; % Regress chemical similarity out from coarse vs fine perceptual similarities
nboot = 1000;  % Number of bootstrap samples. 1000 for testing. 10000 samples for final analyses. 
load_boot = false;
permer = false; % Permutation test
grp_ = false; % Group level analysis
task_struct_ctrl = true; % Correct for task structure
spear = true; % Use Spearman rank corr instead of pearson's

% lesions = [];
lesions.post= [1:2]; % Remove intensity and pleasantness for coarse vs fine perceptual similarities
sz_cntrl = false; % Equal num of voxels chosen.
sz_sam = 70; % Choose a minimum of 70 voxels. Only valid if sz_cntrl = true;

ofile = 'Fullmet_FIR'; % Odor responses
norm_ofile = true;
maskfile =  'ARC3_anatgw.nii'; % Gray matter mask
mask = (spm_read_vols(spm_vol(fullfile(statpath, maskfile)))); % Mask used to construct odor files
mask(isnan(mask))=0;
mask(mask<0.1)=0;
mask = logical(mask);

fmaskfile = 'ARC3_fanatgw3_pos.nii'; 
fmask = (spm_read_vols(spm_vol(fullfile(statpath, fmaskfile)))); % Mask used to examine voxels in RSA
fmask(isnan(fmask))=0;
fmask = logical(fmask); % odor responsive voxels (p<0.001)
fmask_1d = fmask(mask);
check_noodor = false;

if check_noodor
    behav_options.only_detect = 0.2; % Remove odors with low detectability
end

behav_file = fullfile(statpath,sprintf('behav_ratings_NEMO%02d.mat',s));
chem_fname = fullfile(statpath,sprintf('chem_corr_NEMO%02d.mat',s));
anat_names = {'PirF','PirT','AMY','OFC'};
anat_masks = {'rwAPC.nii', 'rwPPC.nii','rwAmygdala.nii','rwofc.nii'};

% Gray matter mask
gmmaskfile =  'ARC3_anatgw.nii'; % Must be the same mask used to construct ofile responses
gmmask = (spm_read_vols(spm_vol(fullfile(statpath, gmmaskfile)))); % Mask used to construct odor files
gmmask(isnan(gmmask))=0;
gmmask(gmmask<0.1)=0;
gmmask = logical(gmmask);

% Model names
masks_set = {};
for ii = 1:length(anat_masks)
    m1 = spm_read_vols(spm_vol(fullfile(maskpath,anat_masks{ii})));
    m1(isnan(m1)) = 0;
    m1(m1<=0)=0;
    m1(m1>0) = 1;
     m1 = and(m1,gmmask);
        m1 = m1(mask);
        masks_set{ii}=m1(fmask_1d);
        fprintf('Anat %s: %04d\n',anat_names{ii},sum(m1(fmask_1d)))
    nvoxS = sum(m1(fmask_1d));
end
% masks_set(isnan(masks_set))=0;
linux_config = false;
warning('off','all')

% Behavior analysis settings
behav_options.normalization = true;

% Voxel_responses
Odor_D = fullfile(statpath,ofile);
load(Odor_D,'odor_responses')
if norm_ofile
    odor_responses = zscore(odor_responses,[],3); % Renormalized version of the responses. Only for RSA.
end
anat_idx = [true(1,length(anat_names))]; % Analyze these ROIs
P = [6 6 6; 5 5 6; 5 4 5; 4 4 6]; % Max FIR times for [PirF, PirT, AMY and OFC];
peak_wins_ = P(:,s);

%% Construction of behavioral similarity matrices
% Behav-data
load(behav_file,'behav')
behav.ratings(isnan(behav.ratings))=0;
behav = analyse_behav(behav,behav_options,lesions);

% Load RSMs for P and C
Behav_RSM_P = corrcoef(behav.ratings');

% For special case of intensity and valence studied separately
if size(behav.ratings,2)==1
    val_vec = vs_normalizer(behav.ratings);
    Behav_RSM_P = val_vec.*val_vec';
end

utl_mask = logical(triu(ones(size(Behav_RSM_P)),1));
Behav_RSM_vals_P_2 = Behav_RSM_P(utl_mask); % Perceptual identity
load(chem_fname)

if dister
    B_dist = pdist(behav.ratings,@maxcorrdist);
    B_dist_mat = squareform(B_dist);
    Behav_RSM_vals_P_ = 1-B_dist_mat(utl_mask); % Perceptual category
    Behav_RSM_vals_C_ = Behav_RSM(utl_mask); % Chemical
    if chem_dister_regress
        a_ = regressmeout(Behav_RSM_vals_P_',Behav_RSM_vals_C_')';
        b_ = regressmeout(Behav_RSM_vals_P_2',Behav_RSM_vals_C_')';
    else
        a_ = Behav_RSM_vals_P_;
        b_ = Behav_RSM_vals_P_2;
    end
else
    if valint %If only studying valence or intensity separately
        assert(isempty(lesions)) % Must not remove additional behavioral descriptors
        val_vec = vs_normalizer(behav.ratings(:,1));
        Behav_RSM_P = val_vec.*val_vec';
        val_vec = vs_normalizer(behav.ratings(:,2));
        Behav_RSM_C = val_vec.*val_vec';
        a_ = Behav_RSM_P(utl_mask);
        b_ = Behav_RSM_C(utl_mask);
    else
        b_ = Behav_RSM(utl_mask); % Chemical
        if chem_dister_regress
            a_ = regressmeout(Behav_RSM_vals_P_2',b_')'; % Perceptual
            b_ = regressmeout(b_',Behav_RSM_vals_P_2')'; % Chemical
        else
            a_ = Behav_RSM_vals_P_2;
        end
    end
end
ab = sqrt(abs(a_.*b_)); % Mutual Valence shared between the similarity matrices

if task_struct_ctrl
    load(fullfile(statpath,'task_struct.mat'))
     task_struct = task_struct(utl_mask);
     task_run = task_run(utl_mask);
end

%% Bootstrap/permutation estimation
% Bootstrap CV around peak

% Create or load bootstrap samples
y = 1:length(a_);
if ~permer % Bootstrap
    [~,bootsam] = bootstrp(nboot,@(x) x,y);
else % Shuffle
    bootsam = repmat(y',1,nboot);
end

corr_voxel_final_t1 = zeros(length(anat_names),3);
corr_voxel_final_p1 = zeros(length(anat_names),3);
corr_voxel_final_p2 = zeros(length(anat_names),3);

p_mat_tt = zeros(length(anat_names),size(bootsam,2)); % Hyperparam for perceptual
c_mat_tt = zeros(length(anat_names),size(bootsam,2));
pc_mat_tt = zeros(length(anat_names),size(bootsam,2));

fprintf('\n')
oc_vals = cell(length(find(anat_idx)),1); % Neural similarity matrices for each ROI
for ii = find(anat_idx)
    fprintf('anat area: %02d\n',ii)
%     mask_ind = find(masks_set{ii}); % Only use for size control
    
    for foldid = 1:size(bootsam,2)
        
        foldind = ismember(y,bootsam(:,foldid));
        Odor_mat = squeeze(odor_responses(fmask_1d,peak_wins_(ii),:));
        odor_vals = Odor_mat(logical(masks_set{ii}),:);
        if permer % Shuffle behavioral data for permutation test
            odor_vals = odor_vals(:,randperm(size(odor_vals,2)));
        end
        if sz_cntrl % Sample with replacement a fixed number of voxels for the neural similarity
            mask_ind_trial = datasample(1:size(odor_vals,1),sz_sam);
            odor_vals = odor_vals(mask_ind_trial,:);
        end
        [r,~] = find(isnan(odor_vals));
        odor_vals(r,:) = [];
        
        % Remove trials with low detectability if needed
        if check_noodor
            detect_odorsID = behav.detect>behav_options.only_detect;
            odor_vals = odor_vals(:,detect_odorsID);
        end
        
        odor_corr = corrcoef(odor_vals);
        odor_corr_vals = odor_corr(utl_mask);
        oc_vals{ii} = odor_vals;
               
        if task_struct_ctrl     
            odor_corr_vals = regressmeout(odor_corr_vals', task_struct')';
            odor_corr_vals = regressmeout(odor_corr_vals', task_run')';
        end
        
        if spear
            p_mat_tt(ii,foldid) = corr(odor_corr_vals(foldind),a_(foldind),'Type','Spearman');
            c_mat_tt(ii,foldid) = corr(odor_corr_vals(foldind),b_(foldind),'Type','Spearman');
            pc_mat_tt(ii,foldid) = corr(odor_corr_vals(foldind),ab(foldind),'Type','Spearman');
        else            
            p_mat_tt(ii,foldid) = fastcorr(odor_corr_vals(foldind),a_(foldind));
            c_mat_tt(ii,foldid) = fastcorr(odor_corr_vals(foldind),b_(foldind));
            pc_mat_tt(ii,foldid) = fastcorr(odor_corr_vals(foldind),ab(foldind));
        end
    end
    % For statistics
    corr_voxel_final_t1(ii,1) = mean(p_mat_tt(ii,:));
    corr_voxel_final_t1(ii,2) = mean(c_mat_tt(ii,:));
    corr_voxel_final_t1(ii,3) = mean(pc_mat_tt(ii,:));
    corr_voxel_final_p1(ii,1) = prctile(p_mat_tt(ii,:),97.5);
    corr_voxel_final_p1(ii,2) = prctile(c_mat_tt(ii,:),97.5);
    corr_voxel_final_p1(ii,3) = prctile(pc_mat_tt(ii,:),97.5);
    corr_voxel_final_p2(ii,1) = prctile(p_mat_tt(ii,:),2.5);
    corr_voxel_final_p2(ii,2) = prctile(c_mat_tt(ii,:),2.5);
    corr_voxel_final_p2(ii,3) = prctile(pc_mat_tt(ii,:),2.5);
end

%% Plotting and statistics

corr_voxel_final_t1 = corr_voxel_final_t1(:,[1 2]);
corr_voxel_final_p1 = corr_voxel_final_p1(:,[1 2]);
corr_voxel_final_p2 = corr_voxel_final_p2(:,[1 2]);
% p_mat_tt = p_mat_tt([1 4:8],:);
% c_mat_tt = c_mat_tt([1 4:8],:);
% anat_names = {'wm2','A1+','PirF','PirT','AMY','OFC'};

bar_mean = (corr_voxel_final_t1);
figure('Position',[0.5 0.5 480 320])
hold on
bar(bar_mean)
hold on
ngroups = size(bar_mean, 1);
nbars = size(bar_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
x_m = [];
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_mean(:,i), corr_voxel_final_p1(:,i)-bar_mean(:,i),corr_voxel_final_p2(:,i)-bar_mean(:,i), 'k.');
    %     errorbar(x, bar_mean(:,i), bar_std(:,i)-bar_mean(:,i),bar_std2(:,i)-bar_mean(:,i), 'k.');
    x_m = [x_m; x];
end
xticks(1:length(anat_names))
xticklabels(anat_names);
savefig(fullfile(respath,'ROI_RSM_r2.fig'))
print(fullfile(respath,'ROI_RSM_r2'),'-dsvg')

M_c = [];
M_inv = [];
for jj = 1:length(anat_names) % Num areas
    t1 = p_mat_tt(jj,:);
    t2 = c_mat_tt(jj,:);
    M_c(jj) = bstrap_hyp(t1,t2);
    M_inv(jj,1) = 2*min(invprctile(t1,0),100-invprctile(t1,0))/100; % Distribution is symmetric. 
    M_inv(jj,2) = 2*min(invprctile(t2,0),100-invprctile(t2,0))/100; 
end
clear odor_responses odor_responses_nn bootsam Odor_mat

% Plotting and statistics
M_bar = [];
M_std = [];
M_std2 = [];
for jj = 1:length(anat_names) %
    t1 = c_mat_tt(jj,:);
    t2 = p_mat_tt(jj,:);
    ds = (t1)-(t2);
    M_bar(:,jj) = ds;
    M_std(jj) = prctile(ds,2.5);
    M_std2(jj) = prctile(ds,97.5);
end

% r2 diff
M_comp_mat = [];
ks = 0;
for jj=1:length(anat_names)
    for jj2 = jj+1:length(anat_names)
        ks= ks+1;
        t1 = M_bar(:,jj);
        t2 = M_bar(:,jj2);
        M_comp_mat(ks) = bstrap_hyp(t2,t1);
    end
end

figure('Position',[0.5 0.5 480 320])
M_bar_ = mean(M_bar,1);
bar(1:length(anat_names),M_bar_)
hold on
errorbar(1:length(anat_names),M_bar_,M_std-M_bar_,M_std2-M_bar_,'.')
% Line art per subject
xticklabels(anat_names)
ylabel('Representational r')
savefig(fullfile(respath,'RSA_var'))
print(fullfile(respath,'RSA_var'),'-dsvg')
save(fullfile(respath,'ROI.mat'))

%% Group analysis
% Combine results across subjects
if grp_
    dirs = {fullfile(root,'RSA_01','S1_pc');
            fullfile(root,'RSA_02','S2_pc');
            fullfile(root,'RSA_03','S3_pc')};
    nS = length(dirs);
    matname = 'ROI.mat';
    anat_names = {'PirF','PirT','AMY','OFC'};
    p_mat_tt = variable_extract(dirs,matname,'p_mat_tt',false);
    p_mat_tt = cat(3,p_mat_tt{:});
    c_mat_tt = variable_extract(dirs,matname,'c_mat_tt',false);
    c_mat_tt = cat(3,c_mat_tt{:});
    
    B_mat(:,1,:,:) = p_mat_tt;
    B_mat(:,2,:,:) = c_mat_tt;
    mean_b_mat_ = squeeze(mean(B_mat,3));
    std_b_mat_ = squeeze(std(B_mat,[],3));
    p_mat = squeeze(mean(mean_b_mat_,1))';
    err_mat = squeeze(std(mean_b_mat_,1))';
    
    %--------------------- Mean subject ratings-------------------------------
    mean_b_mat = mean(mean_b_mat_,3);
    
    % Toggle s.e.m bootstrap bars on off. Combine bootstrap error across
    % subjects.
    
    % Confidence intervals
    M = [];
    M2 = [];
    for jj = 1:size(B_mat,1) % Num areas
        t1 = squeeze(B_mat(jj,1,:,:));
        t2 = squeeze(B_mat(jj,2,:,:));
        t1 = mean(t1,2);
        t2 = mean(t2,2);
        M(jj,1) = prctile(t1,97.5);
        M(jj,2) = prctile(t2,97.5);
        M2(jj,1) = prctile(t1,2.5);
        M2(jj,2) = prctile(t2,2.5);
    end
    bar_mean = mean_b_mat;
    bar_std = M; % 97.5 percentile
    bar_std2 = M2; % 2.5 percentile
    
    % RSA plot
    figure('Position',[0.5 0.5 480 320])
    bar(bar_mean)
    hold on
    ngroups = size(bar_mean, 1);
    nbars = size(bar_mean, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    x_m = [];
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, bar_mean(:,i), bar_std(:,i)-bar_mean(:,i),bar_std2(:,i)-bar_mean(:,i), 'k.');
        x_m = [x_m; x];
    end
    xticks(1:length(anat_names))
    xticklabels(anat_names);
    % Subject data points
    c_s = {'r','g','b'}; % Data dots for subjects
    for ii = 1:length(anat_names) % For bars for perceptual, chemical and combinations
        for jj = 1:nS
            plot(x_m(:,ii),mean_b_mat_(ii,:,jj),c_s{jj})
        end
    end
    ylabel('Representational r')
    savefig('anat_pfm')
    print('anat_pfm','-deps')
    
    % Group test for pairwise P:C
    M_c = [];
    M_eff = [];
    for jj = 1:size(B_mat,1) % Num areas
        t1 = squeeze(B_mat(jj,1,:,:));
        t2 = squeeze(B_mat(jj,2,:,:));
        t1 = mean(t1,2);
        t2 = mean(t2,2);
        M_eff(jj) = mean(t1-t2);
        M_c(jj) = bstrap_hyp(t2,t1);
    end
    save(fullfile('ROI_std2.mat'))
    
    %%% - r2 difference
    M_bar = [];
    M_std = [];
    M_std2 = [];
    M_comp = [];
    for jj = 1:size(B_mat,1) % Num areas
        t1 = squeeze(B_mat(jj,2,:,:));
        t2 = squeeze(B_mat(jj,1,:,:));
        t1 = mean(t1,2);
        t2 = mean(t2,2);
        %     ds = ((t1.^2)-(t2.^2))./((t1.^2+t2.^2));
        % ds = (atanh(t1)-atanh(t2))./(atanh(t1)+atanh(t2));
        ds = (t1.^2)-(t2.^2);
        M_bar(:,jj) = ds;
        M_std(jj) = prctile(ds,2.5);
        M_std2(jj) = prctile(ds,97.5);
        if jj>1
            t1 = M_bar(:,jj);
            t2 = M_bar(:,jj-1);
            t_sq = mean(t1)-mean(t2);
            tdist = (t1-t2)-mean(t1-t2);
            M_comp(jj) = (100-invprctile(tdist,t_sq))/100;
        end
    end
    
    % r2 diff
    M_comp_mat = [];
    ks = 0;
    for jj=1:4
        for jj2 = jj+1:4
            ks= ks+1;
            t1 = M_bar(:,jj);
            t2 = M_bar(:,jj2);
            M_comp_mat(ks) = bstrap_hyp(t2,t1);
        end
    end
    [~,M_comp_mat2] =fdr_benjhoc(M_comp_mat);
    
    figure('Position',[0.5 0.5 480 320])
    M_bar_ = mean(M_bar,1);
    bar(1:length(anat_names),M_bar_)
    hold on
    errorbar(1:length(anat_names),M_bar_,M_std-M_bar_,M_std2-M_bar_,'.')
    % Line art per subject
    subwise_r2 = squeeze(mean_b_mat_(:,2,:).^2-mean_b_mat_(:,1,:).^2);
    plot(subwise_r2)
    xticklabels(anat_names)
    ylabel('Representational r^2')
    savefig('RSA_var')
    print(fullfile('RSA_var'),'-deps')
end