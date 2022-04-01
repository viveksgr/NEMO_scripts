%% Analyze the prediction accuracy map and break it down by sub regions
% Analyze argmaxes (time-basis) for different brain regions for different descriptors.

% Major inputs:
% Output file ('EEM.mat') of EEM_LOOCV.
% Needs fullmat.mat: Matrix of voxels x HRF bases x Odors
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% Gray matter masks of ROIs
% anat_gw.nii: A binary mask of gray matter voxels such that number of
% voxels in fullmat.mat = sum(anat_gw,'all')

s = 1; % Subject
root = 'C:\Data\NEMO';
anat_names = {'PirF','PirT','AMY','OFC'};
anat_masks = {'rwAPC.nii','rwPPC.nii','rwAmygdala.nii', 'rwOFC.nii'};

nanat = length(anat_names);
anatpath = fullfile('C:\Data\NEMO\',sprintf('NEMO_%02d',s),'\imaging\nii\masks');
maskfile =  'anat_gw.nii'; % Gray matter mask
mask = (spm_read_vols(spm_vol(fullfile(anatpath, maskfile)))); % Mask used to construct odor files
mask(isnan(mask))=0;
mask = logical(mask);
fmaskfile = 'f_anat_gw3.nii'; % Functional mask (only olfactory voxels)
fmask = (spm_read_vols(spm_vol(fullfile(anatpath, fmaskfile)))); % Mask used to examine voxels in RSA
fmask(isnan(fmask))=0;
fmask = logical(fmask);
fmask_1d = fmask(mask);

sr_thresh_ = 0.05; % Choose voxels above this p-value
tscorer = true; % Run analyses on tscores of weights instead of betas
kmediods_ = false; % Use Spearmann correlation on kmediods instead of correlation on kmeans

% Model path
statpath = fullfile(root,sprintf('NEMO_%02d',s),'\imaging\1stlevelmodels\FIR_EM');
nodors = 160;

% Model names
masks_set = [];
for ii = 1:length(anat_masks)
    %     masks_set(:,:,:,ii) = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
    m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
    m1(m1<=0)=0;
    m1(m1>0) = 1;
    m1 = m1(mask);
    masks_set(:,ii)=m1(fmask_1d);
end
masks_set(isnan(masks_set))=0;
% all_masks = sum(masks_set,4);

%% Model comp
mask_scores = cell(length(anat_masks),1);
mask_tw = cell(length(anat_masks));

load(fullfile(statpath,'EEM.mat'),'corr_voxel_3d','r_thresh','train_cell','corr_voxel_final','behav')
scorr_voxel_3d = corr_voxel_3d; % Run the analysis on raw voxels
sr_thresh = r_thresh;%t_thresh;

temp_array_ = scorr_voxel_3d(logical(mask));
temp_array_ = temp_array_(fmask_1d);
figure()
hold on
for nn = 1:length(anat_masks)
    temp_array = temp_array_(logical(masks_set(:,nn)));
    mask_scores{nn} = temp_array;
    subplot(2,3,nn)
    histogram(temp_array)
    title(anat_names{nn})
    xlabel('EEM performance')
end
savefig('hists')

gp = masks_set.*[1:nanat];
gp_id = zeros(size(gp,1),1);
for ii = 1:size(gp,1)
    if sum(gp(ii,:))>0
        gp_id(ii) = find(gp(ii,:)>0,1);
    else
        gp_id(ii) = 0;
    end
end

%% Dimensionality analysis (subject-wise)
% Mask flow:
% "mask" to convert 3d to 1d with nvox  = gray matter voxels across all ROIs.
% fmask_1d to find vox with sig sniff activity
% masks_set or mask_anat to further find vox in specific ROIs
% scorr_mask to find sig voxels in EEM

scorr_voxel_3d = corr_voxel_3d; % Run the analysis on raw voxels
sr_thresh = r2t(sr_thresh_,160);%t_thresh;

% Sig voxesl in the anat areas
sig_train = train_cell(fmask_1d,:,:); % Odor_evoked
mask_anat = logical(sum(masks_set,2));
sig_anat_tr = sig_train(mask_anat,:,:); % Restrict to ROIs

% Make Adj matrix
scorr_mask = scorr_voxel_3d(mask);
scorr_mask = scorr_mask(fmask_1d);
scorr_mask = scorr_mask(mask_anat)>sr_thresh;
train_mat = sig_anat_tr(scorr_mask,:,:); % Training wt matrix

vxl_nodes = mean(train_mat,3);
err_ = std(train_mat,[],3);
vxl_nodes = vxl_nodes./(err_);

% Deg box plots across ROIs
masks_set_mod = masks_set(mask_anat,:);
masks_set_em = masks_set_mod(scorr_mask,:);
gp = masks_set_em;
gp_id2 = zeros(size(gp,1),1);
for ii = 1:size(gp,1)
    if sum(gp(ii,:))>0
        gp_id2(ii) = find(gp(ii,:)>0,1);
    else
        gp_id2(ii) = 0;
    end
end

%% 2(a) Complexity estimation
c = categorical(behav.percepts);
n_cutoff = 10000; % Bootstrap 10000 samples
ndescrip = 12; % Max number of components examined. Must be less than  or equal to 18.
size_control = true; % Control for size differences in the ROIs

p_count = zeros(n_cutoff,length(anat_masks)); % Num PC dimensions in each ROI
var1 = zeros(n_cutoff,length(anat_masks),ndescrip); % Scree plot for PCA in all ROI
coeffs_pc = zeros(n_cutoff,length(anat_masks),ndescrip,ndescrip); % Coefficients for PCA

% Compute statistics for Anova1
weights_roi = zeros(n_cutoff,18,length(anat_masks)); % Only for bootstrap
for nn = 1:length(anat_masks)
    if nn>=1
        temp_array = abs(train_mat(gp_id2==nn,:,:));
    end
    fprintf('Num voxels: %02d\n',size(temp_array,1))
    bars_ = mean(temp_array,3);
    if ~size_control
        if size(temp_array,1)>ndescrip
            [~,bootsam] = bootstrp(n_cutoff,@(x) x,1:size(bars_,1));
        else
            nvox = ndescrip+1;
            bootsam = zeros(nvox,n_cutoff);
            for zz = 1:n_cutoff
                bootsam(:,zz) = datasample(1:size(bars_,1),nvox)';
            end
        end
    else
        nvox = 25;
        bootsam = zeros(nvox,n_cutoff);
        for zz = 1:n_cutoff
            bootsam(:,zz) = datasample(1:size(bars_,1),nvox)';
        end
    end
    for idx = 1:n_cutoff
        bars = bars_(bootsam(:,idx),:);
        if rank(bars)<ndescrip
            % Check this manually. Adjust ncomponents, size of voxels else.
            fprintf('Rank deficiency in bootstrap area:%02d, boot:%02d\n',nn,idx)
        end
        temp_store = temp_array(bootsam(:,idx),:,:); %3D weights of descriptors
        weights_roi(idx,:,nn) = mean(mean(temp_store,3)./std(temp_store,[],3));
        
        % Number of PCs
        [coeff_,~,~,~,var] = pca(bars,'NumComponents',ndescrip);
        %         coeffs_pc(idx,nn,:,:) = coeff_;
        var_vec = cumsum(var);
        if length(var_vec)<ndescrip
            v2 = var_vec(end)*ones(1,ndescrip);
            v2(1:length(var_vec)) = var_vec;
            var_vec = v2;
        end
        var1(idx,nn,:) =  var_vec(1:ndescrip);
        p_count(idx,nn) = find_nearest(var1(idx,nn,:),70,true);
    end
end
save('eem_weights.mat')

%% Complexity plots
% Scree line plot
ndescrip = 12;
var1 = var1(:,:,1:ndescrip);
var_mean = squeeze(mean(var1,1));
var_std = squeeze(std(var1,[],1));
c_s = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
%     if ~size_control
figure('Position',[0.5 0.5 300 200])
hold on
for ii=1:size(var_mean,1)
    shaded_plot(1:size(var_mean,2),var_mean(ii,:),var_std(ii,:),c_s(ii,:));
end
ylim([20 100])
xlim([1 12])
savefig('Scree')
print('Scree','-deps')

% % Scree bar plot
var1 = 2-(var1./50);
var_bar = mean(var1,3);
figure('Position',[0.5 0.5 250 150])
hold on
bar(mean(var_bar))
hold on
errorbar(1:4,mean(var_bar),prctile(var_bar,97.5)-mean(var_bar),prctile(var_bar,2.5)-mean(var_bar),'.')
% for ii  = 1:3
%     plot([1:4],nvoxes_s(:,ii))
% end
xticks(1:4)
xticklabels(anat_names);
% xlim([1 4])
savefig('Scree_bar_cntrl')
print('Scree_bar','-deps')

% Statistics
M = [];
for x1 = 1:4
    for x2 = 1:4
        t1 = var_bar(:,x1);
        t2 = var_bar(:,x2);
        M(x1,x2) = bstrap_hyp_2(t1,t2);
    end
end
save('stats.mat','M')


%% 2(b) Functional maps
figure()
set(gcf, 'Position',  [100, 100, 1200, 200])
hold on
rank_order = []; % Rank of descriptors
rank_count = []; % Num voxels in each ROI
bar_mat = {};
mind_mat = [];
for nn = 1:length(anat_masks)
    temp_array = abs(train_mat(gp_id2==nn,:,:));
    % Plotting
    bars_ = mean(temp_array,3);
    err_ = std(temp_array,[],3);
    bars_ = bars_./(err_);
    % For group plots, change s.e.m to C.I.
    errs = std(bars_)./sqrt(size(bars_,1));
    
    % Anova measurements
    [~,~,stats_] = anova1(bars_,[],'off');
    [c_,m] = multcompare(stats_,'Display','off');
    % Check means (m) with number of descriptors that are different.
    c_sq = squareform(c_(:,6));
    c_sq(c_sq<0.05)=nan;
    c_sq(~isnan(c_sq))=0;
    c_sq(isnan(c_sq))=1;
    c_sq = c_sq-eye(size(c_sq));
    [m_sort, m_ind] = sort(m(:,1),'descend'); % Sorted index of significant percepts
    sig_vectors = c_sq(m_ind(end),:);
    sig_vectors = sig_vectors(m_ind);
    
    bars = mean(bars_);
    subplot(1,4,nn)
    hold on
    bar(1:length(behav.percepts),(bars(m_ind)))
    errorbar(1:length(behav.percepts), bars(m_ind),errs(m_ind),'k.')
    xticks(1:length(behav.percepts))
    xticklabels(behav.percepts(m_ind))
    xtickangle(90)
    title(sprintf('%s: Num Comp: %02d, nVox: %02d',anat_names{nn},p_count(nn),size(temp_array,1)))
    yl = ylim;
    %     sig_vectors
    area(1:length(c),4*yl(2)*sig_vectors,'FaceColor','k','FaceAlpha',0.2)
    ylim(yl)
    [~,rank_order(nn,:)] = sort(m_ind);
    rank_count(nn) = size(temp_array,1);
    bar_mat{nn} = bars_;
    mind_mat = [mind_mat m_ind];
end

if size_control
    save('eem_weights_sz.mat')
else
    save('eem_weights.mat')
end

savefig('anat_percepts')
print('anat_percepts','-dpng')

%% Subject-level visualization -------------------------------------------

% Average EM plots
% Average of prediction accuracy
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg1_tscore'};
matname = 'eem_weights.mat';
varname = 'mask_scores';
labs = {'PirF','PirT','AMY','OFC'};
nvoxes = variable_extract(dirs,matname,varname,false);
nvoxes = horzcat(nvoxes{:});
nanat = 4;
nsub = 3;
nboot = 10000;

% % % Bootstrap means and C.I.
% % cd('C:\Data\NEMO\NEMO_all\EM\EM_pfmance')
% bar_ = zeros(1,nanat); % Mean voxel responses
% std_ = zeros(2,nanat); % Std. dev
% for ii = 1:nanat % anatomical areas
%     temp_dist = zeros(1,nboot);
%     for jj = 1:nsub % num subjects
%         nvoxes{ii,jj}(isnan(nvoxes{ii,jj}))=[];
%         [~,bootsam] = bootstrp(nboot,@(x) x,1:length(nvoxes{ii,jj}));
%         temp_dist = temp_dist+mean(nvoxes{ii,jj}(bootsam));
%     end
%     temp_dist = temp_dist/3;
%     bar_(ii) = mean(temp_dist);
%     std_(1,ii) = std(temp_dist);
% end

% % Bootstrap means and C.I.
% % Percent voxels that are significant
% cd('C:\Data\NEMO\NEMO_all\EM\EM_pfmance')
sthresh = variable_extract(dirs,'None.mat','r_thresh',false);
bar_mat = zeros(nanat,nsub);
for ii = 1:nanat % anatomical areas
    for jj = 1:nsub % num subjects
        nvoxes{ii,jj}(isnan(nvoxes{ii,jj}))=[];
        bar_mat(ii,jj) = sum(nvoxes{ii,jj}>sthresh{jj})./length(nvoxes{ii,jj});
    end
end
bar_mat = bar_mat*100;
bar_ = mean(bar_mat,2);
std_ = 1.96*std(bar_mat,[],2)./sqrt(3);

% % cd('C:\Data\NEMO\NEMO_all\EM\EM_pfmance')
bar_mat = zeros(nanat,nsub);
for ii = 1:nanat % anatomical areas
    for jj = 1:nsub % num subjects
        nvoxes{ii,jj}(isnan(nvoxes{ii,jj}))=[];
        bar_mat(ii,jj) = mean(nvoxes{ii,jj});%(nvoxes{ii,jj}>sthresh{jj}));
    end
end
bar_ = mean(bar_mat,2);
std_ = 1.96*std(bar_mat,[],2)./sqrt(3);

figure('Position',[0.5 0.5 350 250])
bar(bar_)
hold on
errorbar(1:nanat,bar_,std_,'k.')
for jj = 1:3 % For bars for perceptual, chemical and combinations
    plot(bar_mat(:,jj),'MarkerSize',10)
end
xticks(1:nanat)
xticklabels(labs)
xtickangle(90)
savefig('mean_em_pfm')
print('mean_em_pfm','-dpng')

%% Average of k_complexity
% Average of prediction accuracy
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg1_tscore'};
matname = 'eem_weights_boot.mat';
varname = 'var1';
anat_names = {'PirF','PirT','AMY','OFC'};
nvoxes = variable_extract(dirs,matname,varname,false);
ndescrip = 12;

% Remove PirF in S3 
nvoxes{3}(:,1,:)=nan*nvoxes{3}(:,1,:);

nvoxes = cat(4,nvoxes{:});
nvoxes = nvoxes(:,:,1:ndescrip,:);
nvoxes_s = squeeze(nanmean(nvoxes));
nvoxes_s = squeeze(nanmean(nvoxes_s,2));
var1 = nanmean(nvoxes,4);

var_mean = squeeze(nanmean(var1,1));
var_std = squeeze(std(var1,[],1));
c_s = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
figure('Position',[0.5 0.5 300 200])
hold on
for ii=1:size(var_mean,1)
    shaded_plot(1:size(var_mean,2),var_mean(ii,:),var_std(ii,:),c_s(ii,:));
end
ylim([20 100])
xlim([1 18])
savefig('Scree_cntrl')
print('Scree_cntrl','-dpng')

% Scree bar plot
var1 = 2-(var1./50);
nvoxes_s = 2-(nvoxes_s./50);
var_bar = mean(var1,3);
figure('Position',[0.5 0.5 250 150])
hold on
bar(mean(var_bar))
hold on
errorbar(1:4,mean(var_bar),prctile(var_bar,97.5)-mean(var_bar),prctile(var_bar,2.5)-mean(var_bar),'.')
for ii  = 1:3
    plot([1:4],nvoxes_s(:,ii))
end
xticks(1:4)
xticklabels(anat_names);
%     ylim([20 100])
%     xlim([1 4])
    
M = [];
for x1 = 1:4
    for x2 = 1:4
        t1 = var_bar(:,x1);
        t2 = var_bar(:,x2);       
        M(x1,x2) = bstrap_hyp(t1,t2);
    end
end
% [~,~,~,~,var_] = pca(behav.ratings);
% yline(100-mean(cumsum(var_)))
savefig('Scree_bar_cntrl')
print('Scree_bar_cntrl','-dpng')
save('stats.mat','M')
% Run graph codes from eem_weightanalysis.m (2)

%% Group plots of EM ROI bars
% Load behav labs manually
% May not work because of wrong eem_weights
notmerged = true;
swise = 3;
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg1_tscore'};
matname = 'eem_weights_boot.mat';
varname = 'weights_roi';
anat_names = {'PirF','PirT','AMY','OFC'};
nvoxes = variable_extract(dirs,matname,varname,false);

idx1 = [1:18 19 19 19];
idx2 = [1:10 19 11:14 19 15 19 16:18]; % Use this is argsort in sub(2 and 3) to match the labels
if ~swise
ndescrip = length(behav_labs);
nan_mats = nan*ones(10000,1,4);
for ii=1:3
    nvoxes{ii} = cat(2,nvoxes{ii},nan_mats);
end
nvoxes{1} = nvoxes{1}(:,idx1,:);
nvoxes{2} = nvoxes{2}(:,idx2,:);
nvoxes{3} = nvoxes{3}(:,idx2,:);
nvox_mat = cat(4,nvoxes{:});
n_ratings = squeeze(nanmean(nvox_mat,4));
else
    ndescrip = 18;
    n_ratings = nvoxes{swise};
end

figure('Position',[0.5 0.5 900 200])
hold on
sig_vectors = zeros(4,ndescrip);
for ii=1:4
    subplot(1,4,ii)
    hold on
    temp_b = squeeze(n_ratings(:,:,ii));
    [~,argsort] = sort(mean(temp_b),'descend');
    if notmerged % Only running this module and not collapsing plots across areas
        temp_b_sorted = temp_b(:,argsort);
    else
        temp_b_sorted = temp_b;
    end
    [~, argmin] = min(mean(temp_b_sorted));
    % Statistical tests
    %     sig_vectors = zeros(1,ndescrip);
    for jj = setdiff(1:ndescrip,argmin)
        sig_vectors(ii,jj) = double((bstrap_hyp(temp_b_sorted(:,jj),temp_b_sorted(:,argmin)))<0.05/nchoosek(18,2));
    end
    
    bar(1:length(behav_labs),mean(temp_b_sorted))
    errorbar(1:length(behav_labs), mean(temp_b_sorted),std(temp_b_sorted),'k.')
    yl = ylim;
    area(1:length(behav_labs),4*yl(2)*sig_vectors(ii,:),'FaceColor','k','FaceAlpha',0.2)
    xticks(1:length(behav_labs))
    xticklabels(behav_labs(argsort))
    xtickangle(90)
    ylim(yl)
%         savefig('bars_grouped_')
    title(anat_names{ii})
end

%% Group plot of EM ROI bars in single panel
% Rating order
% cd('C:\Data\NEMO\NEMO_all\EM\EM_pfmance\functional_profiles')
% load('fprofile_unsort.mat','n_ratings','behav_labs','anat_names','ndescrip','sig_vectors')
% new_lab_order = [1 6 12 2 11 21 15 14 18 4 3 19 5 10 9 7 16 8 13 17 20];
% new_lab_order = [1 6 14 19 12 15 11 2 21 18 4 3 5 10 9 7 16 8 13 17 20];
% new_lab_order = [1 11 19 16 6 7 2 9 21 8 3 15 18 14 12 5 4 17 10 20 13];
% new_lab_order = [1 11 16 6 7 2 9 8 3 15 18 14 12 5 4 17 10 13];

M_sorter = [-sum(sig_vectors)' -mean(mean(n_ratings,1),3)'];
[~,new_lab_order] = sortrows(M_sorter);

figure('Position',[0.5 0.5 400 250])
hold on
% c_s = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
% Color scheme
c_s = [24,15,66;15,59,92;50,80,35;100,41,16]./100;
c_s1 = [57,55,70;57,75,88;67,80,61;100,76,67]./100;
delta = [-3 3 1 -1]*0.05; % Jitter ratings
for ii=1:4
    %     subplot(4,1,ii)
    hold on
    temp_b = squeeze(n_ratings(:,:,ii));
    %     temp_b(:,~logical(sig_vectors(ii,:)))=nan;
    temp_b_sorted = temp_b(:,new_lab_order);
    %     plot(1:length(behav_labs),mean(temp_b_sorted),'Color', c_s(ii,:),'MarkerSize',5)
    errorbar((1:length(behav_labs))+delta(ii), mean(temp_b_sorted),prctile(temp_b_sorted,97.5)-mean(temp_b_sorted)...
        ,prctile(temp_b_sorted,2.5)-mean(temp_b_sorted),'Color', c_s1(ii,:),'CapSize',1)
    
    temp_b(:,~logical(sig_vectors(ii,:)))=nan;
    temp_b_sorted = temp_b(:,new_lab_order);
    %     plot(1:length(behav_labs),mean(temp_b_sorted),'Color', c_s(ii,:),'MarkerSize',5)
    errorbar((1:length(behav_labs))+delta(ii), mean(temp_b_sorted),prctile(temp_b_sorted,97.5)-mean(temp_b_sorted)...
        ,prctile(temp_b_sorted,2.5)-mean(temp_b_sorted),'Color', c_s(ii,:),'CapSize',1,'LineWidth',1)
    %     s.Color(4) = 0.1;
    
    %     yl = ylim;
    %     area(1:length(behav_labs),4*yl(2)*sig_vectors,'FaceColor','k','FaceAlpha',0.2)
    
    xlim([0 21])
    if ii==4
        xticks(1:length(behav_labs))
        xticklabels(behav_labs(new_lab_order))
        xtickangle(90)
    end
    %     ylim(yl)
    %     savefig('bars_grouped_')
    %     title(anat_names{ii})
end
% legend(anat_names)
ylabel('Regression coefficients (s.n.r.)')


