%% Dimensionality analysis of the Encoding Model
% Perform basic visualization and comnpute dimensionality of encoding of
% the encoding model.

% Major inputs:
% EEM.mat = Output file of EEM_LOOCV. Specifies prediction accuracies and
% training weights assigned to voxels.
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% Gray matter masks of ROIs
% anat_gw.nii: A binary mask of gray matter voxels such that number of
% voxels in fullmat.mat = sum(anat_gw,'all')

%% General settings
grp_ = false; % Group level visualization (Combine all subjects)
% Disregard PirF in S3
s = 3; % Subject 

modelname = 'S1_EM_LOOCV';
root = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\NEMO_scripts\Data'; % Path to the Data folder
respath = fullfile(root,sprintf('EM_%02d',s),modelname);
statpath = fullfile(root,'supporting_files',sprintf('NEMO_s%02d',s));
maskpath = statpath;

% % Old repository
% root = 'C:\Data\NEMO';
% statpath = fullfile(root,sprintf('NEMO_%02d',s),'\imaging\1stlevelmodels\FIR_EM\deg_tscore_new2');
% respath = statpath;
% maskpath = fullfile('C:\Data\NEMO\',sprintf('NEMO_%02d',s),'\imaging\nii\masks');

anat_names = {'PirF','PirT','AMY','OFC'};
anat_masks = {'rwAPC.nii', 'rwPPC.nii','rwAmygdala.nii','rwofc.nii'};
nanat = length(anat_names);
maskfile =  'ARC3_anatgw.nii'; % Gray matter mask
mask = (spm_read_vols(spm_vol(fullfile(maskpath, maskfile)))); % Mask used to construct odor files
mask(isnan(mask))=0;
mask = logical(mask);
fmaskfile = 'ARC3_fanatgw3_pos.nii'; % Functional mask (only olfactory voxels)
fmask = (spm_read_vols(spm_vol(fullfile(maskpath, fmaskfile)))); % Mask used to examine voxels in RSA
fmask(isnan(fmask))=0;
fmask = logical(fmask);
fmask_1d = fmask(mask);
sr_thresh_ = 0.05; % Choose voxels above this p-value
nodors = 160;

% Model names
masks_set = [];
for ii = 1:length(anat_masks)
    %     masks_set(:,:,:,ii) = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
    m1 = spm_read_vols(spm_vol(fullfile(maskpath,anat_masks{ii})));
    m1(m1<=0)=0;
    m1(m1>0) = 1;
    m1 = m1(mask);
    masks_set(:,ii)=m1(fmask_1d);
end
masks_set(isnan(masks_set))=0;
% all_masks = sum(masks_set,4);
fname = 'EEM.mat';
%% Model comp
% Visualize subjectwise model performance
mask_scores = cell(length(anat_masks),1);
mask_tw = cell(length(anat_masks));
load(fullfile(respath,fname),'corr_voxel_3d','r_thresh','train_cell','corr_voxel_final','behav')
scorr_voxel_3d = corr_voxel_3d; % Run the analysis on raw voxels

temp_array_ = scorr_voxel_3d(logical(mask));
temp_array_ = temp_array_(fmask_1d);
% Visualize EEM performance
figure()
hold on
for nn = 1:length(anat_masks)
    temp_array = temp_array_(logical(masks_set(:,nn)));
    mask_scores{nn} = temp_array;
    subplot(2,2,nn)
    histogram(temp_array)
    title(anat_names{nn})
    xlabel('EEM performance')
end
savefig(fullfile(respath,'hists'))

gp = masks_set.*[1:nanat];
gp_id = zeros(size(gp,1),1);
for ii = 1:size(gp,1)
    if sum(gp(ii,:))>0
        gp_id(ii) = find(gp(ii,:)>0,1);
    else
        gp_id(ii) = 0;
    end
end

%% Dimensionality analysis (subject-wise) - general settings
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

%% Dimensionality estimation
c = categorical(behav.percepts);
n_cutoff = 10000; % Bootstrap 10000 samples. 1000 for test-run
ndescrip = 12; % Max number of components examined. Must be less than  or equal to 18.
size_control = false; % Control for size differences in the ROIs

p_count = zeros(n_cutoff,length(anat_masks)); % Num PC dimensions in each ROI
var1 = zeros(n_cutoff,length(anat_masks),ndescrip); % Scree plot for PCA in all ROI
coeffs_pc = zeros(n_cutoff,length(anat_masks),ndescrip,ndescrip); % Coefficients for PCA

% Compute statistics for Anova1
weights_roi = zeros(n_cutoff,18,length(anat_masks)); % Only for bootstrap
bar_mat = {};
% if s==3, nn_i = 2; else nn_i = 1; end % Exclude PirF in S3
for nn = 1:length(anat_masks)
    fprintf('Dim estimation for anat area: %02d\n',nn)
    if nn>=1
        temp_array = abs(train_mat(gp_id2==nn,:,:));
        bars_ = mean(temp_array,3);
        err_ = std(temp_array,[],3);
        bars_ = bars_./(err_);
        bar_mat{nn} = bars_;
    end
    
%     fprintf('Num voxels: %02d\n',size(temp_array,1))
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
        nvox = 25; % Default size for dimensionality analysis
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
%         p_count(idx,nn) = find_nearest(var1(idx,nn,:),70,true);
    end
end


%% Dimensionality plots (subjectwise)
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
savefig(fullfile(respath,'Scree'))
% print(fullfile(respath,'Scree'),'-deps')

% % Scree bar plot
var0 = 2-(var1./50); % Renormalize to 0-1. (Area_max - Area)/(Area_max - Area_min) = (1800-18(mean_area))/(1800-900)
var_bar = mean(var0,3);
if s==3; var_bar(:,1) = 0; end % Exclude PirF S3
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
savefig(fullfile(respath,'Scree_bar'))
% print(fullfile(respath,'Scree_bar'),'-deps')

% % Statistics
% M = [];
% for x1 = 1:4
%     for x2 = 1:4
%         t1 = var_bar(:,x1);
%         t2 = var_bar(:,x2);
%         M(x1,x2) = bstrap_hyp(t1,t2);
%     end
% end
save(fullfile(respath,'eem_weights_boot.mat'),'bar_mat','mask_scores','var1','weights_roi')


%% Subject-level visualization 
% Perform group level visualization and combine results across subjects
if grp_
% Average EM performance plots
% Directories where results from encoding models and eem_dimensionality are stored in each
% subject

root = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\NEMO_scripts\Data'; % Path to the Data folder
respath = fullfile(root,'EM_all');
mkdir(respath)
dirs = {fullfile(root,'EM_01','S1_EM_LOOCV');
        fullfile(root,'EM_02','S2_EM_LOOCV');
        fullfile(root,'EM_03','S3_EM_LOOCV')};
matname = 'eem_weights_boot.mat';
varname = 'mask_scores';
labs = {'PirF','PirT','AMY','OFC'};
nvoxes = variable_extract(dirs,matname,varname,false);
nvoxes = horzcat(nvoxes{:});
nanat = 4;
nsub = 3;
nboot = 10000;

sthresh = variable_extract(dirs,'EEM.mat','r_thresh',false);

% % Voxels that are significant
% bar_mat = zeros(nanat,nsub);
% for ii = 1:nanat % anatomical areas
%     for jj = 1:nsub % num subjects
%         nvoxes{ii,jj}(isnan(nvoxes{ii,jj}))=[];
%         bar_mat(ii,jj) = sum(nvoxes{ii,jj}>sthresh{jj})./length(nvoxes{ii,jj});
%     end
% end
% bar_mat = bar_mat*100;
% bar_ = mean(bar_mat,2);
% std_ = 1.96*std(bar_mat,[],2)./sqrt(3);

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
savefig(fullfile(respath,'mean_em_pfm'))

%% Average of dimensionality across subjects
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
xlim([1 12])
savefig(fullfile(respath,'Scree_cntrl'))

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

% Significance testing
M = [];
for x1 = 1:4
    for x2 = 1:4
        t1 = var_bar(:,x1);
        t2 = var_bar(:,x2);       
        M(x1,x2) = bstrap_hyp(t1,t2);
    end
end
savefig(fullfile(respath,'Scree_bar_cntrl'))

%% Group plots of EM ROI bars
% Plot encoding model weights in different areas
% Combine descriptors in all subjects
% May need sorting for visualization
behav_labs = {'Intensity';' Pleasantness';' Fishy';' Burnt';' Sour';' Decayed';...
    ' Musky';' Fruity';' Sweaty';' Cool';' Chemical';'Floral';' Sweet';' Warm';...
    ' Bakery';' Garlic';' Spicy';' Acidic';'Ammonia';'Edible';'Familiar'};

matname = 'eem_weights_boot.mat';
varname = 'weights_roi';
anat_names = {'PirF','PirT','AMY','OFC'};
nvoxes = variable_extract(dirs,matname,varname,false);

idx1 = [1:18 19 19 19];
idx2 = [1:10 19 11:14 19 15 19 16:18]; % Use this is argsort in sub(2 and 3) to match the labels
ndescrip = length(behav_labs);
nan_mats = nan*ones(10000,1,4);
for ii=1:3
    nvoxes{ii} = cat(2,nvoxes{ii},nan_mats);
end
nvoxes{1} = nvoxes{1}(:,idx1,:);
nvoxes{2} = nvoxes{2}(:,idx2,:);
nvoxes{3} = nvoxes{3}(:,idx2,:);

% remove S3 PirF
nvox_mat = cat(4,nvoxes{:});

n_ratings = squeeze(nanmean(nvox_mat,4));
% Statistical significance for training weights in the anatomical areas
sig_vectors = zeros(4,ndescrip);
for ii=1:4
    temp_b = squeeze(n_ratings(:,:,ii));
    [~,argsort] = sort(mean(temp_b),'descend'); 
    temp_b_sorted = temp_b;
    [~, argmin] = min(mean(temp_b_sorted));
    % Statistical tests
    for jj = setdiff(1:ndescrip,argmin)
        sig_vectors(ii,jj) = double((bstrap_hyp(temp_b_sorted(:,jj),temp_b_sorted(:,argmin)))<0.025/nchoosek(18,2));
    end
end
M_sorter = [-sum(sig_vectors)' -mean(mean(n_ratings,1),3)']; % Sort descriptors
[~,new_lab_order] = sortrows(M_sorter);

% EEM weights in anatomical areas
figure('Position',[0.5 0.5 400 250])
hold on
c_s = [24,15,66;15,59,92;50,80,35;100,41,16]./100;
c_s1 = [57,55,70;57,75,88;67,80,61;100,76,67]./100;

delta = [-3 3 1 -1]*0.05; % Jitter ratings
for ii=1:4
    hold on
    temp_b = squeeze(n_ratings(:,:,ii));
    temp_b_sorted = temp_b(:,new_lab_order);
    errorbar((1:length(behav_labs))+delta(ii), mean(temp_b_sorted),prctile(temp_b_sorted,97.5)-mean(temp_b_sorted)...
        ,prctile(temp_b_sorted,2.5)-mean(temp_b_sorted),'Color', c_s1(ii,:),'CapSize',1)
    
    temp_b(:,~logical(sig_vectors(ii,:)))=nan;
    temp_b_sorted = temp_b(:,new_lab_order);
    errorbar((1:length(behav_labs))+delta(ii), mean(temp_b_sorted),prctile(temp_b_sorted,97.5)-mean(temp_b_sorted)...
        ,prctile(temp_b_sorted,2.5)-mean(temp_b_sorted),'Color', c_s(ii,:),'CapSize',1,'LineWidth',1)   
    xlim([0 21])
    if ii==4
        xticks(1:length(behav_labs))
        xticklabels(behav_labs(new_lab_order))
        xtickangle(90)
    end
end
ylabel('Regression coefficients (s.n.r.)')
end

