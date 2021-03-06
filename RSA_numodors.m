%% RSA performance as a function of number of odors
% RSA analysis to study representation of odor category vs identity as a
% function of different number of odors in OFC

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

% options:
% num_od: vector specifying different number of odors

% Major outputs:
% Pred_M: A matrix specifying bootstrapped RSA performance for odor identity in OFC for varying number of odors
% Pred_MC: A matrix specifying bootstrapped RSA performance for odor category in OFC for varying number of odors

% Vivek Sagar (VivekSagar2016@u.northwestern.edu). April 12, 2022

%% General Settings
% Model path
s =1; % Subject
hrf_bases_id = 6; % peak time of HRF in OFC in a given subject
grp = false; % Run group level visualization
num_od = 5:5:160;
statpath = pwd;

ofile = 'fullmat'; % Odor file
norm_ofile = true;
maskfile =  'anat_gw.nii';
mask = (spm_read_vols(spm_vol(fullfile(statpath, maskfile)))); % Mask used to construct odor files
mask(isnan(mask))=0;
% mask(mask<0.1)=0;
mask = logical(mask);
fmaskfile = 'f_anat_gw3.nii';
fmask = (spm_read_vols(spm_vol(fullfile(statpath, fmaskfile)))); % Mask used to examine voxels in RSA
fmask(isnan(fmask))=0;
fmask = logical(fmask);
fmask_1d = fmask(mask);

behav_file = fullfile(statpath,sprintf('behav_ratings_NEMO%02d.mat',s));
behav_filec = fullfile(statpath,sprintf('behav_ratings_NEMO%02dC.mat',s));
anat_names = {'APC','PPC','Amygdala','OFC'};
anat_masks = {'rwAPC.nii', 'rwPPC.nii','rwAmygdala.nii','rwofc.nii'};
anatpath = fullfile('C:\Data\NEMO\',sprintf('NEMO_%02d',s),'\imaging\nii\masks');
nanat = length(anat_names);

% Model names
masks_set = [];
for ii = 1:length(anat_masks)
    m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
    m1(m1<=0)=0;
    m1(m1>0) = 1;
    m1 = m1(mask);
    masks_set(:,ii)=m1(fmask_1d);
    nvoxS(ii) = sum(m1(fmask_1d))./sum(m1);
end
masks_set(isnan(masks_set))=0;
linux_config = false;
warning('off','all')

% Behavior analysis settings
behav_options.normalization = true;
lesions.post = false;

if linux_config
    modelname = '';
    addpath('/home/vsh3681/spm12')
    statpath = fullfile('/home/vsh3681/results', modelname);
    addpath('/home/vsh3681/Scripts/imaging_analysis/common_functions')
    parpool(20)
else
    root = 'C:\Data\NEMO';
    addpath('C:\spm12')
end

% Voxel_responses
Odor_D = fullfile(statpath,ofile);
load(Odor_D,'odor_responses_nn','odor_responses')
odor_responses = odor_responses_nn;

% Behav-data. Perceptual ratings
load(behav_file,'behav')
behav.ratings(isnan(behav.ratings))=0;
behavP = analyse_behav(behav,behav_options,lesions);

%% Boot category vs identity
% Make odor IDs for bootstrap and number analysis

dfs = arrayfun(@(x) nchoosek(x,2),num_od);
nboot = 10000;
num_odors = 160; % Total odors
y = 1:1:num_odors;
oid_mat = cell(1,length(num_od));
for oo = 1:length(num_od)
    oid_temp = zeros(nboot,num_od(oo));
    for zz = 1:nboot
        %         oid_temp(zz,:) = randperm(num_odors,num_od(oo));
        oid_temp(zz,:) = datasample(1:num_odors,num_od(oo));
    end
    oid_mat{oo} = oid_temp;
end

% Predictions
Pred_M = zeros(nanat,nboot,length(num_od)); % Matrix of predictions for odor identity
Pred_MC = zeros(nanat,nboot,length(num_od)); % Matrix of predictions for odor category
test_sig = zeros(nanat,length(num_od));

for oo = 1:length(num_od)
    bootsam = oid_mat{oo}';
    for ii = 4
        Odor_mat = squeeze(odor_responses(fmask_1d,hrf_bases_id,:));
        odor_vals = Odor_mat(logical(masks_set(:,ii)),:);
        [r,~] = find(isnan(odor_vals));
        odor_vals(r,:) = [];
        for foldid = 1:nboot
            foldind = ismember(y,bootsam(:,foldid));
            % RSMs
            Behav_RSM_P = corrcoef(behavP.ratings(foldind,3:end)');
            utl_mask = logical(triu(ones(size(Behav_RSM_P)),1));
            a_ = Behav_RSM_P(utl_mask); % Perceptual RSM
            
            B_dist = pdist(behavP.ratings(foldind,3:end),@maxcorrdist);
            B_dist_mat = squareform(B_dist);
            b_ = 1-B_dist_mat(utl_mask);
            %             Behav_RSM_C = corrcoef(behav.ratings(foldind,:)');
            %             b_ = Behav_RSM_vals_C(utl_mask); % Chemical RSM or Category RSM
            
            odor_vals_temp = odor_vals;
            odor_corr = corrcoef(odor_vals_temp(:,foldind));
            odor_corr_vals = odor_corr(utl_mask); % Neural RSM
            
            Pred_M(ii,foldid,oo) = fastcorr(odor_corr_vals,a_);
            Pred_MC(ii,foldid,oo) = fastcorr(odor_corr_vals,b_);
        end
        test_sig(ii,oo) = bstrap_hyp_2(squeeze(Pred_M(ii,:,oo)),squeeze(Pred_MC(ii,:,oo)));
    end
end

% Final t-score plot for OFC
figure('Position',  [100, 100, 320, 240])
hold on
t1 = tanh(mean(atanh(squeeze(Pred_M(4,:,:))),1));
[~,tsc1] = r2p(t1,dfs);
t2 = tanh(mean(atanh(squeeze(Pred_MC(4,:,:))),1));
[~,tsc2] = r2p(t2,dfs);
plot(num_od,t1.^2-t2.^2,'r')
% plot(num_od,tsc1,'r')
% plot(num_od,tsc2,'b')
yl = ylim;
area(num_od,yl(2)*double(test_sig(4,:)<0.05),'FaceColor','k','FaceAlpha',0.2)
ylabel('Representational r sq')
legend({'identity','category'})
title(anat_names{ii})
savefig('ROI_RSM_odorline_sq.fig')
print(fullfile(statpath,'odorline_ofc'),'-dpng')
save(fullfile(statpath,'odorline.mat'))

%% Average across subjects
% Run previous sections for all subjects
if grp
statpath = pwd;
num_od = 5:5:160;
nanat = 4;
dfs = arrayfun(@(x) nchoosek(x,2),num_od);
nboot = 10000;
num_odors = 160; % Total odors

% Statpaths of various subjects
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\RSA_FIR';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\RSA_FIR'};

matname = 'odorline.mat';
anat_names = {'PirF','PirT','AMY','OFC'};
Pred_M = variable_extract(dirs,matname,'Pred_M',false);
Pred_MC = variable_extract(dirs,matname,'Pred_MC',false);

% Concat matrices across subjects
Pred_M = cat(4,Pred_M{:});
Pred_MC = cat(4,Pred_MC{:});

t1 = squeeze(mean(Pred_M(4,:,:,:),2));
t2 = squeeze(mean(Pred_MC(4,:,:,:),2));
% zsc = (atanh(t1)-atanh(t2))./sqrt(2./(dfs'-3));
for ii = 1:3
    for oo = 1:length(num_od)
        [~,deltar(oo,ii)] = r2p(sqrt(t1(oo,ii).^2-t2(oo,ii).^2),dfs(oo));
    end
end
plot(num_od,t1.^2-t2.^2)

% Concat matrices
Pred_M = mean(Pred_M,4);
Pred_MC = mean(Pred_MC,4);
% Hypothesis testing
test_sig = zeros(nanat,length(num_od));
for oo = 1:length(num_od)
    for ii =1:nanat
        test_sig(ii,oo) = bstrap_hyp(squeeze(Pred_M(ii,:,oo)),squeeze(Pred_MC(ii,:,oo)));
    end
end

% Plot delta t-score
figure('Position',  [100, 100, 320, 240])
hold on
for ii = 4:4
    hold on
    t_vec = test_sig(ii,:); % For t-score
    tsc = tinv(1-t_vec,dfs);
    plot(num_od,tsc,'r')
    hold on
    tsc(test_sig(ii,:)>0.025)=nan;
    plot(num_od,tsc,'b','linewidth',1)
    ylabel('Representational t-value')
    xlabel('No. of odors')
    title(anat_names{ii})
end

savefig('ROI_RSM_odorline_allareas.fig')
print(fullfile(statpath,'odorline_allareas'),'-dpng')
end