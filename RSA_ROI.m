%% General Settings

% Major inputs:
% Needs fullzcored.mat: Matrix of voxels x HRF bases x Odors
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% Gray matter masks of ROIs
% anat_gw.nii: A binary mask of gray matter voxels such that number of
% voxels in fullzscored.mat = sum(anat_gw,'all')


% Model path
s =4; % Subject
modelname = 'classic_chem';

statpath = pwd;
dister = false; % Change correlation to max correlation. For checking odor category vs identity representation.
check_pos = false; % Is correlation driven by pos values?
chem_dister_regress = false;
nboot = 10000;
load_boot = true;

lesions = false; %.post= [1 2];
sz_cntrl = false; % Equal num of voxels chosen.
sz_sam = 35; % Choose a minimum of 35 voxels
load_normals = false; % Load orthogonalized bases rather than calculating them
use_normals = false; % Only for load_normals = false

b_r_t = 'b'; % 'b' for standardized coefficients, t = tscore, r = partial r sq.
ofile = 'full_zscored_11basis_s'; % Odor file
norm_ofile = true;
hrf_ids = [1 1 1 1];
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
chem_fname = 'chem_corr.mat';

behav_file = fullfile(statpath,sprintf('behav_ratings_NEMO%02d.mat',s));
anat_names = {'APC','PPC','Amygdala','OFC'};
anat_masks = {'rwAPC.nii', 'rwPPC.nii','rwAmygdala.nii','rwofc.nii'};
anatpath = fullfile('C:\Data\NEMO\',sprintf('NEMO_%02d',s),'\imaging\nii\masks');

% Model names
masks_set = [];
for ii = 1:length(anat_masks)
    m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
    m1(m1<=0)=0;
    m1(m1>0) = 1;
    m1 = m1(mask);
    masks_set(:,ii)=m1(fmask_1d);
    nvoxS = sum(m1(fmask_1d));
end
masks_set(isnan(masks_set))=0;
linux_config = false;
warning('off','all')

% Behavior analysis settings
chem_or_per = 'P'; % C for chemical and 'P' for perceptual basis
behav_options.normalization = true;
hrf_bases_id = hrf_ids(s); % 16 for NEMO01 and 1 for NEMO02 for canonical hrf

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
if norm_ofile
odor_responses = odor_responses_nn;
end

% Behav-data
load(behav_file,'behav')
behav.ratings(isnan(behav.ratings))=0;
behav = analyse_behav(behav,behav_options,lesions);

% Load RSMs for P and C
Behav_RSM_P = corrcoef(behav.ratings');
utl_mask = logical(triu(ones(size(Behav_RSM_P)),1));
Behav_RSM_vals_P = Behav_RSM_P(utl_mask);
load(chem_fname)

if dister
    B_dist = pdist(behav.ratings,@maxcorrdist);
    B_dist_mat = squareform(B_dist);
    Behav_RSM_vals_P_ = 1-B_dist_mat(utl_mask); 
    
    if check_pos
        Behav_RSM_P= pdist(behav.ratings,@meancorrdist);
        Behav_RSM_P = squareform(Behav_RSM_P);
        Behav_RSM_vals_P_2 = 1-Behav_RSM_P(utl_mask);
    else            
        Behav_RSM_vals_P_2 = Behav_RSM_P(utl_mask);
    end

    Behav_RSM_vals_C_ = Behav_RSM(utl_mask);
    if chem_dister_regress
        Behav_RSM_vals_P = regressmeout(Behav_RSM_vals_P_',Behav_RSM_vals_C_')';
        Behav_RSM_vals_C = regressmeout(Behav_RSM_vals_P_2',Behav_RSM_vals_C_')';
    else
        Behav_RSM_vals_P = Behav_RSM_vals_P_;
        Behav_RSM_vals_C = Behav_RSM_vals_P_2;       
    end
else
    Behav_RSM_vals_C = Behav_RSM(utl_mask);
    if chem_dister_regress
        Behav_RSM_vals_P = regressmeout(Behav_RSM_vals_P',Behav_RSM_vals_C')';
        Behav_RSM_vals_C = regressmeout(Behav_RSM_vals_C',Behav_RSM_vals_P')';
    end
end

% Bootstrap CV around peak
anat_idx = [true(1,4)]; % Analyze these ROIs
peak_wins = 1;
P = [6 6 6 6; 5 5 6 6; 5 4 6 6; 4 4 6 6];
 peak_wins_ = P(:,s);
% peak_wins = 1:3;
% peak_wins_ = [4:6; 4:6; 4:6; 4:6];

% peak_wins_ = [6; 5; 5; 4];
% Define train-test split at x%
if load_normals
    load(fullfile(statpath,'var_partitioner.mat'));
else
    Behav_RSM_vals_P = Behav_RSM_vals_P-mean(Behav_RSM_vals_P);
    Behav_RSM_vals_P = Behav_RSM_vals_P/norm(Behav_RSM_vals_P);
    Behav_RSM_vals_C = Behav_RSM_vals_C-mean(Behav_RSM_vals_C);
    Behav_RSM_vals_C = Behav_RSM_vals_C/norm(Behav_RSM_vals_C);
    if use_normals
        for ii=1:1000
            [a_,b_,ab,t] = var_partitioner(Behav_RSM_vals_P,Behav_RSM_vals_C);
            if min(fastcorr(a_,ab),fastcorr(b_,ab))<r2t(0.05,12720)
                [a_,b_,ab,t] = var_partitioner(Behav_RSM_vals_P,Behav_RSM_vals_C,t,true);
                break
            end
            if ii==1000
                fprintf('None found');
            end
        end
        save('var_partitioner.mat','a_','b_','ab')
    else
        a_ = Behav_RSM_vals_P;
        b_ = Behav_RSM_vals_C;
        ab = sqrt(abs(Behav_RSM_vals_P.*Behav_RSM_vals_C));
    end
end

y = 1:length(Behav_RSM_vals_P);
if load_boot
    load('bootsam.mat')
else
   [~,bootsam] = bootstrp(nboot,@(x) x,y); 
end

% Backup RSA
corr_voxel_final_t1 = zeros(length(anat_names),3);
corr_voxel_final_p1 = zeros(length(anat_names),3);
corr_voxel_final_p2 = zeros(length(anat_names),3);

p_mat_tr = zeros(length(anat_names),size(bootsam,2)); % Hyperparam for perceptual
c_mat_tr = zeros(length(anat_names),size(bootsam,2));
pc_mat_tr = zeros(length(anat_names),size(bootsam,2));

p_mat_tt = zeros(length(anat_names),size(bootsam,2)); % Hyperparam for perceptual
c_mat_tt = zeros(length(anat_names),size(bootsam,2));
pc_mat_tt = zeros(length(anat_names),size(bootsam,2));

fprintf('\n')
oc_vals = cell(length(find(anat_idx)),length(peak_wins));
for ii = find(anat_idx)
    fprintf('anat area: %02d\n',ii)
    mask_ind = find(masks_set(:,ii)); % Only use for size control    
    for foldid = 1:size(bootsam,2)
        foldind = ismember(y,bootsam(:,foldid));
        p_score_tr = zeros(3,length(peak_wins));        
        p_score_tt = zeros(3,length(peak_wins));
        for hrf_bases_id = 1:length(peak_wins)
            if exist('peak_wins_','var')
                Odor_mat = squeeze(odor_responses(fmask_1d,peak_wins_(ii,hrf_bases_id),:));
            else
                Odor_mat = squeeze(odor_responses(fmask_1d,peak_wins(hrf_bases_id),:));
            end
            odor_vals = Odor_mat(logical(masks_set(:,ii)),:);
            
            if sz_cntrl
                mask_ind_trial = datasample(1:size(odor_vals,1),sz_sam);
                odor_vals = odor_vals(mask_ind_trial,:);
            end           
            [r,~] = find(isnan(odor_vals));
            odor_vals(r,:) = [];
            odor_corr = corrcoef(odor_vals);
            odor_corr_vals = odor_corr(utl_mask);
            oc_vals{ii,hrf_bases_id} = odor_vals;
                       
            p_score_tr(1,hrf_bases_id) = fastcorr(odor_corr_vals(~foldind),a_(~foldind));
            p_score_tr(2,hrf_bases_id) = fastcorr(odor_corr_vals(~foldind),b_(~foldind));
            p_score_tr(3,hrf_bases_id) = fastcorr(odor_corr_vals(~foldind),ab(~foldind));
                       
            p_score_tt(1,hrf_bases_id) = fastcorr(odor_corr_vals(foldind),a_(foldind));
            p_score_tt(2,hrf_bases_id) = fastcorr(odor_corr_vals(foldind),b_(foldind));
            p_score_tt(3,hrf_bases_id) = fastcorr(odor_corr_vals(foldind),ab(foldind));
        end
        [~,argmax_p] = max(p_score_tr(1,:));
        [~,argmax_c] = max(p_score_tr(2,:));
        [~,argmax_pc] = max(p_score_tr(3,:));
        
        p_mat_tr(ii,foldid) = argmax_p;
        c_mat_tr(ii,foldid) = argmax_c;
        pc_mat_tr(ii,foldid) = argmax_pc;
              
        p_mat_tt(ii,foldid) = p_score_tt(1,argmax_p);
        c_mat_tt(ii,foldid) = p_score_tt(2,argmax_c);
        pc_mat_tt(ii,foldid)= p_score_tt(3,argmax_pc);        
    end       
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
figure()
set(gcf, 'Position',  [100, 100, 400, 300])
bar(corr_voxel_final_t1)
hold on
edges1 = 0.3;
edges2 = 0.1;
for ii = 1:length(anat_names)
    plot([ii-edges1,ii-edges2-0.05],[corr_voxel_final_p1(ii,1),corr_voxel_final_p1(ii,1)],'k');
    plot([ii-edges2+0.025,ii+edges2-0.025],[corr_voxel_final_p1(ii,2),corr_voxel_final_p1(ii,2)],'k');
    plot([ii+edges2+0.05,ii+edges1],[corr_voxel_final_p1(ii,3),corr_voxel_final_p1(ii,3)],'k');
    if strcmp(b_r_t,'b')
        plot([ii-edges1,ii-edges2-0.05],[corr_voxel_final_p2(ii,1),corr_voxel_final_p2(ii,1)],'k');
        plot([ii-edges2+0.025,ii+edges2-0.025],[corr_voxel_final_p2(ii,2),corr_voxel_final_p2(ii,2)],'k');
        plot([ii+edges2+0.05,ii+edges1],[corr_voxel_final_p2(ii,3),corr_voxel_final_p2(ii,3)],'k');
        plot([ii+(-edges1-edges2-0.05)/2,ii+(-edges1-edges2-0.05)/2],[corr_voxel_final_p1(ii,1),corr_voxel_final_p2(ii,1)],'k');
        plot([ii,ii],[corr_voxel_final_p1(ii,2),corr_voxel_final_p2(ii,2)],'k');
        plot([ii+(edges1+edges2+0.05)/2,ii+(edges1+edges2+0.05)/2],[corr_voxel_final_p1(ii,3),corr_voxel_final_p2(ii,3)],'k');
    end
end
yl = ylim;
ylim([0 yl(2)])
xlim([0 sum(anat_idx)+0.5])
% ylim([0 0.03])
xticklabels(anat_names)
xtickangle(90)
title('Std b, RSA')
mkdir(fullfile(statpath,modelname));
savefig(fullfile(statpath,modelname,'ROI_RSM_r2.fig'))
print(fullfile(statpath,modelname,'ROI_RSM_r2'),'-dpng')


% figure('Position',[0.5 0.5 320 240])
% bar(nvoxS)
% xticklabels(anat_names)
% xtickangle(90)
% title('Sig. odor-evoked voxels')
% savefig('nvoxes')
% print(fullfile(statpath,'nvoxes'),'-dpng')

M_c = [];
for jj = 1:length(anat_names) % Num areas
    t1 = p_mat_tt(jj,:);
    t2 = c_mat_tt(jj,:);
    M_c(jj) = bstrap_hyp(t1,t2);
end
save(fullfile(statpath,modelname,'ROI.mat'))

%% Subject level analysis. P:C ratio and stats
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR\classic_chem';
        'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\RSA_FIR\classic_chem';
        'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\RSA_FIR\classic_chem'};
nS = length(dirs);
matname = 'ROI.mat';
anat_names = {'PirF','PirT','AMY','OFC'};
p_mat_tt = variable_extract(dirs,matname,'p_mat_tt',false);
p_mat_tt = cat(3,p_mat_tt{:});
c_mat_tt = variable_extract(dirs,matname,'c_mat_tt',false);
c_mat_tt = cat(3,c_mat_tt{:});
pc_mat_tt = variable_extract(dirs,matname,'pc_mat_tt',false);
pc_mat_tt = cat(3,pc_mat_tt{:});

B_mat(:,1,:,:) = p_mat_tt;
B_mat(:,2,:,:) = c_mat_tt;
% B_mat(:,3,:,:) = pc_mat_tt;
mean_b_mat_ = squeeze(mean(B_mat,3));
std_b_mat_ = squeeze(std(B_mat,[],3));
p_mat = squeeze(mean(mean_b_mat_,1))';
err_mat = squeeze(std(mean_b_mat_,1))';

% For across model types:
% B_mat = cat(2,B_mat1(:,1,:,:),B_mat2(:,1,:,:));
% % b_mat = cat(2,b_mat1(:,1,:),b_mat2(:,1,:));
% mean_b_mat_ = cat(2,mean_b_mat_1(:,1,:),mean_b_mat_2(:,1,:));
% std_b_mat_ = cat(2,std_b_mat_1(:,1,:),std_b_mat_2(:,1,:));
%-------------------------------------------------------------------------

% bar_mean = p_mat;
% bar_std = err_mat;
% figure()
% bar(bar_mean)
% hold on
% ngroups = size(bar_mean, 1);
% nbars = size(bar_mean, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, bar_mean(:,i), bar_std(:,i), 'k.');
% end
% legend({'Perceptual','Chemical'})
% xticklabels({'S01', 'S02', 'S03'})
% savefig('PC_interaction')
% print(fullfile('PC_interaction'),'-dpng')

%--------------------- Mean subject ratings-------------------------------
mean_b_mat = mean(mean_b_mat_,3);

% Toggle s.e.m bootstrap bars on off. Combine bootstrap error across
% subjects.
M = [];
M2 = [];
for jj = 1:size(B_mat,1) % Num areas
    t1 = squeeze(B_mat(jj,1,:,:));
    t2 = squeeze(B_mat(jj,2,:,:));
%     t3 = squeeze(B_mat(jj,3,:,:));
    t1 = mean(t1,2);
    t2 = mean(t2,2);
%     t3 = mean(t3,2);
%     t1 = t1(:);
%     t2 = t2(:);
    M(jj,1) = prctile(t1,97.5);
    M(jj,2) = prctile(t2,97.5);
%     M(jj,3) = prctile(t3,97.5);
    M2(jj,1) = prctile(t1,2.5);
    M2(jj,2) = prctile(t2,2.5);
%     M2(jj,3) = prctile(t3,2.5);
end

anat_names = {'PirF','PirT','Amy','OFC'};
bar_mean = mean_b_mat;
bar_std = M;
bar_std2 = M2;


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
%     errorbar(x, bar_mean(:,i), bar_std(:,i)-bar_mean(:,i),bar_std2(:,i)-bar_mean(:,i), 'k.');
    x_m = [x_m; x];
end
xticks(1:4)
xticklabels(anat_names);
% legend({'Perceptual','Chemical','Mutual'})
% legend()
% Subject data points
c_s = {'r','g','b'}; % Data dots for subjects
for ii = 1:4 % For bars for perceptual, chemical and combinations
    for jj = 1:nS
        plot(x_m(:,ii),mean_b_mat_(ii,:,jj),c_s{jj})
    end
end
% ylim([0 0.16])
ylabel('Representational r')
savefig('anat_pfm')
print('anat_pfm','-dpng')


% %AIC
% bar_cat = 

% Subject wise tests
% M2 = [];
% for mm = 1:3
%     for jj = 1:size(B_mat,1) % Num areas
%         t1 = squeeze(B_mat(jj,1,:,mm));
%         t2 = squeeze(B_mat(jj,2,:,mm));
%         t_sq = mean(t1)-mean(t2);
%         tdist = (t1-t2)-mean(t1-t2);
%         M2(mm,jj) = (100-invprctile(tdist,t_sq))/100;
%     end
% end
    
% Group test for pairwise P:C 
M_c = [];
for jj = 1:size(B_mat,1) % Num areas
    t1 = squeeze(B_mat(jj,1,:,:));
    t2 = squeeze(B_mat(jj,2,:,:));
    t1 = mean(t1,2);
    t2 = mean(t2,2);
    M_c(jj) = bstrap_hyp_2(t2,t1);
end
save(fullfile('ROI_std2.mat'))

%-----P-C r^2 difference---------------------------------------------------
M_bar = [];
M_std = [];
M_std2 = [];
M_comp = [];
for jj = 1:size(B_mat,1) % Num areas
    t1 = squeeze(B_mat(jj,1,:,:));
    t2 = squeeze(B_mat(jj,2,:,:));
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
bar(1:4,M_bar_)
hold on
errorbar(1:4,M_bar_,M_std-M_bar_,M_std2-M_bar_,'.')
% Line art per subject
subwise_r2 = squeeze(mean_b_mat_(:,1,:).^2-mean_b_mat_(:,2,:).^2);
plot(subwise_r2)
xticklabels(anat_names)
ylabel('Representational r^2')
savefig('RSA_var')
print(fullfile('RSA_var'),'-dpng')
% [~,p]=corr(M_bar',[1:4]','Type','Spearman')
save(fullfile('ROI_std2.mat'))

save('RSA_pvals')

%% Neural RSM
% Run RSA for S=1 and ii (anat idx) = 4, hrf_bases for 6 sec
% Load swampsunset

% % Make colormap
% mat = odor_vals';
% mat2 = odor_corr;
% 
% nidx = 1:10;
% k_sumd = zeros(1,length(nidx));
% for ii = nidx
%     [~,~,sumd]=kmeans(mat,ii,'Distance','correlation');
%     k_sumd(ii) = mean(sumd);
% end

% figure('Position',[0.5 0.5 640 480])
% plot(k_sumd)
for ii = 1:4
mat = oc_vals{ii};
mat2 = corrcoef(mat);
[kidx] = kmeans(mat',4,'distance','correlation');
[~,argsort] =sort(kidx);
figure('Position',[0.5 0.5 640 480])
imagesc(mat2(argsort,argsort))
colormap(swampsunset)
caxis([-0.5 0.7])
xlabel('Odors')
ylabel('Odors')
title(sprintf(anat_names{ii}))
colorbar
colormap(swampsunset)
savefig(fullfile(statpath,anat_names{ii}))
print(fullfile(statpath,anat_names{ii}),'-deps')
end

% figure('Position',[0.5 0.5 640 480])
% m3 = spectral_reorder(1-odor_corr);
% imagesc(1-m3)
% caxis([-0.5 0.7])
% xlabel('Odors')
% ylabel('Odors')
% title(sprintf(anat_names{ii}))
% colorbar
% colormap(swampsunset)
% savefig(fullfile(statpath,anat_names{ii}))
% print(fullfile(statpath,anat_names{ii}),'-dpng')

% hex = ['#0D646B';'#798279';'#C3D1C3';'#E17812';'#FAE0AD'];
% map = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;

% Load RSMs from RSA analysis

%% Lesion Models
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR\RSA_final';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\RSA_FIR\RSA_final';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\RSA_FIR\RSA_final'};
dirs_int = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_lesions\Int_orthog';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\RSA_lesions\Int_orthog';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\RSA_lesions\Int_orthog'};
dirs_pls = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_lesions\Pls_orthog';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\RSA_lesions\Pls_orthog';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\RSA_lesions\Pls_orthog'};

matname = 'ROI.mat';
anat_names = {'PirF','PirT','AMY','OFC'};
p_mat_tt = variable_extract(dirs,matname,'p_mat_tt',false);
p_mat_tt = cat(3,p_mat_tt{:});
p_mat_tt_int = variable_extract(dirs_int,matname,'p_mat_tt',false);
p_mat_tt_int = cat(3,p_mat_tt_int{:});
p_mat_tt_pls = variable_extract(dirs_pls,matname,'p_mat_tt',false);
p_mat_tt_pls = cat(3,p_mat_tt_pls{:});

B_mat(:,1,:,:) = p_mat_tt;
B_mat(:,2,:,:) = p_mat_tt_pls;
B_mat(:,3,:,:) = p_mat_tt_int;
mean_b_mat_ = squeeze(mean(B_mat,3));
std_b_mat_ = squeeze(std(B_mat,[],3));

%--------------------- Mean subject ratings-------------------------------
mean_b_mat = mean(mean_b_mat_,3);
% std_b_mat = rms(std_b_mat_,3);
% std_2 = std(mean_b_mat_,[],3)./sqrt(3);
% std_b_mat = 1.96*sqrt(std_b_mat.^2+std_2.^2);

% Toggle s.e.m bootstrap bars on off. Combine bootstrap error across
% subjects.
M = [];
M2 = [];
for jj = 1:size(B_mat,1) % Num areas
    for kk=1:size(B_mat,2)
        t1 = squeeze(B_mat(jj,kk,:,:));
        t1 = mean(t1,2);
        M(jj,kk) = prctile(t1,97.5);
        M2(jj,kk) = prctile(t1,2.5);
    end
end

anat_names = {'FPC','TPC','Amygdala','OFC'};
bar_mean = mean_b_mat;
bar_std = M;
bar_std2 = M2;

figure()
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
%     errorbar(x, bar_mean(:,i), bar_std(:,i)-bar_mean(:,i),bar_std2(:,i)-bar_mean(:,i), 'k.');
    x_m = [x_m; x];
end
xticks(1:4)
xticklabels(anat_names);
% legend()
% Subject data points
c_s = {'r','g','b'}; % Data dots for subjects
for ii = 1:4 % For bars for perceptual, chemical and combinations
    for jj = 1:3
        plot(x_m(:,ii),mean_b_mat_(ii,:,jj),c_s{jj})
    end
end
savefig('anat_pfm')
print('anat_pfm','-dpng')

% Subject wise tests
% M2 = [];
% for mm = 1:3
%     for jj = 1:size(B_mat,1) % Num areas
%         t1 = squeeze(B_mat(jj,1,:,mm));
%         t2 = squeeze(B_mat(jj,2,:,mm));
%         t_sq = mean(t1)-mean(t2);
%         tdist = (t1-t2)-mean(t1-t2);
%         M2(mm,jj) = (100-invprctile(tdist,t_sq))/100;
%     end
% end
    
% Group test for pairwise P:C 
M = [];
for jj = 1:size(B_mat,1) % Num areas
    t1 = squeeze(B_mat(jj,1,:,:));
    t2 = squeeze(B_mat(jj,3,:,:));
    t1 = mean(t1,2);
    t2 = mean(t2,2);
    t_sq = mean(t1)-mean(t2);
    tdist = (t1-t2)-mean(t1-t2);
    M(jj) = (100-invprctile(tdist,t_sq))/100;
end
save(fullfile('ROI_std2.mat'))

%-----P-C r^2 difference---------------------------------------------------
M_bar = [];
M_std = [];
M_std2 = [];
M_comp = [];
for jj = 1:size(B_mat,1) % Num areas
    t1 = squeeze(B_mat(jj,1,:,:));
    t2 = squeeze(B_mat(jj,2,:,:));
    t1 = mean(t1,2);
    t2 = mean(t2,2);
    M_bar(:,jj) = (t1.^2)-(t2.^2);
    M_std(jj) = prctile(t1.^2-t2.^2,2.5);
    M_std2(jj) = prctile(t1.^2-t2.^2,97.5);
     if jj>2
        t1 = M_bar(:,jj);
        t2 = M_bar(:,jj-2);
        t_sq = mean(t1)-mean(t2);
        tdist = (t1-t2)-mean(t1-t2);
        M_comp(jj) = (100-invprctile(tdist,t_sq))/100;
    end
end
figure('Position',[0.5 0.5 300 200])
M_bar_ = mean(M_bar,1);
bar(1:4,M_bar_)
hold on
errorbar(1:4,M_bar_,M_std-M_bar_,M_std2-M_bar_,'.')
% Line art per subject
subwise_r2 = squeeze(mean_b_mat_(:,1,:).^2-mean_b_mat_(:,2,:).^2);
plot(subwise_r2)
xticklabels(anat_names)
savefig('RSA_var')
% [~,p]=corr(M_bar',[1:4]','Type','Spearman')

save('RSA_pvals')

%% Neural RSM: kmeans cluster
scores = zeros(4,9);
for zz = 1:4
    scores(zz,:) = k_mean_clusterer(oc_vals{zz}');
end
figure('Position',[0.5 0.5 350 250])
plot([2:10],scores')
xlabel('Num clusters')
ylabel('On-Off diagonal')
legend(anat_names)
savefig('kmeans_clusters_inc')
