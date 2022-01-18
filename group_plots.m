
%% Average EM plots
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

%% Model comparision for two EM
% Average of prediction accuracy
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg1_tscore'};
dirs2 = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\PC70';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\PC70';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\PC70'};
matname = 'eem_weights.mat';
varname = 'mask_scores';
labs = {'PirF','PirT','AMY','OFC'};
nvoxes = variable_extract(dirs,matname,varname,true);
nvoxes2 = variable_extract(dirs2,matname,varname,true);
nanat = 4;
nsub = 3;
nboot = 1000;

% Bootstrap means and C.I.
bar_mean = []; % Mean voxel responses
bar_std = []; % 97.5 prctile
bar_std2 = []; % 2.5 prctile 
M_c = zeros(nanat,1);
dists_ = {};
for ii = 1:nanat % anatomical areas
    temp_dist = zeros(1,nboot);
    temp_dist2 = zeros(1,nboot);
    for jj = 1:nsub % num subjects
        nvoxes{ii,jj}(isnan(nvoxes{ii,jj}))=[];
        nvoxes2{ii,jj}(isnan(nvoxes2{ii,jj}))=[];
        [~,bootsam] = bootstrp(nboot,@(x) x,1:length(nvoxes{ii,jj}));
        [~,bootsam2] = bootstrp(nboot,@(x) x,1:length(nvoxes2{ii,jj}));
        temp_dist = temp_dist+mean(nvoxes{ii,jj}(bootsam));
        temp_dist2 = temp_dist2+mean(nvoxes2{ii,jj}(bootsam));
    end
    temp_dist = temp_dist/3;
    temp_dist2 = temp_dist2/3;
    bar_mean(ii,1) = mean(temp_dist);
%     bar_mean(ii,2) = mean(temp_dist2);
    bar_std(ii,1) = prctile(temp_dist,97.5);
%     bar_std(ii,2) = prctile(temp_dist2,97.5);
    bar_std2(ii,1) = prctile(temp_dist,2.5);
%     bar_std2(ii,2) = prctile(temp_dist2,2.5);
    M_c(ii) = bstrap_hyp_2(temp_dist,temp_dist2);
    dists_{ii} = temp_dist;
end

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
bar_mat = cellfun(@mean,nvoxes);
for jj = 1:3 % For bars for perceptual, chemical and combinations
    plot(bar_mat(:,jj),'MarkerSize',10)
end
savefig('EM_comp')
save('EM_comparision.mat')

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

% nvoxes{1}(:,1,:)=nan*nvoxes{1}(:,1,:);
% nvoxes{2}(:,1,:)=nan*nvoxes{2}(:,1,:);
% nvoxes{2}(:,2,:)=nan*nvoxes{2}(:,2,:);
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

%% ROI wise subjectivity
% Order of PCs not maintained during bootstrap. Hence, correlation among
% PCs for a given ROI for a given subject ~= 0. Also, PCs computed on all
% 18 descriptors not on common descriptors
% load('C:\Data\NEMO\NEMO_all\EM\ROI_subjectivity_PC\ROI_settings.mat')

dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg1_tscore'};
matname = 'eem_weights.mat';
anat_names = {'PirF','PirT','AMY','OFC'};
weights_ = variable_extract(dirs,matname,'bar_mat',true);
ndescrip = 18;
nanat = 4;
nS = 3; %length(dirs);
fname = 'ROI_generalizability_nonboot';
diff_mod = false; % Take difference instead of overall mean
diff_sub = false ; % S1_a.S2_a - S1_a.S1_a - S2_a.S2_a. 
% If false, S1_a.S2_a - S1_a.S2_b and so on..., where a and b are areas

% Extrude common descriptors across subjects
inds = true(3,ndescrip);
inds(1,[11 16 18]) = false;
inds([2 3],[16:18]) = false;

% Weight modifier
weights_norm = cell(weights_);
for ii = 1:nanat
    for jj = 1:nS
        weights_norm{jj,ii} = weights_{jj,ii}(:,inds(jj,:));
    end
end

% Correlation of weights per voxel
mean_cons = zeros(1,nanat);
std_cons = zeros(1,nanat);
Ms = cell(1,nanat); % correlation values
Ms2 = cell(1,nanat);
M_corr = cell(1,nanat);
for ii = 1:nanat
    for jj1 = 1:nS
        for jj2 = jj1+1:nS
            M = corrcoef_2(weights_norm{jj1,ii},weights_norm{jj2,ii});
            Ms{ii} = cat(1,Ms{ii},M(:));
            if diff_sub
                M_corr_mat = corrcoef_2(weights_norm{jj1,ii},weights_norm{jj1,ii});
                utl = triu(true(size(M_corr_mat)),1);
                M_corr{ii} = cat(1,M_corr{ii},M_corr_mat(utl));
            else
                for ii_left = setdiff(1:nanat,ii)
                    M_corr_mat = corrcoef_2(weights_norm{jj1,ii},weights_norm{jj1,ii_left});
                    M_corr{ii} = cat(1,M_corr{ii},M_corr_mat(:));
                end
            end
        end
    end
    if diff_mod
        mean_cons(ii) = mean(Ms{ii}.^2)-mean(M_corr{ii}.^2);
        std_cons(ii) = sqrt((var(Ms{ii}.^2)./length(Ms{ii})+var(M_corr{ii}.^2)./length(M_corr{ii}))/2);
    else
        mean_cons(ii) = mean((Ms{ii}.^2));
        std_cons(ii) = std(Ms{ii}.^2)./sqrt(length(Ms{ii}));
    end 
end

figure('Position',  [100, 100, 320, 240])
bar(mean_cons)
hold on
errorbar(mean_cons,1.96*std_cons,'.')
xticks(1:4)
xticklabels(anat_names)
ylabel('Consistency')
% savefig(fname)
% print(fname,'-dpng')
M_test1(1) = ranksum(Ms{1},Ms{2});
M_test1(2) = ranksum(Ms{2},Ms{3});
M_test1(3) = ranksum(Ms{3},Ms{4});

% % Plot correlation matrix of mean-vectors from each ROI across subjects
% num_vecs = zeros(nanat*length(dirs),sum(inds(1,:)));
% labels = {};
% kk = 0;
% for ii = 1:nanat
%     for jj = 1:length(dirs)
%         kk = kk+1;
%         num_vecs(kk,:) = mean(weights_{jj,ii}(:,inds(jj,:)));
%         labels{kk} = sprintf('S%01d: %01s',jj,anat_names{ii});
%     end
% end
% k_C_corr = corrcoef(num_vecs');
% figure('Position',  [100, 100, 400, 300])
% imagesc(k_C_corr.^2)
% xticks(1:nanat*length(dirs))
% xticklabels(labels)
% xtickangle(90)
% yticks(1:nanat*length(dirs))
% yticklabels(labels)
% colorbar
% title('Mean profiles')
% savefig('Meanprofile')
% save('ROI_subjectivity')

%% PCA of OFC
load('C:\Data\NEMO\NEMO_all\EM\ROI_subjectivity_PC\ROI_settings.mat')
load('C:\Data\NEMO\swampsunset.mat')
num_pcs = 4;
ii = 4;
reorder_pc = true;

if reorder_pc
    w_cat =[];
    for jj = 1:length(dirs)
        w_cat= cat(1,w_cat, weights_{jj,ii}(:,inds(jj,:)));
    end
%     [coeffs_cat,~,~,~,vars] = pca(w_cat);
        [coeffs_cat,~,~,~,vars] = pca(weights_{2,ii}(:,inds(2,:)));
end
coeffs_ = zeros(length(dirs),sum(inds(1,:)),sum(inds(1,:)));
for jj = 1:length(dirs)
    if ~reorder_pc
        coeffs_(jj,:,:) = pca(weights_{jj,ii}(:,inds(jj,:)));
    else
        % Gale-Shapely algorithm
        temp_coeff = pca(weights_{jj,ii}(:,inds(jj,:)));
        % Proposer-preference
        rank_order = corrcoef_cos2(temp_coeff(:,1:num_pcs)',coeffs_cat(:,1:num_pcs)');
        rank_id1 = zeros(size(rank_order));
        % Proposee-preference
        rank_id2 = zeros(size(rank_order));
        for zz = 1:length(rank_order)
            [~,temp] = sort(rank_order(zz,:),'descend');
            [~,temp2] = sort(temp);
            [~,rank_id1(zz,:)] = sort(temp2);
            
            % Proposee ranks
            [~,temp] = sort(rank_order(:,zz),'descend');
            [~,temp2] = sort(temp);
            [~,rank_id2(zz,:)] = sort(temp2);
        end
        col_order = galeshapley(num_pcs,rank_id1,rank_id2);
        col_order = [col_order;((num_pcs+1):sum(inds(1,:)))'];
        coeffs_(jj,:,:) =  temp_coeff(:,col_order);
    end
end

% Sort_id
% coeffs_(:,:,[2 3 4]) = coeffs_(:,:,[3 4 2]);

% Arrange the components and group similar PCs for subjects
num_vecs = zeros(num_pcs*length(dirs),sum(inds(1,:)));
labels = {};
kk = 0;
block = []; % Create indices for extracting within and across cluster elements
clustermat = ones(nS);
for pp = 1:num_pcs
    for jj = 1:length(dirs)
        kk = kk+1;
        num_vecs(kk,:) = squeeze(coeffs_(jj,:,pp));
        labels{kk} = sprintf('S%01d, PC%01d',jj,pp);
    end
    block = blkdiag(block,clustermat);
end
block = logical(block);
utl = logical(triu(ones(size(block)),1));
M_on = and(block,utl); % Indices for within PC correlation across subjects
M_off = and(~block,utl); % Indices for across PC correlations

% Remove offdiagonal trivial entries
M_id = []; for nn = 1:size(M_off,1); for mm = 1:size(M_off,2); M_id(nn,mm)=sub2ind(size(M_off),nn,mm) ; end; end
for pp = 1:(num_pcs-1)
    inds_ = diag(M_id,length(dirs)*pp);
    M_off(inds_)=false;
end

% Indices of components
M_on_i = zeros(size(M_on,1),size(M_on,1),num_pcs);
for pp = 1:num_pcs
    M_on_i(length(dirs)*(pp-1)+(1:length(dirs)),length(dirs)*(pp-1)+(1:length(dirs)),pp)...
        = M_on(length(dirs)*(pp-1)+(1:length(dirs)),length(dirs)*(pp-1)+(1:length(dirs)));
end

% For each ROI, plot correlation matrix of PCs across subjects
figure('Position',  [100, 100, 400, 300])
k_C_corr = corrcoef(num_vecs');
k_C_corr = sign(k_C_corr).*(k_C_corr.^2);
% k_C_corr = atanh(k_C_corr);
imagesc(k_C_corr)
xticks(1:3:num_pcs*length(dirs))
xticklabels(labels(1:3:end))
xtickangle(90)
yticks(1:3:num_pcs*length(dirs))
yticklabels(labels(1:3:end))
colorbar
colormap(CustomColormap)
title(sprintf('Correlation PCs: %s',anat_names{ii}))
% savefig('ofc_pc_mat')

% Bar graph of PC-wise subjectivity
m_cell = cell(1,num_pcs);
m_labels = cell(1,num_pcs); % For Anova
for zii=1:num_pcs
    m_cell{zii}=k_C_corr(logical(squeeze(M_on_i(:,:,zii))));
    m_labels{zii} = zii*ones(length(dirs),1);
end
% m_cell{end+1} = k_C_corr(M_off);
% m_labels{end+1} = (num_pcs+1)*ones(sum(M_off(:)),1);
m_cell_ = vertcat(m_cell{:});
m_labels_ = vertcat(m_labels{:});
[~,~,stats] = anova1(m_cell_,m_labels_);

bars = cellfun(@mean,m_cell);
stds = cellfun(@std,m_cell)./sqrt(3);
bar(bars)
hold on
errorbar(bars,stds,'.')

% Tuning curves
figure('Position',  [100, 100, 400, 300])
hold on
coeffs_m = squeeze(mean(coeffs_(:,:,1:num_pcs)));
plot(coeffs_m)
thr = 0.2;
yline(thr)
% yline(-thr)
xticks(1:sum(inds(1,:)))
xticklabels(behav_labs)
xtickangle(90)
%     legend()
title('Tuning curves')
legend({'PC1','PC2','PC3','PC4','PC4'})
% savefig('tuning_curves')

%% ROI wise subjectivity - bootstrapped
% Order of PCs not maintained during bootstrap. Hence, correlation among
% PCs for a given ROI for a given subject ~= 0. Also, PCs computed on all
% 18 descriptors not on common descriptors
% load('C:\Data\NEMO\NEMO_all\EM\ROI_subjectivity_PC\ROI_settings.mat')

dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg1_tscore';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg1_tscore'};
matname = 'eem_weights_boot.mat';
sz_cntrl = true;
anat_names = {'PirF','PirT','AMY','OFC'};
weights_ = variable_extract(dirs,matname,'bar_mat',true);
ndescrip = 18;
nanat = 4;
nS = 3; %length(dirs);
fname = 'ROI_generalizability_boot_corr';
nboot = 1000;
combine_weights = false; % Add APC and PPC
% If false, S1_a.S2_a - S1_a.S2_b and so on..., where a and b are areas

% Extrude common descriptors across subjects
inds = true(3,ndescrip);
inds(1,[11 16 18]) = false;
inds([2 3],[16:18]) = false;

if combine_weights
    weights_comb = cell(nS,nanat);
    for jj = 1:nS
        weights_comb{jj,1} = cat(1,weights_{jj,1},weights_{jj,2});
        weights_comb{jj,2} = weights_{jj,3};
        weights_comb{jj,3} = weights_{jj,4};
    end
else
    weights_comb = weights_;
end

% Weight modifier
weights_norm = cell(size(weights_comb));
for ii = 1:nanat
    for jj = 1:nS
        weights_norm{jj,ii} = weights_comb{jj,ii}(:,inds(jj,:));
    end
end

% Correlation of weights per voxel
mean_cons = zeros(1,nanat);
std_cons = zeros(1,nanat);
std_cons2 = zeros(1,nanat);
Ms = cell(1,nanat); % correlation values
for ii = 1:nanat
    kk = 0;
    for jj1 = 1:nS
        for jj2 = jj1+1:nS
            kk=kk+1;
            M = corrcoef_2(weights_norm{jj1,ii},weights_norm{jj2,ii});   
%             Ms{kk,ii} = M(:);
            Ms{ii} = cat(1,Ms{ii},(M(:)));
        end
    end
    if sz_cntrl
        nsam = 10000;
        Ms_temp = [];
        for data_id = 1:nsam
            Ms_id = datasample(1:length(Ms{ii}),939);
            Ms_temp(data_id)=mean(Ms{ii}(Ms_id));
        end
        Ms{ii} = sign(Ms_temp').*(Ms_temp'.^2);
    else
        Ms{ii} =  bootstrp(nboot, @mean, Ms{ii});
    end
    mean_cons(ii) = mean(Ms{ii},'all');
    std_cons(ii) = prctile(Ms{ii},97.5)-mean_cons(ii);
    std_cons2(ii) = prctile(Ms{ii},2.5)-mean_cons(ii);
end

figure('Position',  [100, 100, 320, 240])
bar(1:nanat,mean_cons)
hold on
errorbar(1:nanat,mean_cons,std_cons,std_cons2,'.')
xticks(1:4)
xticklabels(anat_names)
ylabel('Consistency')
% savefig(fname)
% print(fname,'-dpng')
M_test1(1) = bstrap_hyp_2(Ms{1},Ms{2});
M_test1(2) = bstrap_hyp_2(Ms{1},Ms{3});
M_test1(3) = bstrap_hyp_2(Ms{1},Ms{4});

% % Plot correlation matrix of mean-vectors from each ROI across subjects
% num_vecs = zeros(nanat*length(dirs),sum(inds(1,:)));
% labels = {};
% kk = 0;
% for ii = 1:nanat
%     for jj = 1:length(dirs)
%         kk = kk+1;
%         num_vecs(kk,:) = mean(weights_{jj,ii}(:,inds(jj,:)));
%         labels{kk} = sprintf('S%01d: %01s',jj,anat_names{ii});
%     end
% end
% k_C_corr = corrcoef(num_vecs');
% figure('Position',  [100, 100, 400, 300])
% imagesc(k_C_corr.^2)
% xticks(1:nanat*length(dirs))
% xticklabels(labels)
% xtickangle(90)
% yticks(1:nanat*length(dirs))
% yticklabels(labels)
% colorbar
% title('Mean profiles')
% savefig('Meanprofile')
% save('ROI_subjectivity')

%% Subjectwise HRFs
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\NEMO_fir';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\NEMO04_RSA_FIR'};
matname = 'hrf.mat';
anat_names = {'PirF','PirT','AMY','OFC'};
hrfs = variable_extract(dirs,matname,'hrfs',false);
err_bars = variable_extract(dirs,matname,'err_bars',false);

figure('Position',[0 0 800 200])
hold on
cc = {'r','g','b'};
for ii = 1:4
    subplot(1,4,ii)
    for jj = 1:3
    mt = hrfs{jj}(:,ii)';
    mt = mt-min(mt);
    mt_ = mt./max(mt);
    shaded_plot(1:length(hrfs{jj}(:,ii)),mt_,(err_bars{jj}(:,ii)')./max(mt),cc{jj})
    end
    title(anat_names{ii})
    xticks([0 5 10])
    xlim([1 10])
end

%% Sniff evoked voxels
% nvox sniff: C:\Data\NEMO\NEMO_all\Sniff\Sniff_bar
bars = mean(nvoxes);
std_ = std(nvoxes);
nanat = {'PirF','PirT','AMY','OFC'};
bar(1:4,bars)
hold on
errorbar(1:4,bars,std_/sqrt(3),'.')
plot(1:4,nvoxes,'Marker','.','MarkerSize',10)

%% EM methods simulation
cellsz = 16;

a1 = [rand*ones(cellsz,1); rand*ones(cellsz,1); rand(cellsz,1); rand*ones(cellsz,1)];
figure('Position',[0.2 0.2 60 640])
imagesc(a1)
colormap(swamp)

a2 = [];    
for zz = 1:2*cellsz
a2 = [a2 [rand*ones(cellsz,1); rand*ones(cellsz,1); rand(cellsz,1); rand*ones(cellsz,1)]];
end
figure('Position',[0.2 0.2 120 640])
imagesc(a2)
colormap(swamp)

%% Perceptual generalization
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR\RSA_final';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\NEMO_fir\RSA_final';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\NEMO04_RSA_FIR\RSA_final'};
matname = 'ROI.mat';
anat_names = {'PirF','PirT','AMY','OFC'};
percepts = variable_extract(dirs,matname,'Behav_RSM_vals_P',true);

% Construct subsets of cids for individual subjects for generalization:
load('C:\Data\NEMO\NEMO_all\NEMO_perceptual.mat')
cid1 = ismember(behav(1).cid,behav(2).cid);
cid2 = ismember(behav(2).cid,behav(1).cid);
utl =  triu(true(160),1);
mat1 = false(160);
mat1(cid1,cid1)=true;
mat_id(:,1) = mat1(utl);
mat2 = false(160);
mat2(cid2,cid2)=true;
mat_id(:,2) = mat2(utl);
mat_id(:,3) = mat2(utl);

% % Perceptual correlations
% p_corr = [];
% for ii = 1:3; p_corr(:,ii)=percepts(mat_id(:,ii),ii); end
% p_corrmat = corrcoef(p_corr);
% figure('Position',[0.2 0.2 640 480])
% bar(p_corrmat([2 3 6]))
% yline(r2t(0.025,7750))
% xticklabels({'12','13','23'})
% ylabel('perceptual similarity generalizability')
% savefig('Percept_g')
p_corr = [];
for ii = 1:3; p_corr(:,ii)=percepts(mat_id(:,ii),ii); end
p_corrmat = corrcoef(p_corr);
figure('Position',[0.2 0.2 640 480])
m = p_corrmat([2 3 6]);
bar(mean(m))
hold on
errorbar([1],mean(m),std(m)/sqrt(3))
yline(r2t(0.025,7750))
plot([1 1 1],m,'r.','MarkerSize',15)
xticklabels({'12','13','23'})
ylabel('perceptual similarity generalizability')
savefig('Percept_g')

% Anatomical ROI_corrs
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR\RSA4-7';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\NEMO_fir\RSA4-7';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\NEMO04_RSA_FIR\RSA4-7'};
matname = 'ROI.mat';
oc_vals = variable_extract(dirs,matname,'oc_vals',false);
M = [];
for ii = 1:length(anat_names)
    kk = 0;
    for jj1 = 1:length(oc_vals)
        for jj2 = jj1+1:length(oc_vals)
            kk = kk+1;
            M(ii,kk) = fastcorr(oc_vals{jj1}{ii,3}(mat_id(:,jj1)),oc_vals{jj2}{ii,3}(mat_id(:,jj2)));
        end
    end
end
figure()
bar(mean(M,2))
errorbar(mean(M,2),std(M,[],2)/sqrt(3))
xticks(1:4)
xticklabels(anat_names)
savefig('RSA_gen')



