%% Analyze the prediction accuracy map and break it down by sub regions
% Analyze argmaxes (time-basis) for different brain regions for different descriptors.

% Major inputs:
% Output file ('None.mat') of EEM_LOOCV.
% Needs fullzcored.mat: Matrix of voxels x HRF bases x Odors
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% Gray matter masks of ROIs
% anat_gw.nii: A binary mask of gray matter voxels such that number of
% voxels in fullzscored.mat = sum(anat_gw,'all')


s = 2; % Subject
model_comp = true;
argmax_comp = false;

root = 'C:\Data\NEMO';
anat_names = {'PirF','PirT','AMY','OFC'};
anat_masks = {'rwAPC.nii','rwPPC.nii','rwAmygdala.nii', 'rwOFC.nii'};
% anat_names = {'APC','PPC','Amyg','OFC','Hipp','EC'};
% anat_masks = {'rwAPC.nii','rwPPC.nii','rwAmygdala.nii', 'rwOFC.nii','rwHipp.nii','rwEC.nii'};
nanat = length(anat_names);
anatpath = fullfile('C:\Data\NEMO\',sprintf('NEMO_%02d',s),'\imaging\nii\masks');
maskfile =  'anat_gw.nii';
mask = (spm_read_vols(spm_vol(fullfile(anatpath, maskfile)))); % Mask used to construct odor files
mask(isnan(mask))=0;
mask = logical(mask);
fmaskfile = 'f_anat_gw3.nii';
fmask = (spm_read_vols(spm_vol(fullfile(anatpath, fmaskfile)))); % Mask used to examine voxels in RSA
fmask(isnan(fmask))=0;
fmask = logical(fmask);
fmask_1d = fmask(mask);
draw_fig =true;
kmeans_pos = true;
plot_pos = true; % Even if k-means is constructed on signed weights, use positive value for components
rem_lowdeg = false; % remove voxels with low degree in graph
deg_thresh = 2;
conn_thresh = 0.001; % p-value of significance of connections
sr_thresh_ = 0.1; % Choose voxels above this p-value
tscorer = true; % Run analyses on tscores of weights instead of betas
kmediods_ = false; % Use Spearmann correlation on kmediods instead of correlation on kmeans
kmeans_pos_pregraph = false; % Absolute value of weights before graph construction

% Model path
statpath = fullfile(root,sprintf('NEMO_%02d',s),'\imaging\1stlevelmodels\FIR_EM\PC40');
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

content_sw = 'EM'; %'RSA_p','RSA_c','EM'
%% Model comp
mask_scores = cell(length(anat_masks),1);
mask_tw = cell(length(anat_masks));

load(fullfile(statpath,'None.mat'),'corr_voxel_3d','r_thresh','behav','corr_voxel_final')
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
end
savefig('hists')

% Box Plots
% grouping_ind
% colororder({'k','k'})
gp = masks_set.*[1:nanat];
gp_id = zeros(size(gp,1),1);
for ii = 1:size(gp,1)
    if sum(gp(ii,:))>0
        gp_id(ii) = find(gp(ii,:)>0,1);
    else
        gp_id(ii) = 0;
    end
end
figure('Visible','on')
set(gcf, 'Position',  [100, 100, 400, 300])
hold on
boxplot(temp_array_(gp_id>0),gp_id(gp_id>0),'symbol', '')
xticks(1:nanat)
xticklabels(anat_names)
xtickangle(90)
plot([0.5 6.5],[sr_thresh sr_thresh],'k')
% plot([0.5 6.5],[r2t(0.01,160) r2t(0.01,160)],'k--')
% yticklabels({})
% yyaxis right
% ytickangle(90)
% Sign rank test
yl = ylim;
pvals = splitapply(@signrank,temp_array_(gp_id>0),gp_id(gp_id>0));
paxis = 1:nanat;
ylaxis = yl(2)*ones(size(paxis));
p_1 = pvals<0.001;
p_2 = pvals<0.05/6;
p_3 = pvals<0.05;

plot(paxis(p_1),ylaxis(p_1)-0.05,'r*')
plot(paxis(p_2),ylaxis(p_2)-0.025,'r*')
plot(paxis(p_3),ylaxis(p_3),'r*')
ylim([-inf yl(2)+0.01])
% legend({'p<0.05, whole brain FDR corrected','p<0.01 uncorrected'})
ylabel('Prediction r')
if draw_fig
    savefig('boxp')
    print('boxp','-dpng')
end
save('eem_weights.mat')

%% Graph analysis
% Mask flow:
% "mask" to convert 3d to 1d with nvox  =gray matter voxels across all ROIs.
% fmask_1d to find vox with sig sniff activity
% masks_set or mask_anat to further find vox in specific ROIs
% scorr_mask to find sig voxels in EM

load(fullfile(statpath,'None.mat'),'corr_voxel_3d','r_thresh','train_cell','corr_voxel_final','behav')
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


%     if kmeans_pos
%         train_mat = abs(train_mat);
%     end
    graph_nodes = mean(train_mat,3);
    err_ = std(train_mat,[],3);
    graph_nodes = graph_nodes./(err_);
    if kmeans_pos_pregraph
        graph_nodes = abs(graph_nodes);
    end


if ~kmediods_
    A = pdist(graph_nodes,'correlation');
else
    A = pdist(graph_nodes,'Spearman');
end
A = 1-squareform(A);
A = A-eye(size(A));
A(A<r2t(conn_thresh,size(train_mat,2)))=0;
%   A(A<0.67)=0;
% A = logical(A);

% Create graph ----------------------------------------
G = graph(A);
degs2 = G.degree;

if rem_lowdeg
    graph_nodes(degs2<deg_thresh,:)=[];
    if ~kmediods_
        A = pdist(graph_nodes,'correlation');
    else
        A = pdist(graph_nodes,'Spearman');
    end
    A = 1-squareform(A);
    A = A-eye(size(A));
    A(A<r2t(conn_thresh,size(train_mat,2)))=0;
    %      A(A<0.67)=0;
    %     A = logical(A);
    
    % Create graph ----------------------------------------
    G = graph(A);
    degs = G.degree;
else
    degs = G.degree;
end

figure('Visible','off')
set(gcf, 'Position',  [100, 100, 400, 300])
hold on
subplot(2,1,1)
histogram(degs,10)
ylabel('freq(degree)')
title('Histogram of voxel degrees')
subplot(2,1,2)
hold on
[H, edges, bins] = histcounts(degs,10);
cents = (edges(2:11)+edges(1:10))/2;
% cents(1) = [];
xlog = log(cents);
ylog = log(H);
% ylog(1) = [];
plot(xlog,ylog,'.')
X = [ones(size(xlog')) xlog'];
ylog(isinf(ylog))=nan;
mdl = fitlm(X,ylog');
plot(xlog,mdl.Coefficients{2,1}+xlog*mdl.Coefficients{3,1},'r')
title(sprintf('loglog histogram m: %02f, p: %02f',mdl.Coefficients{3,1},mdl.Coefficients{3,4}))
if draw_fig
    savefig('hist_degs')
    print('hist_degs','-dpng')
end

% Deg box plots across ROIs
masks_set_mod = masks_set(mask_anat,:);
masks_set_em = masks_set_mod(scorr_mask,:);
gp = masks_set_em;

if rem_lowdeg
    gp(degs2<deg_thresh,:) = [];
end

gp_id2 = zeros(size(gp,1),1);
for ii = 1:size(gp,1)
    if sum(gp(ii,:))>0
        gp_id2(ii) = find(gp(ii,:)>0,1);
    else
        gp_id2(ii) = 0;
    end
end

figure('Visible','off')
set(gcf, 'Position',  [100, 100, 400, 300])
hold on
boxplot(degs(gp_id2>0),gp_id2(gp_id2>0))
xticks(1:nanat)
xticklabels(anat_names)
xtickangle(90)
ylabel('ndegs')
if draw_fig
    savefig('boxp_degs')
    print('boxp_edges','-dpng')
end

% Analysis of visual clustering ----------------------------
G2 = G;
figure()
set(gcf, 'Position',  [100, 100, 400, 300])
H = plot(G2,'NodeLabel',{});
if nanat==6
    ccode = {'y','r','k','b','g','c'};
else
    ccode = {'k','r','g','b'};
end
for ii = 1:nanat
    highlight(H,find(gp(:,ii)),'NodeColor',ccode{ii})
end
if draw_fig
    savefig('Net')
    print('Net','-dpng')
end

%% Graph interpretation
% % for running k-means on absolute values uncomment:
if kmeans_pos
    graph_nodes = abs(graph_nodes);
end

% Overlay k-means on graph
nidx = 1:10;
k_sumd = zeros(1,length(nidx));
for ii = nidx
    if ~kmediods_
        [~,~,sumd]=kmeans(graph_nodes,ii,'Distance','correlation');
    else
        [~,~,sumd]=kmedoids(graph_nodes,ii,'Distance','spearman');
    end
    k_sumd(ii) = mean(sumd);
end
figure('Visible','off')
plot(1:length(k_sumd),k_sumd./k_sumd(1))
ylabel('Distortion score')
xlabel('num clusters')
xlim([1 inf])
savefig('kmeans_elbow')
print('kmeans_elbow','-dpng')

nC = 3; % Number of clusters
if ~kmediods_
    [kidx,k_C]=kmeans(graph_nodes,nC,'Distance','correlation','Replicates',100);
else
    [kidx,k_C]=kmedoids(graph_nodes,nC,'Distance','correlation','Replicates',20);
end
target_sort = [2 3 1];
kidx = k_index_sorter(kidx,target_sort );
k_C = k_C(target_sort,:);
figure('Position',  [100, 100, 400, 300])
H = plot(G2,'NodeLabel',{});
ccode = {[1 0 1],[0 1 1],[1 1/2 0]};
% ccode = {'y','r','k','b','g','c'};
for ii = 1:nC
    highlight(H,find(kidx==ii),'NodeColor',ccode{ii})
end
savefig('kmeans_vs_graph')

% Cross tab analysis
[cross_tab_,~,p] = crosstab(kidx,gp_id2);
fprintf('Cross tab p: %.4f',p)
% Normalize columns and plot across subregions
sum_tab = sum(cross_tab_,2);
cross_tab_col = cross_tab_./sum_tab;
figure('Position',  [100, 100, 400, 300])
b = bar([1:nC],cross_tab_col,'stacked');
if nanat==6
    ccode = {'y','r','k','b','g','c'};
else
    ccode = {'k','r','g','b'};
end
for ii = 1:length(unique(gp_id2))
    b(ii).FaceColor = ccode{ii};
end
legend(anat_names)
savefig('crosstab_anal')
print('crosstab_anal_stacked','-dpng')

%% 1.Component labels and error bars
err_bars = zeros(nC,size(behav.ratings,2));
meaner = zeros(size(err_bars));
c = categorical(behav.percepts);
figure('Position',  [100, 100, 900, 200])

if nC==6
    ccode = {'y','r','k','b','g','c'};
else
    ccode = {[1 0 1],[0 1 1],[1 1/2 0]};
end
% nC = 3;

c_mat = [];
m_mat = [];
p_count = [];

for ii = 1:nC
    node_count = graph_nodes(kidx==ii,:);
    [~,~,stats_] = anova1(node_count,[],'off');
    [c_mat(ii,:,:),m_mat(ii,:,:)] = multcompare(stats_);
    
    % Number of PCs
    [~,~,~,~,var] = pca(node_count);
    var = cumsum(var);
    p_count(ii) = find_nearest(var,70,true);
end

% Plot figures.
m_count = {};
for ii = 1:nC
    hold on
    %----------------- Anovas for significance of weights
    c_ = squeeze(c_mat(ii,:,:));
    m = squeeze(m_mat(ii,:,:));
    %     gcf
    %     hold on
    % Check means (m) with number of descriptors that are different.
    c_sq = squareform(c_(:,6));
    c_sq(c_sq<0.05)=nan;
    c_sq(~isnan(c_sq))=0;
    c_sq(isnan(c_sq))=1;
    c_sq = c_sq-eye(size(c_sq));
    
    [m_sort, m_ind] = sort(m(:,1),'descend');
    sig_vectors = c_sq(m_ind(end),:);
    sig_vectors = sig_vectors(m_ind);
    err_bars(ii,:) = 1.96*std(graph_nodes(kidx==ii,:))./sqrt(sum(kidx==ii)); % Standard error across voxels
    
    if nC>3
        subplot(2,nC/2,ii)
    else
        subplot(1,nC,ii)
    end
    hold on
    m_count{ii} = graph_nodes(kidx==ii,:);
    meaner(ii,:) = mean(graph_nodes(kidx==ii,:));
    %     bar(1:length(c), meaner(ii,:)',ccode{ii})
    %     errorbar(1:length(c), meaner(ii,:)',err_bars(ii,:)','k.')
    if plot_pos
        temp_mean = abs(meaner(ii,:)');
        bar(1:length(c),temp_mean(m_ind))
        temp_err = err_bars(ii,:)';
        errorbar(1:length(c), temp_mean(m_ind),temp_err(m_ind),'k.')
        yl = ylim;
        area(1:length(c),4*yl(2)*sig_vectors,'FaceColor','k','FaceAlpha',0.2)
    else
        bar(1:length(c),(meaner(ii,:)'),ccode{ii})
        errorbar(1:length(c), (meaner(ii,:)'),err_bars(ii,:)','k.')
    end
    
    xticks(1:length(c))
    xticklabels(c(m_ind))
    xtickangle(90)
    title(sprintf('cluster: %02d, num PC: %02d',ii, p_count(ii)))
    ylim(yl)
end
if plot_pos
    savefig('kmeans_bars_pos')
    print('kmeans_bars_pos','-dpng')
else
    savefig('kmeans_bars')
    print('kmeans_bars','-dpng')
end

%% 2(a) Complexity estimation
c = categorical(behav.percepts);
n_cutoff = 10000; % Bootstrap 1000 samples
ndescrip = 12; %min(size(behav.ratings,2),50);
size_control = true;
shuffler = false; % Shuffle the order of rows in each column.

% if rem_lowdeg
%     train_mat(degs2<deg_thresh,:,:) = [];
% end

p_count = zeros(n_cutoff,length(anat_masks)); % Num PC dimensions in each ROI
var1 = zeros(n_cutoff,length(anat_masks),ndescrip); % Scree plot for PCA in all ROI
coeffs_pc = zeros(n_cutoff,length(anat_masks),ndescrip,ndescrip); % Coefficients for PCA
% Compute statistics for Anova1
weights_roi = zeros(n_cutoff,18,length(anat_masks)); % Only for bootstrap
for nn = 1:length(anat_masks)
    if nn>=1
        temp_array = abs(train_mat(gp_id2==nn,:,:));
%     else
%         temp_array = abs(train_mat(gp_id2<=nn,:,:));
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
% nvoxes_s = 2-(nvoxes_s./50); % Only for group EM
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

M = [];
for x1 = 1:4
    for x2 = 1:4
        t1 = var_bar(:,x1);
        t2 = var_bar(:,x2);       
        M(x1,x2) = bstrap_hyp(t1,t2);
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
        err2 = std(bars_)./sqrt(size(bars_,1));
%         errs = 1.96*(err2);

    
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

