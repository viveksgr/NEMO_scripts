%% Estimation of subjectity of the encoding model 
% Compute the inter-subject variability of the encoding model across
% different anatomical areas

% Major inputs:
% Output file ('EEM_weights_boot.mat') of EEM_dimensionality.m
% behav.mat: behavioral file with behav.ratings = odors x perceptual bases
% Gray matter masks of ROIs
% anat_gw.nii: A binary mask of gray matter voxels such that number of
% voxels in fullmat.mat = sum(anat_gw,'all')
% swampsunset.m: A custom colormap

%% General settings
% Subject directories with output from EEM_dimensionality.m
root = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\NEMO_scripts\Data'; % Path to the Data folder
respath = fullfile(root,'EM_all');
% mkdir(respath)
dirs = {fullfile(root,'EM_01','S1_EM_LOOCV');
        fullfile(root,'EM_02','S2_EM_LOOCV');
        fullfile(root,'EM_03','S3_EM_LOOCV')};
%     
% dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\FIR_EM\deg_tscore_new2';
%     'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\FIR_EM\deg_tscore_new2';
%     'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\FIR_EM\deg_tscore_new2'};

matname = 'eem_weights_boot.mat';
sz_cntrl = true; % Control for size
anat_names = {'PirF','PirT','AMY','OFC'};
weights_ = variable_extract(dirs,matname,'bar_mat',true);

ndescrip = 18;
nanat = 4;
nS = 3; %length(dirs);
fname = 'ROI_generalizability_boot_corr';
nboot = 10000;

% Extrude common descriptors across subjects
inds = true(3,ndescrip);
inds(1,[11 16 18]) = false;
inds([2 3],[16:18]) = false;

% Alignment of descriptors across subjects
weights_norm = cell(size(weights_));
for ii = 1:nanat
    for jj = 1:nS
        weights_norm{jj,ii} = weights_{jj,ii}(:,inds(jj,:));
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
    
    % Correction for size
    nsam = 10000;
    Ms_temp = [];
    for data_id = 1:nsam
        Ms_id = datasample(1:length(Ms{ii}),2500);
        Ms_temp(data_id)=mean(Ms{ii}(Ms_id));
    end
            Ms{ii} = sign(Ms_temp').*(Ms_temp'.^2);
%     Ms{ii} = atanh(Ms_temp');
    
    mean_cons(ii) = mean(Ms{ii},'all');
    std_cons(ii) = prctile(Ms{ii},97.5)-mean_cons(ii);
    std_cons2(ii) = prctile(Ms{ii},2.5)-mean_cons(ii);
end

figure('Position',  [100, 100, 320, 240])
% bar(1:nanat,mean_cons)
% hold on 
% errorbar(1:nanat,mean_cons,std_cons,std_cons2,'.')
Ms_mat = cell2mat(Ms);
boxplot(Ms_mat);
xticks(1:4)
xticklabels(anat_names)
ylabel('Consistency')
ylim([0 inf])
% savefig(fname)
% print(fname,'-dpng')
% savefig(fname)
% print(fname,'-dpng')
% % Statistical testing
% M_test1(1) = bstrap_hyp(Ms{1},Ms{2});
% M_test1(2) = bstrap_hyp(Ms{2},Ms{3});
% M_test1(3) = bstrap_hyp(Ms{1},Ms{4});
% M_test1(4) = bstrap_hyp(Ms{2},Ms{4});
% M_test1(5) = bstrap_hyp(Ms{3},Ms{4});

%% PCA of OFC
% Only using perceptual descriptors common in all subjects
behav_labs = {' Intensity';' Pleasantness';' Fishy';' Burnt';' Sour';' Decayed';...
    ' Musky';' Fruity';' Sweaty';' Cool';'Floral';' Sweet';' Warm';' Bakery';' Spicy'};
load('C:\Data\NEMO\NEMO_all\EM\ROI_subjectivity_PC\ROI_settings.mat')
load('C:\Data\NEMO\swampsunset.mat')
num_pcs = 4;
ii = 4; % ID for OFC  = 4
% Alignment of PC components using Gale_Shapely
reorder_pc = true;

if reorder_pc
    w_cat =[];
    for jj = 1:length(dirs)
        w_cat= cat(1,w_cat, weights_{jj,ii}(:,inds(jj,:)));
    end
    [coeffs_cat,~,~,~,vars] = pca(w_cat);
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
colormap(swampsunset)
title(sprintf('Correlation PCs: %s',anat_names{ii}))
savefig('ofc_pc_mat')

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
figure('Position',  [100, 100, 400, 300])
bar(1:4,bars)
hold on
errorbar(1:4,bars,stds,'.')
for ii = 1:3
    plot(1:4,m_cell_(ii:3:end))
end
ylim([-0.2 0.8])
savefig('bars_q')
print('bars_q','-dpng')

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
savefig('tuning_curves')





