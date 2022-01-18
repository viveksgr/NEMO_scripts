% NEMO_01 - 
% Load data
mdfa = mdf{:,2:end};
mdf_cid = mdf{:,1};
sdfa = sempruned{:,2:end};
sdf_cid = sempruned{:,1};
mdf_cid_log = ismember(mdf_cid,behav.cid);
sdf_cid_log = ismember(sdf_cid,behav.cid);
percept = sdfa(sdf_cid_log,:);
chemical = mdfa(mdf_cid_log,:);

corr_ = zeros(size(percept,2),size(chemical,2));
for pp = 1:size(percept,2)
    for cc = 1:size(chemical,2)
        temp = corrcoef(percept(:,pp),chemical(:,cc));
        corr_(pp,cc) = temp(2);
    end
end
%        
% % Threshold
% thresh_ = r2t(0.001/size(percept,2),size(percept,1));
% corr_logical = corr_abs>thresh_;
% thresh_col = sum(corr_logical,1)>0;
% 
% [coeff, pca_r,~,~,var] = pca(chemical(:,thresh_col));
% cum_var = cumsum(var);
% ind_pca = find(cum_var>100*pca_threshold);

% Old chem basis, chem_ is the pruned version of chemical
% %% Choose subset of predictive chemical props
% n_val = zeros(1,size(chem_,2));
% for ii = 1:length(n_val)
%     n_val(ii)=length(unique(chem_(:,ii)));
% end
% 
% n_val_vec = n_val>3;
% behav.ratings = chem_(:,n_val_vec);
% 
% percepts = mdfprunednew.Properties.VariableNames(2:end)';
% behav.percepts = percepts(n_val_vec);

%% MI to find descriptors
MI_ = zeros(size(percept,2),size(chemical,2));
for pp = 1:size(percept,2)
    pp
    for cc = 1:size(chemical,2)
        temp = mi(percept(:,pp),chemical(:,cc));
        MI_(pp,cc) = temp;
    end
end

chem_score = max(MI_,[],1);
[~,argsort] = sort(chem_score,'descend');

num_feat = 120;
chem_pruned = chemical(:,argsort(1:num_feat));
behav.ratings = chem_pruned;
percepts = argsort(1:num_feat);
behav.percepts = cell(1,num_feat);
for ii = 1:num_feat
    behav.percepts{ii}=percepts(ii);
end

%% Entropy, PCA
% Chemical descriptor analysis. Choose only the most informative chemical
% descriptors.

group_edges = [1,45,77,150,196,233,281,831,1044,1140,1185,1208,1532,1570,1660,1740,1950,2174,2288,2561,2602,2756,2871,3041,3191,4787,4823,4843,4869];
% If reading directly from dragon data, add 1. CID column removed.
chem_mat = [];
chem_id = [];
for ii = 1:length(group_edges)-1
    group = chemical(:,group_edges(ii):group_edges(ii+1));
    group_id = group_edges(ii):group_edges(ii+1);
    ent = entropy_vector(group);
    [~,argsort] = sort(-ent);
    group = group(:,argsort);
    group_id = group_id(:,argsort)+1; % Dragon indices and data file indices differ by 1. 
    chem_mat = [chem_mat group(:,1:2)];
    chem_id = [chem_id; group_id(1:2)'];
end
behav.ratings = chem_mat;
behav.percepts = mat2cell(chem_id',1,ones(1,56))';

% PCA analysis
chem_pca = [];
pca_c = [];
for ii = 1:length(group_edges)-1
    group = chemical(:,group_edges(ii):group_edges(ii+1));
    [~,sc,~,~,c]=pca(group);
    chem_pca = [chem_pca sc(:,1:2)];
    pca_c = [pca_c c(1:15)];
end
behav.ratings = chem_pca;
behav.percepts = mat2cell(1:56,1,ones(1,56))';

%% Choose chemical descriptors
group_edges = [1,45,77,150,196,233,281,831,1044,1140,1185,1208,1532,1570,1660,1740,1950,2174,2288,2561,2602,2756,2871,3041,3191,4787,4823,4843,4869];
% If reading directly from dragon data, add 1. CID column removed.
X = [];
for ii = 1:length(group_edges)-1
    group = chemical(:,group_edges(ii):group_edges(ii+1));
    [wt, sc] = pca(group);
    end_ind = min(size(group,2),3);
    X = [X sc(:,1:end_ind)];
end

% Take linear model
% fold_ind = crossvalind('Kfold',160,4);
% krange = linspace(0.001,0.3,50);
% wts = zeros(size(X,2),size(percept,2));
% for pp = 1:size(percept,2)
%     pp
%     y = percept(:,pp);
%     lambda = max_lambda(y,X,fold_ind,krange,1);
%     wts(:,pp) = lasso(X,y,'lambda',lambda);
% end
% wts_binary = logical(wts);
% wts_score = sum(wts_binary,2);
% 
% X2 = X(:,wts_score>0);

% Choose chemical to predict percept
X2 = X;

entropy_sc = entropy_vector(X2);
entropy_ind = entropy_sc>5;
X3 = X2(:,entropy_ind);

% To get 70: X1 = 4 components, X2 = [], X3 = entropy (5);
% To get 49: X1 = 3 components, X2 = [], X3 = entropy (5);
% To get 73: X1 = 4 or 5, X2 = true, X3 = 4,5 or 5.5

%% Check if binary descriptors are correlated with intensity and pleasantness
corr = zeros(2,19);
for ii = 2:19
    corr(1,ii) = fastcorr(behav.ratings(:,2),behav.ratings(:,ii));
    corr(2,ii) = fastcorr(behav.ratings(:,3),behav.ratings(:,ii));
end

%% NC and functional groups and mo wt.
func_id = [2603:2756]; %ID of func groups
func_groups = cmdf_log(:,func_id);
func_bin = logical(func_groups);
idx = ~logical(sum(func_bin)); % Useless groups
func_bin_ = func_bin(:,~idx);

mw = exp(cmdf_log(:,3))-1;
nc = exp(cmdf_log(:,25))-1;

