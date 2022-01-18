%% Bar plots of anatomical connectivity in encoding model
% Load anatdist from encoding model k-means file
p(1) = 70125;
p(2) = 110685;
p(3) = 26106;
r_chance = []; for ii = 1:3; r_chance = [r_chance r2t(0.05,p(ii))]; end
figure('Position',  [100, 100, 400, 300])
bar([1:3],[anatdist1 anatdist2 anatdist3]);
hold on
for ii = 1:3
    plot([ii-0.5,ii+0.5],[r_chance(ii),r_chance(ii)],'k');
end
xticks(1:3)
xticklabels({'S01','S02','S03'})
savefig('Anatcluster')
print('Anatcluster','-dpng')

%% Comparision of k-means clustering on graph:
% ii = 0;
% k_C_mat = [];

% Hop to each subject's folder and run this section ---------------------
nS = 3;
nC = 3;
nanat = 6;
statpath = 'C:\Data\NEMO\NEMO_all\EM\kmeans';
anat_or_func = false;
funcpos = false; % Take absolute value of clustering

if anat_or_func
    %--------------------------------------------------------------------
    ii = ii+1;
    load('kmeans.mat','kidx','gp_id2')
    a_c = crosstab(kidx,gp_id2);
          
    if size(a_c )<nanat
       a_c (:,nanat) = 0;
    end
    
    sum_tab = sum(a_c,2);
    a_c = a_c./sum_tab;
    ac_mat(:,:,ii) = a_c;  

    %--------------------------------------------------------------------
    k_C_mat = ac_mat;
    k_C_mat_ = permute(k_C_mat,[1 3 2]);
    k_C_mat_ = reshape(k_C_mat_,[],size(k_C_mat_,3),1);
    k_C_corr = corrcoef(k_C_mat_');    
else  
    %--------------------------------------------------------------------
    ii = ii+1;
    load('kmeans.mat','k_C')
    if ii==1
        k_C(:,[11 17 18]) = [];
    else
        k_C(:,[16:18])=[];
    end
    k_C_mat(:,:,ii) = k_C;
    if funcpos
        k_C_copy = k_C_mat;
        k_C_mat = abs(k_C_mat);
    end
    %--------------------------------------------------------------------
    % k_C_mat(nC+1:end,:,:) = [];
    % Compute correlations
    k_C_mat_ = permute(k_C_mat,[1 3 2]);
    k_C_mat_ = reshape(k_C_mat_,[],size(k_C_mat_,3),1);
    k_C_corr = corrcoef(k_C_mat_');
end

% Clusters are together
mat = zeros(nC, nS);
for ii = 1:numel(mat); mat(ii) =ii; end
sorter = reshape(mat',[],1);
labels = cell(size(sorter));
kk = 0;
for ii = 1:nC
    for jj = 1:nS
        kk = kk+1;
        labels{kk} = sprintf('S%01d C%01d',jj,ii);
    end
end
k_C_corr_ = k_C_corr(sorter,sorter);

figure('Position',  [100, 100, 400, 300])
imagesc(k_C_corr_)
xticks(1:length(sorter))
xticklabels(labels)
xtickangle(90)
yticks(1:length(sorter))
yticklabels(labels)
colorbar
if anat_or_func
    savefig(fullfile(statpath ,'anat correlations'))
    print(fullfile(statpath,'anat correlations'),'-dpng')
else
    savefig(fullfile(statpath ,'Component correlations'))
    print(fullfile(statpath,'Component correlations'),'-dpng')
end

% bar plots for 2X2 matrices
block = [];
clustermat = ones(nS);
for ii = 1:nC
    block = blkdiag(block,clustermat);
end
block = logical(block);
utl = logical(triu(ones(size(block)),1));

M_on = and(block,utl);
M_off = and(~block,utl);
M = k_C_corr_ ;
m_on = M(M_on);
m_off = M(M_off);
fprintf('p_val of difference: %.04f',ranksum(m_on,m_off))

figure('Position',  [100, 100, 400, 300])
bar([1:2],[mean(m_on) mean(m_off)])
hold on
errorbar([1:2],[mean(m_on) mean(m_off)],[std(m_on)./sqrt(length(m_on)) std(m_off)./sqrt(length(m_off))],'k.')
xticks(1:2)
xticklabels({'Within Cluster','Cross-cluster'})

if anat_or_func
    savefig(fullfile(statpath,'anat_bars'))
    print(fullfile(statpath,'anat_bars'),'-dpng')
else
   savefig(fullfile(statpath,'func_bars'))
    print(fullfile(statpath,'func_bars'),'-dpng')
end
save('funcmat.mat')