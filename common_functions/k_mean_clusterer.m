function scores = k_mean_clusterer(mat)
% Compute within cluster - outside cluster correlations (scores) for a
% matrix mat with different number of clusters. Split the rows of mat using
% columns as features.
nmax = 10; % Max 10 clusters
scores = zeros(1,nmax-1);
for ii = 2:nmax
    kidx = kmeans(mat,ii,'distance','correlation');
    scores_on = zeros(1,ii);
    scores_off = zeros(1,nchoosek(ii,2));
    c = 0;
    for c_i = 1:ii
        mat_1 = mat(kidx==c_i,:);
        mat_on = corrcoef(mat_1');
        utl = logical(triu(ones(size(mat_on)),1));
        scores_on(c_i) = mean(atanh(mat_on(utl))); % On diagonal clusters for a given ii.
        for c_j = c_i+1:ii
            c = c+1;
            mat_2 = mat(kidx==c_j,:);
            mat_off = corrcoef_2(mat_1,mat_2);
            scores_off(c) = mean(atanh(mat_off(:))); 
        end
    end
    scores(ii-1) = (mean(scores_on)-mean(scores_off));
end
    