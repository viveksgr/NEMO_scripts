function [entropy_val] = entropy_vector(X)
% Computes the entropy value for each column of X. All cols are
% normalizaed between -1 and 1 for entropy calculation. 
% X = NxM array of M vectors of size N. entropy_val = 1XM vals of entropy
% values for each column.
num_col = size(X,2);
entropy_val = zeros(1,num_col);
p_basis = 0:0.01:1;
for ii = 1:num_col
    vect = X(:,ii);
    vect = vect-min(vect);
    vect = vect/max(vect);
    p_counts = histcounts(vect,p_basis);
    p_ = p_counts/sum(p_counts);
    p_(p_==0) = [];
    entropy_val(ii) = -sum(p_.*log(p_)/log(2));
end

