function [shuffed_vector] = vector_shuffler(vector,niter)
% Takes a 1XN vector and creates niterXN matrix of shuffled copies of the
% vector for permutation tests. Default niter = 1e4
shuffed_vector = zeros(length(vector),niter);
for ii = 1:niter
    shuffed_vector(:,ii) = vector(randperm(length(vector)));
end