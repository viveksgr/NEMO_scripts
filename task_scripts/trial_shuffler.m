function percept_mat = trial_shuffler(lister, order, dtype)

% lister = ids of odor-repeats: 1:nRepeats
% order = pseudorandom order of odor repeats: 1:nRep*nStim
% percept_mat = lister expanded to 1:nRep*nStim such that
% percept_mat(order==odor) = lister. 

if strcmp(dtype,'cell')  
    percept_mat = cell(size(order));
elseif strcmp(dtype, 'mat')
    percept_mat = zeros(size(order));
end

odor_list = unique(order);
ntrials_per_odor = ceil(length(order)/length(odor_list));

if length(lister)~=ntrials_per_odor
    error('Number of perceptual descriptors is not appropriate')
end

for ii = 1:length(odor_list)
    shuffled = randperm(ntrials_per_odor);
    percept_mat(order==odor_list(ii))=lister(shuffled);
end
end

        

 