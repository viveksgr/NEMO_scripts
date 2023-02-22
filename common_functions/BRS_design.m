function [brs_onsets,brs_names] = BRS_design(odor_id,onset_set,num_rep)

if nargin<3
    num_rep = 10; % Number of times an odor is repeated in a single set.
end

assert(length(odor_id)==length(onset_set)); % Odor ID and Odor onsets for each trial
brs_names = cell(1,length(odor_id)/num_rep);
brs_onsets = cell(1,length(odor_id)/num_rep);

odor_keys = unique(odor_id);
for ii = 1:length(odor_keys)
    mask = odor_id ==odor_keys(ii);
    brs_names{ii}=sprintf('%d',odor_keys(ii));
    brs_onsets{ii}=onset_set(mask);
end
end


