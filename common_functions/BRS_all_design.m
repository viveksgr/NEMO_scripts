function [brs_onsets,brs_names] = BRS_all_design(odor_id,onset_set)

assert(length(odor_id)==length(onset_set)); % Odor ID and Odor onsets for each trial

brs_names = cell(1,length(unique(odor_id)));
brs_onsets = cell(1,length(unique(odor_id)));

odor_keys = unique(odor_id);
for ii = 1:length(odor_keys)
    mask = odor_id ==odor_keys(ii);
    brs_names{ii}=sprintf('%d',odor_keys(ii));
    brs_onsets{ii}=onset_set(mask);
end
end


