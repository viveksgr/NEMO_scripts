function [brs_onsets,brs_names] = BRS_design_run(odor_id,onset_set,nsess)
% Valid for NEMO01 only
% odor_id = List of CID of odors
% onset_set = List of time stamps for each trial

assert(length(odor_id)==length(onset_set)); % Odor ID and Odor onsets for each trial

brs_names = cell(1,nsess*length(unique(odor_id)));
brs_onsets = cell(1,nsess*length(unique(odor_id)));

sess_counter_id = repmat([ones(1,4) 2*ones(1,4) 0*ones(1,4)],1,4);
sess_trials = repmat([100,90],4*2*nsess,1);
sess_trials = sess_trials(:);
sess_counter = [];
for tt = 1:length(sess_counter_id)
    sess_counter = [sess_counter; sess_counter_id(tt)*ones(sess_trials(tt),1)];
end

odor_keys = repmat(unique(odor_id),1,nsess)';
odor_keys = odor_keys(:);

for ii = 1:length(odor_keys)
    set_id = mod(ii,3);
    mask_od = odor_id ==odor_keys(ii);
    mask_set = sess_counter ==set_id;
    mask = and(mask_od,mask_set);
    brs_names{ii}=sprintf('%d_%d',odor_keys(ii),set_id);
    brs_onsets{ii}=onset_set(mask);
end
end



