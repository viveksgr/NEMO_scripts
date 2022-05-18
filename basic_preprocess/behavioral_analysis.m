%% Preprocessing of perceptual ratings.
% Major inputs: behavioral ratings stored in the format specified by
% task_scripts>fn_select_odor_olf.m

% data_run = NX13, where N = number of trials in a given block and columns are:
% 1     odor order (1-10)
% 2     CID of odor
% 3     onset of odor trigger from t0
% 4     onset of sniff cue from t0
% 5     percept to be rated
% 6     Detect rating (0 = no smell)
% 7     Mouse button pressed
% 8     Detect RT
% 9     Time at which 9 submitted
% 10    Percept rating
% 11    Percept RT
% 12    Time at which 11 submitted
% 13    Initial starting pt of scale

% Major output = set_ratings. % Perceptual ratings for odors in a given
% session for IDs listed in set_cid.
%% Number of sets, sessions etc.
sn = 1;
set_i = 1;
set_f = 4;
nsets = set_f-set_i+1;
sess_ = [2:4];
nsess = length(sess_);
n_odors_sess = 40;


root = 'C:\Data\NEMO';
% Each cell consists of results from different sets.
% set_detect = cell(1,nsets);
set_detect = [];
set_rt_detect = cell(1,nsets);
set_ratings = cell(nsets,nsess);
set_cid = cell(nsets,nsess);

for set_ = set_i:set_f
    setpath = fullfile(root, sprintf('NEMO_%02d',sn),...
        'behavior','behavior_task',sprintf('set_%02d',set_));
    % Initialize arrays to store sess data for this set.
    data_detect = repmat({cell(n_odors_sess,1)},1,nsess);
%   data_detect = repmat({zeros(n_odors_sess,1)},1,nsess); % For
%   detectability
    data_rt_detect = repmat({zeros(n_odors_sess,1)},1,nsess);
    data_ratings = repmat({zeros(n_odors_sess,18)+NaN},1,nsess);
    for sess = sess_
        sesspath = fullfile(setpath,sprintf('sess_%02d',sess));
        run_files = dir(fullfile(sesspath, sprintf('NEMO*.mat')));
        data_run = cell(1,length(run_files));
        odornames = {};
        run_l = length(run_files);
        for run_ = 1:run_l
            runpath = fullfile(sesspath,run_files(run_).name);
            data_run_raw = load(runpath);
            
            if isfield(data_run_raw.res,'sp_mat')
                s = rmfield(data_run_raw.res,'sp_mat');
                data_run_raw.res = s; % No idea why this doesn't run inplace.
            end
            
            % Extract run data
            data_run{run_} = data_run_raw.res;
            try
                odorname = vertcat(data_run_raw.res.odornames(1:end-1,2));
                odornames = [odornames; odorname];
            catch
                odorname = {};
            end
        end
        % Collapse run data into one array for each session.
        data_sess_raw = cat(1,data_run{:});
        
        % Construct a dictionary to map each odor trial to correct row.
        try
        cid_keys = vertcat(data_sess_raw(:).odor_id);
        set_ratings{set_,sess}=vertcat(data_sess_raw(:).p_mat);    
        catch
            cid_keys = [];
            set_ratings{set_,sess} = [];
        end       
    % Trials from same odor are collated in the same row. 
    % Data from entire set in one cell
%     if ismember(set_,[1,2])
%         set_detect{set_} = nansum(horzcat(data_detect{:}),2)/n_trials12;
%         set_rt_detect{set_} = nansum(horzcat(data_rt_detect{:}),2)/n_trials12;
%     else
%         set_detect{set_} = nansum(horzcat(data_detect{:}),2)/n_trials34;
%         set_rt_detect{set_} = nansum(horzcat(data_rt_detect{:}),2)/n_trials34;
%     end
    set_cid{set_,sess} = cid_keys;
    names{set_,sess} = odornames;
    end
end

%% Test retest reliability for sets
% set_ = []; 
% % Corrections in database for sessions 1 and 3 for NEMO02, NEMO04
% set_cid{3,1}(set_cid{3,1}==1032)=1001;
% set_cid{3,1}(set_cid{3,1}==7122)=7335;
% set_cid{4,1}(set_cid{4,1}==2214)=6501;
% set_cid{4,1}(set_cid{4,1}==8174)=7937;
% set_cid{4,1}(set_cid{4,1}==8697)=556940;
% set_cid{4,1}(set_cid{4,1}==14525)=11569;

sess1_ = sess_(1);
sess2_ = sess_(2);
% toggle_replacement=true;
rat1 = vertcat(set_ratings{:,sess1_});
rat3 = vertcat(set_ratings{:,sess2_});
cid1 = vertcat(set_cid{:,sess1_});
cid3 = vertcat(set_cid{:,sess2_});
[cid,argsort1] = sort(cid1);
rat1 = (rat1(argsort1,:));
[~,argsort3] = sort(cid3);
rat3 = (rat3(argsort3,:));
assert(isequal(sort(cid1),sort(cid3)))

rat1(isnan(rat1))=0;
rat3(isnan(rat3))=0;

% test retest reliability for odors
for ii = 1:length(cid)
    nan_ = isnan(rat3(ii,:));
    corr__ = corrcoef(rat1(ii,~nan_),rat3(ii,~nan_)); 
    if length(corr__)>2
        corr_(ii) = corr__(2);
    end
end

% test retest reliability for descriptors
for ii = 1:size(rat1,2)
    nan_ = isnan(rat3(:,ii));
    corr__ = corrcoef(rat1(~nan_,ii),rat3(~nan_,ii)); 
    corr_2(ii) = corr__(2); 
end

