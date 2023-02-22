%% Preprocessing of perceptual ratings.
% Compilation of behavioral dataset from scanning sessions for S1. Do not
% use this for S2 and S3.

% Major inputs: behavioral ratings stored in the format specified by
% task_scripts>fn_select_odor_olf.m

% data_run = NX13, where N = number of trials in a given block and columns
% are: 1     odor order (1-10) 2     CID of odor 3     onset of odor
% trigger from t0 4     onset of sniff cue from t0 5     percept to be
% rated 6     Detect rating (0 = no smell) 7     Mouse button pressed 8
% Detect RT 9     Time at which 9 submitted 10    Percept rating 11
% Percept RT 12    Time at which 11 submitted 13    Initial starting pt of
% scale

% Major output = rat_scan % Perceptual ratings for odors acquired in the
% scanner for S1.

%% Number of sets, sessions etc.
root = 'C:\Data\NEMO\'; % Or the root directory in which NEMO data is stored.

sn = 1;
set_i = 1;
set_f = 4;
nsets = set_f-set_i+1;
sess_i = 2;
sess_f = 4;
nsess = sess_f;
n_odors_sess = 40;

% Each cell consists of results from different sets.
set_detect = [];
set_rt_detect = cell(1,nsets);
set_ratings = cell(nsets,nsess);
set_cid = cell(1,nsets);
for set_ = set_i:set_f   
    setpath = fullfile(root, sprintf('NEMO_%02d',sn),...
        'behavior',sprintf('set_%02d',set_));   
    % Initialize arrays to store sess data for this set.
    data_detect = repmat({cell(n_odors_sess,1)},1,nsess);
    data_rt_detect = repmat({zeros(n_odors_sess,1)},1,nsess);
    data_ratings = repmat({zeros(n_odors_sess,18)+NaN},1,nsess);
    for sess = sess_i:sess_f
        sesspath = fullfile(setpath,sprintf('sess_%02d',sess));
        run_files = dir(fullfile(sesspath, sprintf('NEMO*.mat')));
        data_run = cell(1,length(run_files));
        for run_ = 1:length(run_files)
            runpath = fullfile(sesspath,run_files(run_).name);
            data_run_raw = load(runpath);
            % Extract run data
            
            data_run{run_} = data_run_raw.res.select_odor{1};
            
        end
        % Collapse run data into one array for each session.
        data_sess_raw = cat(1,data_run{:});
        
        % Construct a dictionary to map each odor trial to correct row.
        if sess==sess_i
            cid_keys = sort(unique(data_sess_raw(:,2)));
            assert(length(cid_keys)==40)
            
            cid_vals = (1:1:length(cid_keys))';
            cid_map = containers.Map(cid_keys, cid_vals);
        end
        
        % Trials from same odor are collated in the same row.
        cid_v = []; % For debugging something       
        for trial = 1:length(data_sess_raw)
            cid_key = data_sess_raw(trial,2);
            if cid_key==220674
                fprintf('Note correction in ID %d%d%d\n',set_,sess,run_)
                cid_key = 7165; % If compound was replaced in experimental set but not database
            end
            cid_val = cid_map(cid_key);
            cid_v = [cid_v cid_val];
            data_detect{sess}{cid_val} = [data_detect{sess}{cid_val} data_sess_raw(trial,6)];
            data_rt_detect{sess}(cid_val) = nansum([data_rt_detect{sess}(cid_val),data_sess_raw(trial,8)]);
            set_ratings{set_,sess}(cid_val,data_sess_raw(trial,5))=data_sess_raw(trial,10);
        end
    end
    set_cid{set_} = cid_keys;
end

%% Compile data across sets
set_ratings{3,1}(:,end+1)=nan;
set_ratings{4,1}(:,end+1)=nan;
set_ratings{3,2}(:,end+1)=nan;
set_ratings{4,2}(:,end+1)=nan;
set_ratings{1,3}(:,end+1:end+4)=nan;
set_ratings{2,3}(:,end+1:end+4)=nan;
rat(:,:,1) = [vertcat(set_ratings{:,1}) nan*ones(160,8)];
rat(:,:,2) = vertcat(set_ratings{:,2});
rat(:,:,3) = vertcat(set_ratings{:,3});
rat(rat==0) = nan;
ratings = nanmean(rat,3);

cid = vertcat(set_cid{:,:});
[~,argsort] = sort(cid);
rat_scan = ratings(argsort,:);

% The missing data is linearly interpolated from the DREAM dataset.

