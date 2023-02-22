%% Preprocessing of perceptual ratings.
% Compilation of behavioral dataset from scanning sessions for S1.

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

%% General
sn = 2;
set_i = 1;
set_f = 4;
nsets = set_f-set_i+1;
sess_i = 2;
sess_f = 4;
nsess = sess_f;
n_odors_sess = 40;

root = 'C:\Data\NEMO\';

% Each cell consists of results from different sets.
% set_detect = cell(1,nsets);
set_detect = [];

set_rt_detect = cell(1,nsets);
set_ratings = cell(nsets,nsess);
set_cid = cell(1,nsets);
for set_ = set_i:set_f
 
        setpath = fullfile(root, sprintf('NEMO_%02d',sn),...
        'behavior','imaging_task',sprintf('set_%02d',set_));
   
    % Initialize arrays to store sess data for this set.
    data_detect = repmat({cell(n_odors_sess,1)},1,nsess);
%     data_detect = repmat({zeros(n_odors_sess,1)},1,nsess); % For
%     detectability
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
                data_run{run_} = data_run_raw.res.select_odor;
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
                cid_key = 7165; % Compound replaced in odor set but not database
            end
            cid_val = cid_map(cid_key);

            cid_v = [cid_v cid_val];
            data_detect{sess}{cid_val} = [data_detect{sess}{cid_val} data_sess_raw(trial,6)];
            if cid_val ==6
                fprintf('%02d,%01d\n',trial,data_sess_raw(trial,6))
            end
            data_rt_detect{sess}(cid_val) = nansum([data_rt_detect{sess}(cid_val),data_sess_raw(trial,8)]);
            set_ratings{set_,sess}(cid_val,data_sess_raw(trial,5))=data_sess_raw(trial,10);
            

        end
    end
    

        set_det = horzcat(data_detect{:});
        set_det_mat = {};
%         for niter = 1:size(set_det,1); set_det_mat{niter} = horzcat(set_det{niter,:}); end
%         set_detect{set_} = cellfun(@(x) nansum(x)/27, set_det_mat);


    set_cid{set_} = cid_keys;
end

% clearvars -except set_ratings
% load('which_desc.mat')
%% Compile sets
sess1_ = 2;
sess2_ = 3;
% toggle_replacement=true;
rat1 = vertcat(set_ratings{:,sess1_});
rat3 = vertcat(set_ratings{:,sess2_});
cid_unsorted = vertcat(set_cid{:});
[cid,argsort] = sort(cid_unsorted);
rat1 = (rat1(argsort,:));
rat3 = (rat3(argsort,:));
rat1 = [rat1 zeros(size(rat1))];
rat_scan = rat1+rat3;

% reorder columns in-line with sequence of descriptors in behavioral task
rat_scan = rat_scan(:,[1 2 5:18 3 4]); 



