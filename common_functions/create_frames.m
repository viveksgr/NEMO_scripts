% Create matrix of frames for the task
% VivekSagar2016@u.northwestern.edu

% Non microtime resolved data in frame_struct
% Microtime data in frame_

statpath = 'C:\Data\NEET\Pilot';
cue_ids = 1:4; %{'A','E','N','c'};
nback_ids = 1:1:9;
n_trials = 8; % Number of trials with each copy of cue
nback_success = 0.1; % Fraction of times n_back matches
no_nback = false; % Create a task with no nback
exp_nitis_type = false; % Is the nback ITI histogram exponential or uniform

% Timing parameters - all parameters in sec
time_.cues_t = [7 7 12 12]; % Map-times for cues for peak reward
time_.micro_t = 1/10; % Screen refresh rate

time_.cue_start = 5; % Minimum ITI for cues. ITIs are periods of 0 rewards.
time_.cue_end = 10; % Maximum ITI for cues
time_.cue_taper = 2; % Total window of reward = 2*tapertime. Trial length = cues_t+taper+iti
time_.cue_shape = 1; % Mean of unnormalized exponential pdf, choose shape parameter less than iti_start for steeper pdf
time_.cue_dur = 1; % Cue duration
time_.cue_pause = 1; % Pause screen for this seconds for reward

% Similiar parameters for nback task
time_.back_start = 0.2;
time_.back_end = 2;
time_.back_shape = 1/4; % Mean ITI before clipping
time_.back_dur = 0.5; % 2-back number duration
time_.back_perfm = 100; % Check the nback performance over a moving window of time_.back_perfm seconds.
%% Create cue_stream
% Cue and nback streams
cue_stream_id = create_cue(1:length(cue_ids),n_trials);
cue_dur = [0 (max(time_.cues_t)+time_.cue_taper)*ones(1,length(cue_stream_id)-1)];%;[0 time_.cues_t(cue_stream_id(1:end-1))+time_.cue_taper]; % Comment out if different trials need to be of different duration.

% Derived paramaters
tn_trials = n_trials*length(cue_ids);
itis_uniform = linspace(time_.cue_start,time_.cue_end,1000);
exp_dist = pdf('exponential',itis_uniform,time_.cue_shape);
exp_dist = exp_dist/sum(exp_dist);
itis = inverse_hist(itis_uniform,exp_dist,tn_trials+1);

cue_start_diff = itis(1+end-1)+cue_dur; % Gaps between cue onset times
frame_struct.cue_start = cumsum(cue_start_diff);
frame_struct.cue_stream = cue_ids(cue_stream_id); % Stream of cues
run_length = (frame_struct.cue_start(end)+time_.cue_taper+3*itis(end)); % Total length of run time in seconds. Last trial runs for 3*iti(1).

%% Create nback stream
min_length = time_.back_dur+time_.back_start; % Minimum length of an nback trial
max_n_trials = ceil(run_length/min_length); % Maximum number of nback trials

% Create ITIs of nback

if exp_nitis_type
    nitis_uniform = linspace(time_.back_start,time_.back_end,1000);
    exp_dist = pdf('exponential',nitis_uniform,time_.back_shape);
    exp_dist = exp_dist/sum(exp_dist);
    nitis = inverse_hist(nitis_uniform,exp_dist,max_n_trials+1);
else
    nitis = time_.back_start+rand(1,max_n_trials+1)*(time_.back_end-time_.back_start);
end

nback_start = cumsum(nitis+time_.back_dur);
nend_id = find(nback_start>run_length,1,'first'); % Number of nback trials

frame_struct.nback_stream = create_nback(nback_ids, nend_id, ceil(nback_success*nend_id));
frame_struct.nback_start = nback_start(1:nend_id);

check_stream = double(~logical(circshift(frame_struct.nback_stream,2)-frame_struct.nback_stream));
frame_struct.check_nback = [0 0 check_stream(3:end)];

%% Reward profiles
reward_kernel_axis = -time_.cue_taper:time_.micro_t:time_.cue_taper;

% Design reward kernel
reward_kernel = pdf('normal',reward_kernel_axis,0,time_.cue_taper/3);
[~,idx_plateau_u] = min(abs(reward_kernel_axis-time_.cue_taper/3));
[~,idx_plateau_l] = min(abs(reward_kernel_axis+time_.cue_taper/3));
reward_kernel(idx_plateau_l:idx_plateau_u) = reward_kernel(idx_plateau_l);
frame_struct.reward_kernel = reward_kernel./max(reward_kernel);
frame_struct.reward_times = frame_struct.cue_start+time_.cues_t(cue_stream_id);

%% Create microsecond frames
% frame_.index = microsecond time (in s)
% frame_.cue = At the given time which cue, if any, is presented.
% frame_.nback = At the given time which nback number, if any is presented.
% frame_.reward = At the given time, what is value of the reward received
% frame_.cue_onset = true for the first frame when cue is presnted
% frame_.nback_onset = true for the first frame when nback is presnted
% frame_.val_response = frame for which nback response is valid

frame_.index = (0:time_.micro_t:run_length+time_.micro_t)';
frame_length = length(frame_.index);

% Cue stream
cue_start = find_nearest(frame_.index,frame_struct.cue_start);
cue_end = find_nearest(frame_.index,frame_struct.cue_start+time_.cue_dur);
frame_.cues_ = zeros(size(frame_.index));
for ii = 1:tn_trials
    size_cell = cue_end(ii)-cue_start(ii)+1;
    frame_.cues_(cue_start(ii):cue_end(ii))=repmat(frame_struct.cue_stream(ii),size_cell,1);
end

% Nback stream
nback_start = find_nearest(frame_.index,frame_struct.nback_start);
nback_end = find_nearest(frame_.index,frame_struct.nback_start+time_.back_dur);
frame_.nback = zeros(size(frame_.index));
for ii = 1:nend_id
    size_cell = nback_end(ii)-nback_start(ii)+1;
    frame_.nback(nback_start(ii):nback_end(ii))=repmat(frame_struct.nback_stream(ii),size_cell,1);
end

% Reward stream
reward_ms = zeros(size(frame_.index));
reward_ = find_nearest(frame_.index,frame_struct.reward_times);
reward_ms(reward_,:) = 1;
frame_.reward = conv(reward_ms,frame_struct.reward_kernel,'same');

% Cue_onset binary
frame_.cue_onset = false(size(frame_.index));
frame_.cue_onset(cue_start) = true;

% nback_onset binary
frame_.nback_onset = false(size(frame_.index));
frame_.nback_onset(nback_start) = true;

% nback_response
% Might sometimes get error if the last entry in check_nback = 1. Just
% rerun the script in that case.
frame_struct.check_nback(end) = false; % last iteration of nback anyway doesn't count.
response_start = nback_start(logical(frame_struct.check_nback));
response_end = nback_start([false logical(frame_struct.check_nback(1:end-1))]);
frame_.val_response = false(size(frame_.index));
for ii = 1:length(response_start)
    size_cell = response_end(ii)-response_start(ii);
    frame_.val_response(response_start(ii):response_end(ii)-1)=true(size_cell,1);
end

if no_nback
    frame_.nback = zeros(size(frame_.nback));
    frame_.nback_onset = false(size(frame_.nback));
    frame_.val_response = false(size(frame_.nback));
end

save(fullfile(statpath,'frames.mat'),'time_','frame_struct','frame_')
%% Test runs
% addpath('C:\Cogent2000v1.32\Toolbox');
% config_display(0, 3, [0.5, 0.5, 0.5], [1,1,1]); % cogent function. 4 indicates 4 buffers
% config_keyboard;
% config_mouse(10);
% start_cogent;
% cell_ = {'a','b','c','e'};
%
% for ii = 1:4
%     clearpict(1)
%     clearkeys
%     settextstyle('Agathodaimon', 200);
%     preparestring(cell_{ii},1,0,0);
%     drawpict(1);
%     wait(5000)
%     readkeys
%     [KeyResp{ii}, KeyTime{ii}, n{ii}] = getkeydown;
% end
% stop_cogent
%
% % a = zeros(1,5);
% % for ii = 1:5
% %     a(ii) = ii;
% %     if ii==4
% %         a(ii-2)
% %     end
% % end