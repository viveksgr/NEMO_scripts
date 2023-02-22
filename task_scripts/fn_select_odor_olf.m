% Modified July 2, 2019
function d = fn_select_odor_olf(res)

sub = res.subj;
clock_int = round(clock);
time_str = sprintf('%02.0f%02.0f',clock_int(4), clock_int(5));
subject = sprintf('Complete_%s_%s', sub, time_str);

%% Variables
odornames = res.odornames;
sniffDur = 2000;                    % Time for which the odor is presented
oddelay = res.delay;                % Delay between stimulus presentation and stimulus reception
pre_cue_dur = 1500;                 % Takes values in (0,2000)
Q1_dur = 2600;                      % Total time spent between the appearance of Ques1 and Ques2
Q1_wait = 800;                      % Total time spent between response of Ques 1 and appearance of Ques2 (<3000)
cueDur = sniffDur+oddelay+Q1_dur;   % Cue duration
nRep = res.nRep;                    % Number of times each odor is presented
nStim = res.nStim;                  % Number of odors presented in each run
ntrials = nRep*nStim;               % Total number of trials
% s = RandStream('mt19937ar','Seed',res.seed);
% RandStream for deterministic randomness

des_id=(1:1:nRep)';
if res.sess==1
    percept_list = {' Intensity';' Pleasantness';' Edible'; ' Familiar';...
    ' Fishy';' Burnt';' Sour'; ' Decayed'; ' Musky'};
    pd_id = (1:1:nRep)';
elseif res.sess ==2
    percept_list = {' Intensity';' Pleasantness';' Edible'; ' Familiar';...
       ' Fishy';' Burnt';' Sour'; ' Decayed'; ' Musky'};
    pd_temp = [1:9]';
    pd_id = pd_temp(1:nRep);
else
    percept_list = {' Fruity'; ' Sweaty'; ' Cool';' Floral';' Sweet'...
        ' Warm';' Bakery-like'; ' Spicy';' Ammonia'};
    pd_temp = [10:18]';
    pd_id = pd_temp(1:nRep);
end
%% Initialize

% Results table
d = zeros(ntrials,12);
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

% pseudo-random stimulus sequence
[~, order] = sort(rand(nStim,nRep));
order = order(:)';

d(:,1) = order;

% Odor IDs for each trial
CID_list = cell2mat(odornames(1:end-1,1));
d(:,2) = CID_list(order);

% Perceptual descriptor rated for each trial
QID = trial_shuffler(des_id, order, 'mat');
percept_mat = percept_list(QID);
d(:,5) = pd_id(QID);

% Randomize itis (uniform)
% itis = randi([8000,10000],[1,ntrials-1]);

% Randomize itis (exponential)
itis = trial_shuffler(res.pd, order, 'mat');
% pd2 = (11100+exprnd(1000,nRep,1))/1000;


% Initialize Olfactometer
if res.run_odors
    daq = OlfConfigDaq;
    OlfFlowControl(daq, odornames{end,4}(1), odornames{end,4}(2)); % set MFCs
    OlfOpenLines(odornames{end,3}, daq, 16);
end
% open clean air lines, send trigger

%% START
settextstyle('Arial', 30);
setforecolour(1,1,1);

% wait for triggerpulse
% if run in scanner, mode must be 1

preparestring('Please wait.', 1, 0, -50);
drawpict(1);
waitkeydown(inf,59);
clearpict(1);
drawpict(2)

wait(2000) % Small delay of 2 seconds
t0 = time;
if res.run_odors
    OlfOpenLines(odornames{end,3},daq,14);     % 10bit signal in LabChart
    wait(50);
    OlfOpenLines(odornames{end,3},daq,16);     % return to 16bit
end

% Ctrigger
Ctrigger = cumsum([t0+500 sniffDur + itis]);

for t =1:ntrials
    % Release Odor
    waituntil(Ctrigger(t));
    if res.run_odors
        OlfFlowControl(daq, odornames{d(t,1),4}(1), odornames{d(t,1),4}(2));
        OlfOpenLines(odornames{d(t,1),3}, daq, d(t,1));
    end
    tOdorTrigOns = time;
    d(t,3) = tOdorTrigOns-t0;   
    wait(oddelay-pre_cue_dur);
    
    % Fixation cross
    setforecolour(1,1,1);
    settextstyle('Arial',40);
    preparestring('+',2,0,0);
    drawpict(2);
    wait(pre_cue_dur);
    
    % Present sniff cue: "+" changes color
    setforecolour(0,0,1);
    settextstyle('Arial',40);
    preparestring('+',3,0,0);
    drawpict(3);
    clearpict(2);
    
    tsniffCueOns = drawpict(3);
    d(t,4) = tsniffCueOns-t0;
    
    % Set odor to clean
    if res.run_odors
        OlfFlowControl(daq,odornames{end,4}(1),odornames{end,4}(2));
        OlfOpenLines(odornames{end,3},daq,16);
    end
    
    wait(sniffDur);
    
    % Perceptual ratings: different for each trial
    [detect, rt, sw] = detectability(order(t),Q1_dur-Q1_wait,res.mode);
    d(t,6) = detect;
    d(t,7) = sw;
    d(t,8) = rt;
    if or(detect==0,res.sess==4) % Empty trials for session 4 incorporated
        d(t,10) = nan;
        d(t,11) = nan;
    else
        wait(Q1_wait);
        [d(t,10), d(t,11), d(t,13)] = fn_scale_horizontal(percept_mat{t},4000);
    end
    
    drawpict(4);
    
    % Record response time
    d(t,9) = d(t,3)+cueDur+d(t,8);
    d(t,12) = d(t,3)+cueDur+d(t,11);
end
wait(2000);
clearpict(4);

%% RESULTS
info_d = {...
    '1, Odor Order (1-10)';...
    '2, N/A in this experiment';...
    '3, Onset of odor trigger from t0';...
    '4, Onset of sniff cue from t0';...
    '5, Perceptual rating';...
    '6, Rating RT (ms)';...
    '7, Time of rating submitted from t0';...
    };

%% SAVE
if ~exist('res', 'dir')
    mkdir('res');
end

save(fullfile('res', sprintf('%s.mat', subject)));
