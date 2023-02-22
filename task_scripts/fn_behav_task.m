function [res] = fn_behav_task(res)

sub = res.subj;
clock_int = round(clock);
time_str = sprintf('%02.0f%02.0f',clock_int(4), clock_int(5));
subject = sprintf('Complete_%s_%s', sub, time_str);

%% Timing variables
odornames = res.odornames;
% Time for which the odor is presented

options.delay = res.delay;
options.run_odors = res.run_odors;
options.sniffDur = res.sniffDur;

percept_list = {' Intensity';' Pleasantness';' Fishy';' Burnt';...
    ' Sour'; ' Decayed'; ' Musky';' Fruity'; ' Sweaty'; ' Cool';' Floral';...
    ' Sweet'; ' Warm';' Bakery-like'; ' Spicy';' Ammonia'; ' Edible'; ' Familiar'};
percept_list = percept_list(1:res.nRep);

% Results table (res):
% .odor_id =      CID of odor
% .odor_onsets =      time trial onset of odor
% .p_      = perceptual ratings
% .rt = reaction time

% Initialize all variables that are saved
res.p_mat = zeros(res.nStim,length(percept_list)); % percept_ratings. Acquired if and only if odors are detectable.
res.rt_mat = zeros(size(res.p_mat)); % reaction time of each rating
res.rt_det = zeros(res.nStim,2); % reaction time of the detection task
res.sp_mat = zeros(size(res.p_mat)); % Starting pt of scale
res.odor_onsets = zeros(res.nStim,1); % time onset of each odor trial

% pseudo-random stimulus sequence
res.odor_id = cell2mat(odornames(1:end-1,1));
[~,arg_percept] = sort(rand(length(percept_list),res.nStim));
res.arg_percept = arg_percept; % Order of perceptual descriptors

% Initialize Olfactometer
if res.run_odors
    daq = OlfConfigDaq;
    OlfFlowControl(daq, odornames{end,4}(1), odornames{end,4}(2)); % set MFCs
    OlfOpenLines(odornames{end,3}, daq, 16);
end

%% START
settextstyle('Arial', 40);
setforecolour(1,1,1);

% Wait to start the task
if res.chunk_number~=res.start_run
    preparestring('Please ring the bell.', 1, 0, -50);
    drawpict(1);
    waitkeydown(inf,59); % 59=enter
end
t0 = time;
clearpict(1);
res.odor_counts = zeros(1,res.nStim);

preparestring('You will be presented with 10 odors.', 1, 0, 75);
preparestring('Please describe the odor in terms of perceptual ratings.', 1, 0, 25);
preparestring('Left click to continue', 1, 0, -110);
drawpict(1);
waitmouse(4);
clearpict(1)

str = sprintf('Press left mouse button to smell the odor again.');
preparestring(str, 1, 0, 25);
preparestring('You do not have to memorize the ratings.', 1, 0, -25);
preparestring('Left click to continue', 1, 0, -125);
drawpict(1);
waitmouse(4);
clearpict(1)

preparestring('Odor will be presented after a countdown of 3.', 1, 0, 25);
preparestring('Smell the odor when the blue fixation appears.', 1, 0, -25);
preparestring('Use the mouse to make responses.', 1, 0, -75);
preparestring('Left click to start', 1, 0, -125);
drawpict(1);
waitmouse(4);
clearpict(1)

for tt =1:res.nStim
    counter = 0;
    setforecolour(1,1,1);
    settextstyle('Arial', 40);
    preparestring('Left click to smell new odor', 1, 0, -50);
    drawpict(1);
    waitmouse(4);
    drawpict(2);
    
    res.odor_onsets(tt) = time;
    odorname = odornames([tt size(odornames,1)],:);
    odor_present_backup(odorname,options); % Present the odor with a countdown
    wait(1000)
    % Perceptual ratings: different for each trial
    [detect, rt_temp, ~] = detectability(odornames{tt,1},60000,res.mode); % Maximum time allowed = 1 min;
    wait(1000)
    res.rt_det(tt,1) = rt_temp;
    
    if detect==0 % No detection
        clearpict(1)
        drawpict(2)
        wait(1000)
        odor_present_backup(odorname,options); % Present the odor with a countdown
        counter = counter+1;
        wait(1000)
        [detect, rt_temp, ~] = detectability(odornames{tt,1},60000,res.mode); % Maximum time allowed = 1 min;
        res.rt_det(tt,2) = rt_temp;
    end
    
    if detect==0 % No detection even second time
        res.p_mat(:,tt) = NaN;
        res.rt_mat(:,tt) = NaN;
    else % Odor detected or no response
        pp = 1;
        while (pp <= (length(percept_list)))
            
            pp_id = arg_percept(pp,tt); % Perceptual ratings are shuffled.
            switcher = task_switcher(); % Smell the odor or rate the odor
            
            if switcher==1 % rate the odor
                [res.p_mat(tt,pp_id), res.rt_mat(tt,pp_id), res.sp_mat(tt,pp_id)] = fn_scale_horizontal(percept_list{pp_id},60000,false);
                pp = pp+1;
                clearpict(1)
                drawpict(2)
                wait(1000)
            else
                counter = counter+1;
                odor_present_backup(odorname,options); % Smell the odor again (even if no response)
            end
        end
    end
    res.odor_counts(tt) = counter;
    drawpict(4);
end
wait(5000);
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
