function [res] = fn_behav_training(res)

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

percept_list = {'Intensity'; 'Pleasantness'; 'Fishy'; 'Burnt'; 'Sour'; 'Decayed'; 'Musky'; 'Fruity'; 'Sweaty'; 'Cool'; 'Floral'; ...
    'Sweet'; 'Warm'; 'Bakery'; 'Spicy'; 'Ammonia'};
% percept_list = percept_list(1:res.nRep);

quality_list = {'Low','High'}; 

% Initialize all variables that are saved.
res.det =zeros(res.nStim,1); % reaction time of each rating
res.rt_det = zeros(res.nStim,1);        % reaction time of the detection task
res.rating = zeros(res.nStim,2);
res.rt_rating = zeros(res.nStim,2);
res.odor_onsets = zeros(res.nStim,2);   % time onset of each odor trial

% pseudo-random stimulus sequence
res.odor_id = cell2mat(odornames(1:end-1,1));

% Initialize Olfactometer
if res.run_odors
    daq = OlfConfigDaq;
    OlfFlowControl(daq, odornames{end,4}(1), odornames{end,4}(2)); % set MFCs
    OlfOpenLines(odornames{end,3}, daq, 16);
end

%% START
settextstyle('Arial', 30);
setforecolour(1,1,1);
% Wait to start the task

preparestring('You will be presented with 5 odors.', 1, 0, 75);
percept_id = 2*res.chunk_number+1;
preparestring('The odors have been usually described by high', 1, 0, 25);
str = sprintf('or low values of "%s" descriptor.',percept_list{percept_id});
preparestring(str, 1, 0, -10);
preparestring('Left click to continue', 1, 0, -110);
drawpict(1);
waitmouse(4);
clearpict(1)

str = sprintf('Smell the odors and observe how "%s" descriptor is rated.',percept_list{percept_id});
preparestring(str, 1, 0, 75);
str = sprintf('Try to interpret what "%s" descriptor means for you.',percept_list{percept_id});
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

% Initialize Olfactometer
if res.run_odors
    daq = OlfConfigDaq;
    OlfFlowControl(daq, odornames{end,4}(1), odornames{end,4}(2)); % set MFCs
    OlfOpenLines(odornames{end,3}, daq, 16);
end

for tt =1:res.nStim/2
    setforecolour(1,1,1);
    settextstyle('Arial', 40);
    str = sprintf('Odor%d: %s rating of %s',tt, quality_list{res.ratings(tt)+1}, percept_list{percept_id});
    preparestring(str, 1, 0, 25);
    preparestring('Left click to smell', 1, 0, -75);
    drawpict(1);
    waitmouse(4);
    drawpict(2);
    
    res.odor_onsets(tt,1) = time;
    odorname = odornames([tt size(odornames,1)],:)
    odor_present_backup(odorname,options); % Present the odor with a countdown
    wait(1000)
    % Perceptual ratings: different for each trial
%     [detect, rt_temp, ~] = detectability(tt,60000,res.mode); % Maximum time allowed = 1 min;
%     wait(1000)
%     res.det(tt) = detect;
%     res.rt_det(tt) = rt_temp;
    
    [percept_rating,rt] = fn_scale_horizontal(percept_list{percept_id},60000,false);
    res.rating(tt,1) = percept_rating;
    res.rt_rating(tt,1) = rt;
end

for tt =1:res.nStim/2
    setforecolour(1,1,1);
    settextstyle('Arial', 40);
    str = sprintf('Odor%d: %s rating of %s',tt, quality_list{res.ratings(tt)+1}, percept_list{percept_id});
    preparestring(str, 1, 0, 25);
    preparestring('Left click to smell', 1, 0, -75);
    drawpict(1);
    waitmouse(4);
    drawpict(2);
    
    res.odor_onsets(tt,2) = time;
    odorname = odornames([tt size(odornames,1)],:)
    odor_present_backup(odorname,options); % Present the odor with a countdown
    wait(1000)
    % Perceptual ratings: different for each trial
%     [detect, rt_temp, ~] = detectability(tt,60000,res.mode); % Maximum time allowed = 1 min;
%     wait(1000)
%     res.det(tt) = detect;
%     res.rt_det(tt) = rt_temp;
    
    [percept_rating,rt] = fn_scale_horizontal(percept_list{percept_id},60000,false);
    res.rating(tt,2) = percept_rating;
    res.rt_rating(tt,2) = rt;
end

settextstyle('Arial', 40);
preparestring('You will be presented with 5 odors.', 1, 0, 75);
percept_id = 2*res.chunk_number+2;
preparestring('The odors have been usually described by high', 1, 0, 25);
str = sprintf('or low values of "%s" descriptor.',percept_list{percept_id});
preparestring(str, 1, 0, -10);
preparestring('Left click to continue', 1, 0, -110);
drawpict(1);
waitmouse(4);
clearpict(1)

str = sprintf('Smell the odors and observe how "%s" descriptor is rated.',percept_list{percept_id});
preparestring(str, 1, 0, 75);
str = sprintf('Try to interpret what "%s" descriptor means for you.',percept_list{percept_id});
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

for tt =res.nStim/2+1:res.nStim
    setforecolour(1,1,1);
    settextstyle('Arial', 40);
    str = sprintf('Odor%d: %s rating of %s',tt, quality_list{res.ratings(tt)+1}, percept_list{percept_id});
    preparestring(str, 1, 0, 25);
    preparestring('Left click to smell', 1, 0, -75);
    drawpict(1);
    waitmouse(4);
    drawpict(2);
    
    res.odor_onsets(tt,1) = time;
    odorname = odornames([tt size(odornames,1)],:)
    odor_present_backup(odorname,options); % Present the odor with a countdown
    wait(1000)
    % Perceptual ratings: different for each trial
%     [detect, rt_temp, ~] = detectability(tt,60000,res.mode); % Maximum time allowed = 1 min;
%     wait(1000)
%     res.det(tt) = detect;
%     res.rt_det(tt) = rt_temp;
    
    [percept_rating,rt] = fn_scale_horizontal(percept_list{percept_id},60000,false);
    res.rating(tt,1) = percept_rating;
    res.rt_rating(tt,1) = rt;
end

for tt =res.nStim/2+1:res.nStim
    setforecolour(1,1,1);
    settextstyle('Arial', 40);
    str = sprintf('Odor%d: %s rating of %s',tt-5, quality_list{res.ratings(tt)+1}, percept_list{percept_id});
    preparestring(str, 1, 0, 25);
    preparestring('Left click to smell', 1, 0, -75);
    drawpict(1);
    waitmouse(4);
    drawpict(2);
    
    res.odor_onsets(tt,2) = time;
    odorname = odornames([tt size(odornames,1)],:)
    odor_present_backup(odorname,options); % Present the odor with a countdown
    wait(1000)
    % Perceptual ratings: different for each trial
%     [detect, rt_temp, ~] = detectability(tt,60000,res.mode); % Maximum time allowed = 1 min;
%     wait(1000)
%     res.det(tt) = detect;
%     res.rt_det(tt) = rt_temp;
    
    [percept_rating,rt] = fn_scale_horizontal(percept_list{percept_id},60000,false);
    res.rating(tt,2) = percept_rating;
    res.rt_rating(tt,2) = rt;
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

settextstyle('Arial', 30);
setforecolour(1,1,1);
preparestring('Please ring the bell.', 1, 0, -50);
drawpict(1);
waitkeydown(inf,59); % 59=enter
clearpict(1);

settextstyle('Arial', 30);
setforecolour(1,1,1);
preparestring('Please wait.', 1, 0, -50);
drawpict(1);
waitkeydown(inf,59); % 59=enter
clearpict(1);
