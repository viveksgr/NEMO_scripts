function sw = task_switcher()
% Ask the user if they want to sniff the odor again or rate the odor on the
% perceptual task.
% sw = 1 for rating the odor
% sw = 0 for resmelling the odor
% sw = -1 for no response

t_total = 60000; % Give a minute to respond to the prompt
t_start = time;
readmouse;
map = getmousemap;

box_locations = [90,-90]; % Where to put the choice boxes
results = {1, 0}; % Results for the switcher task 
sw = -1; % for the possibility that the subject doesn't respond

% % If task is shuffled
% coin_toss = randperm(2); % Shuffle mapping of responses to stimuli.
% box_locations = [60,-60]; % Where to put the choice boxes
% box_locations = box_locations(coin_toss);
% results = {1, 0}; % Results
% results = results(coin_toss);
% sw = -1;

while 1
    clearpict(1);
       
    setforecolour(1,1,1);
    settextstyle('Arial',30);
    preparestring('Rate the odor',1, box_locations(1), 0);
    preparestring('Smell again',1, box_locations(2), 0);
    

    
    drawpict(1);
    readmouse;

    if (time-t_start) < t_total
        if getmouse(map.Button1) == 128 % Left button
            sw = results{2};
            break
        elseif getmouse(map.Button2) == 128
            sw = results{1};
            break
        end
    else
        sw = -1;
        rating = NaN;
        rt = NaN;
        break;
    end
end

if sw > -1
    setforecolour(0.7,0.7,0.7);
    settextstyle('Arial',30);
    preparestring('Rate the odor',1, box_locations(1), 0);
    preparestring('Smell again',1, box_locations(2), 0);
    drawpict(1);
    wait(500)
end

clearpict(1)
end
        
