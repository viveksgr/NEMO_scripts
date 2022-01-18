function [rating, rt, sw] = detectability(n, t_total,mode)

if nargin<4
    cid = 0;
end

% n = index of each trial. Can be used to uniquely identify each run of
% this function.
% t_total = maximum time before this goes away
% mode = toggle debugging on/off
% Ask if the odor is detectable
% sw = 0 (right click) or 1 (left click)

t_start = time;
readmouse;
map = getmousemap;

coin_toss = randperm(2); % Shuffle mapping of responses to stimuli.
box_locations = [60,-60]; % Where to put the choice boxes
box_locations = box_locations(coin_toss);
results = {1, 0}; % Results
results = results(coin_toss);
sw = -1;

while 1
    setforecolour(1,1,1);
    clearpict(1);
    
    if mode==99
        cue_text = strcat(sprintf('Smell anything?_%d',n));
    else
        cue_text = strcat('Smell anything?');
    end
    
    settextstyle('Arial',30);
    preparestring(cue_text,1, 0, 60);
    
    setforecolour(1,1,1);
    preparestring('Yes',1, box_locations(1), 0);
    preparestring('No',1, box_locations(2), 0);
    
    setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    if mode==1
        preparestring('+',1,0,0);
    end
    
    drawpict(1);
    readmouse;

    if (time-t_start) < t_total
        if getmouse(map.Button1) == 128
            rating = results{2};
            sw = 1;
            break
        elseif getmouse(map.Button2) == 128
            rating = results{1};
            sw = 0;
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
    rt = time - t_start;
    setforecolour(0.7,0.7,0.7);
    settextstyle('Arial',30);
    preparestring('Yes',1, box_locations(1), 0);
    preparestring('No',1, box_locations(2), 0);
    drawpict(1);
    if mode==1
        wait(t_total-rt)
    end
end

clearpict(1)
end