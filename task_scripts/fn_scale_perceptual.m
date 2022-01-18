function [rating, rt] = fn_scale_perceptual(percept, t_time,waitmode)
% Wait for the scale to disappear after
if nargin<3
    waitmode = true;
end

% Rating scale for the given percept.
coin_toss = randperm(2,1);
seq = [1:5;5:-1:1];
seq_this_trial = seq(coin_toss,:);

scalevalue = [-1, -0.5, 0, 0.5, 1];
max_y = 120;
phys_y_range = [-max_y, max_y];
markers = scalevalue*max_y; % points on scale with labels

if strcmp(percept,' Intensity')
    labels = {'Very weak'; ''; ''; ''; 'Very strong'};
elseif strcmp(percept, ' Pleasantness')
    labels = {'Dislike'; ''; 'Neutral'; ''; 'Like very much'};
else
    labels = {'Not at all'; ''; 'Somewhat'; ''; 'Extremely'};
end

labels = labels(seq_this_trial);
scale_x_offset = 0; % rating scale offset on x axis

if length(markers)~=length(labels)
    error('Labels and markers need to have the same size')
end

% Get mouse map
readmouse;
map = getmousemap;

yPos = 0;
setforecolour(1,1,1);

offScale = 38;
mOffset = 50+ones(length(labels),1)*offScale;
tRatOns = time;
while 1
    clearpict(1);
    
    % Prepare rating scale
    cgpencol(1,1,1);
    cgpenwid(3);
    cgdraw(scale_x_offset, phys_y_range(1), scale_x_offset, phys_y_range(2));
    
    % prepare markers and labels
    for i=1:length(labels)
        cgfont('Arial', 23);
        cgtext(labels{i}, scale_x_offset+mOffset(i), markers(i));
        cgfont('Arial', 30);
        preparestring('--',1, scale_x_offset, markers(i));
    end
    
    % prepare curser
    setforecolour(0,0,0);
    cgfont('Arial', 60);
    preparestring('--',1, scale_x_offset, yPos);
    setforecolour(1,1,1);
    
    % add title at top
    setforecolour(1,1,1);
    cgfont('Arial', 40);
    title_text = percept;
    preparestring(title_text,1, 0, 200);
       
    % present scale
    drawpict(1);
    
    % read responses and update position
    readmouse;
    yPos = yPos - sum(getmouse(map.Y));
    % xPos = round(xPos/stepsize)*stepsize; % get a non-contineous rating for Likert-scales
    
    % set bounds to position
    if yPos < phys_y_range(1)
        yPos = phys_y_range(1);
    elseif yPos > phys_y_range(2)
        yPos = phys_y_range(2);
    end
    
    % prepare curser - After selection
    setforecolour(0.3,0.3,0.3);
    cgfont('Arial', 60);
    preparestring('--',1, scale_x_offset, yPos);
    setforecolour(1,1,1);
     
    % check if ratings are submitted
    if getmouse(map.Button1) == 128
        rating = ((-1)^(coin_toss+1))*(yPos/(max_y/10)); % rating [range(1):range(end)];
        rt = time - tRatOns; % RT of rating [ms]
        drawpict(1);
        if waitmode
            wait(t_time-rt);
        end
        break;
    elseif (time-tRatOns) > t_time % too slow
        rating = NaN; % rating [range(1):range(end)];
        rt = NaN; % RT of rating [ms]
        break;
    end
    
end
clearpict(1)
end