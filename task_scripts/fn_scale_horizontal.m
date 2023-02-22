function [rating, rt,sp] = fn_scale_horizontal(percept, t_time,waitmode)
% Wait for the scale to disappear after
if nargin<3
    waitmode = true;
end

% Rating scale for the given percept.
% coin_toss = 1;
coin_toss = randperm(2,1);
seq = [1:5;5:-1:1];
seq_this_trial = seq(coin_toss,:);

scalevalue = [-1, -0.5, 0, 0.5, 1];
max_x = 120;
phys_x_range = [-max_x, max_x];
ythick = 12; % Thickness of triangle
markers = scalevalue*max_x; % points on scale with labels

if strcmp(percept,' Intensity')
    labels = {'Weak'; ''; ''; ''; 'Strong'};
elseif strcmp(percept, ' Pleasantness')
    labels = {'Dislike'; ''; 'Neutral'; ''; 'Like'};
else
    labels = {'Not at all'; ''; 'Somewhat'; ''; 'Extremely'};
end

labels = labels(seq_this_trial);
scale_y_offset = -2; % rating scale offset on x axis

if length(markers)~=length(labels)
    error('Labels and markers need to have the same size')
end

% Get mouse map
readmouse;
map = getmousemap;

% if strcmp(percept,' Intensity')
%     xPos = 0;
% elseif strcmp(percept, ' Pleasantness')
%     xPos = 0;
% else
%     xPos = ((-1)^(coin_toss))*120;
% end

sp_scale = -120:10:120;
sp = sp_scale(randperm(25,1));
xPos = sp;

setforecolour(1,1,1);

offScale = 38;
mOffset = 10+ones(length(labels),1)*offScale;
tRatOns = time;
while 1
    clearpict(1);
    
    % Prepare rating scale
    cgpencol(1,1,1);
    cgpenwid(3);
    cgdraw(phys_x_range(1), scale_y_offset, phys_x_range(2), scale_y_offset);
    
    x_triangle = ([-max_x max_x max_x])*((-1)^(coin_toss+1));
    cgpolygon(x_triangle,[scale_y_offset scale_y_offset+ythick scale_y_offset-ythick])
    
    % prepare markers and labels
    for i=1:length(labels)
        cgfont('Arial', 20);
        cgtext(labels{i}, markers(i), -scale_y_offset-mOffset(i));
        cgfont('Arial', 40);
        preparestring('|',1, markers(i), 0);
    end
    
    % prepare curser
    setforecolour(0,0,0);
    cgfont('Arial', 60);
    if strcmp(percept,' Intensity')
        preparestring('|',1, xPos,0);
    elseif strcmp(percept, ' Pleasantness')
        preparestring('|',1, xPos,0);
    else
        preparestring('|',1, xPos,0);
    end
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
    xPos = xPos + sum(getmouse(map.X));
    
    % set bounds to position
    if xPos < phys_x_range(1)
        xPos = phys_x_range(1);
    elseif xPos > phys_x_range(2)
        xPos = phys_x_range(2);
    end
    
    % prepare curser - After selection
    setforecolour(0.3,0.3,0.3);
    cgfont('Arial', 60);
    preparestring('|',1, xPos,0);
    setforecolour(1,1,1);
     
    % check if ratings are submitted
    if getmouse(map.Button1) == 128
        rating = ((-1)^(coin_toss+1))*(xPos/(max_x/10)); % rating [range(1):range(end)];
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