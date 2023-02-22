function p_val = bstrap_hyp(t1,t2)
% Assuming t1 and t2 have equal number of bootsamples. 
% std(t1) is the error in the sampling dist of the mean for var 1 across
% all bootsamples
% tdist is symmetric.

% t-stat
if isequal(t1,t2)
    p_val = 1;
else
t_sq = (mean(t1)-mean(t2))./sqrt(var(t1)+var(t2));

% Null distribution
grp1 = t1-mean(t1)+mean(t1+t2); % For group 1
grp2 = t2-mean(t2)+mean(t1+t2); % For group 2
tdist = (grp1-grp2)./sqrt(var(grp1)+var(grp2));

% one-tailed p-value t2>t1
p_val = invprctile(tdist,t_sq)/100; % one-sided

% Convert to two tailed
p_val2 = 1-(invprctile(tdist,t_sq)/100); 
p_val = 2*min(p_val,p_val2);
end


