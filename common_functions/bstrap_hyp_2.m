function p_val = bstrap_hyp_2(t1,t2)

if isequal(t1,t2)
    p_val = 1;
else
    t_sq = mean(t1)-mean(t2);
    tdist = (t1-t2)-mean(t1-t2);
    p_val = invprctile(tdist,t_sq)/100; % one-sided
    p_val2 = invprctile(tdist,-t_sq)/100;
    p_val = 2*min(p_val,p_val2);
end
