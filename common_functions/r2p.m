function [p,t_sq] = r2p(r,df)
% One tailed p_val for a given r (only use in EM). Accepts r and df as vectors, unless p is needed.
t_sq = (r.^2/(1-r.^2)).*(df-2);
p = 1-tcdf(sqrt(t_sq),df); % Not vectorized.
end