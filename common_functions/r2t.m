function [r,t_st] = r2t(p,df)
% Value of pearson r needed to ensure p value with df-2 degree of
% freedom. Should be technically named p2r, but too much of a hassle to
% change name now in all functions.
t_st = tinv(p,df);
r = 1/sqrt(1+((df-2)/t_st^2));
end