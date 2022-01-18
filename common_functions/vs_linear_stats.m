function [r2val,p_val,t_val,b] = vs_linear_stats(y,x1,x2)
% Perform linear regression with y against two independent variables and
% find the t_values of betas.
% Reshape them as col vectors
% t_val = t values of weights
% p_val = p_val of weights
% r2_val = partial r^2 of weights 

y = y/std(y);
x1 = x1/std(x1);
x2 = x2/std(x2);

X = [ones(size(x1)) x1 x2];

[b, bint] = regress(y,X); 

berror = (bint(:,2)-bint(:,1))/3.92; % 95% CI is 1.96 std error

t_val(1) = b(2)/berror(2);
t_val(2) = b(3)/berror(3);
p_val(1) = 1-tcdf(t_val(1),length(x1)-1);
p_val(2) = 1-tcdf(t_val(2),length(x1)-1);

r2val = (t_val.^2)./(t_val.^2+length(x1)-1);