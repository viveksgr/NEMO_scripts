function coeff_correct = coch_orc(X,y,n_odors,interc)

% Computes X\y assuming that the error is serially correlated.
% Do not provide intercept in X
% November 2, 2018 @VivekSagar


if nargin <4
    interc = true;   
end

if interc
     X = [X ones(size(X,1),1)]; % Just one intercept
end

coeff_ =  X\y;
residue = y - X*coeff_;  % Residual

r = zscore(residue(1:end-1))'*zscore(residue(2:end))/length(y); % Coefficient for cochrane orcutt

y_t1 = y(1:end-1); % Series at t-1
y_t = y(2:end); % Series at t
y_star = y_t-r*y_t1;

x_t1 = X(1:end-1,:); % Regressors at t-1
x_t = X(2:end,:); % Regressors at t
x_star = x_t-r*x_t1;

% The value of the intercept is inaccurate.Technically, the constant
% regressor should be removed and then added again, but we don't need it
% anyway.

coeff_correct = x_star\y_star;
coeff_correct = coeff_correct(1:n_odors);
end




