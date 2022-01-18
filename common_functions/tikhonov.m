function weight = tikhonov(X,y,lambda)

% Performs generalized ridge regression.
% Lambda = 2X1. Regularization parameters such that first two regressors
% are regularized with lambda(1) and rest with lambda(2)

num_ = 2;

alpha_mat = lambda(2)*eye(size(X,2));
inds = 1:num_; 
inds_ = sub2ind(size(alpha_mat),inds,inds);
alpha_mat(inds_) = lambda(1);

weight = (X'*X+alpha_mat'*alpha_mat)\(X'*y);