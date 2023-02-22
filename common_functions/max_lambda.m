function [lambda] = max_lambda(y,X,fold_ind,krange,method)
% Figure out the regularization hyperparameter in lasso or ridge regression
% y = Nx1 labels
% X = NxM predictors
% fold_ind = NX1 vector of cross validation indices created using function
% crossvalind
% krange = range of hyperparameter values to be tried
% method = 1 for L1 lasso, 2 for L2 ridge.
% lambda  = 1x1 value of the optimized hyperparameter for L1 and L2. 2X1
% values of parameters for Elastic Net

num_folds = length(unique(fold_ind));

switch method
    case 1 % Lasso
        scores = zeros(length(krange),1);
        for fold = 1:num_folds
            test_mask = fold_ind==fold;
            X_train = X(~test_mask,:);
            X_test = X(test_mask,:);
            y_train = y(~test_mask);
            y_test = repmat(y(test_mask),1,length(krange));
            wts = lasso(X_train,y_train,'lambda',krange);
            y_pred = X_test*wts;
            scores = scores+iter_corr(y_test',y_pred');
            [~,argmax]=max(scores);
            lambda = krange(argmax);
        end
        
    case 2 % Ridge
        scores = zeros(length(krange),1);
        for fold = 1:num_folds
            test_mask = fold_ind==fold;
            X_train = X(~test_mask,:);
            X_test = X(test_mask,:);
            y_train = y(~test_mask);
            y_test = repmat(y(test_mask),1,length(krange));
            wts = ridge(y_train,X_train,krange);
            y_pred = X_test*wts;
            scores = scores+iter_corr(y_test',y_pred');
            [~,argmax]=max(scores);
            lambda = krange(argmax);
        end
        
    case 1.5 % Elastic Net
        k2range = [0.25 0.5 0.75];
        scores = zeros(length(krange)*length(k2range),1);
        for fold = 1:num_folds
            test_mask = fold_ind==fold;
            X_train = X(~test_mask,:);
            X_test = X(test_mask,:);
            y_train = y(~test_mask);
            y_test = repmat(y(test_mask),1,length(krange)*length(k2range));
            wts_mat = zeros(size(X,2),length(krange)*length(k2range));
            for ii = 1:length(k2range)
                wts_mat(:,(ii-1)*length(krange)+1:ii*length(krange)) = lasso(X_train,y_train,'lambda',krange,'alpha',k2range(ii));
            end
            y_pred = X_test*wts_mat;
            scores = scores+iter_corr(y_test',y_pred');
            [~,argmax]=max(scores);
            [argmax1, argmax2] = ind2sub([length(krange) length(k2range)],argmax);
            lambda(1) = krange(argmax1);
            lambda(2) = k2range(argmax2);
        end
        
        
    case 2.5 % Generalized ridge
        k2range = krange;
        scores = zeros(length(krange)*length(k2range),1);
        for fold = 1:num_folds
            test_mask = fold_ind==fold;
            X_train = X(~test_mask,:);
            X_test = X(test_mask,:);
            y_train = y(~test_mask);
            y_test = repmat(y(test_mask),1,length(krange)*length(k2range));
            wts_mat = zeros(size(X,2),length(krange)*length(k2range));
            kk = 1;
            for ii = 1:length(k2range)
                for jj = 1:length(krange)
                    wts_mat(:,kk) = tikhonov(X_train,y_train,[krange(ii) k2range(ii)]);
                    kk = kk+1;
                end
            end
            y_pred = X_test*wts_mat;
            scores = scores+iter_corr(y_test',y_pred');
            [~,argmax]=max(scores);
            [argmax1, argmax2] = ind2sub([length(krange) length(k2range)],argmax);
            lambda(1) = krange(argmax1);
            lambda(2) = k2range(argmax2);
        end
end
end

