function X2 = var_partioner_mat(X)
% Symmetrically orthog the first column from the rest in X
X2 = X;
for ii = 2:size(X,2)
%     if abs(fastcorr(X(:,1),X(:,ii)))>r2t(0.05,length(X(:,ii)))
        [~,X2(:,ii),t] = var_partitioner(X(:,1),X(:,ii));
%     end
end