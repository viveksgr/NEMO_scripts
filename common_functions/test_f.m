function r_mat = test_f(p_mean_)
r_mat = arrayfun(@(x) r2t(x,160),p_mean_);

% For all folds
% r_mat = zeros(size(p_mean_));
% for ii = 1:size(p_mean_,1)
%     for jj = 1:size(p_mean_,2)
%         for kk = 1:size(p_mean_,3)
%             r_mat(ii,jj,kk) = r2t(p_mean_(ii,jj,kk),160);
%         end
%     end
% end

end