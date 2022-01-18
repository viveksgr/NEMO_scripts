function pred_acc = decoder_mat(mat,shuff)
if nargin <2
    shuff = false;
end
pred_acc = zeros(1,size(mat,2));

if ~shuff
    for ii = 1:size(mat,2)
        pred_acc(ii) = mat(ii,ii)==max(mat(:,ii));
%         pred_acc(ii) = (pred_acc(ii)+ mat(ii,ii)==max(mat(ii,:)))/2;
        %     pred_acc(ii) = atanh(mat(ii,ii))>mean(atanh(mat([1:ii-1 ii+1:end],ii)));
    end
else
    ii_s = randperm(size(mat,2),1);
    for ii = 1:size(mat,2)
        pred_acc(ii) = mat(ii_s,ii)==max(mat(:,ii));
    end
end
