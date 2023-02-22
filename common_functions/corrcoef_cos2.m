function z = corrcoef_cos2(x,y)
% x = NxD. y = NxD.
% z = NxN
% Only upper triangle provided

for ii = 1:size(x,1)
    for jj = 1:size(x,1)
        z(ii,jj) = x(ii,:)*y(jj,:)';
    end
end




