function z = corrcoef_2(x,y)
% x = N1xD. y = N2xD.
% z = N1xN2 matrix of correlations of rows in x and y

z = zeros(size(x,1),size(y,1));
for ii = 1:size(x,1)
    for jj = 1:size(y,1)
        z(ii,jj) = fastcorr(x(ii,:),y(jj,:));
    end
end
