function D2 = maxcorrdist(ZI,ZJ)
% Custom distance (not similiarity) between 1xN vector ZI and mXN vector ZJ. D2 = mX1
% Find the coordinate with minimum distance and take the product of that
% index from two vectors
% Make each column zero-mean
ZI = bsxfun(@minus,ZI,mean(ZI,2));
ZJ = bsxfun(@minus,ZJ,mean(ZJ,2));

% % L2 normalize each column
% if min(min(ZI(:)),min(ZJ(:)))<0
%     error('Feed non negative vectors');
% end

ZI = bsxfun(@times,ZI,1./sqrt(sum(ZI.^2,2)));
ZJ = bsxfun(@times,ZJ,1./sqrt(sum(ZJ.^2,2)));

M = (ZI.*ZJ);
colind = ZI>0;
colind2 = logical(sum(ZJ>0,1));
% M(:,~and(colind,colind2)) = 0;
% Take the dot product of the columns and then sum
D2=1-max(M,[],2);
