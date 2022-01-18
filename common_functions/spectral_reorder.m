function [DO i_sort]=spectral_reorder(D)

% reorders the distance matrix D, so that off diagonal 
% elements are on average as close as possible to the diagonal.
% Following the algorythm described in Johansen-Berg eet al, PNAS.
% D has to be symmetric matrix with positiv entries only and should
% be full rank.

if min(D(:)<0)
    error('D must have positive elements only');
end
% if any(D~=D')
%     error('D must be a symmetric matrix');
% end

sz=size(D);
Q=-D;
Q(find(eye(sz(1))))=0;
Q(find(eye(sz(1))))=sum(-Q,2);

T=zeros(size(Q));
T(find(eye(sz(1))))=1./sqrt(sum(D,2));

[V,E]=eig(T*Q*T);

[v,i_sort]=sort(V(:,2));

DO=D(:,i_sort);
DO=DO(i_sort,:);