function [Y] = vs_normalizer(X)
% Column wise range normalization.
Y = zeros(size(X));
for ii = 1:size(X,2)
    vec = X(:,ii);
    vec = (vec-min(vec))/(max(vec)-min(vec));
    vec = 2*vec-1;
% %     vec = vec/2+0.5;
%     vec = 100.^vec;
%     vec = log(2+vec);
    Y(:,ii)=vec;
end
end

