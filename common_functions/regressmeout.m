function [outputM] = regressmeout(M1,M2)
% Regress matrix of row vectors of M2 from M1

outputM = zeros(size(M1));
for vidx = 1:size(M1,1)
%   vec2 = zscore(M1(vidx,:)');
%   vec1 = zscore(M2(vidx,:)');
% %   outputM(vidx,:) = (vec2-((vec1'*vec2)./norm(vec1)^2)*vec1)*std(M1(vidx,:)')+mean(M1(vidx,:)');
% 
%   temp = var_partitioner(vec2,vec1);
%   outputM(vidx,:) = temp;

    vec2 = (M1(vidx,:)');
    vec1 = (M2(vidx,:)');
    bs = [ones(size(vec1)) vec1]\vec2;
    outputM(vidx,:) = vec2-[ones(size(vec1)) vec1]*bs;
end