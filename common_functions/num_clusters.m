function num_clusters(X,max_num)
% Finds the number of k-means clusters using elbow method. Find the elbow
% in the plot generated.

warning('off','all')

if nargin<2
    max_num = 500;
end

cent_score = zeros(1,max_num); % Minimum num of clusters is 2.
h = waitbar(0,'Initializing waitbar...');
for ii = 1:max_num
    [~,~,sumd] = kmeans(X,ii);
    cent_score(ii) = sum(sumd(:));
    perc = ii*100/max_num;
    waitbar(perc/100,h,sprintf('%.2f%% Completed',perc))
end
close(h)