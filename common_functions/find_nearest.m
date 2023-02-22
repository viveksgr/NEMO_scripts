function closest_id = find_nearest(target_, query_,order_sw)
% Find the ID of elements in target that are closest to elements in query. 
% Column vector target NX1, query = MX1, closest_id = MX1
% VivekSagar2016@u.northwestern.edu
% Order switch = true if target_ (cue time) < query_ (button press). Only works for sorted arrays.


if nargin<3
    order_sw = false;
end

closest_id = zeros(size(query_));
for ii = 1:length(query_)
    if ~order_sw
        [~,closest_id(ii)] = min(abs(target_-query_(ii)));
    else
        listsort = target_-query_(ii);
        closest_id(ii)=find(listsort<0,1,'last');
        if closest_id(ii)<0
            warning('Negative Index! Check the code')
        end
    end
end
