function [p_thresh, p_correct] = fdr_benjhoc(p_vals, q)
% Employs Benjamin Hochberg FDR correction to adjust p_vals with FDR<q. Default q = 0.05
% p_vals = Nx1, column vector of uncorrected p_values
% q = 1X1, FDR rate, 0.05 default
% p_thresh = 1X1, uncorrected p_value such that indices with
% p_vals<=p_thresh have null hypothesis rejected
% p_correct = NX1, adjusted p_values
%-------------------------------------------------------------------------
% Vivek Sagar, VivekSagar2016@u.northwestern.edu, May 2019
% Kahnt Lab, Northwestern University
% Adapted from Author: David M. Groppe

if nargin<2
    q = 0.05; % Default q_value in FDR
end

if size(p_vals,2)>1
    p_vals = p_vals(:); % Reshape into a column vector
end

m = length(p_vals); % Number of tests
% Compute p_thresh --------------------------------------------------------
[sorted_p, argsort] = sort(p_vals);
% Diagonal line for thresh
k_line = ((1:m)*q/m)';
p_thresh_ind = find(sorted_p<=k_line,1,'last');
p_thresh = sorted_p(p_thresh_ind);

% Adjust p_values----------------------------------------------------------
if nargout>1
    [~, unsort_ids] = sort(argsort);
    raw_adjusted =m*sorted_p./(1:m)';
    p_correct = zeros(size(p_vals));
    [wtd_p_sorted, wtd_p_sindex] = sort(raw_adjusted);
    nextfill = 1;
    for k = 1:m
        if wtd_p_sindex(k)>=nextfill
            p_correct(nextfill:wtd_p_sindex(k))=wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill>m
                break;
            end
        end
    end
    p_correct=p_correct(unsort_ids);
end
