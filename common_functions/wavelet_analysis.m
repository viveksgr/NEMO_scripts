function [wv_ratings]=wavelet_analysis(ratings,num_level)

% For each column in ratings, apply wavelet transformation to give
% num_level x size(ratings,2) number of columns. 
% Columns are range normalized automatically. 


% Load behavior data
if nargin<2
    num_level = 3;
end

% Inputs to the wavelet function
ratings1 = repmat(ratings,1,num_level);

shift_template = -floor((num_level)/2):1:floor((num_level)/2);
shift_template = vs_normalizer((shift_template-median(shift_template))');
% shift_template = shift_template(2:end-1);

shift_row = repmat(shift_template',size(ratings,2),1); % Make sure that shift_template is column transposed to a row. 
shift_row = shift_row(:)';
shift_ratings = repmat(shift_row,size(ratings,1),1);
input_mat = ratings1-shift_ratings;

wv_ratings = arrayfun(@(x) wvlt(x),input_mat);
end

function [ricker] = wvlt(x,a)
if nargin <2
    a=0.3; % Variance in the function
end

ricker = (1-(x/a)^2)*exp(-(x/a)^2);
% ricker = exp(-(x/a)^2); % Test raw gaussian
end