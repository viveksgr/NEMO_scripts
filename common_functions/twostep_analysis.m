function [odor_pred_lr,training_lr] = twostep_analysis(behav_ratings,odor_responses,test_mask,mask,prog_bar)

% GLM predictions for odors given in test_mask using voxel_responses and
% perceptual ratings. N rows in all matrices should be same = number of
% odors. Toggle prog_bar to show progress bar - default yes.
% Training doesn't support lasso right now.
% Vivek Sagar, VivekSagar2016@u.northwestern.edu - February 6, 2019
% ------------------------------------------------------------------------
if nargin<5
    prog_bar = true; % Show progress bar.
end

% Training data
train_mask = ~test_mask;
behav_tr = behav_ratings(train_mask,:); % Training behavior
odor_resp_tr = odor_responses(:,:,:,train_mask);

% Test data
behav_tt = behav_ratings(test_mask,:);
odor_resp_tt = odor_responses(:,:,:,test_mask);

% GLM creation and prediction
dims = size(odor_resp_tr);
training_lr = zeros([dims(1:3),size(behav_tr,2)]);
odor_pred_lr = zeros(size(odor_resp_tt));

if prog_bar
    fprintf('GLM Progress:\n');
    fprintf(['\n' repmat('.',1,dims(1)) '\n\n']);
end

for ii = 1:dims(1)
    if prog_bar
        fprintf('\b|\n');
    end
    for jj = 1:dims(2)
        for kk = 1:dims(3)
            if mask(ii,jj,kk)~=0
                voxel_resp = squeeze(odor_resp_tr(ii,jj,kk,:));
                training_lr(ii,jj,kk,:) = (voxel_resp\behav_tr)'; % Model training
                odor_pred_lr(ii,jj,kk,:) = behav_tt*squeeze(training_lr(ii,jj,kk,:));                
            end
        end
    end
end
end


