function odor_evoked = BRS_analysis(voxel_act,DM_mat,n_odors, basis_fourier)

% BRS analysis extracts odor-specific responses from a voxel time series using cochrane-orcutt algorithm.
% voxel_act = Tx1 array of voxel time series (concatenated across sessions)
% DM_mat = SPM.mat with the design matrix of shape Tx((FxM)+N) where F is the number of fourier basis functions,
% M  = number of unique odors and N = the number of nuisance regressors. 
% n_odors = M = number of unique odors.
% n_iter = number of iterations in the cochrane-orcutt algorithm.

% -----------------------------------------------------------------------
% Vivek Sagar, Kahnt Lab

if nargin<4
    basis_fourier = true;
end

n_iter = 5; % Default value of n_iter

DM = sparse(DM_mat.SPM.xX.X); % Sparse design matrix to save memory

if basis_fourier
    num_fourier = 2*DM_mat.SPM.xBF.order+1; % Number of fourier basis functions
else
    num_fourier = DM_mat.SPM.xBF.order; % Number of gamma bases functions
end

h = ones(1,n_odors); 
c = zeros(1,n_odors*num_fourier); 


for iter = 1:n_iter
    extend_h_1 = repmat(h,[num_fourier 1]); 
    extend_h_2 = extend_h_1(:); 
    extend_h = [extend_h_2' ones(1,(size(DM,2)-length(extend_h_2)))]; % Pad with 1s to match h's dimension with DM
    
    DM = DM.*repmat(extend_h,[length(DM),1]); % DM multiplied with h
    
    master_b_1 = mldivide(DM,voxel_act); % c and b obtained through least square fit, b = nuisance regressor weights
    
    c = master_b_1(1:length(c))'; 
    b = master_b_1(length(c)+1:end)'; 
    
    voxel_act_pred = DM*master_b_1; % Estimate r for cochrane-orcutt
    residue = voxel_act - voxel_act_pred; 
    r = circshift(residue,1)\residue;
    
    voxel_act = voxel_act-r*circshift(voxel_act,1); %AR(1) correction
    DM = DM-r.*circshift(DM,1); 
    
    % DM needs to be transformed for the second part of the algorithm. Add
    % all fourier components in the design matrix:
    
    DM_t_b = eye(length(b)); % Leave nuisance regressors as it is. 
    c_mat = reshape(c,num_fourier,[]); % DM transformer for odor evoked act    
    DM_t_c = [];      
    for ii = 1:size(c_mat,2)
        DM_t_c = blkdiag(DM_t_c,c_mat(:,ii));
    end  
    
    DM_transformer = blkdiag(DM_t_c,DM_t_b);
    DM_2_iter = sparse(DM*DM_transformer); % Transformed DM
    
    master_b_2 = DM_2_iter\voxel_act; % All h's and b's
    
    h = master_b_2(1:length(h)); % h updated
    
    voxel_act_pred = DM_2_iter*master_b_2; % Estimate r for cochrane-orcutt
    residue = voxel_act - voxel_act_pred;
    r = circshift(residue,1)\residue; 
    
    voxel_act = voxel_act-r*circshift(voxel_act,1); %AR(1) correction 
    DM = DM-r.*circshift(DM,1); 
end

odor_evoked = zscore(h);
end

