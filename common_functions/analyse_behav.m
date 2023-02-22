function [behav_data, detect_odorsID] = analyse_behav(behav,options,lesions)

% behav  = struct of behavior data
% behav.ratings = NxM behavior data: N = odors and M = perceptual or chemical features
% behav.cid = Nx1 array of odor pubchem IDs
% behav.percepts = Mx1 cell array of perceptual descriptors' descriptions
% behav.detect = Nx1 detection probability of odors
%
% Options include (in the order of operations):
% options.only_detect = (float) Use only odors that can be detected. Specify threshold. Warning = this option changes the size of rows of behav.ratings
% options.pca_anal = PCA(float, between 0 and 1) = cutoff needs to be specified
% options.polybase = (int) Polynomial order of features
% options.wavelet = (int) Convolve with morlet wavelets. Specify the number of basis functions
% options.perrun = (int) Repeat matrix for perrun analysis
% options.orthog = (any) Orthogonalize the array
% options.normalization = (any) Normalize the data between -1 and 1
% options.intercept = Intercept(any data type/logical) = append ones in the beginning
% options.binarize = binarize the vectors with threshold provided by
% options.binarize. Or use 'medsplit' for median split.
%
% lesions.pre.id = (int) Specify the ID of the descriptor to be removed before analysis
% lesions.pre.regress_out = (any data type/logical) Regress out the descriptor from the rest of the descriptors
% lesions.post = (int) Specify the ID of the descriptor to be removed after analysis
% behav_data = (NxM', N' = number of odors chosen, M' = new descriptors.) Ratings data post-behavior analysis
% detect_odorsID = (logical) odors used in analysis
% Vivek Sagar, December 6, 2018
% VivekSagar2016@u.northwestern.edu
% -------------------------------------------------------------------------

% If no descriptors need to be removed
if nargin<3
    lesions = struct();
end

% Choose only the odors that can be detected
if isfield(options,'only_detect')
    detect_odorsID = behav.detect>options.only_detect;
    behav.cid = behav.cid(detect_odorsID);
    behav.ratings = behav.ratings(detect_odorsID,:);
end

% Remove or regress out a descriptor before analysis
if isfield(lesions,'pre')
    if isfield(lesions.pre,'regress_out')
        if lesions.pre.regress_out
            behav.ratings(:,[1,lesions.pre.id])=behav.ratings(:,[lesions.pre.id,1]);
            behav.ratings = var_partioner_mat(behav.ratings);
            behav.ratings(:,1)=[];
        end
    else
        behav.ratings(:,lesions.pre.id)=[];
    end
end

% Do PCA on descriptors
if isfield(options,'pca_anal')
    [~, pca_r,~,~,var] = pca(behav.ratings);
    cum_var = cumsum(var);
    ind_pca = find(cum_var>100*options.pca_anal);
    behav.ratings = pca_r(:,1:ind_pca(1));
end

% Use wavelet kernels
if isfield(options,'wavelet')
    behav.ratings = wavelet_analysis(behav.ratings,options.wavelet);
end

% Expand behav.ratings for perrun analysis
if isfield(options,'perrun')
    n_odors = size(behav.ratings,1);
    ind2 = repmat(1:n_odors,options.perrun,1);
    behav.ratings = behav.ratings(ind2(:),:);
    behav.cid = behav.cid(ind2(:));
end

% Orthogonalization
if isfield(options,'orthog')
    if options.orthog
        behav.ratings = orth(behav.ratings);
    end
end

% Normalization
if isfield(options,'normalization')
    if options.normalization
%        behav.ratings = vs_normalizer(behav.ratings);
        behav.ratings = zscore(behav.ratings);
    end
end

% Regress out - orthogonalized
if isfield(options,'regress_orth') % Vectors must be z-scored
    fprintf('Running orthog, check if z-scored')
    behav.ratings(:,[1,options.regress_orth])=behav.ratings(:,[options.regress_orth,1]);
    behav_mat = zeros(size(behav.ratings,1),size(behav.ratings,2)-1);
    vec1 = behav.ratings(:,1);
    for ii = 2:size(behav.ratings,2)
        vec2 = behav.ratings(:,ii);
        behav_mat(:,ii-1)= vec2-((vec1'*vec2)./norm(vec1)^2)*vec1;
    end
    behav.ratings = [vec1 behav_mat];
    behav.ratings(:,[1,options.regress_orth])=behav.ratings(:,[options.regress_orth,1]);
end

% Append powers of descriptors
if isfield(options,'polybase')
    behav_ratings_powers = [];
    for poly = 1:options.polybase
        behav_ratings_powers = [behav_ratings_powers (behav.ratings).^poly];
    end
    behav.ratings = behav_ratings_powers;
end

% Fourier
if isfield(options, 'fourier')
    behav.ratings = pi*behav.ratings;
    ratings = zeros(size(behav.ratings,1),size(behav.ratings,2)*options.fourier*2);
    ll = 1;
    for ii = 1:size(behav.ratings,2)
        for kk = 1:options.fourier
            ratings(:,ll) = sin(kk*behav.ratings(:,ii));
            ratings(:,ll+1) = cos(kk*behav.ratings(:,ii));
            ll = ll+2;
        end
    end
    behav.ratings = ratings;
end

% Regress out - linear regress
if isfield(options,'regress')
    behav.ratings(:,[1,options.regress])=behav.ratings(:,[options.regress,1]);
    behav_mat = zeros(size(behav.ratings,1),size(behav.ratings,2)-1);
    vec1 = behav.ratings(:,1);
    for ii = 2:size(behav.ratings,2)
        vec2 = behav.ratings(:,ii);
        b = [ones(size(vec2)) vec2]\vec1;
        behav_mat(:,ii-1) = vec2-vec1*b(2)-b(1);
    end
    behav.ratings = [vec1 behav_mat];
end

% binarize the regressors or not
if isfield(options,'binarize')
    if strcmp(options.binarize,'medsplit')
        ratings_ = behav.ratings;
        meds = median(ratings_);
        ratings_bin = [];
        for ii = 1:length(meds)
            ratings_bin(:,ii) = double(ratings_(:,ii)>meds(ii));
        end
        behav.ratings=ratings_bin;
    else
        ratings_ = behav.ratings(:,options.binarize(2):end);
        ratings_bin = double(ratings_>options.binarize(1));
        behav.ratings(:,options.binarize(2):end)=ratings_bin;
    end
    
end

% Ternary expansion
if isfield(options,'ternarize')
    ratings_ = vs_normalizer(behav.ratings(:,5:end));
    ratings_bin = [];
    for ii = 1:size(ratings_,2)
        t1 = prctile(ratings_(:,ii),100/3);
        t2 = prctile(ratings_(:,ii),200/3);
        ratings_bin(:,ii) = double(ratings_(:,ii)>=t2)-double(ratings_(:,ii)<t1);
    end
    behav.ratings(:,3:18)=ratings_bin;
end

% Log scale
if isfield(options,'log_scale')
    ratings_ = vs_normalizer(behav.ratings(:,3:18));
    ratings_bin = (ratings_./abs(ratings_)).*log(abs(ratings_));
    behav.ratings(:,3:18)=ratings_bin;
end

% Lesion out descriptors post analysis
if isfield(lesions,'post')
    behav.ratings(:,lesions.post)=[];
end

% Append intercept
if isfield(options,'intercept')
    if options.intercept
        behav.ratings = [ones(size(behav.ratings,1),1) behav.ratings];
    end
end

% Remove descriptors of low values
if isfield(options,'rem_low')
    if options.rem_low
        nm = (sum(behav.ratings<=0,1)<160);
        behav.ratings = behav.ratings(:,nm);
        behav.percepts = behav.percepts(nm);
    end
end

behav_data = behav;
end
