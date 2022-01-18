% Create training set of odors
% C:\Data\NEMO\NEMO_01\behavior\dream_data_21descriptors

load('training_odors.mat')
cid = behav.cid;
lia = ismember(sempruned{:,'CID'},cid);
ratings_ = sempruned(lia,:);
odornames = ratings_.IdentifierOdor;

ratings_(:,1:3) = [];
columns_ = [1 2 6 11 14 5 9 21 13 3 8 12 18 19]; % Need to be mentioned manually
ratings_ = ratings_(:,columns_);
ratings_ = ratings_{:,:};
ratings_ = vs_normalizer(ratings_);

n_high = 5;
n_low = 5;
highs_ = zeros(n_high,length(columns_));
high_vals = zeros(n_high,length(columns_));
lows_ = zeros(n_low,length(columns_));
low_vals = zeros(n_low,length(columns_));

% Create chunks of descriptors
for ii = 1:length(columns_)
    [rating, argsort] = sort(ratings_(:,ii));
    cid_percept = cid(argsort);
    names_percept = odornames(argsort);
    
    high_vals(:,ii) = rating(end-n_high+1:end);
    highs_(:,ii) = cid_percept(end-n_high+1:end);
    low_vals(:,ii) = rating(1:n_low);
    lows_(:,ii) = cid_percept(1:n_low);
    
    struct_.CID = [lows_(:,ii); highs_(:,ii)];
    struct_.Odor = [names_percept(1:n_low); names_percept(end-n_high+1:end)];
    struct_.Valves = [(1:5)'; (8:12)'];
    struct_.Air = [7*ones(5,1); 6*ones(5,1)];
    ratings = [zeros(5,1); ones(5,1)];
    
    % Shuffle the rows
    if ii>7
        statpath = 'C:\Data\NEMO\NEMO_02\behavior\set_00\sess_02';
        kk = ii-8; % file_id
    else
        statpath = 'C:\Data\NEMO\NEMO_02\behavior\set_00\sess_01';    
        kk = ii-1; % file_id
    end
    argsort = randperm(10);
    ratings = ratings(argsort);
    filename = fullfile(statpath,sprintf('ratings%d.mat',kk));
    save(filename,'ratings')
    
    struct_.ratings = ratings;
    T = struct2table(struct_);    
    T{:,1} = T{argsort,1};
    T{:,2} = T{argsort,2};

    chunkname = fullfile(statpath,sprintf('chunk%d.csv',kk));
    writetable(T,chunkname,'Delimiter','\t')  
end

% 
% % Check non-zero sempruned.
% semtable = sempruned1{lia,:};
% semtable = vs_normalizer(semtable);
% semtable_sn = sign(semtable);
% semtable_sn(semtable_sn<1) = 0;
% num_odors = sum(semtable_sn);
% descr = sempruned1.Properties.VariableNames;