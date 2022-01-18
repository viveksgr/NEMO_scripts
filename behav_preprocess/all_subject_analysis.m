%% Preprocess behavior data

% Number of sets, sessions etc.
sn = 1;
set_i = 1;
set_f = 4;
nsets = set_f-set_i+1;
sess_ = [2 4];
nsess = length(sess_);
n_odors_sess = 40;
% n_trials12 = 30;
% n_trials34 = 27;

root = 'C:\Data\NEMO';
% Each cell consists of results from different sets.
% set_detect = cell(1,nsets);
set_detect = [];
set_rt_detect = cell(1,nsets);
set_ratings = cell(nsets,nsess);
set_cid = cell(nsets,nsess);

for set_ = set_i:set_f
    setpath = fullfile(root, sprintf('NEMO_%02d',sn),...
        'behavior','behavior_task',sprintf('set_%02d',set_));
    % Initialize arrays to store sess data for this set.
    data_detect = repmat({cell(n_odors_sess,1)},1,nsess);
%   data_detect = repmat({zeros(n_odors_sess,1)},1,nsess); % For
%   detectability
    data_rt_detect = repmat({zeros(n_odors_sess,1)},1,nsess);
    data_ratings = repmat({zeros(n_odors_sess,18)+NaN},1,nsess);
    for sess = sess_
        sesspath = fullfile(setpath,sprintf('sess_%02d',sess));
        run_files = dir(fullfile(sesspath, sprintf('NEMO*.mat')));
        data_run = cell(1,length(run_files));
        odornames = {};
        run_l = length(run_files);
        for run_ = 1:run_l
            runpath = fullfile(sesspath,run_files(run_).name);
            data_run_raw = load(runpath);
            
            if isfield(data_run_raw.res,'sp_mat')
                s = rmfield(data_run_raw.res,'sp_mat');
                data_run_raw.res = s; % No idea why this doesn't run inplace.
            end
            
            % Extract run data
            data_run{run_} = data_run_raw.res;
            try
                odorname = vertcat(data_run_raw.res.odornames(1:end-1,2));
                odornames = [odornames; odorname];
            catch
                odorname = {};
            end
        end
        % Collapse run data into one array for each session.
        data_sess_raw = cat(1,data_run{:});
        
        % Construct a dictionary to map each odor trial to correct row.
        try
        cid_keys = vertcat(data_sess_raw(:).odor_id);
        set_ratings{set_,sess}=vertcat(data_sess_raw(:).p_mat);    
        catch
            cid_keys = [];
            set_ratings{set_,sess} = [];
        end       
    % Trials from same odor are collated in the same row. 
    % Data from entire set in one cell
%     if ismember(set_,[1,2])
%         set_detect{set_} = nansum(horzcat(data_detect{:}),2)/n_trials12;
%         set_rt_detect{set_} = nansum(horzcat(data_rt_detect{:}),2)/n_trials12;
%     else
%         set_detect{set_} = nansum(horzcat(data_detect{:}),2)/n_trials34;
%         set_rt_detect{set_} = nansum(horzcat(data_rt_detect{:}),2)/n_trials34;
%     end
    set_cid{set_,sess} = cid_keys;
    names{set_,sess} = odornames;
    end
end
% clearvars -except set_ratings
% load('which_desc.mat')

%% Test retest reliability for sets
% set_ = []; 
% Mistake in database for sessions 1 and 3 for NEMO02, NEMO03, NEMO04
% 7122-7335
% 1032-1001
% 2214-6501
% 8174-7937
% 8697-556940
% 14525-11569

set_cid{3,1}(set_cid{3,1}==1032)=1001;
set_cid{3,1}(set_cid{3,1}==7122)=7335;
set_cid{4,1}(set_cid{4,1}==2214)=6501;
set_cid{4,1}(set_cid{4,1}==8174)=7937;
set_cid{4,1}(set_cid{4,1}==8697)=556940;
set_cid{4,1}(set_cid{4,1}==14525)=11569;

sess1_ = sess_(1);
sess2_ = sess_(2);
% toggle_replacement=true;
rat1 = vertcat(set_ratings{:,sess1_});
rat3 = vertcat(set_ratings{:,sess2_});
cid1 = vertcat(set_cid{:,sess1_});
cid3 = vertcat(set_cid{:,sess2_});
[cid,argsort1] = sort(cid1);
rat1 = (rat1(argsort1,:));
[~,argsort3] = sort(cid3);
rat3 = (rat3(argsort3,:));
assert(isequal(sort(cid1),sort(cid3)))

rat1(isnan(rat1))=0;
rat3(isnan(rat3))=0;
% if toggle_replacement
% end

% test retest
for ii = 1:length(cid)
    nan_ = isnan(rat3(ii,:));
    corr__ = corrcoef(rat1(ii,~nan_),rat3(ii,~nan_)); 
    if length(corr__)>2
        corr_(ii) = corr__(2);
    end
end

% test retest
for ii = 1:size(rat1,2)
    nan_ = isnan(rat3(:,ii));
    corr__ = corrcoef(rat1(~nan_,ii),rat3(~nan_,ii)); 
    corr_2(ii) = corr__(2); 
end

bar(1:length(percepts),corr_2)
xticks(1:18)
xticklabels(percepts)
xtickangle(90)
hold on
plot([0 xl(2)],[r2t(0.05,160) r2t(0.05,160)],'r')

%% Make figures from all subjects 
behav1_labels = behav.percepts;
behav2_labels = behav_2.percepts;

% Manually change the indices, so they match across subjects.
behav_new = behav2_labels;

% quartiles
med1 = median(behav.ratings);
med2 = median(behav_2.ratings);
q25_1 = prctile(behav.ratings,25);
q25_2 = prctile(behav_2.ratings,25);
q75_1 = prctile(behav.ratings,75);
q75_2 = prctile(behav_2.ratings,75);

% Polar
theta = deg2rad(linspace(0,360,22));
figure()
polarplot(theta,[med1 med1(1)],'b')
set(gca,'Thetatick',rad2deg(theta))
set(gca,'ThetaTickLabel',behav_new)
set(gca,'Rtick',[])
rlim([0 1])
hold on
polarplot(theta,[med2 med2(1)],'r')
% polarplot(theta,[med3 med3(1)],'b.')
% polarplot(theta,[q25_2 q25_2(1)],'b.')
% polarplot(theta,[q75_2 q75_2(1)],'b.')
% polarplot(theta+deg2rad(5),[med1 med1(1)],'g.')
% polarplot(theta+deg2rad(5),[q25_1 q25_1(1)],'g.')
% polarplot(theta+deg2rad(5),[q75_1 q75_1(1)],'g.')
% title('r 01, b 02')
savefig(fullfile(pwd,'polarscorr_s3.fig'))

%% RSA correlations:
utl = logical(triu(ones(size(Behav_RSMC)),1));
Behav_RSM = corrcoef(behav.ratings');
bdist = (1-Behav_RSM)/2;
mat = spectral_reorder(bdist);
imagesc((1 - 2*mat))
colorbar
colormap('jet')

rsac = Behav_RSMC(utl);
rsap = Behav_RSM(utl);
fastcorr(rsac,rsap)

densityplot(rsac,rsap)
hold on
X = [ones(size(rsap)) rsac];
b = X\rsap;
y = X*b;
plot(rsac,y,'k')
title(sprintf('corr%f',fastcorr(rsac,rsap)))
save('chemcorr')
save('chemcorr.mat')

%% Behavioral Bar plots
% Load NEMO_perceptual from NEMO_all
nS = 3;
behav_lab1 = behav(1).percepts;
behav_lab2 = behav(2).percepts;
for ii=1:nS
    behav(ii).rel(end+1)=0;
end

idx2 = [1:10 19 11:14 19 15 19 16:18]; % Use this is argsort in sub(2 and 3) to match the labels
idx1 = [1:18 19 19 19];
behav_labs = {behav_lab1{:} behav_lab2{end-2:end}};

% Reordering.
% ord = [8 18; 10 19; 5 17; 12 11; 9 5; 14 2];
% ord2 = [11 16; 19 13; 14 17; 9 8];
% behav_labs = mat_shuffler(behav_labs,ord);
% behav_labs = mat_shuffler(behav_labs,ord2);
rels = zeros(nS,length(behav_labs));
for ii = 1:nS
    if ii==1
        rels(ii,:)=behav(ii).rel(idx1);
    else
        rels(ii,:)=behav(ii).rel(idx2);
    end
%     rels(ii,:) = mat_shuffler(rels(ii,:),ord);
%     rels(ii,:) = mat_shuffler(rels(ii,:),ord2);
end

% Theta plots
theta = deg2rad(linspace(0,360,22));
cs = {'r.','g.','b.'};
figure('Position',[0.5 0.5 400 300])
polarplot(theta,[rels(1,:) rels(1,1)],cs{1})
hold on
for ii = 2:nS
    polarplot(theta,[rels(ii,:) rels(ii,1)],cs{ii},'MarkerSize',5)
end
set(gca,'Thetatick',rad2deg(theta))
set(gca,'ThetaTickLabel',behav_labs)
set(gca,'Rtick',[0.13 0.5 1])
rlim([0 0.7])
savefig(fullfile(pwd,'polarscorr.fig'))

% Bar plots
rels2 = rels;
rels2(rels2==0)=nan;
rels2_m = tanh(nanmean(atanh(rels2)));
% [~,argsort] = sort(rels2_m,'descend');
% Remove bars with 0 ratings manually
bar(1:length(behav_labs),rels2_m)
hold on
plot(1:length(behav_labs),rels2,'.','MarkerSize',5)
yline(0.13)
xticks(1:length(behav_labs))
xtickangle(90)
xticklabels(behav_labs)

% Plotting odor profiles of ethyl valerate and cyclopentanethiol ---------- 
med1 = behav(2).ratings(82,:);
med2 = behav(2).ratings(106,:);
med1(end+1) = 0;
med2(end+1) = 0;
med1 = med1(idx2);
med2 = med2(idx2);
med1 = mat_shuffler(med1,ord);
med2 = mat_shuffler(med2,ord);
med1 = mat_shuffler(med1,ord2);
med2 = mat_shuffler(med2,ord2);
behav_new = behav_labs;
theta = deg2rad(linspace(0,360,22));
% figure()
polarplot(theta,[med1 med1(1)],'b')
set(gca,'Thetatick',rad2deg(theta))
set(gca,'ThetaTickLabel',behav_new)
% set(gca,'Rtick',[])
rlim([0 1])
hold on
polarplot(theta,[med2 med2(1)],'r')
savefig(fullfile(pwd,'polarscorr_full.fig'))
% % Histograms
% xedges = -0.2:0.1:1;
% figure('Position',[0.5 0.5 400 300])
% histogram(NEMO02,xedges,'FaceAlpha',0.3,'FaceColor','g')
% hold on
% histogram(NEMO04,xedges,'FaceAlpha',0.3,'FaceColor','b')

% % Box plots of ratings
% figure('Position',[0.5 0.5 300 200])
% ii = 3;
% boxplot(behav(ii).ratings)
% xticks(1:18)
% xticklabels(behav(ii).percepts)
% xtickangle(90)

%% Rating discriminability
load('NEMO_perceptual.mat')
nS = 3;
bwidth = 0.25;
edges = 0:bwidth:7;
edgec = edges(1:end-1)+(bwidth/2);
mz = zeros(length(edges)-1,nS);
for zz = 1:nS
    rat_temp = behav(zz).ratings;
    rat_temp = zscore(rat_temp);
    
    minz = zeros(length(rat_temp));
    for ii = 1:length(rat_temp)
        for jj = ii:length(rat_temp)
            des_diff = abs(rat_temp(ii,:)-rat_temp(jj,:));
            minz(ii,jj) = max(des_diff);
        end
    end
    utl = logical(triu(ones(size(minz)),1));
    minz_vec = minz(utl);
    mz(:,zz) = histcounts(minz_vec,edges);   
end
mz2 = mean(mz,2);
bar(edgec,mz2,'hist')

minz(minz==0)=nan;
histogram(minz(:))
minz(isnan(minz))=[];
sum(minz>2)/length(minz);

%% RSMs for 1 subject

% Make colormap
hex = ['#0D646B';'#798279';'#C3D1C3';'#E17812';'#FAE0AD'];
map = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;

% Load RSMs from RSA analysis
mat = NEMO01_chemical.ratings;
mat2 = Behav_RSM_C;

nidx = 1:10;
k_sumd = zeros(1,length(nidx));
for ii = nidx
    [~,~,sumd]=kmeans(mat,ii,'Distance','correlation');
    k_sumd(ii) = mean(sumd);
end

[kidx] = kmeans(mat,4,'distance','correlation');
mat2 = corrcoef(mat');
[~,argsort] =sort(kidx);
figure('Position',[0.5 0.5 320 320])
imagesc(mat2(argsort,argsort))
colormap(swampsunset)
% colorbar
xticks([1 80 160])
yticks([1 80 160])
caxis([-1 1])
savefig('chems_mat')

%% Identity vs category methods
% Load behavlabs. Remove Intensity and pleasantness.
behav(1).ratings = vs_normalizer(behav(1).ratings);
behav_lab1 = behav(1).percepts;
behav_lab2 = behav(2).percepts;
behav_labs = {behav_lab1{3:end} behav_lab2{end-2:end}};

B_dist = 1-pdist(behav(1).ratings,@maxcorrdist);
B_dist_mat = squareform(B_dist);
[m_sort, argsort] = sort(B_dist_mat(:),'descend');
[rr, cc] = ind2sub(size(B_dist_mat),argsort(1));

% Reduce odor ratings to include relevant descriptors
% idx2 = [1:10 19 11:14 19 15 19 16:18];
idx2 = [1:18 19 19 19];
temp = [ones(18,1); 0];
mask1d = temp(idx2);
o1 = rr;
o2 = cc;
rat1 = behav(1).ratings(o1,:);
rat1 = unmasker(rat1,mask1d);
rat1 = rat1(3:end);
rat2 = behav(1).ratings(o2,:);
rat2 = unmasker(rat2,mask1d);
rat2 = rat2(3:end);
rat3 = rat1.*rat2;

% Reorder
% ord = [7 16; 16 8; 16 13; 16 14];
% behav_labs = mat_shuffler(behav_labs,ord);
% rat1= mat_shuffler(rat1,ord);
% rat2= mat_shuffler(rat2,ord);
% rat3= mat_shuffler(rat3,ord);
% ord = [9 8];
% ord = [1 2; 4 5; 11 13; 12 13];
ord = [6 7];
behav_labs = mat_shuffler(behav_labs,ord);
rat1= mat_shuffler(rat1,ord);
rat2= mat_shuffler(rat2,ord);
rat3= mat_shuffler(rat3,ord);

% Theta plots
theta = deg2rad(linspace(0,360,20));
delta = deg2rad(1.5);
cs = {'r','g','b'};
figure('Position',[0.5 0.5 400 300])
polarplot(theta-delta,[rat1; rat1(1)],cs{1})
hold on
polarplot(theta+delta,[rat2; rat2(1)],cs{2})
polarplot(theta,[rat3; rat3(1)],cs{3})
set(gca,'Thetatick',rad2deg(theta))
set(gca,'ThetaTickLabel',behav_labs)
set(gca,'Rtick',[])
% rlim([0 0.7])
savefig(fullfile(pwd,'polarscorr.fig'))

%% Perceptual vs chemical collinearity
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\RSA_FIR\RSA_final';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\RSA_FIR\RSA_final';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\RSA_FIR\RSA_final'};
s_names = {'S1','S2','S3'};
behavP = variable_extract(dirs,'ROI.mat','Behav_RSM_vals_P',true);
behavC = variable_extract(dirs,'ROI.mat','Behav_RSM_vals_C',true);
nboot = 10000;

% mbar = zeros(1,3);
% mstd = zeros(1,3);
% mstd2 = zeros(1,3);
% bsam = bootstrp(nboot,@(x) x, 1:length(behavP));
% for ii = 1:length(s_names)
%     temp = zeros(1,nboot);
%     for ss = 1:nboot
%         pvec = behavP((bsam(ss,:)),ii);
%         cvec = behavC((bsam(ss,:)),ii);
%         temp(ss) = fastcorr(pvec,cvec);
%     end
%     mbar(ii) = mean(temp);
%     mstd(ii) = prctile(temp,97.5);
%     mstd2(ii) = prctile(temp,2.5);
% end

mdist = zeros(1,nboot);
bsam = bootstrp(nboot,@(x) x, 1:length(behavP));
temp = zeros(length(dirs),nboot);
for ii = 1:length(dirs)
    for ss = 1:nboot
        pvec = behavP((bsam(ss,:)),ii);
        cvec = behavC((bsam(ss,:)),ii);
        temp(ii,ss) = fastcorr(pvec,cvec);
    end
end
temp = mean(temp);
m_temp = mean(temp);
s_temp1 = prctile(temp,97.5);
s_temp2 = prctile(temp,2.5);
p_val = (100-invprctile(temp-m_temp,t_sq))/100; 
% figure('Position',[0 0 800 200])
% bar(1:3,mbar)
% hold on
% errorbar(1:3,mbar,mstd-mbar,mstd2-mbar,'.')
% savefig('pvsc_collinearity')
figure('Position',[0 0 800 200])
bar(1,mean(mbar))
hold on
errorbar(1,mean(mbar),std(mbar)/sqrt(3),'.')
plot([1 1 1],mbar,'r.','MarkerSize',15)
savefig('pvsc_collinearity')
    