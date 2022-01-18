% Make nuisance regressors and align breathing and behavioral data
%% Basic Settings
linux_config = false; % Run on windows or quest
if linux_config
    root = '/projects/p30489/Data/NEMO'; % Quest Version
    addpath('/home/vsh3681/spm12')
else
    root = 'C:\Data\NEMO'; % Windows Version
end
subjects = [4];
disp = false; % Display plots
max_badvol = 22; % How many bad volumes should be eliminated. 
slice_diff_ = true; % Analyze slice differences
spiro_sw = true; 
% Add data from spirometer as well as breathing belt. By default breathing belt is on channel 4.
% If spirometer is connected on line                    

set_1 = 1; % Initial set
nsets = 4; % Final set
sess_1 = 2; % Initial session
nsess = 4; % Final session
nruns = {[1:4] [1:4] [1:4] [1:4]}; % Choose runs in each set
repTime = 1.4; % TR

s = subjects;
sn_ = sprintf('NEMO_%02d', s); 
sn = sprintf('NEMO%02d', s);

respath = fullfile(root,   sn_); % Subject path

NR = [];  %NR Nuisanace Regressors
nTR = [];
%% Get motion data
% Loads realignment parameters and plots translations and rotations.
for set_ =set_1:nsets
    for sess =sess_1:nsess
        
        for r=nruns{set_}
            datapath = fullfile(respath, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d', r));
            filename = fullfile(datapath, 'rp_*.txt');
            n = dir(filename);                          %rp file contains 6 columns represnting positions
            mp = load(fullfile(datapath, n(1).name));  % mp motion parameters
            nTR{set_, sess,r}  = size(mp,1);
            
            % add diff, squared mp, and squared diff
            NR{set_, sess,r} = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];
            
            if disp
                figure((s*1000)+(set_*100)+(sess*10)+(r));
                x= 1:nTR{set_, sess,r};
                subplot(4,2,1); plot(x, NR{set_,sess,r}(:,1:3)); xlabel('scans'); ylabel('mm'); subplot(4,2,2); plot(x,NR{set_,sess,r}(:,4:6)); xlabel('scans'); ylabel('deg');
                subplot(4,2,3); plot(x, NR{set_,sess,r}(:,7:9)); xlabel('scans'); ylabel('diff(mm)'); subplot(4,2,4); plot(x,NR{set_,sess,r}(:,10:12)); xlabel('scans'); ylabel('diff(deg)');
                subplot(4,2,5); plot(x, NR{set_,sess,r}(:,13:15)); xlabel('scans'); ylabel('mm^2');subplot(4,2,6);  plot(x,NR{set_,sess,r}(:,16:18)); xlabel('scans'); ylabel('deg^2');
                subplot(4,2,7); plot(x, NR{set_,sess,r}(:,19:21)); xlabel('scans'); ylabel('diff(mm)^2'); subplot(4,2,8);  plot(x,NR{set_,sess,r}(:,22:24)); xlabel('scans'); ylabel('diff(deg)^2');
            end
        end
    end
end

%% Analyze physiological data
adjust_mat = [];
for set_ = set_1:nsets
    for sess =sess_1:nsess
        for r=nruns{set_}
            
            filename = fullfile(respath, 'breathing','imaging_task',  sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('%s_set%02d_sess%02d_run%02d.mat', sn, set_, sess, r));
            load(filename)        
            % grab the data
            ben = data(datastart(1):dataend(1));     % when the first blip appeared. The first variable in labchart
            trigger = data(datastart(5):dataend(5));  % The 5th variable in labchart. The one with multiple blips for triggers
            
            % determine first trigger
            tr = find(round(trigger./max(trigger))==1);
            tr = tr(1);
            
            % determine t0 adjustment
            start = find(round(ben*10)<=42); % Adjust this threshold manually
            start = start(1);
            adjust = (start-tr)/samplerate(1);
            adjust_mat = [adjust_mat adjust];
            % sprintf('adjustment for subject %s, set %01d,  session %01d, run %01d: %1.4f sec', sn, set_, sess, r, adjust)
            save(fullfile(respath, 'breathing', 'imaging_task', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('time_adjust_%s_set_%02d_sess_%02d_run_%02d.mat',sn, set_, sess, r)), 'adjust');
            
%            % grab breathing data
            R = data(datastart(4):dataend(4));   % respiration variable from labchart,4th variable in labchart
%             
%           % take task data
            end_time = tr+ floor(nTR{set_,sess,r}*repTime*samplerate(1));
            if length(R)>=end_time
                R = R(tr:end_time);
            else
                R_y = length(R)+1:end_time;
                R_val = R(end)*ones(1,length(R_y));
                R = [R(tr:end) R_val];
            end
%             
%           % smooth
            R = smoothdata(R,'movmean',250)';
            
%             % high pass filter
            K.RT = 1/samplerate(1);
            K.row = ones(length(R),1);
            K.HParam = 50;
            R = spm_filter(K,R);
            R = zscore(R);
            if spiro_sw
                R = cumsum(R);
            end
            
            % down sample to scan resolution
            ssR = R(700:1400:end-700);
% %             
%             if disp
%                 figure((s*100)+(set_*10)+(sess));
%                 subplot(2,2,r)
%                 plot(ssR); xlabel('scans'), ylabel('breathing')
%             end
            
%             % add breathing trace
            NR{set_,sess,r} = [ssR, NR{set_,sess,r}];
            
            if spiro_sw
                if length(datastart)>5
                    spiro_ = data(datastart(6):dataend(6));
                    if length(spiro_)>=end_time
                        spiro_ = spiro_(tr:end_time);
                    else
                        spiro_ = [spiro_(tr:end) R_val];
                    end
                    
                    spiro_smooth = smoothdata(spiro_,'movmean',250)';
                    dr = spm_filter(K,spiro_smooth);
                    dr = zscore(dr);                    
                else
                    dr = [0; diff(R)];
                    dr = zscore(dr);                    
                end
                sdr = dr(700:1400:end-700);
                NR{set_,sess,r} = [sdr.^2, sdr, cumsum(sdr), NR{set_,sess,r}];
            end
            
        end
    end
    
end

%% Analyze slice difference
if slice_diff_
    subdir = fullfile(respath, 'imaging', 'nii');
    voldir = fullfile(subdir,  sprintf('set_%02d', set_1), sprintf('sess_%02d',sess_1), sprintf('run_%02d', 1) );
    filenames = dir(fullfile(voldir,'f*.nii'));
    v = spm_read_vols(spm_vol(fullfile(voldir,filenames(1).name)));
    nslices = size(v,3);
    idx1 = 1:2:nslices;   %odd slices
    idx2 = 2:2:nslices;   %even slices
    slicediff= [];
    slicevar =[];
    rr = 0;
    for set_ = set_1:nsets
        for sess =sess_1:nsess
            for r=nruns{set_}
                rr
                rr = rr +1;
                voldir = fullfile(subdir,  sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d', r) );
                filenames = dir(fullfile(voldir,'f*.nii'));
                
                for vol = 1:length(filenames)
                    v = spm_read_vols(spm_vol(fullfile(voldir,filenames(vol).name)));
                    odd_slice = nanmean(nanmean(nanmean(v(:,:,idx1),3)));
                    even_slice = nanmean(nanmean(nanmean(v(:,:,idx2),3)));
                    slicediff{rr}(vol) = (odd_slice - even_slice);
                    slicevar{rr}(vol) = var(nanmean(squeeze(nanmean(v))));
                end             
            end
        end
    end
    
    
    %% Get mean and std across all runs
    % start_l = []; for ii = set_1:length(nruns); start_l = [start_l length(nruns{ii})]; end % Length of runs per set
    % start_id = cumsum([1 start_l(1:end-1)]); % Start Id of run for each set
    %
    % tmp = []; for ii = set_1:nsets; tmp{ii} = horzcat(slicediff{start_id(ii-set_1+1):start_id(ii-set_1+1)+start_l(ii-set_1+1)-1}); end % Setwise summation
    % tmp1 = []; for ii = set_1:nsets; tmp1{ii} = horzcat(slicevar{start_id(ii-set_1+1):start_id(ii-set_1+1)+start_l(ii-set_1+1)-1}); end
    % tmp = horzcat(slicediff{:});
    % tmp1 = horzcat(slicevar{:});
    % aztmp = [max(abs((tmp-mean(tmp))./std(tmp))), max(abs((tmp1-mean(tmp1))./std(tmp1)))]; % calcuate the absolute of the z score
    bv_counter = 0;
    rcntr = 0;
    for set_ = set_1:nsets
        for sess =sess_1:nsess
            for r=nruns{set_}
                rcntr = rcntr+1;
                slicediff{rcntr} = abs((slicediff{rcntr}-mean(slicediff{rcntr}))./std(slicediff{rcntr}));
                slicevar{rcntr} = abs((slicevar{rcntr}-mean(slicevar{rcntr}))./std(slicevar{rcntr}));
                
                % set criterion for above mean
                crit = [5,5]; % criterion for slicediff and slicevar, eyeball it looking at the peaks
                
                if disp
                    figure(rcntr);
                    sh = subplot(2, 1, 1);
                    plot(slicediff{rcntr});
                    title(sprintf('slice diff, run %01d', rcntr));
                    hold on, plot(find(slicediff{rcntr}> crit(1)), slicediff{rcntr}(slicediff{rcntr}> crit(1)), '*r');
                    
                    sh = subplot(2, 1, 2);
                    plot(slicevar{rcntr});
                    title(sprintf('slice var, run %01d', rcntr));
                    hold on, plot(find(slicevar{rcntr}> crit(2)), slicevar{rcntr}(slicevar{rcntr}> crit(2)), '*r');
                end
                
                % include regressors for bad volumes
                badvol = unique([find(slicediff{rcntr}> crit(1)), find(slicevar{rcntr}> crit(2))]);
                ctrr1(s,rcntr) = length(badvol); % count bad volumes using crit
                bv = zeros(nTR{set_,sess,r} ,max_badvol);
                for i = 1:ctrr1(s,rcntr)
                    bv_counter = bv_counter+1;
                    bv(badvol(i),bv_counter) = 1;
                end
                
                % add slicediff and diff of slicediff and bad volume regressors
                NR{set_,sess,r} = [NR{set_,sess,r},...
                    slicediff{rcntr}', slicevar{rcntr}',...
                    [0; diff(slicediff{rcntr})'], [0; diff(slicevar{rcntr})'],...
                    [0; diff(slicediff{rcntr})'].^2, [0; diff(slicevar{rcntr})'].^2,...
                    bv];
            end
        end
    end
    max_badvol = sum(ctrr1(s,:));
end

%% Write Files
for set_ = set_1:nsets
    for sess =sess_1:nsess
        for r=nruns{set_}
            datapath = fullfile(respath, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d', r));
            % save nusiance regressors as txt
            fname = fullfile(datapath, sprintf('nuisance_regresssors_%s_set_%02d_sess_%02d_run_%02d.txt', sn, set_, sess, r));
            dlmwrite(fname, zscore(NR{set_, sess,r}), 'delimiter', ' ', 'precision', '%.8f');
        end
    end
end

% % Merge breathing traces if accidentally paused recording.
% if ~linux_config
%     % Break into two matrices
%     M1 = [];
%     M3 = [];
%     for ii = 1:6
%         M1 = [M1; data(datastart(ii,1):dataend(ii,1))];
%         M3 = [M3; data(datastart(ii,2):dataend(ii,2))];
%     end
%     % ------- Manually change breathing trace start point for M3
%     M3(4,1:100) = M3(4,100);
%     % Create appropriate filler
%     % First create time series - this needs to be done manually
%     T_series = (0:1:size(M1,2)-1)./samplerate(1);
%     spike_thresh = 1; % Manually choose this threshold
%     t_lastspike_s1 = find(M1(5,:)>1,1,'last'); % Last TR trigger in series1
%     t_firstspike_s2 = find(M3(5,:)>1,1); % First TR trigger in series2
%
%     % Manually check that 2TRs have elapsed during pause.
%     % Solve this: (length(M1)-t_lastspike_s1)+time_gap+t_firstspike_s2 =
%     % 2800 (2*TR --- choose this manually)
%     time_gap = 2800-(length(M1)-t_lastspike_s1)-t_firstspike_s2;
%     T_series2 = (1:1:time_gap)/1000+T_series(end); % Timeseries across gap
%     T_series3 = (1:1:length(M3))/1000+T_series2(end); % Timeseries of second recording
%
%     % Interpolated data
%     M2 = [];
%     for ii=[1:3 5]
%         M2(ii,:) = (M1(ii,end)+M3(ii,1))/2*(ones(1,length(T_series2)));
%     end
%
% %     % Try interpolation
% %     % Interpolated breathing data
% %     M2(4,:)= interp1([T_series T_series3],[M1(4,:) M3(4,:)],T_series2,'pchip');
%
%     % Take a standard breathing trial - check visually
%     breathing_clip = M3(4,1385:3952);
%     b_start = find_nearest(breathing_clip,M1(4,end));
%     b_end = find_nearest(breathing_clip(1001:end),M3(4,1))+1000;
%     breathing_ = breathing_clip(b_start:b_end);
%     M2(4,:)=interp1(linspace(1,100,length(breathing_)),breathing_,linspace(1,100,length(T_series2)));
%
%     % M_plot example
%     M_before = M1(4,end-10000:end);
%     M_after = M3(4,1:10000);
%     M = [M_before M2(4,:) M_after];
%     plot(M)
%     hold on
%     plot([10000 10000],[-0.03 0.06],'r')
%     plot([10000+length(M2) 10000+length(M2)],[-0.03 0.06],'r')
%
%     % Interp for spiro data
%     M2(6,:)= interp1([T_series T_series3],[M1(6,:) M3(6,:)],T_series2,'pchip');
%     M_before = M1(6,end-10000:end);
%     M_after = M3(6,1:10000);
%     M = [M_before M2(6,:) M_after];
%     plot(M)
%     hold on
%     plot([10000 10000],[-0.05 0.02],'r')
%     plot([10000+length(M2) 10000+length(M2)],[-0.05 0.02],'r')
%
%     % Join the 3 matrices
%     M = [M1 M2 M3];
%     M_sharp = M';
%     data = M_sharp(:)';
%     datastart = [];
%     dataend = [];
%
%     for ii = 1:6
%         datastart = [datastart; (ii-1)*length(M)+1];
%         dataend = [dataend; ii*length(M)]
%     end
% end
%
