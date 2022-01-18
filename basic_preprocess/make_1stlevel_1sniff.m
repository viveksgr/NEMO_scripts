unction make_1stlevel_1sniff

%dbstop if error
modelname = '1sniffCue_set24_sniff';
root = 'C:\Data\NEMO\';
addpath('C:\spm12')
spm fmri
subjects = [4];
nsubs = length(subjects);

% More options
detection_con = true; % Contrast between detected and non-detected odors
parametric_intensity = false; % Parametric modulation with intensity

n_vol = 876;

set_i = 1;
nsets = 4;
nsess = 4;
sess_1 = 2;
nruns = {[1:4] [1:4] [1:4] [1:4]};

for s = subjects
    sn = sprintf('NEMO_%02d', s);
    matlabbatch = [];    
    statpath = fullfile(root, sn, 'imaging', '1stlevelmodels',modelname); %Sessname folder, from sleepdepr matrix go to first and second dim
    mkdir(statpath)
    % Added the script for the model below this spm script and the function is called here.
    % input variables in bracket
    %% 3odor_rating MODEL
    matlabbatch{1}.spm.stats.fmri_spec.dir = {statpath};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.4;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% Selecting rf images for model analysis
    filename = [];
    r=0;
    for set_ =set_i:nsets
        for sess = sess_1:nsess
            for rr =nruns{set_}
                make_conditions(s, root, statpath, modelname, set_, sess, rr)  % make_conditions is the funtion file we created for the stats model script.
                
                r=r+1;                                                          % add counter to move from one run to another run within the for loop
                datapath = fullfile(root,sn,'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));             
                n = dir(fullfile(datapath,sprintf('srf*.nii')));
               
                for i= 1:length(n)
                    fname = fullfile(datapath, n(i).name);
                    filename{i,r}= sprintf('%s,1',fname);
                end
                
                nc = dir(fullfile(statpath,sprintf('conditions_%s_%s_set_%02d_sess_%02d_run_%02d.mat',modelname,sn,set_,sess,rr)));
                condfile = fullfile(statpath,nc(1).name);
                loader = load(condfile);
                onset_set(:,r)=loader.onsets{1};
                
                detect_set(:,r) = loader.detect;
                cid_set(:,r) = loader.cid;
                
                
                n = dir(fullfile(datapath,sprintf('nuisance_regresssors_NEMO%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess,rr)));
                nuisanceregfile1 = fullfile(datapath,n(1).name);
                loader_2 = load(nuisanceregfile1);
                noiser{r} = loader_2;
                
                
            end
        end
    end
end

file_list = cat(1,filename(:));
r_col = (n_vol*1.4)*(0:1:size(onset_set,2)-1);
onset_set = onset_set+r_col;
 
if detection_con
    det_ = detect_set(:);
    det_(isnan(det_)) = 0;
    det_ = logical(det_);   
    onset_ = onset_set(:);
    onsets{1} = onset_(det_);
    onsets{2} = onset_(~det_);
    durations{1} = 0;
    durations{2} = 0;
    names{1} = 'sniff';
    names{2} = 'nosniff';    
else
    onsets{1} =  onset_set(:);
    durations{1} = 0;
    names{1} = 'sniffCue';
    if parametric_intensity
        cids = cid_set(:);
        
        load(fullfile(root, sprintf('NEMO_%02d',s), 'behavior', 'imaging_task','ratings_.mat'))
        cid_keys = cid;
        cid_vals = (1:1:length(cid_keys))';
        cid_map = containers.Map(cid_keys, cid_vals);
        
        pmod = struct('name',{''},'param',{},'poly',{});
        pmod(1).name{1}  = 'Intensity';
        
        idx_ = [];
        for ii = 1:length(cids)
            idx_ = [idx_; cid_map(cids(ii))];
        end
        
        pmod(1).param{1} = ratings_(idx_,1);
        pmod(1).poly{1}  = 1;
    end
end

save(fullfile(statpath, sprintf('conditions_%s_%s.mat',modelname, sn)), 'onsets','durations', 'names');
nc = dir(fullfile(statpath,sprintf('conditions_%s_%s.mat',modelname,sn)));
condfile = fullfile(statpath,nc(1).name);

noise_profiles = cat(1,noiser{:});
dlmwrite(fullfile(statpath, sprintf('nusiance_%s_%s.txt',modelname, sn)),noise_profiles)
n = dir(fullfile(statpath,sprintf('nusiance_%s_%s.txt',modelname, sn)));
nuisanceregfile1 = fullfile(statpath,n(1).name);

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = file_list;   
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {condfile};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nuisanceregfile1};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

path = fullfile(root, sn, 'jobs');
mkdir(path);
fname = fullfile(path, sprintf('job_%s_%s.mat',modelname,  sn));
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);

% Concatenation of basis functions
scans = n_vol*ones(1,length(horzcat(nruns{set_i:nsets})));
model_loc = fullfile(statpath,'SPM.mat');
spm_fmri_concatenate(model_loc,scans)

% estimate
matlabbatch = [];
fname = fullfile(statpath, 'SPM.mat');
matlabbatch{1}.spm.stats.fmri_est.spmmat =  {fname};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);


end

% function file for the model with multiple inputs
function make_conditions(s, root, statpath, modelname, set_, sess, rr)

respath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', 'imaging_task',sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
filename = fullfile(respath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr));
load(filename)
breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', 'imaging_task', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
load(fullfile(breathpath, sprintf('time_adjust_NEMO%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
% 1     odor order (1-10)
% 2     CID of odor
% 3     onset of odor trigger from t0
% 4     onset of sniff cue from t0
% 5     percept to be rated
% 6     Detect rating (0 = no smell)
% 7     Mouse button pressed
% 8     Detect RT
% 9     Time at which 8 submitted from t0
% 10    Percept rating
% 11    Percept RT
% 12    Time at which 11 submitted from t0
% if res.sess==1
%     percept_list = {' Intensity';' Pleasantness';' Fishy';' Burnt';' Sour';' Decayed';' Musky';' Fruity'; ' Sweaty';' Cool'};
% elseif res.sess ==2
%     percept_list = {' Intensity';' Pleasantness';' Chemical-like';' Flowery';' Sweet';' Warm';' Bakery-like';' Garlic-like'; ' Spicey';' Acidic'};
% else
%     percept_list = {' Intensity';' Pleasantness';' Fishy';' Burnt';' Sour';' Decayed';' Chemical-like';' Flowery';' Sweet';' Warm'};
% end

d = res.select_odor;

% adjust times
d(:,[3,4,9,12]) = d(:,[3,4,9,12])/1000 + adjust;

onsets{1} = d(:,4);
detect = d(:,6);
cid = d(:,2);
if set_==3
    cid(cid==1032)=1001;
end
save(fullfile(statpath, sprintf('conditions_%s_NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat',modelname, s, set_, sess, rr)), 'onsets', 'detect','cid');
end
