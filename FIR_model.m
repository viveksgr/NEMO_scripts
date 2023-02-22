%% Voxel-wise estimate of odor responses
% FIR Model to extract scalar response per voxel per odor for a range of
% times.

% Specify:
% datapath: Path to preprocessed T2 images and nuisance regressors
% behavpath: Path to behavioral dataset
% breathpath: Path to breathing data
% statpath: Gray matter mask and save location of betamaps

% Output
% name_chunk: NxTxO matrix for N voxels in the Gray matter mask, T times
% and O odors. For canonical HRF, T =1.

% Vivek Sagar (VivekSagar2016@u.northwestern.edu). April 9, 2022

%% Settings
s = 2; % Subject
n_volumes = 876; % Number of volumes per session
hrf_basestype = true; % Bases: HRF (true) or FIR (false)
nfir_comp = 1; % Set this == 1 for canonical HRF.
nfir_l = 11; % Choose number of FIR components for FIR
set_i = 1; % Initial set
set_f = 4; % Final set. Use set_i == set_f for setwise data
sess_i = 1; % Initial session (for all sets)
sess_f = 3; % Final session (for all sets). Use sess_i == sess_f for sesswise data
modelname = 'NEMO02_can_sess03';
nruns = {[1:4],[1:4],[1:4],[1:4]}; % Run indices used in each set
name_chunk = 'fullmat.mat';
nodors = 160; % Number of odors
exec_model = true; % Compute bases. False, if only compiing precomputed bases.
sniff_contrast = true; % Create a contrast for sniff-evoked activity
compile_model = true; % Compile bases functions in a 3D matrix

odor_shift = false; % Sanity check to make sure no perceptual information before odor delivery

% Initial params
root = 'C:\Data\NEMO'; % Windows Version
addpath('C:\spm12')
statpath = fullfile(root, sn, 'imaging', '1stlevelmodels', modelname);
fullfpath = statpath;

if hrf_basestype
    n_bases = length(horzcat(nruns{:}))*10*(set_f-set_i+1);
else
    n_bases = length(horzcat(nruns{:}))*10*nfir_comp;
end
n_sets = set_f-set_i+1;
sn = sprintf('NEMO_%02d', s);
nruns_ = length(horzcat(nruns{:}))*(sess_f-sess_i+1);

% spm defaults fmri
spm_jobman initcfg
spm_get_defaults('cmdline',true)

%% Load/Create SPM.mat
if exec_model
    r = 0;
    onset_set = cell(nruns_,1);
    odor_id = cell(nruns_,1);
    filename = cell(n_volumes,nruns_);
    noise_files = cell(nruns_,1);
    for set_ = set_i:set_f
        for sess = sess_i:sess_f
            for rr=nruns{set_}
                r= r+1;
                
                if s==1
                    behavpath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(behavpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr)));
                    breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(breathpath, sprintf('time_adjust_NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
                    datapath = fullfile(root, sn, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));
                    ns = dir(fullfile(datapath,sprintf('nusiance_regresssors_NEMO_%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess,rr)));
                    d = res.select_odor{1};
                else
                    behavpath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', 'imaging_task',sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(behavpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr)));
                    breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', 'imaging_task', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(breathpath, sprintf('time_adjust_NEMO%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
                    datapath = fullfile(root, sn, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));
                    ns = dir(fullfile(datapath,sprintf('nuisance_regresssors_NEMO%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess,rr)));
                    d = res.select_odor;
                end
                % Col.  Description
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
                if ~odor_shift
                    onset_set{r} = d(:,4)/1000 + adjust; % Condition onset
                else
                    onset_set{r} = d(:,4)/1000 + adjust-odor_shift; % Condition onset
                end
                odor_id{r} = d(:,2); % IDs of odor presented
                
                %             datapath = fullfile(root, sn, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));
                n = dir(fullfile(datapath,'srf*.nii'));
                for i= 1:length(n)
                    fname = fullfile(datapath, n(i).name);
                    filename{i,r}= sprintf('%s,1',fname); % Datafiles
                end
                
                noise_files{r} =  load(fullfile(datapath,ns(1).name));
            end
        end
    end
    
    %% Vertical concatenation of all runs
    file_list = cat(1,filename(:)); % All *.nii images to be used in training
    file_V = spm_vol(file_list);
    file_V = cell2mat(file_V);
    save(fullfile(statpath, sprintf('names_%s_%s.mat',modelname, sn)), 'file_V');
    
    r_col = (n_volumes*1.4)*(0:1:length(onset_set)-1); % Adjust onset times for vertical concatenation
    onset_times = cell(size(onset_set));
    for ii = 1:length(r_col)
        onset_times{ii} = onset_set{ii}+r_col(ii);
    end
    
        [onsets,names] = BRS_design_run(vertcat(odor_id{:}),vertcat(onset_times{:}),sess_f-sess_i+1);

    durations = repmat({0},1,length(onsets));
    
    save(fullfile(statpath, sprintf('conditions_%s_%s.mat',modelname, sn)), 'onsets','durations', 'names');
    nc = dir(fullfile(statpath,sprintf('conditions_%s_%s.mat',modelname,sn)));
    condfile = fullfile(statpath,nc(1).name);
    
    % Noise profiles
    noise_profiles = cat(1,noise_files{:});
    dlmwrite(fullfile(statpath, sprintf('nuisance_%s_%s.txt',modelname, sn)),noise_profiles)
    n = dir(fullfile(statpath,sprintf('nuisance_%s_%s.txt',modelname, sn)));
    nuisanceregfile = fullfile(statpath,n(1).name);
    
    %% Config. SPM-GLM
    % Model specification
    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_spec.dir = {statpath};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.4;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    if hrf_basestype
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    else
        matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = nfir_l;
        matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order = nfir_comp;
    end
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    % matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(statpath,'anat_gw.nii')};
    m = spm_read_vols(spm_vol(fullfile(statpath,'anat_gw.nii')));
    m_cell{1} = m;
    matlabbatch{1}.spm.stats.fmri_spec.mask = m_cell;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = file_list;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {condfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nuisanceregfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    spm_jobman('run', matlabbatch);
    
    % Concatenation of basis functions
    scans = n_volumes*ones(1,nruns_);
    model_loc = fullfile(statpath,'SPM.mat');
    spm_fmri_concatenate(model_loc,scans)
    fprintf('Specification complete...')
    
    matlabbatch = [];
    fname = fullfile(statpath, 'SPM.mat');
    matlabbatch{1}.spm.stats.fmri_est.spmmat =  {fname};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch);
end

%% Create sniff contrast
if sniff_contrast
    statpath = fullfpath;
    matlabbatch = [];
    contrast_vec = ones(1,nodors);
    
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(statpath,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'SniffT';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights =  contrast_vec;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    
    spm_jobman('run', matlabbatch)
end


%% Compile model
if compile_model
    statpath = fullfpath;   
    file_list = dir(fullfile(statpath,sprintf('beta*.nii')));   
    mask = spm_read_vols(spm_vol(fullfile(statpath,'anat_gw.nii')));
    mask(isnan(mask))=0;
    mask = logical(mask);
    
    filename = cell(size(file_list));
    for i= 1:length(file_list)
        fname = fullfile(statpath, file_list(i).name);
        filename{i,1}= sprintf('%s,1',fname); % Datafiles
    end
    file_V = spm_vol(filename);
    file_V = cell2mat(file_V);
    file_V = file_V(1:n_bases,:);
    
    [Y, ~] = spm_read_vols(file_V(1));
    dims = size(Y);
    dims2 = dims(2);
    dims3 = dims(3);
    
    odor_responses = zeros([dims n_bases]);
    odor_responses2 = zeros([dims n_bases]);
    for ii = 1:dims(1)
        fprintf('|')
        %         fprintf(sprintf('%d\n',size(odor_responses,4)))
        for jj = 1:dims2
            for kk = 1:dims3
                if mask(ii,jj,kk)~=0
                    odor_resp = spm_get_data(file_V,[ii;jj;kk]);
                    odor_responses(ii,jj,kk,:) = odor_resp; % Compile beta maps
                end
            end
        end
    end
    
    % Change 4D to 2D array
    mask_1d = logical(reshape(mask,1,[]));
    voxel_ids = find(mask_1d);
    two_D_ = zeros(size(voxel_ids,1),n_bases); % Two dimensional representation of odor amplitude
    two_D_nn = zeros(size(voxel_ids,1),n_bases); % Two dimensional representation of odor amplitude
    for ii = 1:length(voxel_ids)
        [xx,yy,zz]  = ind2sub(size(mask),voxel_ids(ii)); % Index of the mask
        two_D_(ii,:) = odor_responses(xx,yy,zz,:);
    end
    
    % Add third dimension (of bases function used)
    odor_responses = zeros(size(two_D_,1),nfir_comp,ceil(size(two_D_,2)/nfir_comp));
    if hrf_basestype
        odor_responses(:,1,:) = two_D_;
    else
        r_count = 0:(n_bases/nfir_comp)-1;
        for ii = 1:length(voxel_ids)
            for jj = 1:nfir_comp
                r_vec = r_count*nfir_comp+jj;
                odor_responses(ii,jj,:)=two_D_(ii,r_vec);
            end
        end
    end
    save(fullfile(statpath,name_chunk)) 
end