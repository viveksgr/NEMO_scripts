% Compute temporal signal to noise ratio of the bold response 
% Inputs: preprocessed T2 images from make_preprocessing_job.

s = 1; % Subject ID
n_volumes = 856; % Number of volumes in each run
root = 'C:\Data\NEMO'; % Windows Version
statpath = pwd;

% Filelist
nset_i = 1;
nset_f = 4;
nsess_i = 1;
nsess_f = 3;
nruns = {[1:4] [1:4] [1:4] [1:4]}; 

anat_names = {'APC','PPC','Amygdala','OFC'};
anat_masks = {'rwAPC.nii', 'rwPPC.nii','rwAmygdala.nii','rwofc.nii'};
anatpath = fullfile(root,sprintf('NEMO_%02d',s),'\imaging\nii\masks');
maskfile =  'anat_gw.nii';
mask = (spm_read_vols(spm_vol(fullfile(anatpath, maskfile)))); % Mask used to construct odor files
mask(isnan(mask))=0;
mask = logical(mask);
fmaskfile = 'f_anat_gw3.nii';
fmask = (spm_read_vols(spm_vol(fullfile(anatpath, fmaskfile)))); % Mask used to examine voxels in RSA
fmask(isnan(fmask))=0;
fmask = logical(fmask);
sn_ =  sprintf('NEMO_%02d',s);
datapath = fullfile(root, sn_,  'imaging', 'nii'); % sn:subject's name

% Filelist
rcntr = 0;
filename = [];
for set_=nset_i:nset_f
    for sess = nsess_i:nsess_f
        for r = nruns{set_}
            rcntr = rcntr+1;
            % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
            path_ = fullfile(datapath, sprintf('set_%02d', set_),...
                sprintf('sess_%02d', sess), sprintf('run_%02d', r));
            n = dir( fullfile(path_, sprintf('srf*.nii')));
            
            for i=1:length(n)
                % Different runs in different cells, add ",1" for spm
                filename{rcntr}{i,1} = fullfile(path_, n(i).name);
            end
        end
    end
end
filename = vertcat(filename{:});
Res_V = spm_vol(filename);
[~, XYZmm] = spm_read_vols(Res_V{1});
XYZvx = round(Res_V{1}.mat\[XYZmm; ones(1,size(XYZmm,2))]);
Res_Vol = cell2mat(Res_V);

% Loop over masks
m_std = cell(1,length(anat_names));
for mm = 1:length(anat_names)
    fprintf('Area: %02d\n',mm)
    % Mask
    m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{mm})));
    m1(isnan(m1))=0;
    m1(m1<=0)=0;
    m1(m1>0) = 1;
    m1 = and(logical(m1),mask);
    m1 = and(m1,fmask);

    % Load data in batches
    anatmask1D = m1(:); 
    vxl = XYZvx(:,anatmask1D);
    binsize = 100;
    q_ = floor(size(vxl,2)/binsize);
    r_ = rem(size(vxl,2),binsize);
    nbins = q_+1;
    
    mmcell = cell(nbins,1);   
    for nn = 1:nbins
        fprintf('Bin: %02d\n',nn)
        if nn<nbins
            bins = (nn-1)*binsize+1:nn*binsize;
        else
            bins = (nn-1)*binsize+1:size(vxl,2);
        end
        voxel_act = spm_get_data(Res_Vol,vxl(:,bins))'; 
        mmcell{nn} = mean(voxel_act,2)./std(voxel_act,[],2);
    end
    m_std{mm} = vertcat(mmcell{:});
end
mkdir(sprintf('subj_%02d',s))
save(fullfile(sprintf('subj_%02d',s),'tsnr.mat'),'m_std')

sz = cellfun(@length,m_std);
g = {}; for gg = 1:length(anat_names); g{gg} = gg*ones(sz(gg),1); end
boxplot(vertcat(m_std{:}),vertcat(g{:}))
ylabel('t-snr')
xticklabels(anat_names)
savefig(fullfile(sprintf('subj_%02d',s),'tsnr.fig'))

% %% Group plot
% dirs = {'C:\Data\NEMO\NEMO_all\t-snr\subj_01';
%         'C:\Data\NEMO\NEMO_all\t-snr\subj_02';
%         'C:\Data\NEMO\NEMO_all\t-snr\subj_04'};
% nS = length(dirs);
% matname = 'tsnr.mat';
% anat_names = {'PirF','PirT','AMY','OFC'};
% mstd = variable_extract(dirs,matname,'m_std',true);
% 
% % bootstrap
% % nboot = 10000;
% % func = @(x) bootstrp(nboot, @mean, x);
% % mean_b = cellfun(func,mstd,'UniformOutput',false);
% % mean_b_cat = {}; for ii = 1:length(anat_names); mean_b_cat{ii} = mean(horzcat(mean_b{:,ii}),2); end
% % mean_b = vertcat(mean_b{:,:});
% meaner =cellfun(@mean,mstd);
% 
% mean_ = cellfun(@mean,mean_b_cat);
% std_ = cellfun(@std,mean_b_cat);
% cs = {'r','g','b'};
% figure('Position',[0.5 0.5 300 200])
% bar(1:4,mean_)
% hold on
% errorbar(1:4,mean_,std_,'.')
% for s=1:3; plot(1:4,meaner(s,:),cs{s}); end
% xticklabels(anat_names)
% p = anova1(meaner);