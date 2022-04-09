%% Pattern discriminability of odors
% Major inputs:
% fullmat_sess<sess_number>_old: A sesswise estimate of odor responses 
% of shape voxels x hrf bases x odors. Estimated from FIR_model.m
% maskfile: mask of gray matter voxels in the whole brain
% anat_masks: masks of ROIs in native space of subjects.
% fmaskfile: functional mask of voxels with significant odor evoked
% activity (obtained from FIR_model.m, sniff_contrast = true)

% Major output:
% rsa_P1: A matrix specifying bootstrapped estimates of pattern
% discriminability of odors for all areas and averaged across subejcts and sessions.

% Vivek Sagar (VivekSagar2016@u.northwestern.edu). April 9, 2022
%% Settings
vsn = 0;
fmasker = true; % Only select files with significant odor evoked activity.
anat_names = {'APC','PPC','Amygdala','OFC'};
anat_masks = {'rwAPC.nii', 'rwPPC.nii','rwAmygdala.nii','rwofc.nii'};
maskfile =  'anat_gw.nii';
fmaskfile = 'f_anat_gw3.nii';
sess_l = cat(3,nchoosek([1 2 3],2),nchoosek([2 3 4],2),nchoosek([2 3 4],2),nchoosek([2 3 4],2));
dirs = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\Mnemonic';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\Mnemonic';
    '';'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\Mnemonic'};

nperm = 10000;
nodor = 160;
[~,bootsam] = bootstrp(nperm,@(x) x,1:nodor);
utl_mask_diag = logical(eye(nodor));
% load(fullfile(statpath,'fir_cv.mat'))
anat_idx = [true(1,length(anat_names))]; % Analyze these ROIs
rsa_sc = zeros(length(anat_names),nperm);
utl_mask = logical(eye(nodor));

fprintf('\n')
% Subject - index
for s = [1 2 4] % Subject
    fprintf('Subject: %02d\n',s)
    statpath = dirs{s};
    % Gray Matter, Functional and Anatomical Masks
    mask = (spm_read_vols(spm_vol(fullfile(statpath, maskfile)))); % Mask used to construct odor files
    mask(isnan(mask))=0;
    mask = logical(mask);
    fmask = (spm_read_vols(spm_vol(fullfile(statpath, fmaskfile)))); % Mask used to examine voxels in RSA
    fmask(isnan(fmask))=0;
    if fmasker
        fmask = logical(fmask); % Only choose voxels with significant odor evoked activity
    else
        fmask = logical(fmask+0.1); % Select all voxels
    end
    fmask_1d = fmask(mask);
    anatpath = fullfile('C:\Data\NEMO\',sprintf('NEMO_%02d',s),'\imaging\nii\masks');
    
    % Model names
    masks_set = [];
    for ii = 1:length(anat_masks)
        m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
        m1(m1<=0)=0;
        m1(m1>0) = 1;
        m1 = m1(mask);
        masks_set(:,ii)=m1(fmask_1d);
    end
   
    masks_set(isnan(masks_set))=0;
    linux_config = false;
    warning('off','all')
    
    rsa_P1 = zeros(length(anat_names),nperm,3);
    for perms = 1:nperm
        boot_vec = bootsam(:,perms);
        if mod(perms,25)==0
            fprintf('perms: %03d\n',perms)
        end
        % Session - index
        for zz_i = 1:3
            sess1 = squeeze(sess_l(zz_i,1,s));
            sess2 = squeeze(sess_l(zz_i,2,s));
            
            % Voxel_responses
            S1 = load(fullfile(statpath,sprintf('fullmat_sess%02d_old.mat',sess1)),'odor_responses');
            S1_odors = squeeze(S1.odor_responses(:,1,:));
            S1_odors = S1_odors(:,boot_vec);
            S2 = load(fullfile(statpath,sprintf('fullmat_sess%02d_old.mat',sess2)),'odor_responses');
            S2_odors = squeeze(S2.odor_responses(:,1,:));
            S2_odors = S2_odors(:,boot_vec);
            %% Decoding analysis
            % Permutation test and index
                        
            for ii = find(anat_idx)
                %     fprintf('Anat: %02d\n',ii)
                % Session-1
                S1_omat = squeeze(S1_odors(fmask_1d,:));
                S1_omat_vals =  S1_omat(logical(masks_set(:,ii)),:);
                [r1,~] = find(isnan(S1_omat_vals));
                %     S1_omat_vals = zscore(S1_omat_vals,0,2); % Std across all rows (voxels)
                
                % Session-2
                S2_omat = squeeze(S2_odors(fmask_1d,:));
                S2_omat_vals =  S2_omat(logical(masks_set(:,ii)),:);
                [r2,~] = find(isnan(S2_omat_vals));
                %     S2_omat_vals = zscore(S2_omat_vals,0,2);
                
                r = unique([r1; r2]);
                S1_omat_vals(r,:) = [];
                S2_omat_vals(r,:) = [];
                
                % Main RSA
                cross_M = corrcoef_2(S1_omat_vals',S2_omat_vals');
                ondiag = cross_M(utl_mask);
                offdiag = cross_M(~utl_mask);
                %     rsa_P1(ii) =  mean(sign(ondiag).*(ondiag.^2))-mean(sign(offdiag).*(offdiag.^2));
                rsa_P1(ii,perms,zz_i) =  mean(atanh(ondiag))-mean(tanh(offdiag));
            end            
        end
    end
     rsa_P1 = mean(rsa_P1,3);       
     save(fullfile(statpath,sprintf('r-vsn%01d',vsn)),'rsa_P1')
end

%% Average across voxels
dirs2 = {'C:\Data\NEMO\NEMO_01\imaging\1stlevelmodels\Mnemonic';
    'C:\Data\NEMO\NEMO_02\imaging\1stlevelmodels\Mnemonic';
    'C:\Data\NEMO\NEMO_04\imaging\1stlevelmodels\Mnemonic'};
vsn = 0;
P1s = variable_extract(dirs2,'r-vsn0.mat','rsa_P1',false);
rsa_P1 = cat(3,P1s{:});
rsa_P1 = mean(rsa_P1,3);

figure('Position',  [100, 100, 480, 320])
hold on
bar(mean(rsa_P1,2))
errorbar(1:length(anat_names),mean(rsa_P1,2)',prctile(rsa_P1',99.5)-mean(rsa_P1,2)',prctile(rsa_P1',0.5)-mean(rsa_P1,2)','.')
hold on
c = {'r','g','b'};
for jj = 1:3
    plot(1:length(anat_names),mean(P1s{jj},2),c{jj});
end
xticks([1:length(anat_names)])
xticklabels(anat_names)
ylabel('r (on-off)')
save(sprintf('r-vsn%01d',vsn),'rsa_P1','rsa_sc')
savefig(sprintf('r-vsn%01d',vsn))
print(sprintf('r-vsn%01d',vsn),'-dpng')