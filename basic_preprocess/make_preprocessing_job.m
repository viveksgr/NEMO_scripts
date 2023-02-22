% Realignment, coregistration and smoothing of functional images. 
% Needs a T1 reference image.
% Folder structure is specified as follows.
% root: parent directory for each subject
% root>imaging>nii: all nii images
% root>imaging>nii>set[number]>sess[number]>run[number]
% In this dataset, there are 4 sets of 40 odors each. Each set has 3
% sessions, and each session has 4 runs. Each run has data corresponding to
% 10 odors or 90 trials. 

%% Basic settings
linux_config = false;
if linux_config
    root = '/projects/p30489/Data/NEMO'; % Quest Version
    addpath('/home/vsh3681/spm12')
else
    root = 'C:\Data\NEMO'; % Windows Version
%   addpath('C:\spm12')
end
subjects = [1];
nsubs = length(subjects);

nset_i = 1;
nset_f = 4;
nsess_i = 2;
nsess_f = 4;
nsess = nsess_f-nsess_i+1;
nsets = nset_f-nset_i+1;
nruns = {[1:4] [1:4] [1:4] [1:4]}; % Across sets for all sess

for s = subjects
    matlabbatch = [];
    sn_ =  sprintf('NEMO_%02d',s);
    datapath = fullfile(root, sn_,  'imaging', 'nii'); % sn:subject's name
    
   
    %% REALIGNMENT OF WHOLE BRAIN
    filename = cell(1,nsets*nsess);
    ss=0;

    for set_=nset_i:nset_f
        for sess = nsess_i:nsess_f
            ss=ss+1;
            % path for whole brain images.
            path_ = fullfile(datapath,sprintf('set_%02d', set_), ...
                   sprintf('sess_%02d', sess), 'wb'); 
            % *.nii will give all the files.
            n = dir(fullfile(path_, sprintf('fNEMO*.nii')));         
            if set_==nset_i && sess==nsess_i
                wbpath = path_;
                % Filename of the first file in the directory.
                wbfile = n(1).name;                                  
            end
            for i=1:length(n)
                fname = fullfile(path_, n(i).name);
                % SPM job file has format (----.nii,1)
                filename{ss}{i,1} = sprintf('%s,1', fname);          
            end
        end
    end
    
    matlabbatch{1}.spm.spatial.realign.estwrite.data = filename;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
     
    %% Functional scan realignment
    % rnctr is run counter
    rcntr = 0;                                                          
    filename = [];
    for set_=nset_i:nset_f
        for sess = nsess_i:nsess_f
            for r = nruns{set_}
                rcntr = rcntr+1;
                % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
                path_ = fullfile(datapath, sprintf('set_%02d', set_),...
                       sprintf('sess_%02d', sess), sprintf('run_%02d', r));
                n = dir( fullfile(path_, sprintf('fNEMO*.nii')));
                
                if rcntr==1
                    funcpath = path_;
                    funcfile = n(1).name;
                end
                
                for i=1:length(n)
                    fname = fullfile(path_, n(i).name);
                    % Different runs in different cells, add ",1" for spm
                    filename{rcntr}{i,1} = sprintf('%s,1', fname);    
                end
            end
        end
    end
    % Corregistration needs files from different runs in one big array.                           
    fimages = cat(1,filename{:});
        
    matlabbatch{2}.spm.spatial.realign.estwrite.data = filename;      
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    
    %% COREGISTRATION
    % Anatomical and whole brain coreg
    
    % define path for anatomical image
    anatpath = fullfile(datapath, 'anat');    
    n = dir(fullfile(anatpath, sprintf('sNEMO*.nii')));
    anatfile = n(1).name;                           
    fname = fullfile(anatpath, anatfile);
    matlabbatch{3}.spm.spatial.coreg.estimate.ref = {fname};
    fname = fullfile(wbpath, sprintf('mean%s', wbfile));                
    matlabbatch{3}.spm.spatial.coreg.estimate.source = {fname};
    matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % Whole brain and functional scan coreg
    fname = fullfile(wbpath, sprintf('mean%s', wbfile));
    matlabbatch{4}.spm.spatial.coreg.estimate.ref = {fname};
    fname = fullfile(funcpath, sprintf('mean%s', funcfile));
    matlabbatch{4}.spm.spatial.coreg.estimate.source = {fname};
    matlabbatch{4}.spm.spatial.coreg.estimate.other = fimages;
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol =...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    
    %% Normalize
    fname = fullfile(anatpath, anatfile);
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.vol = {fname};
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample =  {fname};
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spm('dir'),'tpm','TPM.nii')};
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{5}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{5}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
    matlabbatch{5}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    
    %     % Normalise write
    %     fname = fullfile(anatpath, sprintf('y_%s', anatfile));
    %     matlabbatch{6}.spm.spatial.normalise.write.subj.def = {fname};
    %     matlabbatch{6}.spm.spatial.normalise.write.subj.resample = fimages;
    %     matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    %     matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    %     matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
    
    %%REALIGN-RESLICE FUNCTIONAL IMAGES
    matlabbatch{6}.spm.spatial.realign.write.data = fimages;
    matlabbatch{6}.spm.spatial.realign.write.roptions.which = [2 1];
    matlabbatch{6}.spm.spatial.realign.write.roptions.interp = 4;
    matlabbatch{6}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{6}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{6}.spm.spatial.realign.write.roptions.prefix = 'r';
    
    
    %% Spatial smoothing
    ii=0;
    filename = [];
    for set_ = nset_i:nset_f
        for sess = nsess_i:nsess_f
            for r = nruns{set_}
                path_ = fullfile(datapath, sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d', r));
                n = dir(fullfile(path_, sprintf('fNEMO*.nii')));
                for i=1:length(n)
                    ii=ii+1;
                    fname = fullfile(path_, sprintf('r%s',n(i).name));
                    filename{ii,1} = sprintf('%s,1', fname);
                end
            end
        end
    end
    matlabbatch{7}.spm.spatial.smooth.data = filename;
    matlabbatch{7}.spm.spatial.smooth.fwhm = [2 2 2];
    matlabbatch{7}.spm.spatial.smooth.dtype = 0;
    matlabbatch{7}.spm.spatial.smooth.im = 0;
    matlabbatch{7}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run', matlabbatch);    % SPM jobman will run the batch
    
end

% %% Average T1 images
% 
% % clearvars -except anatpath
% n = dir(fullfile(anatpath, sprintf('sNEMO*.nii')));
% % Order it manually
% anatfile = n(1).name;
% fname = fullfile(anatpath, anatfile);
% matlabbatch = cell(length(n)-1,1);
% for ii = 2:length(n)
% %     fprintf('Counter: %d',ii)
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.ref = {fname};
%     fname2 = fullfile(anatpath, n(ii).name);
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.source = {fname2};
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.other = {''};
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
%     matlabbatch{ii-1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
% end
% spm_jobman('run', matlabbatch); 