# NEMO_scripts

A repository for the neural encoding models of olfaction (NEMO) project. To run, clone the repository and modify the "root" variable in the scripts to the working directory. 

Dependencies:
1. [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
2. COGENT2000 

Workflow:
1. task_scripts
Script (NEMO_outline.m) and additional functions to run the perceptual rating task. COGENT required.
2. basic_preprocess:
Preprocessing of .nii images (make_preprocessing_job.m), computation of nuisance parameters (make_nuisance_regressors.m), generation of masks (make_inverse_normalize_mask.m), and computation of t-snr of files (t_snr.m). Scripts to analyze the perceptual properties of the odors (behavioral_analysis.m).
3. FIR_model:
FIR model to extract a scalar odor response for a given voxel in a subject at an instant of time.
4. Mnemonic:
Decoding analysis to compare similarity in session-wise response of an odor (compared to similarity between two different odors).
5. RSA_ROI:
Representational Similarity Analysis (for chemical vs perceptual and category vs identity)
6. EEM_LOOCV:
The encoding model to predict odor responses using perceptual descriptors
7. EEM_dimensionality:
Performance of EM for different ROI, dimensionality estimates.
8. EEM_subjectivity:
Computation of subjectivity of encoding

common_functions: library of custom functions used during all parts of the analysis.

Data description:
Files to replicate RSA analysis, encoding model estimation, perceptual weight profiles, dimensionality and subjectivity analyses. 
For each subject, the following files are provided:
1. Masks and nii images: 
    a. sNEMO0*_defaced.nii = T1 structural image
    b. rw<masks>nii = masks for PirF (rw.APC,nii), PirT (rw.PPC.nii), AMY (rwAmygdala.nii), OFC (rwOFC.nii)
    c. ARC3_anatgw.nii = Gray matter mask. 
    d. ARC3_fanatgw3_pos.nii = Functional mask based on odor responses
2. Perceptual and chemical ratings:
    a. behav_ratings_NEMO0*.mat = structure of behavioral ratings specifying odor ids, odor ratings and names of the columns in ratings.
    b. chem_corr_NEMO0*.mat = correlation matrix of molecular ratings
3. Odor responses: Fullmet_FIR.mat. Output from FIR_model.m consisting of odor responses in shape nvoxels X FIR components X number of odors. 
4. Supporting files for the encoding model:
    a. cv_fold_new.mat: index specifying which odors belong to which fold of crossvalidation
    b. argmaxes2.mat: for each voxel, specifies the HRF bin chosen in each fold of crossvalidation
5. behavioral data:
    a. xlsx files specifying which odors were used to train the subjects 2 and 3 on the perceptual rating task
    b. odor_names: odor ids, odor names and concentration of odors. Odor IDs can be used to compare with other studies or search for odor properties in Pubchem.

Scripts may not be compatible for linux and/or MacOS. 
Additional useful sets/scripts not in the repository:
1. [Computing inverse percentile](https://www.mathworks.com/matlabcentral/fileexchange/41131-inverse-percentiles-of-a-sample)
2. [Stable-marriage algorithm](https://www.mathworks.com/matlabcentral/fileexchange/44262-gale-shapley-stable-marriage-algorithm)


Raw Data:
Raw Data is available upon request at: https://zenodo.org/record/7636722
The recommended way to store the data is to specify a root directory and construct a separate subdirectory for each subject. This subject subdirectory should contain further subdirectories for breathing, behavior and imaging data, further organized into sets, sessions and runs. 