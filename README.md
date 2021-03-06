# NEMO_scripts

A repository for the neural encoding models of olfaction (NEMO) project. 

Dependencies:
1. SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

Workflow:
1. task_scripts
Script (NEMO_outline.m) and additional functions to run the perceptual rating task. Requires COGENT2000. 
2. basic_preprocess:
Preprocessing of .nii images (make_preprocessing_job.m), computation of nuisance parameters (make_nuisance_regressors.m), generation of masks (make_inverse_normalize_mask.m), and computation of t-snr of files (t_snr.m). Scripts to analyze the perceptual properties of the odors (behavioral_analysis.m).
3. FIR_model:
FIR model to extract a scalar odor response for a given voxel in a subject at an instant of time.
4. Mnemonic:
Decoding analysis to compare similarity in session-wise response of an odor (compared to similarity between two different odors).
5. RSA_ROI:
Representational Similarity Analysis (for chemical vs perceptual and category vs identity)
6. RSA_numodors:
Performance of RSA as a function of number of odors
7. EEM_LOOCV:
The encoding model to predict odor responses using perceptual descriptors
8. EEM_dimensionality:
Performance of EM for different ROI, dimensionality estimates.
9. EEM_subjectivity:
Computation of subjectivity of encoding

common_functions: library of custom functions used during all parts of the analysis.
