# NEMO_scripts

A repository for the neural encoding models of olfaction (NEMO) project. Documentation is under progress.

Dependencies:
1. SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

Workflow:
1. Behav_process 
Scripts to analyze the perceptual and chemical properties of the odors.
2. Basic preprocess
Preprocessing of .nii images, generation of masks and sniff contrast, and computation of neural RSM and t-snr.
3. FIR_model
FIR model to extract a scalar odor response for a given voxel in a subject at an instant of time.
4. Mnemonic
Decoding analysis to compare similarity in session-wise response of an odor (compared to similarity between two different odors).
5. RSA_ROI
Representational Similarity Analysis (for chemical vs perceptual and category vs identity)
6. RSA_numodors
Performance of RSA as a function of number of odors
7. EEM_LOOCV
The encoding model to predict odor responses using perceptual descriptors
8. group_plots
Performance of EM for different ROI, dimensionality estimates and subjectivity analysis.
