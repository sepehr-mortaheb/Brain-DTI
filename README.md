# Brain Alterations after Spaceflight (Brain-DTI)  

<p align="center">
<img src="img.jpg" alt="" height="400"/>
</p>

## Introduction
With this project, brain MRI scans are acquired in astronauts from the European Space Agency (ESA) and cosmonauts from the Russian Space Agency (Roscosmos). 
MRI scans are acquired before, shortly after, and 6-7 months after a mission to the International Space Station (ISS) to study spaceflight's immediate and long-term effects 
on the brain. With each scan, anatomical T1 and T2-weighted images, diffusion MRI, resting-state functional MRI, and task-based functional MRI are acquired. 
These complementary MRI scans allow us to study brain structural and functional changes after a mission to the ISS. The goal is to detect sites of neuroplasticity 
and to study microgravity-induced fluid shifts.

This repository contains the codes used to analyze different data modalities.

## fMRI Preprocessing Pipeline
The fMRI spatial preprocessing pipeline has been developed based on the [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and [CAT12](https://neuro-jena.github.io/cat/). 
This pipeline consists of the following steps: 

- Slice Timing Correction 
- Realignment
- Segmentation of the structural data
- Coregistration of the functional data into the T1 image.
- Normalization to the MNI space
  - It uses the forward transformation parameters of the segmentation step, to normalize the functional data into the MNI space. 
- Smoothing 

### Practical Info: 

- In this pipeline, the **CAT12** has been used for the segmentation as it gives more precise results than the **SPM** segmentation. So, make sure that you have downloaded the [CAT12](https://neuro-jena.github.io/cat/) package and have put it in the **SPM toolbox** folder.
- In this pipeline, the **ART** toolbox has been used for outlier volume detection. So, make sure that you have downloaded the [Artifact Detection Tools (ART)](https://www.nitrc.org/projects/artifact_detect/) and have added its directory to the Matlab paths.
- As the number of sessions is highly variable across subjects, we create session folders as much as the maximum number of available sessions for one subject,
but they can be empty if that session has not been conducted for the related subject.  
- The subjects usually have two anatomical T1 scans, one related to the DWI acquisition session and the other to the fMRI acquisition session. The pipeline by default takes
the T1 image of the fMRI session. But if the fMRI-related T1 scan is unavailable, the pipeline looks for the T1 image of the DWI scan to proceed with the preprocessing.     
- The pipeline considers that your data is organized in the BIDS format as follows:
  ```
  Data_dir -->
              sub-XXX
              sub-YYY
              .
              .
              .
              sub-ZZZ -->
                         ses-xxx
                         ses-yyy
                         .
                         .
                         .
                         ses-zzz -->
                                    anat -->
                                            sub-ZZZ_ses-zzz_acq-ge_run-dti_T1w.json
                                            sub-ZZZ_ses-zzz_acq-ge_run-dti_T1w.nii
                                            sub-ZZZ_ses-zzz_acq-ge_run-fmri_T1w.json
                                            sub-ZZZ_ses-zzz_acq-ge_run-fmri_T1w.nii
                                    func -->
                                            sub-ZZZ_ses-zzz_acq-ge_task-rest_bold.json
                                            sub-ZZZ_ses-zzz_acq-ge_task-rest_bold.nii
                                            sub-ZZZ_ses-zzz_acq-ge_task-navigation_bold.json
                                            sub-ZZZ_ses-zzz_acq-ge_task-navigation_bold.nii
                                            sub-ZZZ_ses-zzz_acq-ge_task-tennis_bold.json
                                            sub-ZZZ_ses-zzz_acq-ge_task-tennis_bold.nii
  ```
- To run the pipeline, open the `preprocessing.m` file, set the requested directories and parameters, and run the file.
- If you need to change the hyperparameters of different preprocessing steps, open the `func_PreprocBatch.m` file and set the parameters with your values accordingly.
