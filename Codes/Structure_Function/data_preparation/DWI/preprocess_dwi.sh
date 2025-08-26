#!/bin/bash
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Define your data and preprocessing directories
data_dir=/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/bids/ROS
preproc_dir=/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/preprocessed/ROS

# Loop over subjects
for subj_path in ${data_dir}/sub-*; do
    subid=$(basename "$subj_path" | cut -d- -f2)
    echo "Processing subject: $subid"

  # Loop over sessions for each subject
    for ses_path in ${subj_path}/ses-*; do
        sesid=$(basename "$ses_path" | cut -d- -f2)
        echo "  Session: $sesid"

        # Check if session folder is empty
        if [ -z "$(ls -A ${ses_path})" ]; then
            echo "  Warning: Session folder ${ses_path} is empty. Skipping."
            continue
        fi

        dwi_dir=${ses_path}/dwi
    
        # Check for dwi files
        if [ -z "$(ls -A ${dwi_dir})" ]; then
            echo "  Warning: No dwi data found in ${dwi_dir}. Skipping."
            continue
        fi

        # Determine ACQID (scanner type) - fixed per participant
        if ls ${dwi_dir}/sub-${subid}_ses-${sesid}_acq-ge_dir-AP_dwi.nii; then
            acqid="ge"
        elif ls ${dwi_dir}/sub-${subid}_ses-${sesid}_acq-prisma_dir-AP_dwi.nii; then
            acqid="prisma"
        else
            acqid="siemens"
        fi
    
        # Define input files directories 
        dwi_in=${dwi_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_dir-AP_dwi.nii
        bval_in=${dwi_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_dir-AP_dwi.bval
        bvec_in=${dwi_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_dir-AP_dwi.bvec
        dwi_pe_in=${dwi_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_dir-PA_dwi.nii
        
        # Define output directories
        result_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/dwi
        mkdir -p ${result_dir}

        anat_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/anat

        # Convert to mif format
        echo "    Converting to mif format ..."
        dwi=${result_dir}/dwi.mif
        mrconvert ${dwi_in} -fslgrad ${bvec_in} ${bval_in} ${dwi} -force


        # Denoising
        echo "     Denoising ..."
        dwidenoise ${dwi} -noise ${result_dir}/noise.mif ${result_dir}/dwi_den.mif -force

        # Gibbs Ringing Artefact Removal 
        echo "     Gibbs Ringing Artefact Removel ..."
            # For the main dwi data
        mrdegibbs ${result_dir}/dwi_den.mif ${result_dir}/dwi_den_unr.mif -force
            # For the pe dwi data 
        mrconvert ${dwi_pe_in} -coord 3 0 - | mrdegibbs - ${result_dir}/b0_pa_unr.mif -force

        # Motion and Distortion Correction 
        echo "     Motion and Distortion Correction ..."
            # pairing b0 images 
        dwiextract ${result_dir}/dwi_den_unr.mif -bzero - | mrconvert - -coord 3 0 ${result_dir}/b0_ap_unr.mif -force
        mrcat ${result_dir}/b0_ap_unr.mif ${result_dir}/b0_pa_unr.mif -axis 3 ${result_dir}/b0s_paired.mif -force 
            # Correction
        dwifslpreproc ${result_dir}/dwi_den_unr.mif ${result_dir}/dwi_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi ${result_dir}/b0s_paired.mif -eddy_options " --slm=linear" -force 

        # Bias Field Correction
        echo "     Bias Field Correction ..."
        dwibiascorrect ants ${result_dir}/dwi_den_unr_preproc.mif ${result_dir}/dwi_den_unr_preproc_bc.mif -bias ${result_dir}/bias.mif -force 

        # Coregister to T1 space 
        echo "     Coregister dwi to the T1 space ..."
            # Extract b=0 volumes and calculate the mean.
            # Export to NIFTI format for compatibility with FSL.
        dwiextract ${result_dir}/dwi_den_unr_preproc_bc.mif - -bzero | mrmath - mean ${result_dir}/mean_b0_preproc.nii -axis 3 -force 
            # Correct for bias field in the T1w image:
        T1_bc=$(ls ${anat_dir}/sub-*.nii)
        echo ${T1_bc}
            # Perform linear registration with 6 degrees of freedom:
        flirt -in ${result_dir}/mean_b0_preproc.nii -ref ${T1_bc} -dof 6 -cost normmi -omat ${result_dir}/diff2struct_fsl.mat
            # Convert the resulting linear transformation matrix from FSL to MRtrix format:
        transformconvert ${result_dir}/diff2struct_fsl.mat ${result_dir}/mean_b0_preproc.nii ${T1_bc} flirt_import ${result_dir}/diff2struct_mrtrix.txt -force
            # Apply linear transformation to header of diffusion-weighted image:
        mrtransform ${result_dir}/dwi_den_unr_preproc_bc.mif -linear ${result_dir}/diff2struct_mrtrix.txt ${result_dir}/dwi_den_unr_preproc_bc_coreg.mif -reorient_fod 0 -force

        # Brain Mask Estimation
        echo "     Brain Mask Estimation ..."
        dwi2mask ${result_dir}/dwi_den_unr_preproc_bc_coreg.mif ${result_dir}/dwi_mask.mif -force
  done
done

echo "End of Preprocessing of All Subjects"