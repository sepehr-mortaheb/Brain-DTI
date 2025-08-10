#!/bin/bash
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Define your data and preprocessing directories
preproc_dir=/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/preprocessed/ROS
atlas_dir=/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/MNI

# Loop over subjects
for subj_path in ${preproc_dir}/sub-*; do
    subid=$(basename "$subj_path" | cut -d- -f2)
    echo "Processing subject: $subid"

  # Loop over sessions for each subject
    for ses_path in ${subj_path}/ses-*; do
        sesid=$(basename "$ses_path" | cut -d- -f2)
        echo "  Session: $sesid"
    
        # Check if dwi folder exist
        if [ ! -d "${ses_path}/dwi" ]; then
            echo "  Warning: No dwi data found in ${ses_path}. Skipping."
            continue
        fi

        # Define output directories
        anat_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/anat
        result_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/dwi

        # Determine ACQID (scanner type) - fixed per participant
        if ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-ge_*_T1w_bc.nii; then
            echo " GE Scanner Found!"
            continue
        fi
         #   acqid="ge"
        #else
            acqid="prisma"
        #fi

        # Determine RUNID based on ACQID preference
        #if [ ${acqid}=="ge" ]; then
            #  For ge scanner, prefer dti over fmri
         #   if ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-ge_run-dti_T1w_bc.nii; then
           #     runid="dti"
          #  elif ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-ge_run-fmri_T1w_bc.nii; then
            #    runid="fmri"
           # else
            #    echo "  Warning: No dti or fmri run found for ge in ${anat_dir}. Skipping."
             #   continue
            #fi
       # else
            # For siemens scanner, prefer highch over lowch
            if ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-prisma_run-highch_T1w_bc.nii; then
                runid="highch"
            elif ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-prisma_run-lowch_T1w_bc.nii; then
                runid="lowch"
            else
                echo "  Warning: No highch or lowch run found for siemens in ${anat_dir}. Skipping."
                continue
            fi
        #fi

        infile_anat=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.nii

        # Anatomical FSL Preprocessing 
        echo "   Anatomical FSL Preprocessing..."
        fsl_anat --noseg --nosubcortseg -i ${infile_anat}

        # Registering Atlases to the subject space 
        echo "   Registering SCH100 Atlas to the subject space..."
        applywarp --ref=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.anat/T1_biascorr_brain.nii.gz --in=${atlas_dir}/SCH100.nii --warp=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.anat/MNI_to_T1_nonlin_field.nii.gz --out=${anat_dir}/SCH100_subj.nii --interp=nn 
        echo "   Registering SCH400 Atlas to the subject space..."
        applywarp --ref=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.anat/T1_biascorr_brain.nii.gz --in=${atlas_dir}/SCH400.nii --warp=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.anat/MNI_to_T1_nonlin_field.nii.gz --out=${anat_dir}/SCH400_subj.nii --interp=nn 
        echo "   Registering AAL Atlas to the subject space..."
        applywarp --ref=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.anat/T1_biascorr_brain.nii.gz --in=${atlas_dir}/AAL.nii --warp=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.anat/MNI_to_T1_nonlin_field.nii.gz --out=${anat_dir}/AAL_subj.nii --interp=nn 


        # SC Estimation
        echo "   SC Estimation based on SCH100 Atlas..."
        tck2connectome -tck_weights_in ${result_dir}/sift2_weights.txt -symmetric -zero_diagonal -scale_invnodevol ${result_dir}/tracks_10m.tck ${anat_dir}/SCH100_subj.nii.gz ${result_dir}/SC_SCH100.csv
        echo "   SC Estimation based on SCH400 Atlas..."
        tck2connectome -tck_weights_in ${result_dir}/sift2_weights.txt -symmetric -zero_diagonal -scale_invnodevol ${result_dir}/tracks_10m.tck ${anat_dir}/SCH400_subj.nii.gz ${result_dir}/SC_SCH400.csv
        echo "   SC Estimation based on AAL Atlas..."
        tck2connectome -tck_weights_in ${result_dir}/sift2_weights.txt -symmetric -zero_diagonal -scale_invnodevol ${result_dir}/tracks_10m.tck ${anat_dir}/AAL_subj.nii.gz ${result_dir}/SC_AAL.csv
    done
done

echo "End of Tractogram Generation of All Subjects"
