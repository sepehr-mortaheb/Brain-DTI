#!/bin/bash
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Define your data and preprocessing directories
preproc_dir=/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/preprocessed/ROS

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
        result_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/dwi

        # Response Function Estimation
        echo "    Response Function Estimation ..."
        dwi2response dhollander ${result_dir}/dwi_den_unr_preproc_bc_coreg.mif ${result_dir}/wm.txt ${result_dir}/gm.txt ${result_dir}/csf.txt -voxels ${result_dir}/voxels.mif -force


        # FOD Estimation 
        echo "    FOD Estimation ..."
        dwi2fod msmt_csd  ${result_dir}/dwi_den_unr_preproc_bc_coreg.mif -mask ${result_dir}/dwi_mask.mif ${result_dir}/wm.txt ${result_dir}/wmfod.mif ${result_dir}/gm.txt ${result_dir}/gm.mif ${result_dir}/csf.txt ${result_dir}/csf.mif -force

        # Intensity Normalization
        echo "    Intensity Normalization ..."
        mtnormalise ${result_dir}/wmfod.mif ${result_dir}/wmfod_norm.mif ${result_dir}/gm.mif ${result_dir}/gm_norm.mif ${result_dir}/csf.mif ${result_dir}/csf_norm.mif -mask ${result_dir}/dwi_mask.mif -check_factors ${result_dir}/check_factors.txt -check_norm ${result_dir}/check_norm.mif -check_mask ${result_dir}/check_mask.mif   -force
            
  done
done

echo "End of FOD Estimation of All Subjects"