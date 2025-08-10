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
        
        # Define output directories
        result_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/anat
	
	mrconvert ${result_dir}/5tt.mif -coord 3 0 ${result_dir}/gm_cortical.nii -force
	mrconvert ${result_dir}/5tt.mif -coord 3 1 ${result_dir}/gm_subcortical.nii -force
	mrconvert ${result_dir}/5tt.mif -coord 3 2 ${result_dir}/wm.nii -force
	mrconvert ${result_dir}/5tt.mif -coord 3 3 ${result_dir}/csf.nii -force
	
	mrcalc ${result_dir}/gm_cortical.nii ${result_dir}/gm_subcortical.nii -add ${result_dir}/.nii -force       
  done
done

echo "The End!"
