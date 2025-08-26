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
        fs_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/anat/freesurfer/sub-${subid}_ses-${sesid}
        anat_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/anat
        result_dir=${preproc_dir}/sub-${subid}/ses-${sesid}/dwi

        # 5ttgen Segmentation
        echo "   5ttgen Segmentation..."
        5ttgen hsvs ${fs_dir} ${anat_dir}/5tt.mif

        # Tractogram Generation 
        echo "   Tractogram Generation..."
        tckgen -algorithm ifod2 -act ${anat_dir}/5tt.mif -backtrack -seed_dynamic ${result_dir}/wmfod_norm.mif -select 10m ${result_dir}/wmfod_norm.mif ${result_dir}/tracks_10m.tck -force

        # Tractogram Filtering
        echo "   Tractogram Filtering..."
        tcksift2 -act ${anat_dir}/5tt.mif ${result_dir}/tracks_10m.tck ${result_dir}/wmfod_norm.mif ${result_dir}/sift2_weights.txt -out_mu ${result_dir}/sift2_mu.txt -force 
    done
done

echo "End of Tractogram Generation of All Subjects"