#!/bin/bash
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Define your data and preprocessing directories
data_dir=/data/project/cosmonauts/Data/bids_ros_test
preproc_dir=/data/project/cosmonauts/Data/preprocessed

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

    anat_dir=${ses_path}/anat
    # Check for T1w files
    if [ -z "$(ls -A ${anat_dir})" ]; then
      echo "  Warning: No T1w data found in ${anat_dir}. Skipping."
      continue
    fi

    # Determine ACQID (scanner type) - fixed per participant
    if ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-ge_*_T1w.nii; then
      acqid="ge"
    else
      acqid="siemens"
    fi

    # Determine RUNID based on ACQID preference
    if [ ${acqid}=="ge" ]; then
      # For ge scanner, prefer dti over fmri
      if ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-ge_run-dti_T1w.nii; then
        runid="dti"
      elif ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-ge_run-fmri_T1w.nii; then
        runid="fmri"
      else
        echo "  Warning: No dti or fmri run found for ge in ${anat_dir}. Skipping."
        continue
      fi
    else
      # For siemens scanner, prefer highch over lowch
      if ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-siemens_run-highch_T1w.nii; then
        runid="highch"
      elif ls ${anat_dir}/sub-${subid}_ses-${sesid}_acq-siemens_run-lowch_T1w.nii; then
        runid="lowch"
      else
        echo "  Warning: No highch or lowch run found for siemens in ${anat_dir}. Skipping."
        continue
      fi
    fi

    infile=${anat_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w.nii
    echo "    Using file: ${infile}"

    # Define output directories
    result_dir=${preproc_dir}/sub-${subid}/ses-${sesid}
    bias_dir=${result_dir}/anat
    fs_dir=${result_dir}/anat/freesurfer
    mkdir -p ${bias_dir} ${fs_dir}

    # N4 bias field correction
    echo "    Running N4 bias correction..."
    N4BiasFieldCorrection -d 3 -i ${infile} -o ${bias_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.nii

    # FreeSurfer recon-all
    echo "    Running FreeSurfer recon-all..."
    fs_subj=sub-${subid}_ses-${sesid}
    recon-all -sd ${fs_dir} -subject ${fs_subj} -i ${bias_dir}/sub-${subid}_ses-${sesid}_acq-${acqid}_run-${runid}_T1w_bc.nii -all
    echo "    FreeSurfer output saved in ${fs_dir}/${fs_subj}"
  done
done

echo "End of Preprocessing of All Subjects"