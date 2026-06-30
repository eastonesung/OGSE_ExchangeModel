#! /bin/bash

############ Run Yixin's dcm2bids first ##############
# source ~/init_conda
# cd /autofs/cluster/connectome2/Bay8_C2/bids/code/preprocessing_dwi
# python run_v1.py
# Save it under {project_dir}/bids

# Load Data
# Find dicom path using findsession
subj_list=('05-02' '06' '07' '08' '09' '10' '11' '12' '13' '14' '15')
scan_list=('011_2025_06_05' '06_2025_04_14' '07_2025_04_30' '08_2025_05_28' '09_2025_05_30' '10_2025_06_04')
# for subnum in ${subj_list[@]}; do
for si in {3..3}; do
subj_label=sub-ogse0${subj_list[si]}
scan_id=OGSE0${scan_list[si]}
prot='C2'
# subj_label=sub-ogse005-02
# scan_id=OGSE_PhantomTest_20251112_water
script_dir=/autofs/space/allegro_003/users/ds1574/OGSE/code/preprocessing
projectDir=/autofs/space/allegro_003/users/ds1574/OGSE
cd $script_dir

do_dcm2bids=0
do_fullpipeline=0
do_dtidki=0
do_addsynthseg=0
do_ants=0
do_reconall=1
do_amp2sh=0

############ Convert DCM to BIDS #################
if [ ${do_dcm2bids} = 1 ]; then
   ## get DICOM path from the Bourget
   dicomPath=$( findsession ${scan_id} command | grep -v 'siemens/MR-' | sed -n 's/.*PATH   :  //p' )
   echo "DICOM path    :   ${dicomPath}"
   if [ ! -d ${projectDir}/nii/${scan_id} ]; then
      mkdir -p ${projectDir}/nii/${scan_id}
   else
      echo "Directory 'nii' already exists."
      rm -rf ${projectDir}/nii/${scan_id}
      mkdir -p ${projectDir}/nii/${scan_id}
   fi

   # Convert DCM to Nifti
   /usr/pubsw/packages/dcm2niix/v1.0.20230411/dcm2niix -f %n_%d -t y -o ${projectDir}/nii/${scan_id} ${dicomPath}

#    # Reorganize Nifti files in BIDS format
#    matlab -nodesktop -nojvm -nodisplay -r "addpath('${script_dir}'),nii2bids_ogse('${projectDir}','${scan_id}','${subj_label}'),quit"
fi

if [ ${do_fullpipeline} = 1 ]; then
    set -e
    isDenoise=1
    isDegibbs=1
    isTopup=1
    isEddy=1
    isGradNonLinear=1
    isInterpEddyGNC=1
    isSynthSeg=1
    json_file=${projectDir}/bids/${subj_label}/dwi/${subj_label}_acq-30Hz_N1_dwi.json
    sh ${script_dir}/s01_run_preprocessing_ogse_C2.sh $subj_label $json_file $prot $isDenoise $isDegibbs $isTopup $isEddy $isGradNonLinear $isInterpEddyGNC $isSynthSeg
fi

################# DTI & DKI analysis ################
if [ ${do_dtidki} = 1 ]; then
    sh ${script_dir}/s02_dti_dki_analysis.sh $subj_label
fi


# ############## ANTS coregistration #######################
if [ ${do_ants} = 1 ]; then
    do_modelparam=1
    sh ${script_dir}/s04_run_antsRegistration.sh ${projectDir} ${subj_label} ${do_modelparam}
fi

if [ ${do_reconall} = 1 ]; then
    echo "Running recon-all"
    freesurfer_dir=${projectDir}/bids/derivatives/freesurfer
    if [ ! -d ${freesurfer_dir} ]; then
        mkdir -p ${freesurfer_dir}
    else
        echo "Directory 'freesurfer' already exists."
    fi
    input_nii=${projectDir}/bids/${subj_label}/anat/${scan_id}_t1_mprage_noellp_1mm_grappa2.nii.gz
    echo ${input_nii}
    recon-all -i ${input_nii} -s ${subj_label} -sd ${freesurfer_dir} -all
fi

if [ ${do_amp2sh} = 1 ]; then
    # output structure
    sh ${script_dir}/s05_run_amp2sh.sh ${projectDir} $subj_label

fi

############ Additional Synthseg ########
if [ ${do_addsynthseg} = 1 ]; then
    # Define the directory and prefix variables
    synthseg_dir=${projectDir}/bids/derivatives/synthseg/${subj_label}
    preprocessed_dir=${projectDir}/bids/derivatives/preprocessed/${subj_label}
    brainmask=${preprocessed_dir}/${subj_label}_dwi_brain_mask_gncorr_${prot}.nii.gz
    mkdir -p ${synthseg_dir}/all_labels
    # Loop through all the labels
    label_list=(7 8 10 11 12 13 17 18 26 46 47 49 50 51 52 53 54 58)
    for label in "${label_list[@]}"; do
        mri_extract_label "${synthseg_dir}/${subj_label}_acq-30Hz_N1_dwi_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz" "${label}" "${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}.nii.gz"
        mrgrid ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}.nii.gz regrid -interp nearest -template ${brainmask} \
                ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz -force
        fslmaths ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz -bin ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz
    done

    for label in {1001..1035}; do
        mri_extract_label "${synthseg_dir}/${subj_label}_acq-30Hz_N1_dwi_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz" "${label}" "${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}.nii.gz"
        mrgrid ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}.nii.gz regrid -interp nearest -template ${brainmask} \
                ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz -force
        fslmaths ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz -bin ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz
    done

    for label in {2001..2035}; do
        mri_extract_label "${synthseg_dir}/${subj_label}_acq-30Hz_N1_dwi_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz" "${label}" "${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}.nii.gz"
        mrgrid ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}.nii.gz regrid -interp nearest -template ${brainmask} \
                ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz -force
        fslmaths ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz -bin ${synthseg_dir}/all_labels/${subj_label}_synthseg_${label}_res.nii.gz
    done
fi

if [ ${do_jhuatlasreg} = 1 ]; then
    antsRegistrationSyN.sh -d 3 \
    -f subject_FA.nii.gz \
    -m $FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz \
    -o MNI2sub_ -t s
fi

done