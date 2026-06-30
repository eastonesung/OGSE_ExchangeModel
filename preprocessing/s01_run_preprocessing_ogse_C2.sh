#! #! /bin/bash
#
# This is a shell script DWI preprocessing
#
# Dependencies: (1)MRtrix (2)FSL (3)Freesurfer
#
# Creator: Kwok-Shing Chan @ MGH
# kchan2@mgh.harvard.edu
#
# Date created: 12 October 2023
# Date edit: 27 Oct 2024 by Dongsuk Sung (dsung2@mgh.harvard.edu)
############################################################

# Subject label
subj_label=$1
json_file=$2
prot=$3

isDenoise=$4
isDegibbs=$5
isTopup=$6
isEddy=$7
isGradNonLinear=$8
isInterpEddyGNC=$9
isSynthSeg=${10}


# output structure
proj_dir=/autofs/space/allegro_003/users/ds1574/OGSE
# depth #1
bids_dir=${proj_dir}/bids
derivatives_dir=${bids_dir}/derivatives
code_dir=${proj_dir}/code/

# nifti data dir
dwi_dir=${bids_dir}/${subj_label}/dwi

# dwi_acq=('PA' '90Hz_N3' '30Hz_N1' '60Hz_N2' '60Hz_N1' '40Hz_N1' '50Hz_N1')
# dwi_acq=('PA' '90Hz_N3' '30Hz_N1' '60Hz_N2' '60Hz_N1' '40Hz_N1' '50Hz_N1' '70Hz_N2' '80Hz_N2' '80Hz_N3')
###########################################################################
## Step 0: Acquire dwi_acq (prefix)
# Initialize arrays to store the values
combined=()
dwi_acq=('PA')

# Iterate through the files in the directory
for file in "${dwi_dir}"/*_acq-*Hz_N*_dwi.nii.gz; do
    # Extract the values using regular expressions
    if [[ $file =~ _acq-([0-9]{1,3})Hz_N([0-9]{1,3})_dwi.nii.gz ]]; then
        freq_value="${BASH_REMATCH[1]}"
        N_value="${BASH_REMATCH[2]}"
        combined+=("${freq_value}:${N_value}")
    fi
done

# Sort the combined array
IFS=$'\n' sorted=($(sort -t: -k1,1n <<<"${combined[*]}"))
unset IFS

# Separate the sorted values back into b_values and dir_values
for pair in "${sorted[@]}"; do
    IFS=':' read -r freq N <<< "$pair"
    dwi_acq+=("${freq}Hz_N${N}")
done

# Print the dwi_acq list
echo "dwi_acq: ${dwi_acq[@]}"
########################################################################

# Denoising using MRtrix
denoise_dir=${derivatives_dir}/denoise/${subj_label}

# FSL related
topup_dir=${derivatives_dir}/topup/${subj_label}
eddy_dir=${derivatives_dir}/eddy/${subj_label}

# Freesurfer related
synthstrip_dir=${derivatives_dir}/synthstrip/${subj_label}
synthseg_dir=${derivatives_dir}/synthseg/${subj_label}

# gradient nonlinearity correction
gradnonlin_dir=${derivatives_dir}/gnc/${subj_label}
applywarp_dir=${derivatives_dir}/applywarp/${subj_label}

# final preprocessed directory
preprocessed_dir=${derivatives_dir}/preprocessed/${subj_label}

mkdir -p ${denoise_dir}
mkdir -p ${topup_dir}
mkdir -p ${eddy_dir}
mkdir -p ${synthstrip_dir}
mkdir -p ${synthseg_dir}
mkdir -p ${gnc_dir}
mkdir -p ${applywarp_dir}
mkdir -p ${preprocessed_dir}

# add tools
script_dir=/autofs/space/virtuoso_002/users/ds1574/tools/dwi/utils/
export PATH=${PATH}:${script_dir}

## static output filename
acq_param_txt=${topup_dir}/acquisition_parameters_${prot}.txt
topup_output_prefix=${topup_dir}/dwi_topup_res_${prot}
topup_hifi_nii=${topup_dir}/dwi_hifi_b0_${prot}.nii.gz
topup_hifi_Tmean_nii=${topup_dir}/dwi_hifi_b0_${prot}_Tmean.nii.gz
dwi_brain_nii=${synthstrip_dir}/${subj_label}_dwi_${prot}_brain.nii.gz
dwi_brain_mask_nii=${synthstrip_dir}/${subj_label}_dwi_${prot}_brain_mask.nii.gz
topup_input_nii=${topup_dir}/dwi_${prot}_b0_all.nii.gz
eddy_template=${eddy_dir}/target_dwi.nii.gz

set -e

log_file=${dwi_dir}/process_history.log

########################################
## Step 1: Rician floor corrected MPPCA

if [ ${isDenoise} = 1 ]; then

    echo "Step 1: Performing Rician-floor-corrected MPPCA..."

    for acq in ${dwi_acq[@]}; do

        dwi_prefix=${subj_label}_acq-${acq}_dwi
        nthreads=16

        input_nii=${dwi_dir}/${dwi_prefix}.nii.gz
        output_nii=${denoise_dir}/${dwi_prefix}_denoise.nii.gz
        sigma_nii=${denoise_dir}/${dwi_prefix}_denoise_sigma.nii.gz

        echo "Input data: ${input_nii}" | tee -a "$log_file"
        echo "Output data: ${output_nii}" | tee -a "$log_file"
        echo "Noise data: ${sigma_nii}" | tee -a "$log_file"
        # util_dir=/autofs/space/allegro_003/users/ds1574/scripts
        # source $util_dir/rician_correct_mppca.sh -i ${input_nii} -o ${output_nii} -n 10 -c 16
        # sh rician_correct_mppca.sh -i ${input_nii} -o ${output_nii} -n 10 -c ${nthreads}

        # pre-iteration processing
        dwidenoise -noise ${sigma_nii} ${input_nii} ${output_nii} -nthreads ${nthreads} -force

        # # For each b-shell
        # # Extract low and high b-value data
        # if [ "$acq" = ${dwi_acq[0]} ]; then
        #     echo "PA doesn't have bval"
        # else
        #     lb=1000
        #     hb=2000
        #     bval=${dwi_dir}/${dwi_prefix}.bval
        #     bvec=${dwi_dir}/${dwi_prefix}.bvec
        #     dwi_lb=${denoise_dir}/${dwi_prefix}_b${lb}.nii.gz
        #     dwi_hb=${denoise_dir}/${dwi_prefix}_b${hb}.nii.gz
        #     dwiextract ${input_nii} -shells $lb -fslgrad $bvec $bval ${dwi_lb} -force
        #     dwiextract ${input_nii} -shells $hb -fslgrad $bvec $bval ${dwi_hb} -force
        #     output_lb=${denoise_dir}/${dwi_prefix}_denoise_b${lb}.nii.gz
        #     output_hb=${denoise_dir}/${dwi_prefix}_denoise_b${hb}.nii.gz
        #     sigma_lb=${denoise_dir}/${dwi_prefix}_denoise_b${lb}_sigma.nii.gz
        #     sigma_hb=${denoise_dir}/${dwi_prefix}_denoise_b${hb}_sigma.nii.gz
        #     sh rician_correct_mppca.sh -i ${dwi_lb} -o ${output_lb} -n 10 -c 16
        #     sh rician_correct_mppca.sh -i ${dwi_hb} -o ${output_hb} -n 10 -c 16
        #     # dwidenoise -noise ${sigma_lb} ${dwi_lb} ${output_lb} -nthreads 16 -force
        #     # dwidenoise -noise ${sigma_hb} ${dwi_hb} ${output_hb} -nthreads 16 -force
        # fi
    done

    echo "Step 1 is completed."
else
    echo -e "\e[1;32mWarning! Step 1. DENOISE step skipped!\e[0m"
fi


############################################
## Step 2: Gibbs rings removal on b0 image

if [ ${isDegibbs} = 1 ]; then
    echo "Step 2: Performing Gibbs ringing removal"

    for acq in ${dwi_acq[@]}; do    
        echo $acq
        dwi_prefix=${subj_label}_acq-${acq}_dwi

        input_nii=${denoise_dir}/${dwi_prefix}_denoise.nii.gz
        output_nii=${denoise_dir}/${dwi_prefix}_denoise_degibbs.nii.gz

        # remove NaN entries before mrdegibbs
        fslmaths ${input_nii} -nan ${input_nii}

        # degibbs
        mrdegibbs ${input_nii} ${output_nii} -force
    done
    echo "Step 2 is completed." | tee -a "$log_file"
else
    echo -e "\e[1;32mWarning! Step 2. DEGIBBS step skipped!\e[0m"
fi

############################################
# ## Step 3: topup
# number of volumes per dataset for topup
Nvols_topup=1 # Nvols=$(fslnvols "${input_PA_nii}")


if [ ${isTopup} = 1 ]; then

    echo "Step 3: Performing topup to correct susceptibility distortion..."

    # Extract b=0 from the OG waveform acquired closes to b0PA (60Hz-N2)
    acq_label=${dwi_acq[5]}
    dwi_prefix=${subj_label}_acq-${acq_label}_dwi
    echo "dwi_prefix: "${dwi_prefix} | tee -a "$log_file"
    input_nii=${denoise_dir}/${dwi_prefix}_denoise_degibbs.nii.gz
    output_nii=${denoise_dir}/${dwi_prefix}_denoise_b0_degibbs.nii.gz
    bval_txt=${dwi_dir}/${dwi_prefix}.bval
    bvec_txt=${dwi_dir}/${dwi_prefix}.bvec

    echo "Extracting b=0 images from" ${input_nii} | tee -a "$log_file"
    dwiextract -bzero -fslgrad ${bvec_txt} ${bval_txt} ${input_nii} ${output_nii} -force 


    # use the two successive acquisitions
    input_PA_nii=${denoise_dir}/${subj_label}_acq-${dwi_acq[0]}_dwi_denoise_degibbs.nii.gz
    input_AP_nii=${denoise_dir}/${subj_label}_acq-${dwi_acq[5]}_dwi_denoise_b0_degibbs.nii.gz

    ## 3.1 combine PA-AP volumes
    # only use a few good quality b0 volumes
    echo "Combine bliped-up bliped-down images"
    echo "Combine bliped-up bliped-down images" | tee -a "$log_file"
    echo "inputPA: ${input_PA_nii}, inputAP: ${input_AP_nii}, topup_input: ${topup_input_nii}, Nvols_topup:${Nvols_topup}" | tee -a "$log_file"
    script_path=/autofs/space/allegro_003/users/ds1574/scripts
    # matlab -nodesktop -nojvm -nodisplay -r "addpath('${script_path}'),concat_matched_blippedupdown_SSIM_DS('${input_PA_nii}','${input_AP_nii}','${topup_input_nii}',${Nvols_topup}),quit"
    matlab -nodesktop -nojvm -nodisplay -r "addpath('${script_dir}'),concat_matched_blippedupdown_SSIM('${input_PA_nii}','${input_AP_nii}','${topup_input_nii}',${Nvols_topup}),quit"

    ## 3.2 create acquisition_parameter.txt
    # Setting acquisition parameter for topup, may need manual operation
    # 1/-1: PE direction; total_readout_time: time you have had collected all k-space lines
    # 1st dim: LR; 2nd dim: AP; 3rd dim: HF
    # Total_readout=EffectiveEchoSpacing*(ReconMatrixPE-1), eddy has a hard threshold of readout time <=200ms
    # total_readout_time=0.015190049   # 0.00021*(110-1)
    total_readout_time=`grep "TotalReadoutTime" $json_file | awk '{print substr($2, 1, length($2)-1)}'`
    echo "Total Readout Time = "$total_readout_time

    echo "3.2. create acquisition_parameter.txt" | tee -a "$log_file"
    echo -n "" > ${acq_param_txt}
    for ((n=1;n<=$Nvols_topup;n++)); do echo "0 -1 0 $total_readout_time" >> ${acq_param_txt}; done
    for ((n=1;n<=$Nvols_topup;n++)); do echo "0 1 0 $total_readout_time" >> ${acq_param_txt}; done

    echo "Performing topup..."
    echo "Performing topup..." | tee -a "$log_file"
    ## 3.3 topup main
    # if it's a small FOV phantom data change --config=b02b0_1.cnf instead of --config=b02b0.cnf
    topup --verbose \
        --imain=${topup_input_nii} \
        --datain=${acq_param_txt} \
        --config=b02b0_1.cnf \
        --out=${topup_output_prefix} \
        --iout=${topup_hifi_nii} 

    ## 3.4 Brain extraction
    echo "Extracting brain mask on ${topup_hifi_Tmean_nii}"
    echo "Extracting brain mask on ${topup_hifi_Tmean_nii}" | tee -a "$log_file"
    fslmaths ${topup_hifi_nii} -Tmean ${topup_hifi_Tmean_nii}
    mri_synthstrip -i ${topup_hifi_Tmean_nii} -m ${dwi_brain_mask_nii} -o ${dwi_brain_nii}

    echo "Step 3 is completed." | tee -a "$log_file"
else
    echo -e "\e[1;32mWarning! Step 3. TOPUP step skipped!\e[0m"
fi

############################################
## Step 4: eddy

if [ ${isEddy} = 1 ]; then

    echo "Step 4: Performing eddy current distortion correction..."

    ## get target volume for eddy from AP
    input_nii=${denoise_dir}/${subj_label}_acq-${dwi_acq[1]}_dwi_denoise_degibbs.nii.gz
    fslroi ${input_nii} ${eddy_template} 0 1

    bet ${topup_hifi_Tmean_nii} ${dwi_brain_nii} -f 0.5 -m
    fslmaths ${dwi_brain_mask_nii} -dilM -dilM -bin ${eddy_dir}/brain_mask_4eddy.nii.gz
    
    echo "Done brain masking..."

    ## prepare for eddy
    for acq in ${dwi_acq[@]:1}; do

        acq_label=${acq}
        dwi_prefix=${subj_label}_acq-${acq_label}_dwi

        # write slspec text file and export a clean bvals text file
        json_txt=${dwi_dir}/${dwi_prefix}.json
        slspec_txt=${eddy_dir}/${dwi_prefix}_slspec.txt
        bval=${dwi_dir}/${dwi_prefix}.bval
        output_prefix=${eddy_dir}/${dwi_prefix}_corr

        echo "Export slspec file from JSON and clean up bval..."

        matlab -nodesktop -nojvm -nodisplay -nosplash -r "addpath('${script_dir}'),export_slspec_from_json('${json_txt}','${slspec_txt}'),cleanup_bvals_v2('${bval}','${output_prefix}'),quit"

        echo "Save slspec file as ${slspec_txt}"
        echo "Save cleaned up bval file as ${eddy_dir}/${dwi_prefix}_corr_bvals.txt"

        # concat NIFTI data for eddy
        echo "Concatenating NIFTI data with template for eddy..."
        input_nii=${dwi_dir}/${dwi_prefix}.nii.gz
        input_tmp_nii=${eddy_dir}/${dwi_prefix}_4eddy.nii.gz
        fslmerge -t ${input_tmp_nii} ${eddy_template} ${input_nii}

        # input_tmp_nii=${eddy_dir}sub-NEXIC2HC011_C2D13ms_dwi_denoise_degibbs_eddy_tmp.eddy_outlier_free_data.nii.gz

        # concat bval and bvec
        echo "Concatenating bval and bvec for eddy..."
        sed '1,3s/^/0 /' ${dwi_dir}/${dwi_prefix}.bvec > ${eddy_dir}/${dwi_prefix}_bvecs_tmp
        sed '1s/^/0 /' ${eddy_dir}/${dwi_prefix}_corr_bvals.txt > ${eddy_dir}/${dwi_prefix}_bvals_tmp

        nb=$(fslnvols ${input_tmp_nii})
        for ((n=1; n<=$nb; n+=1)); do indx="$indx $((${Nvols_topup}+1))"; done
        index_txt=${eddy_dir}/${dwi_prefix}_index.txt
        echo $indx > ${index_txt}
        unset indx

        # *******masked out non-brain tissue to mitigate fat rings*******
        fslmaths ${input_tmp_nii} -mas ${eddy_dir}/brain_mask_4eddy ${input_tmp_nii}

        # eddy main
        echo "Performing eddy..."
        input_nii=${input_tmp_nii}
        bvec=${eddy_dir}/${dwi_prefix}_bvecs_tmp
        bval=${eddy_dir}/${dwi_prefix}_bvals_tmp
        output=${eddy_dir}/${dwi_prefix}_denoise_degibbs_eddy_tmp
        echo "Input NIFTI data  : ${input_nii}"
        echo "Input bvec        : ${bvec}"
        echo "Input bval        : ${bval}"
        echo "Output data       : ${output}"
        eddy_cuda10.2 --imain=${input_nii} \
            --mask=${dwi_brain_mask_nii} \
            --acqp=${acq_param_txt} \
            --index=${index_txt} \
            --bvecs=${bvec} \
            --bvals=${bval} \
            --topup=${topup_output_prefix} \
            --out=${output} \
            --niter=8 \
            --repol \
            --fwhm=10,0,0,0,0,0,0,0 \
            --data_is_shelled \
            --flm=cubic \
            --mporder=6 \
            --slspec=${slspec_txt} \
            --s2v_niter=5 \
            --s2v_lambda=1 \
            --s2v_interp=trilinear \
            --cnr_maps --residuals \
            --verbose 

        echo "Applying eddy..."
        input_nii=${denoise_dir}/${dwi_prefix}_denoise_degibbs.nii.gz
        input_tmp_nii=${eddy_dir}/${dwi_prefix}_4eddy_add.nii.gz
        fslmerge -t ${input_tmp_nii} ${eddy_template} ${input_nii}
        input_nii=${input_tmp_nii}
        eddy_cuda10.2 --imain=${input_nii} \
            --mask=${dwi_brain_mask_nii} \
            --acqp=${acq_param_txt} \
            --index=${index_txt} \
            --bvecs=${bvec} \
            --bvals=${bval} \
            --topup=${topup_output_prefix} \
            --out=${output} \
            --niter=0 \
            --repol \
            --data_is_shelled \
            --flm=cubic \
            --mporder=6 \
            --slspec=${slspec_txt} \
            --s2v_niter=0 \
            --s2v_lambda=1 \
            --s2v_interp=trilinear \
            --init=${output}.eddy_parameters \
            --init_s2v=${output}.eddy_movement_over_time \
            --verbose --dfields

        echo "Removing template from temporary eddy result..."
        fslroi ${output} ${eddy_dir}/${dwi_prefix}_denoise_degibbs_eddy 1 $((${nb}-1))
        sed 's/^...//' ${output}.eddy_rotated_bvecs > ${eddy_dir}/${dwi_prefix}_denoise_degibbs_eddy.eddy_rotated_bvecs
        rm ${output}.nii.gz \
            ${input_tmp_nii} \
            ${eddy_dir}/${dwi_prefix}_bvecs_tmp \
            ${eddy_dir}/${dwi_prefix}_bvals_tmp

        echo "The final output of eddy is ${eddy_dir}/${dwi_prefix}_denoise_degibbs_eddy.nii.gz"

    done

    echo "Step 4 is completed."
else
    echo -e "\e[1;32mWarning! Step 4. EDDY step skipped!\e[0m"
fi

############################################
## Step 5: gradient nonlinearity correction

if [ ${isGradNonLinear} = 1 ]; then

    echo "Step 5: Performing gradient non-linearity correction..."
    gnc_script_dir='/space/scheherazade/2/users/qfan/tools/preproc/hcps_diff_prep_v2'
    for acq in ${dwi_acq[@]:1}; do
        dwi_prefix=${subj_label}_acq-${acq}_dwi

        input_nii=${eddy_dir}/${dwi_prefix}_denoise_degibbs_eddy.nii.gz
        # input_nii=${eddy_dir}/${dwi_prefix}_4eddy.nii.gz
        output_nii=${gradnonlin_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr.nii.gz
        ${gnc_script_dir}/hcps_gnc_c2.sh -i ${input_nii} -o ${output_nii} -interp spline

        input_nii=${denoise_dir}/${dwi_prefix}_denoise_sigma.nii.gz
        output_nii=${gradnonlin_dir}/${dwi_prefix}_denoise_sigma_gncorr.nii.gz
        ${gnc_script_dir}/hcps_gnc_c2.sh -i ${input_nii} -o ${output_nii} -interp spline

        # ### GNC for b=1000 and b=2000 separately ###
        # lb=1000
        # hb=2000
        # sigma_lb=${denoise_dir}/${dwi_prefix}_denoise_b${lb}_sigma.nii.gz
        # sigma_hb=${denoise_dir}/${dwi_prefix}_denoise_b${hb}_sigma.nii.gz
        # output_lb=${gradnonlin_dir}/${dwi_prefix}_denoise_b${lb}_sigma_gncorr.nii.gz
        # output_hb=${gradnonlin_dir}/${dwi_prefix}_denoise_b${hb}_sigma_gncorr.nii.gz
        # /space/scheherazade/2/users/qfan/tools/preproc/hcps_diff_prep_v2/hcps_gnc_c2.sh -i ${sigma_lb} \
        #     -o ${output_lb} -interp trilinear
        # /space/scheherazade/2/users/qfan/tools/preproc/hcps_diff_prep_v2/hcps_gnc_c2.sh -i ${sigma_hb} \
        #     -o ${output_hb} -interp trilinear

    done
    input_nii=${dwi_brain_mask_nii}
    output_nii=${gradnonlin_dir}/${subj_label}_dwi_brain_mask_gncorr.nii.gz
    /space/scheherazade/2/users/qfan/tools/preproc/hcps_diff_prep_v2/hcps_gnc_c2.sh -i ${input_nii} \
        -o ${output_nii} -interp spline

    echo "Step 5 is completed."
else
    echo -e "\e[1;32mWarning! Step 5. GNC step skipped!\e[0m"
fi

############################################
## Step 6: Interpolate results from EDDY and GNC
if [ ${isInterpEddyGNC} = 1 ]; then

    echo "Step 6: Apply warp combining EDDY and GNC..."
    util_dir=/autofs/space/virtuoso_002/users/ds1574/tools/dwi/utils
    cd $util_dir
    ## prepare for eddy
    for acq in ${dwi_acq[@]:1}; do
        dwi_prefix=${subj_label}_acq-${acq}_dwi
        echo ${dwi_prefix}
        python3 interpolation_eddy_gnc_helper.py "${subj_label}" "${derivatives_dir}" "${dwi_prefix}"
        # echo "Removing template from temporary eddy result..."
        # tmp_nii=${applywarp_dir}/${dwi_prefix}_denoise_degibbs_eddy_gnc_warp_tmp.nii.gz
        # nb=$(fslnvols ${tmp_nii})
        # fslroi ${tmp_nii} ${applywarp_dir}/${dwi_prefix}_denoise_degibbs_eddy_gnc_warp.nii.gz 1 $((${nb}-1))
    done

    
    cd '/autofs/space/allegro_003/users/ds1574/OGSE/code/preprocessing'
else
    echo -e "\e[1;32mWarning! Step 6. APPLYWARP step skipped!\e[0m"
fi

############################################
## Step 7: Quick segmentation on DWI for GM/WM masks

if [ ${isSynthSeg} = 1 ]; then

    echo "Step 7: Simple tissue classification on DWI..."
    mrtrix_bin=/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/mrtrix/3.0.3/bin

    dwi_prefix=${subj_label}_acq-${dwi_acq[1]}_dwi

    ## extract b=0
    bvec_txt=${dwi_dir}/${dwi_prefix}.bvec
    bval_txt=${eddy_dir}/${dwi_prefix}_corr_bvals.txt
    input_nii=${applywarp_dir}/${dwi_prefix}_denoise_degibbs_eddy_gnc_warp.nii.gz
    # input_nii=${gradnonlin_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr.nii.gz
    output_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_b0.nii.gz
    dwiextract -bzero -fslgrad ${bvec_txt} ${bval_txt} ${input_nii} ${output_nii} -force

    ## compute Tmean
    input_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_b0.nii.gz
    output_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_b0_Tmean.nii.gz
    fslmaths ${input_nii} -Tmean ${output_nii} 
    rm ${input_nii}

    ## acquire new brain mask
    b0_Tmean_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_b0_Tmean.nii.gz
    dwi_brain_mask=${synthseg_dir}/${subj_label}_dwi_${prot}_brain_mask_gncorr.nii.gz
    dwi_brain=${synthseg_dir}/${subj_label}_dwi_${prot}_brain_gncorr.nii.gz
    mri_synthstrip -i ${b0_Tmean_nii} -m ${dwi_brain_mask} -o ${dwi_brain}

    ## perform synthseg
    input_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_b0_Tmean.nii.gz
    output_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz
    resample_nii=${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_SynthSeg_Resample.nii.gz
    mri_synthseg --i ${input_nii} --o ${output_nii} --robust --parc --resample ${resample_nii}

    # Extract WM mask 
    mri_extract_label ${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz 2 41 ${synthseg_dir}/${subj_label}_wm_mask.nii.gz
    mrgrid ${synthseg_dir}/${subj_label}_wm_mask.nii.gz regrid -interp nearest -template ${dwi_brain_mask} \
            ${synthseg_dir}/${subj_label}_wm_mask_res.nii.gz -force
    fslmaths ${synthseg_dir}/${subj_label}_wm_mask_res.nii.gz -bin ${synthseg_dir}/${subj_label}_wm_mask_res.nii.gz

    ## Extract cortex mask
    fslmaths ${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz -thr 100 ${synthseg_dir}/${subj_label}_cortex_mask.nii.gz
    mrgrid ${synthseg_dir}/${subj_label}_cortex_mask.nii.gz regrid -interp nearest -template ${dwi_brain_mask} \
             ${synthseg_dir}/${subj_label}_cortex_mask_res.nii.gz -force
    fslmaths ${synthseg_dir}/${subj_label}_cortex_mask_res.nii.gz -bin ${synthseg_dir}/${subj_label}_cortex_mask_res.nii.gz

    ## Sub-Nuclei
    mri_extract_label ${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz 10 11 12 13 49 50 51 52 ${synthseg_dir}/${subj_label}_subnuclei_mask.nii.gz
    mrgrid ${synthseg_dir}/${subj_label}_subnuclei_mask.nii.gz regrid -interp nearest -template ${dwi_brain_mask} \
            ${synthseg_dir}/${subj_label}_subnuclei_mask_res.nii.gz -force
    fslmaths ${synthseg_dir}/${subj_label}_subnuclei_mask_res.nii.gz -bin ${synthseg_dir}/${subj_label}_subnuclei_mask_res.nii.gz

    ## CSF mask
    mri_extract_label ${synthseg_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz 4 5 14 15 43 44 24 ${synthseg_dir}/${subj_label}_csf_mask.nii.gz
    mrgrid ${synthseg_dir}/${subj_label}_csf_mask.nii.gz regrid -interp nearest -template ${dwi_brain_mask} \
            ${synthseg_dir}/${subj_label}_csf_mask_res.nii.gz -force
    fslmaths ${synthseg_dir}/${subj_label}_csf_mask_res.nii.gz -bin ${synthseg_dir}/${subj_label}_csf_mask_res.nii.gz

    echo "Step 7 is completed."
else
    echo -e "\e[1;32mWarning! Step 7. SYNTHSEG step skipped!\e[0m"
fi

############################################
## Final step
echo "Copying all essential files to ${preprocessed_dir}..."

for acq in ${dwi_acq[@]:1}; do

    acq_label=${acq}
    dwi_prefix=${subj_label}_acq-${acq_label}_dwi

    cp ${applywarp_dir}/${dwi_prefix}_denoise_degibbs_eddy_gnc_warp.nii.gz ${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gnc_warp.nii.gz
    cp ${gradnonlin_dir}/${dwi_prefix}_denoise_sigma_gncorr.nii.gz ${preprocessed_dir}/${dwi_prefix}_denoise_sigma_gncorr.nii.gz
    cp ${eddy_dir}/${dwi_prefix}_corr_bvals.txt ${preprocessed_dir}/${dwi_prefix}_corr_bvals.txt
    cp ${eddy_dir}/${dwi_prefix}_denoise_degibbs_eddy.eddy_rotated_bvecs ${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy.eddy_rotated_bvecs
    cp ${dwi_dir}/${dwi_prefix}.bvec ${preprocessed_dir}/${dwi_prefix}.bvec
done
cp ${synthseg_dir}/${subj_label}_dwi_${prot}_brain_mask_gncorr.nii.gz ${preprocessed_dir}/${subj_label}_dwi_brain_mask_gncorr_${prot}.nii.gz
cp ${synthseg_dir}/${subj_label}_acq-${dwi_acq[1]}_dwi_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz ${preprocessed_dir}/${subj_label}_acq-${dwi_acq[1]}_dwi_denoise_degibbs_eddy_gncorr_SynthSeg.nii.gz

echo "The preprocessed DWI data is now available in ${preprocessed_dir}"