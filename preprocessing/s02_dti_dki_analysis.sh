#! /bin/bash

# Subject label
subj_label=$1

# output structure
proj_dir=/autofs/space/allegro_003/users/ds1574/OGSE
# depth #1
bids_dir=${proj_dir}/bids
# nifti data dir
dwi_dir=${bids_dir}/${subj_label}/dwi
derivatives_dir=${bids_dir}/derivatives
preprocessed_dir=${derivatives_dir}/preprocessed/${subj_label}

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



if [ ! -d ${derivatives_dir}/mrtrix_tensor/${subj_label} ]; then
  mkdir -p ${derivatives_dir}/mrtrix_tensor/${subj_label}
else
  echo "Directory 'mrtrix_tensor' already exists."
fi
output_dir=${derivatives_dir}/mrtrix_tensor/${subj_label}

for acq_label in ${dwi_acq[@]:1}; do

  dwi_prefix=${subj_label}_acq-${acq_label}_dwi

  # Convert to mif
  bvec=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy.eddy_rotated_bvecs
  bval=${preprocessed_dir}/${dwi_prefix}_corr_bvals.txt
  dwi=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr.nii.gz
  mask=${preprocessed_dir}/${subj_label}_dwi_brain_mask_gncorr_C2.nii.gz
  mrconvert -fslgrad ${bvec} ${bval} ${dwi} ${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr.mif -force
  input=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr.mif
  

  # DKI
  /autofs/space/linen_001/users/Yixin/mrtrix3-dev/bin/dwi2tensor ${input} ${output_dir}/${dwi_prefix}_tensor.mif -dkt ${output_dir}/${dwi_prefix}_dkt.mif -constrain -force
  /autofs/space/linen_001/users/Yixin/mrtrix3-dev/bin/tensor2metric ${output_dir}/${dwi_prefix}_tensor.mif -dkt ${output_dir}/${dwi_prefix}_dkt.mif -mask ${mask} \
  -mk ${output_dir}/${dwi_prefix}_mk.nii.gz -ak ${output_dir}/${dwi_prefix}_ak.nii.gz -rk ${output_dir}/${dwi_prefix}_rk.nii.gz \
  -fa ${output_dir}/${dwi_prefix}_fa.nii.gz -vector ${output_dir}/${dwi_prefix}_fa_col.nii.gz \
  -ad ${output_dir}/${dwi_prefix}_ad.nii.gz -rd ${output_dir}/${dwi_prefix}_rd.nii.gz -adc ${output_dir}/${dwi_prefix}_adc.nii.gz \
  -cl ${output_dir}/${dwi_prefix}_cl.nii.gz -cp ${output_dir}/${dwi_prefix}_cp.nii.gz -cs ${output_dir}/${dwi_prefix}_cs.nii.gz -force  

  # Extract only low b-value data
  column_count=$(awk '{print NF; exit}' "$bval")
  if [ "$column_count" -eq 1 ]; then
      blow=$(sort -n "$bval" | uniq | head -2 | tail -1) # Adjust head -# to choose the desired b-value
      bhigh=$(sort -n "$bval" | uniq | tail -2 | tail -1) # Adjust head -# to choose the desired b-value
  else
      blow=$(tr ' ' '\n' < "$bval" | sort -n | uniq | head -2 | tail -1) # Adjust head -# to choose the desired b-value
      bhigh=$(tr ' ' '\n' < "$bval" | sort -n | uniq | tail -2 | tail -1) # Adjust head -# to choose the desired b-value
  fi
  echo 'blow='$blow
  echo 'bhigh='$bhigh
  input_lowb=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_lowb.mif
  input_highb=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_highb.mif
  dwiextract ${input} -shells 0,$blow ${input_lowb} -force
  dwiextract ${input} -shells 0,$bhigh ${input_highb} -force

  # DTI (low b shell)
  /autofs/space/linen_001/users/Yixin/mrtrix3-dev/bin/dwi2tensor ${input_lowb} ${output_dir}/${subj_label}_tensor_lowb.mif -force
  /autofs/space/linen_001/users/Yixin/mrtrix3-dev/bin/tensor2metric ${output_dir}/${subj_label}_tensor_lowb.mif -mask ${mask} \
  -fa ${output_dir}/${subj_label}_fa_lowb.nii.gz -vector ${output_dir}/${subj_label}_fa_col_lowb.nii.gz -ad ${output_dir}/${subj_label}_ad_lowb.nii.gz \
  -rd ${output_dir}/${subj_label}_rd_lowb.nii.gz -adc ${output_dir}/${subj_label}_adc_lowb.nii.gz \
  -cl ${output_dir}/${subj_label}_cl_lowb.nii.gz -cp ${output_dir}/${subj_label}_cp_lowb.nii.gz -cs ${output_dir}/${subj_label}_cs_lowb.nii.gz \
  -force  

  # DTI (high b shell)
  /autofs/space/linen_001/users/Yixin/mrtrix3-dev/bin/dwi2tensor ${input_highb} ${output_dir}/${subj_label}_tensor_highb.mif -force
  /autofs/space/linen_001/users/Yixin/mrtrix3-dev/bin/tensor2metric ${output_dir}/${subj_label}_tensor_highb.mif -mask ${mask} \
  -fa ${output_dir}/${subj_label}_fa_highb.nii.gz -vector ${output_dir}/${subj_label}_fa_col_highb.nii.gz -ad ${output_dir}/${subj_label}_ad_highb.nii.gz \
  -rd ${output_dir}/${subj_label}_rd_highb.nii.gz -adc ${output_dir}/${subj_label}_adc_highb.nii.gz \
  -cl ${output_dir}/${subj_label}_cl_highb.nii.gz -cp ${output_dir}/${subj_label}_cp_highb.nii.gz -cs ${output_dir}/${subj_label}_cs_highb.nii.gz \
  -force 

  # # DKI (single shell)
  # input_lowb_nii=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_lowb.nii.gz
  # input_highb_nii=${preprocessed_dir}/${dwi_prefix}_denoise_degibbs_eddy_gncorr_highb.nii.gz
  # lowb_bvec=${output_dir}/${dwi_prefix}_lowbdata.bvecs
  # highb_bvec=${output_dir}/${dwi_prefix}_highbdata.becs
  # lowb_bval=${output_dir}/${dwi_prefix}_lowbdata.bvals
  # highb_bval=${output_dir}/${dwi_prefix}_highbdata.bvals
  # dwiextract ${input} -shells 0,$blow ${input_lowb_nii} -fslgrad ${bvec} ${bval} -export_grad_fsl \
  #           ${lowb_bvec} ${lowb_bval} -force
  # dwiextract ${input} -shells 0,$bhigh ${input_highb_nii} -fslgrad ${bvec} ${bval} -export_grad_fsl \
  #           ${highb_bvec} ${highb_bval} -force

  # dtifit -k ${dwi} -o ${output_dir}/${dwi_prefix}_dkifit -m ${mask} -r ${bvec} -b ${bval} --sse --save_tensor -w --kurt
  # dtifit -k ${input_lowb_nii} -o ${output_dir}/${dwi_prefix}_dkifit_lowb -m ${mask} -r ${lowb_bvec} -b ${lowb_bval} --sse --save_tensor -w --kurt
  # dtifit -k ${input_highb_nii} -o ${output_dir}/${dwi_prefix}_dkifit_highb -m ${mask} -r ${highb_bvec} -b ${highb_bval} --sse --save_tensor -w --kurt

done