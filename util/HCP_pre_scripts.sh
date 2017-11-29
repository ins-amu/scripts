#!/usr/bin/env bash

# check mandatory variables
if [ -z "$PRD" ]; then
  echo "PRD path missing"
  exit 1
fi

if [ -z "$SUBJ_ID" ]; then
  echo "SUBJ_ID path missing"
  exit 1
fi

# unzip data dowloaded from HCP connectomeDB
unzip "$PRD"/"$SUBJ_ID"_3T_Diffusion_preproc.zip -d "$PRD"
unzip "$PRD"/"$SUBJ_ID"_3T_Structural_preproc.zip -d "$PRD"
unzip "$PRD"/"$SUBJ_ID"_3T_Structural_preproc_extended.zip -d "$PRD"

## prepare surface files
mkdir -p "$PRD"/data/T1 "$PRD"/surface
mkdir -p "$FS"/"$SUBJ_ID"

touch "$PRD"/data/T1/T1.nii.gz

mv "$PRD"/"$SUBJ_ID"/T1w/"$SUBJ_ID"/* "$FS"/"$SUBJ_ID"/

##  prepare connectivity files

mkdir -p "$PRD"/connectivity

# TODO: FS files or HCP files?

cp "$PRD"/"$SUBJ_ID"/T1w/aparc+aseg.nii.gz "$PRD"/connectivity/aparc+aseg.nii.gz
cp "$PRD"/"$SUBJ_ID"/T1w/aparc+aseg.nii.gz "$PRD"/connectivity/aparc+aseg_reorient.nii.gz
cp "$PRD"/"$SUBJ_ID"/T1w/aparc+aseg.nii.gz "$PRD"/connectivity/aparcaseg_2_diff.nii.gz
cp "$PRD"/"$SUBJ_ID"/T1w/T1w_acpc_dc_restore_brain.nii.gz "$PRD"/connectivity/brain.nii.gz
cp "$PRD"/"$SUBJ_ID"/T1w/T1w_acpc_dc_restore_brain.nii.gz "$PRD"/connectivity/brain_2_diff.nii.gz
cp "$PRD"/"$SUBJ_ID"/T1w/Diffusion/bvecs "$PRD"/connectivity/bvecs
cp "$PRD"/"$SUBJ_ID"/T1w/Diffusion/bvals "$PRD"/connectivity/bvals

# Crop HCP data to avoid RAM issues
gunzip "$PRD"/"$SUBJ_ID"/T1w/Diffusion/data.nii.gz
mrcrop "$PRD"/"$SUBJ_ID"/T1w/Diffusion/data.nii \
       "$PRD"/connectivity/data_crop.nii.gz \
       -mask "$PRD"/"$SUBJ_ID"/T1w/Diffusion/nodif_brain_mask.nii.gz -force
mrcrop "$PRD"/"$SUBJ_ID"/T1w/Diffusion/nodif_brain_mask.nii.gz \
       "$PRD"/connectivity/nodif_brain_mask_crop.nii.gz \
       -mask "$PRD"/"$SUBJ_ID"/T1w/Diffusion/nodif_brain_mask.nii.gz -force


mrconvert "$PRD"/connectivity/data_crop.nii.gz \
          "$PRD"/connectivity/predwi.mif \
          -fslgrad "$PRD"/connectivity/bvecs "$PRD"/connectivity/bvals \
          -datatype float32 -force
mrconvert "$PRD"/connectivity/data_crop.nii.gz \
          "$PRD"/connectivity/predwi_denoised.mif \
          -fslgrad "$PRD"/connectivity/bvecs "$PRD"/connectivity/bvals \
          -datatype float32 -force
mrconvert "$PRD"/connectivity/data_crop.nii.gz \
          "$PRD"/connectivity/predwi_denoised_preproc.mif \
          -fslgrad "$PRD"/connectivity/bvecs "$PRD"/connectivity/bvals \
          -datatype float32 -force
mrconvert "$PRD"/connectivity/data_crop.nii.gz \
          "$PRD"/connectivity/predwi_denoised_preproc_bias.mif \
          -fslgrad "$PRD"/connectivity/bvecs "$PRD"/connectivity/bvals \
          -datatype float32 -force
mrconvert "$PRD"/connectivity/data_crop.nii.gz \
          "$PRD"/connectivity/predwi_denoised_preproc_bias.mif \
          -fslgrad "$PRD"/connectivity/bvecs "$PRD"/connectivity/bvals \
          -datatype float32 -force
mrconvert "$PRD"/connectivity/nodif_brain_mask_crop.nii.gz \
          "$PRD"/connectivity/mask_native.mif -datatype float32 \
          -force

rm -r "$PRD"/"$SUBJ_ID"/