#!/usr/bin/env bash

unzip "$PRD"/"$SUBJ_ID"_3T_Diffusion_preproc.zip
unzip "$PRD"/"$SUBJ_ID"_3T_Structural_preproc.zip

touch "$PRD"/connectivity/predwi.mif
touch "$PRD"/connectivity/predwi_denoised.mif
touch "$PRD"/connectivity/predwi_denoised_preproc.mif
touch "$PRD"/connectivity/predwi_denoised_preproc_bias.mif
touch "$PRD"/connectivity/mask_native.mif
touch "$PRD"/connectivity/lowb.nii.gz
touch "$PRD"/connectivity/brain.nii.gz
touch "$PRD"/connectivity/aparc+aseg.nii.gz
touch "$PRD"/connectivity/aparc+aseg_reorient.nii.gz

cp "$PRD"/"$SUBJ_ID"/T1w/aparc+aseg.nii.g aparcaseg_2_diff.nii.gz
cp "$PRD"/"$SUBJ_ID"/T1w/T1w_acpc_dc_restore_brain.nii.gz brain_2_diff.nii.gz

mrconvert "$PRD"/"$SUBJ_ID"/T1w/Diffusion/data.nii.gz dwi.mif -fslgrad bvecs bvals -datatype float32 -stride 0,0,0,1
mrconvert "$PRD"/"$SUBJ_ID"/T1w/Diffusion/nodif_brain_mask.nii.gz  mask.mif -datatype float32 -stride 0,0,0,1
