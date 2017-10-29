#!/usr/bin/env bash

# Main surface script for SCRIPTS
# Launch with:
# bash path_to_scripts/main_surface.sh -c path_to_config/config.sh
# if you want to write the output instead of a log file, use:
# bash path_to_scripts/main_surface.sh -c path_to_config/config.sh > logfile.txt

# style guide:
# https://google.github.io/styleguide/shell.xml

# TODO: add explicits echo for what to check in the figures
# TODO: add the missing mrviews
# TOCHECK: test/retest

#### Checks and preset variables

# import and check config
while getopts "c:eqf" opt; do
  case $opt in
    c)
      CONFIG=$OPTARG
      if [ ! -f "$CONFIG" -a "$CONFIG" != "test" ];then
        echo "config file "$CONFIG" unexistent" >&2
        exit 1
      elif [ $CONFIG = "test" ]; then
        echo "test mode"
      else
        echo "Using config file $CONFIG." >&2
        source "$CONFIG"
      fi
      ;;
    e) 
      set -e 
      ;;
    q)
      QUIET="yes"
      ;;
    f)
      FORCE="yes"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    esac
done

if [ -z "$CONFIG" ]; then
  echo "You must provide a config file."
  exit 1
fi

# check mandatory variables
if [ -z "$PRD" ]; then
  echo "PRD path missing"
  exit 1
fi

if [ -z "$SUBJ_ID" ]; then
  echo "SUBJ_ID path missing"
  exit 1
fi

if [ -z "$MATLAB" ]; then
  echo "Matlab path missing"
  exit 1
fi

if [ -z "$SUBJECTS_DIR" ]; then
  echo "you have to set the SUBJECTS_DIR environnement variable for FreeSurfer" >> "$PRD"/log_processing_parameters.txt
  exit 1
else
    export FS="$SUBJECTS_DIR"
fi

# set default parameters if not set in config file
echo "##### $( date ) #####" | tee -a "$PRD"/log_processing_parameters.txt

if [ -z "$FORCE" ] || [ "$FORCE" != "no" -a "$FORCE" != "yes" ]; then
  echo "set FORCE parameter to no" | tee -a "$PRD"/log_processing_parameters.txt
  FORCE="no"
else
  echo "FORCE parameter is "$FORCE"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$QUIET" ] || [ "$QUIET" != "no" -a "$QUIET" != "yes" ]; then
  echo "set QUIET parameter to no" | tee -a "$PRD"/log_processing_parameters.txt
  export QUIET="no"
else
  echo "QUIET parameter is "$QUIET"" | tee -a "$PRD"/log_processing_parameters.txt
  # TODO: finish quiet
  export MRTRIX_QUIET=1
fi

if [ -z "$FSL" ] || [ "$FSL" != "fsl5.0" ]; then
  echo "set FSL parameter to empty" | tee -a "$PRD"/log_processing_parameters.txt
  FSL=""
else
  echo "FSL parameter is "$FSL"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$HCP" ] || [ "$HCP" != "no" -a "$HCP" != "yes" ]; then
  echo "set HCP parameter to no" | tee -a "$PRD"/log_processing_parameters.txt
  HCP="no"
else
  echo "HCP parameter is "$HCP"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$CHECK" ] || [ "$CHECK" != "no" -a "$CHECK" != "yes" -a "$CHECK" != "force" ]; then
  echo "set CHECK parameter to no"| tee -a "$PRD"/log_processing_parameters.txt
  export CHECK="no"
else
  echo "CHECK parameter is "$CHECK""| tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$REGISTRATION" ]; then
  echo "set REGISTRATION parameter to regular"| tee -a "$PRD"/log_processing_parameters.txt
  REGISTRATION="regular"
else
  echo "REGISTRATION parameter is "$REGISTRATION""| tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$REGION_MAPPING_CORR" ]; then
  echo "set REGION_MAPPING_CORR parameter to 0.42"| tee -a "$PRD"/log_processing_parameters.txt
  export REGION_MAPPING_CORR=0.42
else
  echo "REGION_MAPPING_CORR parameter is "$REGION_MAPPING_CORR""| tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$NUMBER_TRACKS" ] || ! [[ "$NUMBER_TRACKS" =~ ^[0-9]+$ ]]; then
  echo "set NUMBER_TRACKS parameter to 10.000.000" | tee -a "$PRD"/log_processing_parameters.txt
  NUMBER_TRACKS=10000000
else
  echo "NUMBER_TRACKS parameter is "$NUMBER_TRACKS"" | tee -a "$PRD"/log_processing_parameters.txt
fi

# TODO: check if list of integers
if [ -z "$K_LIST" ]; then
  echo "set K_LIST parameter to empty" | tee -a "$PRD"/log_processing_parameters.txt
  K_LIST=""
else
  echo "K_LIST parameter is "$K_LIST"" | tee -a "$PRD"/log_processing_parameters.txt
fi


if [ -z "$PARCEL" ] || [ "$PARCEL" != "desikan" -a "$PARCEL" != "destrieux" -a "$PARCEL" != "HCP-MMP" ]; then
  echo "set PARCEL parameter to desikan" | tee -a "$PRD"/log_processing_parameters.txt
  PARCEL="desikan"
else
  echo "PARCEL parameter is "$PARCEL"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$TOPUP" ] || [ "$TOPUP" != "no" -a "$TOPUP" != "eddy_correct" ]; then
  echo "set TOPUP parameter to no" | tee -a "$PRD"/log_processing_parameters.txt
  TOPUP="no"
else
  echo "TOPUP parameter is "$TOPUP"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$ACT" ] || [ "$ACT" != "no" -a "$ACT" != "yes" ]; then
  echo "set ACT parameter to yes" | tee -a "$PRD"/log_processing_parameters.txt
  ACT="yes"
else
  echo "ACT parameter is "$ACT"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$SIFT" ] || [ "$SIFT" != "no" -a "$SIFT" != "sift"  -a "$SIFT" != "sift2" ]; then
  echo "set SIFT parameter to sift2" | tee -a "$PRD"/log_processing_parameters.txt
  SIFT="sift2"
else
  echo "SIFT parameter is "$SIFT"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$SIFT_MULTIPLIER" ] || ! [[ "$NUMBER_TRACKS" =~ ^[0-9]+$ ]]; then
  echo "set SIFT_MULTIPLIER parameter to 10" | tee -a "$PRD"/log_processing_parameters.txt
  SIFT_MULTIPLIER=10
else
  echo "SIFT_MULTIPLIER parameter is "$SIFT_MULTIPLIER"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$SEED" ] || [ "$SEED" != "gmwmi" -a "$SEED" != "dynamic" ]; then
  echo "set SEED parameter to dynamic" | tee -a "$PRD"/log_processing_parameters.txt
  SEED="dynamic"
else
  echo "SEED parameter is "$SEED"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$ASEG" ] || [ "$ASEG" != "fs" -a "$ASEG" != "fsl" ]; then
  echo "set ASEG parameter to fsl" | tee -a "$PRD"/log_processing_parameters.txt
  ASEG="fsl"
else
  echo "ASEG parameter is "$ASEG"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z  "$NB_THREADS" ] || ! [[ "$NUMBER_TRACKS" =~ ^[0-9]+$ ]]; then
  if [ -f ~/.mrtrix.conf ]; then
    number_threads_mrtrix_conf=$(grep 'NumberOfThreads' ~/.mrtrix.conf | cut -f 2 -d " ")
    if [ -n "$number_threads_mrtrix_conf" ]; then 
      echo "set number of threads to \
"$number_threads_mrtrix_conf" according to ~/.mrtrix.conf file" | tee -a "$PRD"/log_processing_parameters.txt
      NB_THREADS="$number_threads_mrtrix_conf"
    else
      echo "set number of threads to 1" | tee -a "$PRD"/log_processing_parameters.txt
      NB_THREADS=1
    fi
  else 
    echo "set number of threads to 1" | tee -a "$PRD"/log_processing_parameters.txt
    NB_THREADS=1
  fi
else
echo "number of threads is "$NB_THREADS"" | tee -a "$PRD"/log_processing_parameters.txt
fi

view_step=0

######## HCP pre_scripts
if [ "$HCP" = "yes" ]; then
  if [ ! -f "$PRD"/connectivity/mask_native.mif ]; then
    bash util/HCP_pre_scripts.sh
  fi
fi


######### build cortical surface and region mapping
if [ ! -f "$PRD"/data/T1/T1.nii.gz ]; then
  echo "generating T1 from DICOM"
  mrconvert $PRD/data/T1/ $PRD/data/T1/T1.nii.gz -nthreads "$NB_THREADS"
fi

###################### freesurfer
if [ ! -d "$FS"/"$SUBJ_ID" ] ; then
  echo "running recon-all of freesurfer"
  recon-all -i $PRD/data/T1/T1.nii.gz -s $SUBJ_ID -all
fi

###################################### left hemisphere
# export pial into text file
mkdir -p "$PRD"/surface
if [ ! -f "$PRD"/surface/lh.pial.asc ]; then
  echo "importing left pial surface from freesurfer"
  mris_convert "$FS"/"$SUBJ_ID"/surf/lh.pial "$PRD"/surface/lh.pial.asc
  # take care of the c_(ras) shift which is not done by FS (thks FS!)
  mris_info "$FS"/"$SUBJ_ID"/surf/lh.pial >& "$PRD"/surface/lhinfo.txt
fi

# triangles and vertices high
if [ ! -f "$PRD"/surface/lh_vertices_high.txt ]; then
  echo "extracting left vertices and triangles"
  python util/extract_high.py lh
fi

# decimation using remesher
if [ ! -f $PRD/surface/lh_vertices_low.txt ]; then
  echo "left decimation using remesher"
  # -> to mesh
  python util/txt2off.py $PRD/surface/lh_vertices_high.txt $PRD/surface/lh_triangles_high.txt $PRD/surface/lh_high.off
  #  decimation
  ./remesher/cmdremesher/cmdremesher $PRD/surface/lh_high.off $PRD/surface/lh_low.off
  # export to list vertices triangles
  python util/off2txt.py $PRD/surface/lh_low.off $PRD/surface/lh_vertices_low.txt $PRD/surface/lh_triangles_low.txt
fi

# create the left region mapping
if [ ! -f "$PRD"/surface/lh_region_mapping_low_not_corrected.txt ]; then
  echo "generating the left region mapping on the decimated surface"
  python util/region_mapping.py lh
fi

# correct
if [ ! -f "$PRD"/surface/lh_region_mapping_low.txt ]; then
    echo "correct the left region mapping"
    python util/correct_region_mapping.py lh
    echo "check left region mapping"
    python util/check_region_mapping.py lh
fi

###################################### right hemisphere
# export pial into text file
if [ ! -f "$PRD"/surface/rh.pial.asc ]; then
  echo "importing right pial surface from freesurfer"
  mris_convert $FS/$SUBJ_ID/surf/rh.pial $PRD/surface/rh.pial.asc
  # take care of the c_(ras) shift which is not done by FS (thks FS!)
  mris_info $FS/$SUBJ_ID/surf/rh.pial >& $PRD/surface/rhinfo.txt
fi

# triangles and vertices high
if [ ! -f "$PRD"/surface/rh_vertices_high.txt ]; then
  echo "extracting right vertices and triangles"
  python util/extract_high.py rh
fi

# decimation using brainvisa
if [ ! -f "$PRD"/surface/rh_vertices_low.txt ]; then
  echo "right decimation using remesher"
  # -> to mesh
  python util/txt2off.py $PRD/surface/rh_vertices_high.txt $PRD/surface/rh_triangles_high.txt $PRD/surface/rh_high.off
  #  decimation
  ./remesher/cmdremesher/cmdremesher $PRD/surface/rh_high.off $PRD/surface/rh_low.off
  # export to list vertices triangles
  python util/off2txt.py $PRD/surface/rh_low.off $PRD/surface/rh_vertices_low.txt $PRD/surface/rh_triangles_low.txt
fi

# create the right region mapping
if [ ! -f "$PRD"/surface/rh_region_mapping_low_not_corrected.txt ]; then
  echo "generating the right region mapping on the decimated surface"
  python util/region_mapping.py rh
fi

# correct
if [ ! -f "$PRD"/surface/rh_region_mapping_low.txt ]; then
  echo " correct the right region mapping"
  python util/correct_region_mapping.py rh
  echo "check right region mapping"
  python util/check_region_mapping.py rh
fi
###################################### both hemisphere
# prepare final directory
mkdir -p $PRD/$SUBJ_ID
mkdir -p $PRD/$SUBJ_ID/surface

# reunify both region_mapping, vertices and triangles
if [ ! -f "$PRD"/"$SUBJ_ID"/surface/region_mapping.txt ]; then
  echo "reunify both region mappings"
  python util/reunify_both_regions.py
fi

# zip to put in final format
pushd . > /dev/null
cd $PRD/$SUBJ_ID/surface > /dev/null
zip $PRD/$SUBJ_ID/surface.zip vertices.txt triangles.txt -q
cp region_mapping.txt ..
popd > /dev/null

########################### subcortical surfaces
# extract subcortical surfaces 
if [ ! -f "$PRD"/surface/subcortical/aseg_058_vert.txt ]; then
  echo "generating subcortical surfaces"
  ./util/aseg2srf -s $SUBJ_ID
  mkdir -p $PRD/surface/subcortical
  cp $FS/$SUBJ_ID/ascii/* $PRD/surface/subcortical
  python util/list_subcortical.py
fi

########################## build connectivity using mrtrix 3
mkdir -p $PRD/connectivity
mkdir -p $PRD/$SUBJ_ID/connectivity


## preprocessing
# See: http://mrtrix.readthedocs.io/en/0.3.16/workflows/DWI_preprocessing_for_quantitative_analysis.html

# handle encoding scheme
if [ ! -f "$PRD"/connectivity/predwi.mif ]; then 
  view_step=1
  select_images="n"
  i_im=1
  echo "generate dwi mif file"
  echo "if asked, please select a series of images by typing a number"
  mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi_"$i_im".mif \
            -export_pe_table $PRD/connectivity/pe_table \
            -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init \
            -datatype float32 -stride 0,0,0,1 -force -nthreads "$NB_THREADS"  
  cp $PRD/connectivity/predwi_1.mif $PRD/connectivity/predwi.mif
  if [ "$FORCE" = "no" ]; then
    echo "Do you want to add another image serie (different phase encoding)? [y, n]"
    read select_images
    while [ "$select_images" != "y" ] && [ "$select_images" != "n" ]; do
      echo " please answer y or n"
      read select_images
    done
    while [ "$select_images" == "y" ]; do
      i_im=$(($i_im + 1))
      mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi_"$i_im".mif \
                -export_pe_table $PRD/connectivity/pe_table \
                -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init \
                -datatype float32 -stride 0,0,0,1 -force -nthreads "$NB_THREADS"
      mrcat $PRD/connectivity/predwi.mif $PRD/connectivity/predwi_"$i_im".mif \
            $PRD/connectivity/predwi.mif -axis 3 -nthreads "$NB_THREADS" -force
      echo "Do you want to add another image serie (different phase encoding)? [y, n]"
      read select_images
      while [ "$select_images" != "y" ] && [ "$select_images" != "n" ]; do
        echo " please answer y or n"
        read select_images
      done
    done
  fi
  mrinfo $PRD/connectivity/predwi.mif \
        -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init \
        -export_pe_table $PRD/connectivity/pe_table -force 
fi
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  view_step=0
  echo "check predwi_*.mif files"
  mrview $PRD/connectivity/predwi_*.mif
fi

# denoising the volumes
if [ ! -f "$PRD"/connectivity/predwi_denoised.mif ]; then
  # denoising the combined-directions file is preferable to denoising \
  # predwi1 and 2 separately because of a higher no of volumes
  # see: https://github.com/MRtrix3/mrtrix3/issues/747
  echo "denoising dwi data"
  view_step=1
  dwidenoise $PRD/connectivity/predwi.mif \
             $PRD/connectivity/predwi_denoised.mif \
             -noise $PRD/connectivity/noise.mif -force -nthreads "$NB_THREADS"
  if [ ! -f $PRD/connectivity/noise_res.mif ]; then
    # calculate residuals noise
    mrcalc $PRD/connectivity/predwi.mif \
           $PRD/connectivity/predwi_denoised.mif \
           -subtract $PRD/connectivity/noise_res.mif -nthreads "$NB_THREADS"
  fi
fi
# check noise file: lack of anatomy is a marker of accuracy
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  # noise.mif can also be used for SNR calculation
  echo "check noise/predwi_*_denoised.mif files"
  echo "lack of anatomy in noise_res is a marker of accuracy"
  view_step=0
  mrview $PRD/connectivity/predwi.mif \
         $PRD/connectivity/predwi_denoised.mif \
         $PRD/connectivity/noise.mif \
         $PRD/connectivity/noise_res.mif  
fi

# topup/eddy corrections
if [ ! -f "$PRD"/connectivity/predwi_denoised_preproc.mif ]; then
  view_step=1
  if [ "$TOPUP" = "eddy_correct" ]; then
    # eddy and maybe topup corrections depending of the encoding scheme
    # TODO: removed repol option for now, as it is not in current FSL debian release
    echo "apply eddy and maybe topup if reverse phase-encoding scheme"
    dwipreproc $PRD/connectivity/predwi_denoised.mif \
               $PRD/connectivity/predwi_denoised_preproc.mif \
               -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final \
               -rpe_header -cuda -force -nthreads "$NB_THREADS"    
  else # no topup/eddy
    echo "no topup/eddy applied"
    mrconvert $PRD/connectivity/predwi_denoised.mif \
              $PRD/connectivity/predwi_denoised_preproc.mif \
              -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final \
              -force -nthreads "$NB_THREADS"
  fi
fi
# check preproc files
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ]  && [ -n "$DISPLAY" ]; then
  echo "check preprocessed mif file (no topup/no eddy)"
  view_step=0
  mrview $PRD/connectivity/predwi.mif \
         $PRD/connectivity/predwi_denoised.mif \
         $PRD/connectivity/predwi_denoised_preproc.mif 
fi

# TOCHECK: Masking step before or after biascorrect?
# Native-resolution mask creation
if [ ! -f "$PRD"/connectivity/mask_native.mif ]; then
  echo "create dwi mask"
  view_step=1
  dwi2mask $PRD/connectivity/predwi_denoised_preproc.mif \
           $PRD/connectivity/mask_native.mif -nthreads "$NB_THREADS"
fi
# check mask file
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ]  && [ -n "$DISPLAY" ]; then
  echo "check native mask mif file"
  view_step=0
  mrview "$PRD"/connectivity/predwi_denoised_preproc.mif \
         -overlay.load $PRD/connectivity/mask_native.mif \
         -overlay.opacity 0.5
fi

# Bias field correction
if [ ! -f "$PRD"/connectivity/predwi_denoised_preproc_bias.mif ]; then
  # ANTS seems better than FSL
  # see http://mrtrix.readthedocs.io/en/0.3.16/workflows/DWI_preprocessing_for_quantitative_analysis.html
  if [ -n "$ANTSPATH" ]; then
    echo "bias correct using ANTS"
    dwibiascorrect $PRD/connectivity/predwi_denoised_preproc.mif \
                   $PRD/connectivity/predwi_denoised_preproc_bias.mif \
                   -mask $PRD/connectivity/mask_native.mif \
                   -bias $PRD/connectivity/B1_bias.mif -ants -force \
                   -nthreads "$NB_THREADS"
  else
    echo "bias correct using FSL"
    dwibiascorrect $PRD/connectivity/predwi_denoised_preproc.mif \
                   $PRD/connectivity/predwi_denoised_preproc_bias.mif \
                   -mask $PRD/connectivity/mask_native.mif \
                   -bias $PRD/connectivity/B1_bias.mif -fsl -force \
                   -nthreads "$NB_THREADS"
  fi
fi
# check bias field correction
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check bias field correction"
  mrview $PRD/connectivity/predwi.mif \
         $PRD/connectivity/predwi_denoised_preproc.mif \
         $PRD/connectivity/predwi_denoised_preproc_bias.mif 
fi


# TOCHECK: why not upsampling to vox=1.25 as recommended in mrtrix?
# upsampling and reorienting a la fsl
# reorienting from mrtrix to FSL, RAS to LAS, means -stride -1,+2,+3,+4
# see: http://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html
# upsampling (Dyrby TB. Neuroimage. 2014 Dec;103:202-13.) can help registration
# with structural and is common with mrtrix3 fixel analysis pipeline
# see: http://community.mrtrix.org/t/upsampling-dwi-vs-tckgen-defaults/998/2
if [ ! -f "$PRD"/connectivity/dwi.mif ]; then
  echo "upsample dwi"
  mrresize $PRD/connectivity/predwi_denoised_preproc_bias.mif - -scale 2 -force | \
  mrconvert - -datatype float32 -stride -1,+2,+3,+4 $PRD/connectivity/dwi.mif -force 
fi

if [ ! -f "$PRD"/connectivity/mask.mif ]; then
  # for dwi2fod step, a permissive, dilated mask can be used to minimize
  # streamline premature termination, see BIDS protocol: 
  # https://github.com/BIDS-Apps/MRtrix3_connectome/blob/master/run.py
  echo "upsample mask"
  view_step=1
  mrresize $PRD/connectivity/mask_native.mif - -scale 2 -force | \
  mrconvert - $PRD/connectivity/mask.mif -datatype bit -stride -1,+2,+3 \
            -force -nthreads "$NB_THREADS"
  maskfilter $PRD/connectivity/mask.mif dilate \
             $PRD/connectivity/mask_dilated.mif -npass 2 -force \
             -nthreads "$NB_THREADS" 
fi
# check upsampled files
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check upsampled mif files"
  view_step=0
  mrview $PRD/connectivity/dwi.mif \
         -overlay.load $PRD/connectivity/mask.mif \
         -overlay.load $PRD/connectivity/mask_dilated.mif \
         -overlay.opacity 0.5 -norealign 
fi


## FLIRT registration
# a comparison of registration methods is available in:
# Ou Y, et al. IEEE Trans Med Imaging. 2014 Oct;33(10):2039-65
# FLIRT -bbr (Kerstin's protocol):
# http://community.mrtrix.org/t/registration-of-structural-and-diffusion-weighted-data/203/8


# low b extraction to FSL
if [ ! -f "$PRD"/connectivity/lowb.nii.gz ]; then
  view_step=1
  if [ "$REGISTRATION" = "regular" ] || [ "$REGISTRATION" = "boundary" ]; then
    echo "extracting b0 vols for registration"
    dwiextract $PRD/connectivity/dwi.mif $PRD/connectivity/lowb.mif \
               -bzero -force -nthreads "$NB_THREADS" 
    # stride from mrtrix to FSL, RAS to LAS
    # see: http://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html
    mrconvert $PRD/connectivity/lowb.mif $PRD/connectivity/lowb.nii.gz \
              -stride -1,+2,+3,+4 -force -nthreads "$NB_THREADS" 
    # for visualization 
    mrmath  $PRD/connectivity/lowb.mif mean $PRD/connectivity/meanlowb.mif \
            -axis 3 -force -nthreads "$NB_THREADS"
  elif [ "$REGISTRATION" = "pseudo" ]; then
    # lowb-pseudo brain for pseudo registration
    echo "generate lowb-pseudo vols for pseudo registration"
    # Generate transform image (dwi) for pseudo registration method: 
    # see: Bhushan C, et al. Neuroimage. 2015 Jul 15;115:269-8
    # used in: https://github.com/BIDS-Apps/MRtrix3_connectome/blob/master/run.py
    dwiextract $PRD/connectivity/dwi.mif -bzero - \
    | mrmath - mean - -axis 3 \
    | mrcalc 1 - -divide $PRD/connectivity/mask.mif -multiply - \
    | mrconvert - - -stride -1,+2,+3 \
    | mrhistmatch - $PRD/connectivity/brain.nii.gz \
                  $PRD/connectivity/lowb.nii.gz
  fi
fi
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check lowb image"
  view_step=0
  mrview $PRD/connectivity/lowb.mif \
         -overlay.load $PRD/connectivity/dwi.mif \
         -overlay.opacity 1. -norealign
fi

# generating FSl brain.mgz
if [ ! -f "$PRD"/connectivity/brain.nii.gz ]; then
  # brain.mgz seems to be superior to diff to T1
  # as the main problem for registration is the wmgm interface that we want to
  # remove and BET stripping is unfortunate in many situations, 
  # and FS pial eddited volumes already present
  # stride from FS to FSL: RAS to LAS
  # see: http://www.grahamwideman.com/gw/brain/fs/coords/fscoords.htm
  # we could do
  # mrconvert $FS/$SUBJ_ID/mri/brain.mgz $PRD/connectivity/brain.nii.gz \
  #           -datatype float32 -stride -1,+2,+3,+4 -force -nthreads "$NB_THREADS" 
  # instead we use the pure brain from aparc+aseg:
    echo "generating masked brain in FSL orientation"
  mri_binarize --i $FS/$SUBJ_ID/mri/aparc+aseg.mgz \
               --o $FS/$SUBJ_ID/mri/aparc+aseg_mask.mgz --min 0.5 --dilate 1 
  mri_mask $FS/$SUBJ_ID/mri/brain.mgz $FS/$SUBJ_ID/mri/aparc+aseg_mask.mgz \
           $FS/$SUBJ_ID/mri/brain_masked.mgz
  mrconvert $FS/$SUBJ_ID/mri/brain_masked.mgz $PRD/connectivity/brain.nii.gz \
            -force -datatype float32 -stride -1,+2,+3
fi



# aparc+aseg to FSL
if [ ! -f "$PRD"/connectivity/aparc+aseg.nii.gz ]; then
  echo "generating FSL orientation for aparc+aseg"
  # stride from FS to FSL: RAS to LAS
  if [ $PARCEL = "desikan" ]; then
    mrconvert $FS/$SUBJ_ID/mri/aparc+aseg.mgz \
            $PRD/connectivity/aparc+aseg.nii.gz -stride -1,+2,+3 -force \
            -nthreads "$NB_THREADS" 
  elif [ $PARCEL = "destrieux" ]; then
    mrconvert $FS/$SUBJ_ID/mri/aparc.a2009s+aseg.mgz \
            $PRD/connectivity/aparc+aseg.nii.gz -stride -1,+2,+3 -force \
            -nthreads "$NB_THREADS" 
  fi
fi

# check orientations
if [ ! -f $PRD/connectivity/aparc+aseg_reorient.nii.gz ]; then
  echo "reorienting the region parcellation"
  view_step=1
  "$FSL"fslreorient2std $PRD/connectivity/aparc+aseg.nii.gz \
                  $PRD/connectivity/aparc+aseg_reorient.nii.gz
fi
# check parcellation to brain.mgz
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  # TODO: mrview discrete colour scheme?
  echo "check parcellation"
  echo "if it's correct, just close the window." 
  echo "Otherwise... well, it should be correct anyway"
  view_step=0
  mrview $PRD/connectivity/brain.nii.gz \
         -overlay.load $PRD/connectivity/aparc+aseg_reorient.nii.gz \
         -overlay.opacity 0.5 -norealign
fi

# aparcaseg to diff by inverse transform
if [ ! -f "$PRD"/connectivity/aparcaseg_2_diff.nii.gz ]; then
  view_step=1
  if [ "$REGISTRATION" = "regular" ] || [ "$REGISTRATION" = "pseudo" ]; then
    # 6 dof; see:
    # http://web.mit.edu/fsl_v5.0.8/fsl/doc/wiki/FLIRT(2f)FAQ.html#What_cost_function_and.2BAC8-or_degrees_of_freedom_.28DOF.29_should_I_use_in_FLIRT.3F
    echo "register aparc+aseg to diff"
    "$FSL"flirt -in $PRD/connectivity/lowb.nii.gz \
                -out $PRD/connectivity/lowb_2_struct.nii.gz \
                -ref $PRD/connectivity/brain.nii.gz \
                -omat $PRD/connectivity/diffusion_2_struct.mat -dof 6 \
                -searchrx -180 180 -searchry -180 180 -searchrz -180 180 \
                -cost mutualinfo
  elif [ "$REGISTRATION" = "boundary" ]; then
    echo "register aparc+aseg to diff using bbr cost function in FLIRT"
    # as per http://community.mrtrix.org/t/registration-of-structural-and-diffusion-weighted-data/203/8
    "$FSL"fast -N -o $PRD/connectivity/brain_fast $PRD/connectivity/brain.nii.gz 
    "$FSL"fslmaths $PRD/connectivity/brain_fast_pve_2.nii.gz -thr 0.5 \
                   -bin $PRD/connectivity/brain_fast_wmmask.nii.gz
    # first flirt to get an init transform mat
    "$FSL"flirt -in $PRD/connectivity/lowb.nii.gz \
                -ref $PRD/connectivity/brain.nii.gz \
                -omat $PRD/connectivity/flirt_bbr_tmp.mat -dof 6 
    # flirt using bbr cost
    "$FSL"flirt -in $PRD/connectivity/lowb.nii.gz \
                -out $PRD/connectivity/lowb_2_struct.nii.gz \
                -ref $PRD/connectivity/brain.nii.gz \
                -omat $PRD/connectivity/diffusion_2_struct.mat \
                -wmseg $PRD/connectivity/brain_fast_wmmask.nii.gz \
                -init $PRD/connectivity/flirt_bbr_tmp.mat \
                -schedule $FSLDIR/etc/flirtsch/bbr.sch -dof 6 \
                -searchrx -180 180 -searchry -180 180 -searchrz -180 180 \
                -cost bbr 
  fi
  transformconvert $PRD/connectivity/diffusion_2_struct.mat \
                   $PRD/connectivity/lowb.nii.gz \
                   $PRD/connectivity/brain.nii.gz \
                   flirt_import $PRD/connectivity/diffusion_2_struct_mrtrix.txt \
                   -force 
  mrtransform $PRD/connectivity/aparc+aseg_reorient.nii.gz \
              $PRD/connectivity/aparcaseg_2_diff.nii.gz \
              -linear $PRD/connectivity/diffusion_2_struct_mrtrix.txt \
              -inverse -datatype uint32 -force -nthreads "$NB_THREADS"
fi

# brain to diff by inverse transform
if [ ! -f "$PRD"/connectivity/brain_2_diff.nii.gz ]; then
  echo "register brain to diff"
  view_step=1
  mrtransform $PRD/connectivity/brain.nii.gz \
              $PRD/connectivity/brain_2_diff.nii.gz \
              -linear $PRD/connectivity/diffusion_2_struct_mrtrix.txt \
              -inverse -force 
fi
# check brain and parcellation to diff
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check parcellation registration to diffusion space"
  echo "if it's correct, just close the window."
  echo "Otherwise you will have to do the registration by hand"
  view_step=0
  mrview $PRD/connectivity/brain_2_diff.nii.gz \
         $PRD/connectivity/lowb.nii.gz \
         -overlay.load $PRD/connectivity/aparcaseg_2_diff.nii.gz \
         -overlay.opacity 0.5 -norealign
fi



# prepare file for act
if [ "$ACT" = "yes" ] && [ ! -f $PRD/connectivity/act.mif ]; then
  echo "prepare files for act"
  view_step=1
  5ttgen fsl $PRD/connectivity/brain_2_diff.nii.gz $PRD/connectivity/act.mif \
         -premasked -force  -nthreads "$NB_THREADS"
  5tt2vis $PRD/connectivity/act.mif $PRD/connectivity/act_vis.mif -force \
        -nthreads "$NB_THREADS"
fi
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
    echo "check tissue segmented image"
    view_step=0
    mrview $PRD/connectivity/act_vis.mif -colourmap 4
fi

# Response function estimation
# Check if multi or single shell
shells=$(mrinfo -shells $PRD/connectivity/dwi.mif)
echo "shell b values are $shells"
nshells=($shells)
no_shells=${#nshells[@]}
echo "no of shells are $no_shells"

if [ "$no_shells" -gt 2 ]; then
# Multishell
  if [ ! -f "$PRD"/connectivity/response_wm.txt ]; then
    view_step=1
    if [ "$ACT" = "yes" ]; then 
      echo "estimating response using msmt algorithm"
      dwi2response msmt_5tt $PRD/connectivity/dwi.mif \
                   $PRD/connectivity/act.mif \
                   $PRD/connectivity/response_wm.txt \
                   $PRD/connectivity/response_gm.txt \
                   $PRD/connectivity/response_csf.txt \
                   -voxels $PRD/connectivity/RF_voxels.mif \
                   -mask $PRD/connectivity/mask.mif -force \
                   -nthreads "$NB_THREADS"
    else
      echo "estimating response using dhollander algorithm"
      dwi2response dhollander $PRD/connectivity/dwi.mif \
                   $PRD/connectivity/response_wm.txt \
                   $PRD/connectivity/response_gm.txt \
                   $PRD/connectivity/response_csf.txt \
                   -voxels $PRD/connectivity/RF_voxels.mif \
                   -mask $PRD/connectivity/mask.mif -force \
                   -nthreads "$NB_THREADS"
    fi
  fi
else
# Single shell only
  if [ ! -f "$PRD"/connectivity/response_wm.txt ]; then
    echo "estimating response using dhollander algorithm"
    view_step=1
    # dwi2response tournier $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt -force -voxels $PRD/connectivity/RF_voxels.mif -mask $PRD/connectivity/mask.mif
    dwi2response dhollander $PRD/connectivity/dwi.mif \
                 $PRD/connectivity/response_wm.txt \
                 $PRD/connectivity/response_gm.txt \
                 $PRD/connectivity/response_csf.txt \
                 -voxels $PRD/connectivity/RF_voxels.mif \
                 -mask $PRD/connectivity/mask.mif -force \
                 -nthreads "$NB_THREADS"
  fi
fi
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check ODF image"
  view_step=0
  mrview $PRD/connectivity/meanlowb.mif \
         -overlay.load $PRD/connectivity/RF_voxels.mif \
         -overlay.opacity 0.5
fi

# Fibre orientation distribution estimation
if [ ! -f "$PRD"/connectivity/wm_CSD.mif ]; then
  # Both for multishell and single shell since we use dhollander in the 
  # single shell case
  # see: http://community.mrtrix.org/t/wm-odf-and-response-function-with-dhollander-option---single-shell-versus-multi-shell/572/4
  echo "calculating fod on multishell or single shell data"
  view_step=1
  dwi2fod msmt_csd $PRD/connectivity/dwi.mif \
          $PRD/connectivity/response_wm.txt \
          $PRD/connectivity/wm_CSD.mif \
          $PRD/connectivity/response_gm.txt \
          $PRD/connectivity/gm_CSD.mif \
          $PRD/connectivity/response_csf.txt \
          $PRD/connectivity/csf_CSD.mif \
          -mask $PRD/connectivity/mask_dilated.mif -force \
          -nthreads "$NB_THREADS"
fi
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check ODF image"
  view_step=0
  mrconvert $PRD/connectivity/wm_CSD.mif - -coord 3 0 \
  -nthreads "$NB_THREADS" -force \
  | mrcat $PRD/connectivity/csf_CSD.mif \
          $PRD/connectivity/gm_CSD.mif - \
          $PRD/connectivity/tissueRGB.mif -axis 3 -nthreads "$NB_THREADS" \
          -force
  mrview $PRD/connectivity/tissueRGB.mif \
         -odf.load_sh $PRD/connectivity/wm_CSD.mif 
fi


# tractography
if [ ! -f "$PRD"/connectivity/whole_brain.tck ]; then
  if [ "$SIFT" = "sift" ]; then
    # temporarily change number of tracks for sift
    number_tracks=$(($NUMBER_TRACKS*$SIFT_MULTIPLIER))
  fi
  native_voxelsize=$(mrinfo $PRD/connectivity/mask_native.mif -vox \
                   | cut -f 1 -d " " | xargs printf "%.3f")
  stepsize=$( bc -l <<< "scale=2; "$native_voxelsize"/2" )
  echo "stepsize parameter for tckgen is $stepsize"
  angle=$( bc -l <<< "scale=2; 90*"$stepsize"/"$native_voxelsize"" )
  echo "angle parameter for tckgen is $angle"
  if [ "$ACT" = "yes" ]; then
    # when using msmt_csd in conjunction with ACT, the cutoff threshold
    # can be reduced to 0.06
    # see: https://github.com/MRtrix3/mrtrix3/blob/master/docs/quantitative_structural_connectivity/ismrm_hcp_tutorial.rst#connectome-generation
    echo "generating tracks using act"
    if [ "$SEED" = "gmwmi" ]; then
      echo "seeding from gmwmi" 
      5tt2gmwmi $PRD/connectivity/act.mif \
                $PRD/connectivity/gmwmi_mask.mif -force \
                -nthreads "$NB_THREADS"
      # TODO: min length check andreas paper
      tckgen $PRD/connectivity/wm_CSD.mif \
             $PRD/connectivity/whole_brain.tck \
             -seed_gmwmi $PRD/connectivity/gmwmi_mask.mif \
             -act $PRD/connectivity/act.mif -select "$NUMBER_TRACKS" \
             -seed_unidirectional -crop_at_gmwmi -backtrack \
             -minlength 4 -maxlength 250 -step "$stepsize" -angle "$angle" \
             -cutoff 0.06 -force -nthreads "$NB_THREADS"
    elif [ "$SEED" = "dynamic" ]; then
       # -dynamic seeding may work slightly better than gmwmi, 
       # see Smith RE Neuroimage. 2015 Oct 1;119:338-51.
      echo "seeding dynamically"   
      tckgen $PRD/connectivity/wm_CSD.mif \
             $PRD/connectivity/whole_brain.tck \
             -seed_dynamic $PRD/connectivity/wm_CSD.mif \
             -act $PRD/connectivity/act.mif -select "$NUMBER_TRACKS" \
             -crop_at_gmwmi -backtrack -minlength 4 -maxlength 250 \
             -step "$stepsize" -angle "$angle" -cutoff 0.06 -force \
             -nthreads "$NB_THREADS"
    fi  
  else
    echo "generating tracks without using act" 
    echo "seeding dynamically" 
    tckgen $PRD/connectivity/wm_CSD.mif \
           $PRD/connectivity/whole_brain.tck \
           -seed_dynamic $PRD/connectivity/wm_CSD.mif \
           -mask $PRD/connectivity/mask.mif -select "$NUMBER_TRACKS" \
           -maxlength 250  -step "$stepsize"  -angle "$angle" -cutoff 0.1 \
           -force -nthreads "$NB_THREADS"
  fi
fi

# postprocessing
if [ ! -e "$PRD"/connectivity/whole_brain_post.tck ]; then
  echo "$PRD"/connectivity/whole_brain_post.tck
  if [ "$SIFT" = "sift" ]; then
    echo "using sift"
    number_tracks=$(($NUMBER_TRACKS/$SIFT_MULTIPLIER))
    if [ "$ACT" = "yes" ]; then
        echo "trimming tracks using sift/act" 
        tcksift $PRD/connectivity/whole_brain.tck \s
                $PRD/connectivity/wm_CSD.mif \
                $PRD/connectivity/whole_brain_post.tck \
                -act $PRD/connectivity/act.mif \
                -out_mu $PRD/connectivity/mu.txt \
                -term_number $NUMBER_TRACKS -fd_scale_gm -force \
                -nthreads "$NB_THREADS"
    else
        echo "trimming tracks using sift/without act" 
        tcksift $PRD/connectivity/whole_brain.tck \
                $PRD/connectivity/wm_CSD.mif \
                $PRD/connectivity/whole_brain_post.tck \
                -out_mu $PRD/connectivity/mu.txt \
                -term_number $NUMBER_TRACKS -force \
                -nthreads "$NB_THREADS"
    fi
  elif [ "$SIFT" = "sift2" ]; then 
    echo "running sift2"
    ln -s $PRD/connectivity/whole_brain.tck $PRD/connectivity/whole_brain_post.tck
    if [ "$ACT" = "yes" ]; then
      echo "using act" 
      tcksift2 $PRD/connectivity/whole_brain.tck \
               $PRD/connectivity/wm_CSD.mif \
               $PRD/connectivity/streamline_weights.csv\
               -act $PRD/connectivity/act.mif \
               -out_mu $PRD/connectivity/mu.txt \
               -out_coeffs $PRD/connectivity/streamline_coeffs.csv \
               -fd_scale_gm -force -nthreads "$NB_THREADS"
    else
      tcksift2 $PRD/connectivity/whole_brain.tck \
               $PRD/connectivity/wm_CSD.mif \
               $PRD/connectivity/streamline_weights.csv \
               -out_mu $PRD/connectivity/mu.txt \
               -out_coeffs $PRD/connectivity/streamline_coeffs.csv \
               -force -nthreads "$NB_THREADS"
    fi
  else 
    echo "not using sift2"
    ln -s $PRD/connectivity/whole_brain.tck \
          $PRD/connectivity/whole_brain_post.tck
  fi
fi

## now compute connectivity and length matrix
if [ ! -f "$PRD"/connectivity/aparcaseg_2_diff_"$ASEG".mif ]; then
  echo "compute parcellation labels"
  python util/compute_luts.py
  labelconvert $PRD/connectivity/aparcaseg_2_diff.nii.gz \
               $PRD/connectivity/lut_in.txt $PRD/connectivity/lut_out.txt\
               $PRD/connectivity/aparcaseg_2_diff_fs.mif \
               -force -nthreads "$NB_THREADS"
  echo "$ASEG"
  if [ "$ASEG" = "fsl" ]; then
    # FS derived subcortical parcellation is too variable and prone to 
    # errors => labelsgmfix) was generated, 
    # see Smith RE Neuroimage. 2015 Jan 1;104:253-65.
    # TODO: check effect on region mapping
    # TODO; -sgm_amyg_hipp option to consider
    echo "fix FS subcortical labels to generate FSL labels"
    labelsgmfix $PRD/connectivity/aparcaseg_2_diff_fs.mif \
                $PRD/connectivity/brain_2_diff.nii.gz $PRD/connectivity/lut_out.txt \
                $PRD/connectivity/aparcaseg_2_diff_fsl.mif -premasked \
                -force -nthreads "$NB_THREADS"   
  fi 
fi

if [ ! -f "$PRD"/connectivity/weights.csv ]; then
  echo "compute connectivity matrix weights"
  if [ "$SIFT" = "sift2" ]; then
    # -tck_weights_in flag only needed for sift2 but not for sift/no processing
    tck2connectome $PRD/connectivity/whole_brain_post.tck \
                   $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif \
                   $PRD/connectivity/weights.csv -assignment_radial_search 2 \
                   -out_assignments $PRD/connectivity/edges_2_nodes.csv \
                   -tck_weights_in $PRD/connectivity/streamline_weights.csv \
                   -force -nthreads "$NB_THREADS"
  else
    tck2connectome $PRD/connectivity/whole_brain_post.tck \
                   $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif \
                   $PRD/connectivity/weights.csv -assignment_radial_search 2 \
                   -out_assignments $PRD/connectivity/edges_2_nodes.csv \
                   -force -nthreads "$NB_THREADS"
  fi
fi

if [ ! -f "$PRD"/connectivity/tract_lengths.csv ]; then
  echo "compute connectivity matrix edge lengths"
  view_step=1
  # mean length result: weight by the length, then average
  # see: http://community.mrtrix.org/t/tck2connectome-edge-statistic-sift2-questions/1059/2 
  # Not applying sift2, as here the mean is \
  # sum(streamline length * streamline weight)/no streamlines, does not make sense
  tck2connectome $PRD/connectivity/whole_brain_post.tck \
                 $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif \
                 $PRD/connectivity/tract_lengths.csv \
                 -assignment_radial_search 2 -zero_diagonal -scale_length \
                 -stat_edge mean -force -nthreads "$NB_THREADS"
fi

# view connectome
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "view connectome edges as lines or streamlines"
  if [ ! -f $PRD/connectivity/exemplars.tck ]; then
    if [ "$SIFT" = "sift2" ]; then
        connectome2tck $PRD/connectivity/whole_brain_post.tck \
                       $PRD/connectivity/edges_2_nodes.csv \
                       $PRD/connectivity/exemplars.tck \
                       -exemplars $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif \
                       -tck_weights_in $PRD/connectivity/streamline_weights.csv \
                       -files single -nthreads "$NB_THREADS"
    else 
        connectome2tck $PRD/connectivity/whole_brain_post.tck \
                       $PRD/connectivity/edges_2_nodes.csv \
                       $PRD/connectivity/exemplars.tck \
                       -exemplars $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif \
                       -files single -nthreads "$NB_THREADS"
    fi
  fi
  # TOCHECK: in mrview, load the lut table ($PRD/connectivity/lut_out.txt) for node correspondence, 
  # and exemplars.tck if wanting to see edges as streamlines 
  mrview $PRD/connectivity/aparcaseg_2_diff_$ASEG.mif \
         -connectome.init $PRD/connectivity/aparcaseg_2_diff_$ASEG.mif \
         -connectome.load $PRD/connectivity/weights.csv 
fi

# view tractogram and tdi
if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "view tractogram and tdi image"
  view_step=0
  if [ ! -f $PRD/connectivity/whole_brain_post_decimated.tck ]; then
    # $(( number_tracks/100)) this follows some recommendations by 
    # JD Tournier to avoid mrview to be to slow
    # (for visualization no more than 100-200K streamlines)
    # => min(100k, number_tracks/100)
    if [ "$SIFT" = "sift2" ]; then
        tckedit $PRD/connectivity/whole_brain_post.tck \
                $PRD/connectivity/whole_brain_post_decimated.tck \
                -tck_weights_in $PRD/connectivity/streamline_weights.csv \
                -number $(($NUMBER_TRACKS<100000?$NUMBER_TRACKS:100000)) \
                -minweight 1 -force -nthreads "$NB_THREADS"
    else 
        tckedit $PRD/connectivity/whole_brain_post.tck \
                $PRD/connectivity/whole_brain_post_decimated.tck \
                -number $(($NUMBER_TRACKS<100000?$NUMBER_TRACKS:100000)) \
                -force -nthreads "$NB_THREADS"
    fi  
  fi
  if [ ! -f $PRD/connectivity/whole_brain_post_tdi.mif ]; then
      if [ "$SIFT" = "sift2" ]; then
          tckmap $PRD/connectivity/whole_brain_post.tck \
                 $PRD/connectivity/whole_brain_post_tdi.mif \
                 -tck_weights_in $PRD/connectivity/streamline_weights.csv \
                 -dec -vox 1 -force -nthreads "$NB_THREADS"
      else
          tckmap $PRD/connectivity/whole_brain_post.tck \
                 $PRD/connectivity/whole_brain_post_tdi.mif \
                 -dec -vox 1 -force -nthreads "$NB_THREADS"
      fi 
  fi
  mrview $PRD/connectivity/aparcaseg_2_diff_$ASEG.mif \
         -overlay.load $PRD/connectivity/whole_brain_post_tdi.mif \
         -overlay.opacity 0.5 -overlay.interpolation_off \
         -tractography.load $PRD/connectivity/whole_brain_post_decimated.tck 
fi


# Compute other files
# we do not compute hemisphere
# subcortical is already done

# compute centers, areas, cortical and orientations
if [ ! -f "$PRD"/"$SUBJ_ID"/connectivity/cortical.txt ]; then
  echo "generate useful files for TVB"
  python util/compute_connectivity_files.py
fi

# zip to put in final format
pushd . > /dev/null
cd $PRD/$SUBJ_ID/connectivity > /dev/null
zip $PRD/$SUBJ_ID/connectivity.zip areas.txt average_orientations.txt \
  weights.txt tract_lengths.txt cortical.txt centres.txt -q
popd > /dev/null 


################### subparcellations
# compute sub parcellations connectivity if asked
if [ -n "$K_LIST" ]; then
  for K in $K_LIST; do
    export curr_K=$(( 2**K ))
    mkdir -p $PRD/$SUBJ_ID/connectivity_"$curr_K"
    if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii.gz ]; then
      echo "compute subparcellations for $curr_K"
      $MATLAB -r "run subparcel.m; quit;" -nodesktop -nodisplay 
      gzip $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii
    fi
    if [ ! -f $PRD/$SUBJ_ID/region_mapping_"$curr_K".txt ]; then
      echo "generate region mapping for subparcellation "$curr_K""
      python util/region_mapping_other_parcellations.py
    fi
    if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif ]; then
      mrconvert $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii.gz \
                   $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif \
                   -datatype float32 -force
    fi
    if [ ! -f $PRD/connectivity/weights_"$curr_K".csv ]; then
      echo "compute connectivity matrix using act for subparcellation "$curr_K""
      if [ "$SIFT" = "sift2" ]; then
      # -tck_weights_in flag only needed for sift2 but not for sift/no processing
      tck2connectome $PRD/connectivity/whole_brain_post.tck \
                     $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif \
                     $PRD/connectivity/weights_"$curr_K".csv -assignment_radial_search 2 \
                     -out_assignments $PRD/connectivity/edges_2_nodes_"$curr_K".csv \
                     -tck_weights_in $PRD/connectivity/streamline_weights.csv \
                     -force -nthreads "$NB_THREADS"
      else
      tck2connectome $PRD/connectivity/whole_brain_post.tck \
                     $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif \
                     $PRD/connectivity/weights.csv -assignment_radial_search 2 \
                     -out_assignments $PRD/connectivity/edges_2_nodes_"$curr_K".csv \
                     -force -nthreads "$NB_THREADS"
      fi
    fi
    if [ ! -f "$PRD"/connectivity/tract_lengths_"$curr_K".csv ]; then
      echo "compute connectivity matrix edge lengths subparcellation "$curr_K""
      view_step=1
      # mean length result: weight by the length, then average
      # see: http://community.mrtrix.org/t/tck2connectome-edge-statistic-sift2-questions/1059/2 
      # Not applying sift2, as here the mean is \
      # sum(streamline length * streamline weight)/no streamlines, does not make sense
      tck2connectome $PRD/connectivity/whole_brain_post.tck \
                     $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif \
                     $PRD/connectivity/tract_lengths_"$curr_K".csv \
                     -assignment_radial_search 2 -zero_diagonal -scale_length \
                     -stat_edge mean -force -nthreads "$NB_THREADS"
      #fi
    fi
    if [ ! -f "$PRD"/"$SUBJ_ID"/connectivity_"$curr_K"/weights.txt ]; then
      echo "generate files for TVB subparcellation "$curr_K""
      python util/compute_connectivity_sub.py $PRD/connectivity/weights_"$curr_K".csv $PRD/connectivity/tract_lengths_"$curr_K".csv $PRD/$SUBJ_ID/connectivity_"$curr_K"/weights.txt $PRD/$SUBJ_ID/connectivity_"$curr_K"/tract_lengths.txt
    fi
    pushd . > /dev/null
    cd $PRD/$SUBJ_ID/connectivity_"$curr_K" > /dev/null
    zip $PRD/$SUBJ_ID/connectivity_"$curr_K".zip weights.txt tract_lengths.txt centres.txt average_orientations.txt -q 
    popd > /dev/null
  done
fi

######################## compute MEG and EEG forward projection matrices
# make BEM surfaces
if [ ! -h "$FS"/"$SUBJ_ID"/bem/inner_skull.surf ]; then
  echo "generating bem surfaces"
  mne_watershed_bem --subject ${SUBJ_ID} --overwrite
  ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_inner_skull_surface \
        ${FS}/${SUBJ_ID}/bem/inner_skull.surf
  ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_outer_skin_surface  \
        ${FS}/${SUBJ_ID}/bem/outer_skin.surf
  ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_outer_skull_surface \
        ${FS}/${SUBJ_ID}/bem/outer_skull.surf
fi

# export to ascii
if [ ! -f "$FS"/"$SUBJ_ID"/bem/inner_skull.asc ]; then
  echo "importing bem surface from freesurfer"
  mris_convert $FS/$SUBJ_ID/bem/inner_skull.surf $FS/$SUBJ_ID/bem/inner_skull.asc
  mris_convert $FS/$SUBJ_ID/bem/outer_skull.surf $FS/$SUBJ_ID/bem/outer_skull.asc
  mris_convert $FS/$SUBJ_ID/bem/outer_skin.surf $FS/$SUBJ_ID/bem/outer_skin.asc
fi

# triangles and vertices bem
if [ ! -f $PRD/$SUBJ_ID/surface/inner_skull_vertices.txt ]; then
  echo "extracting bem vertices and triangles"
  python util/extract_bem.py inner_skull 
  python util/extract_bem.py outer_skull 
  python util/extract_bem.py outer_skin 
fi

if [ ! -f "$FS"/"$SUBJ_ID"/bem/"$SUBJ_ID"-head.fif ]; then
  echo "generating head bem"
  view_step=1
  mkheadsurf -s $SUBJ_ID
  mne_surf2bem --surf ${FS}/${SUBJ_ID}/surf/lh.seghead --id 4 --check \
               --fif ${FS}/${SUBJ_ID}/bem/${SUBJ_ID}-head.fif --overwrite
fi

if [ "$view_step" = 1 -a "$CHECK" = "yes" ] || [ "$CHECK" = "force" ] && [ -n "$DISPLAY" ]; then
  echo "check bem surfaces"
  view_step=0
  # TODO: use mrview instead
  freeview -v ${FS}/${SUBJ_ID}/mri/T1.mgz \
           -f ${FS}/${SUBJ_ID}/bem/inner_skull.surf:color=yellow:edgecolor=yellow \
           ${FS}/${SUBJ_ID}/bem/outer_skull.surf:color=blue:edgecolor=blue \
           ${FS}/${SUBJ_ID}/bem/outer_skin.surf:color=red:edgecolor=red
fi

# Setup BEM
if [ ! -f "$FS"/"$SUBJ_ID"/bem/*-bem.fif ]; then
  worked=0
  outershift=0
  while [ "$worked" == 0 ]; do
    echo "try generate forward model with 0 shift"
    worked=1
    mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 \
                            --outershift $outershift || worked=0 
    if [ "$worked" == 0 ]; then
      echo "try generate foward model with 1 shift"
      worked=1
      mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 --outershift 1 \
      || worked=0 
    fi
    if [ "$worked" == 0 ] && [ "$CHECK" = "yes" ]; then
      echo 'you can try using a different shifting value for outer skull, 
           please enter a value in mm'
      read outershift;
      echo $outershift
    elif [ "$worked" == 0 ]; then
      echo "bem did not worked"
      worked=1
    elif [ "$worked" == 1 ]; then
      echo "success!"
    fi
  done
fi

exit
