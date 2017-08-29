#!/usr/bin/env bash

# Main surface script for SCRIPTS
# Launch with:
# bash path_to_scripts/main_surface.sh -c path_to_config/config.sh
# if you want to write the output instead of a log file, use:
# bash path_to_scripts/main_surface.sh -c path_to_config/config.sh > logfile.txt

# style guide:
# https://google.github.io/styleguide/shell.xml

# TODO: add explicits echo for what to check in the figures

#### Checks and preset variables

# import and check config
while getopts ":c:" opt; do
  case $opt in
  c)
    CONFIG=$OPTARG
    echo "use config file $CONFIG" >&2
    if [ ! -f $CONFIG ]
    then
    echo "config file unexistent" >&2
    exit 1
    fi
    source "$CONFIG"
    ;;
  \?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done

if [ -z "$CONFIG" ]; then
  echo "you must provide a config file"
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

if [ -z "$MATLAB" ] && [ -z "$MCR"]; then
  echo "matlab or MCR path missing"
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

if [ -z "$FSL" ] || [ "$FSL" != "fsl5.0" ]; then
  echo "set FSL parameter to empty" | tee -a "$PRD"/log_processing_parameters.txt
  FSL=""
else
  echo "FSL parameter is "$FSL"" | tee -a "$PRD"/log_processing_parameters.txt
fi

if [ -z "$CHECK" ] || [ "$CHECK" != "no" -a "$CHECK" != "yes" ]; then
  echo "set CHECK parameter to no"| tee -a "$PRD"/log_processing_parameters.txt
  export CHECK="no"
else
  echo "CHECK parameter is "$CHECK""| tee -a "$PRD"/log_processing_parameters.txt
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

if [ -z "$TOPUP" ] || [ "$TOPUP" != "no" -a "$TOPUP" != "reversed"  -a "$TOPUP" != "eddy_correct" ]; then
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
"$number_threads_mrtrix_conf" according to .mrtrix.conf file" | tee -a "$PRD"/log_processing_parameters.txt
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

######### build cortical surface and region mapping
if [ ! -f $PRD/data/T1/T1.nii.gz ]
then
  echo "generating T1 from DICOM"
  mrconvert $PRD/data/T1/ $PRD/data/T1/T1.nii.gz -nthreads "$NB_THREADS"
fi

###################### freesurfer
if [ ! -d $FS/$SUBJ_ID ] 
then
  echo "running recon-all of freesurfer"
  recon-all -i $PRD/data/T1/T1.nii.gz -s $SUBJ_ID -all
fi

###################################### left hemisphere
# export pial into text file
mkdir -p $PRD/surface
if [ ! -f $PRD/surface/lh.pial.asc ]; then
  echo "importing left pial surface from freesurfer"
  mris_convert $FS/$SUBJ_ID/surf/lh.pial $PRD/surface/lh.pial.asc
  # take care of the c_(ras) shift which is not done by FS (thks FS!)
  mris_info $FS/$SUBJ_ID/surf/lh.pial >& $PRD/surface/lhinfo.txt
fi

# triangles and vertices high
if [ ! -f $PRD/surface/lh_vertices_high.txt ]; then
  echo "extracting left vertices and triangles"
  python extract_high.py lh
fi

# decimation using brainvisa
if [ ! -f $PRD/surface/lh_vertices_low.txt ]; then
  echo "left decimation using remesher"
  # -> to mesh
  python txt2off.py $PRD/surface/lh_vertices_high.txt $PRD/surface/lh_triangles_high.txt $PRD/surface/lh_high.off
  #  decimation
  ./remesher/cmdremesher/cmdremesher $PRD/surface/lh_high.off $PRD/surface/lh_low.off
  # export to list vertices triangles
  python off2txt.py $PRD/surface/lh_low.off $PRD/surface/lh_vertices_low.txt $PRD/surface/lh_triangles_low.txt
fi

# create left the region mapping
if [ ! -f $PRD/surface/lh_region_mapping_low_not_corrected.txt ]; then
  echo "generating the left region mapping on the decimated surface"
  if [ -n "$MATLAB" ]; then
      $MATLAB -r "rl='lh';run region_mapping.m; quit;" -nodesktop -nodisplay
  else
      sh region_mapping/distrib/run_region_mapping.sh $MCR
  fi
fi

# correct
if [ ! -f $PRD/surface/lh_region_mapping_low.txt ]; then
    echo "correct the left region mapping"
    python correct_region_mapping.py lh
    echo "check left region mapping"
    python check_region_mapping.py lh
fi

###################################### right hemisphere
# export pial into text file
if [ ! -f $PRD/surface/rh.pial.asc ]; then
  echo "importing right pial surface from freesurfer"
  mris_convert $FS/$SUBJ_ID/surf/rh.pial $PRD/surface/rh.pial.asc
  # take care of the c_(ras) shift which is not done by FS (thks FS!)
  mris_info $FS/$SUBJ_ID/surf/rh.pial >& $PRD/surface/rhinfo.txt
fi

# triangles and vertices high
if [ ! -f $PRD/surface/rh_vertices_high.txt ]; then
  echo "extracting right vertices and triangles"
  python extract_high.py rh
fi

# decimation using brainvisa
if [ ! -f $PRD/surface/rh_vertices_low.txt ]; then
  echo "right decimation using remesher"
  # -> to mesh
  python txt2off.py $PRD/surface/rh_vertices_high.txt $PRD/surface/rh_triangles_high.txt $PRD/surface/rh_high.off
  #  decimation
  ./remesher/cmdremesher/cmdremesher $PRD/surface/rh_high.off $PRD/surface/rh_low.off
  # export to list vertices triangles
  python off2txt.py $PRD/surface/rh_low.off $PRD/surface/rh_vertices_low.txt $PRD/surface/rh_triangles_low.txt
fi

if [ ! -f $PRD/surface/rh_region_mapping_low_not_corrected.txt ]; then
  echo "generating the right region mapping on the decimated surface"
  # create left the region mapping
  if [ -n "$MATLAB" ]; then
      $MATLAB -r "rl='rh'; run region_mapping.m; quit;" -nodesktop -nodisplay
  else
      sh region_mapping/distrib/run_region_mapping.sh $MCR
  fi
fi

# correct
if [ ! -f $PRD/surface/rh_region_mapping_low.txt ]; then
  echo " correct the right region mapping"
  python correct_region_mapping.py rh
  echo "check right region mapping"
  python check_region_mapping.py rh
fi
###################################### both hemisphere
# prepare final directory
mkdir -p $PRD/$SUBJ_ID
mkdir -p $PRD/$SUBJ_ID/surface

# reunify both region_mapping, vertices and triangles
if [ ! -f $PRD/$SUBJ_ID/surface/region_mapping.txt ]; then
  echo "reunify both region mappings"
  python reunify_both_regions.py
fi

# zip to put in final format
pushd . > /dev/null
cd $PRD/$SUBJ_ID/surface > /dev/null
zip $PRD/$SUBJ_ID/surface.zip vertices.txt triangles.txt -q
cp region_mapping.txt ..
popd > /dev/null

########################### subcortical surfaces
# extract subcortical surfaces 
if [ ! -f $PRD/surface/subcortical/aseg_058_vert.txt ]; then
  echo "generating subcortical surfaces"
  ./aseg2srf -s $SUBJ_ID
  mkdir -p $PRD/surface/subcortical
  cp $FS/$SUBJ_ID/ascii/* $PRD/surface/subcortical
  python list_subcortical.py
fi

########################## build connectivity using mrtrix 3
mkdir -p $PRD/connectivity
mkdir -p $PRD/$SUBJ_ID/connectivity


## preprocessing
# See: http://mrtrix.readthedocs.io/en/0.3.16/workflows/DWI_preprocessing_for_quantitative_analysis.html

# TODO: add HCP data

# if single acquisition  with reversed directions
function mrchoose () {
  choice=$1
  shift
  $@ << EOF
$choice
EOF
}

# TODO detect phase encoding automatically
# handle encoding scheme
if [ ! -f $PRD/connectivity/predwi.mif ]; then 
  if [ "$TOPUP" = "reversed" ]; then
    echo "generate dwi mif file for use with reversed phase encoding"
    echo "(use of fsl topup)"
    # strides are arranged to make volume data contiguous in memory for
    # each voxel
    # float 32 to make data access faster in subsequent commands
    mrchoose 0 mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi_1.mif \
                         -datatype float32 -stride 0,0,0,1 -force \
                         -nthreads "$NB_THREADS"
    mrchoose 1 mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi_2.mif \
                         -datatype float32 -stride 0,0,0,1 -force \
                         -nthreads "$NB_THREADS"
    # check mif files
    if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
        echo "check predwi_*.mif files"
        mrview $PRD/connectivity/predwi_1.mif $PRD/connectivity/predwi_2.mif
    fi
    # recombining PE dir files 
    mrcat $PRD/connectivity/predwi_1.mif $PRD/connectivity/predwi_2.mif \
          $PRD/connectivity/predwi.mif -axis 3 -nthreads "$NB_THREADS"
    mrinfo $PRD/connectivity/predwi.mif \
           -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init \
           -export_pe_table $PRD/connectivity/pe_table -force 
  else
    echo "generate dwi mif file for use without topup (fsl)"
    mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi.mif \
              -export_pe_table $PRD/connectivity/pe_table \
              -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init \
              -datatype float32 -stride 0,0,0,1 -force -nthreads "$NB_THREADS"
  fi
  # check mif file
  if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
      echo "check predwi_*.mif file"
      mrview $PRD/connectivity/predwi.mif 
  fi
fi



# denoising the volumes
if [ ! -f $PRD/connectivity/predwi_denoised.mif ]; then
  # denoising the combined-directions file is preferable to denoising \
  # predwi1 and 2 separately because of a higher no of volumes
  # see: https://github.com/MRtrix3/mrtrix3/issues/747
  echo "denoising dwi data"
  dwidenoise $PRD/connectivity/predwi.mif \
             $PRD/connectivity/predwi_denoised.mif \
             -noise $PRD/connectivity/noise.mif -force -nthreads "$NB_THREADS"
  if [ ! -f $PRD/connectivity/noise_res.mif ]; then
    # calculate residuals noise
    mrcalc $PRD/connectivity/predwi.mif \
           $PRD/connectivity/predwi_denoised.mif \
           -subtract $PRD/connectivity/noise_res.mif -nthreads "$NB_THREADS"
    # check noise file: lack of anatomy is a marker of accuracy
    if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
      # noise.mif can also be used for SNR calculation
      echo "check noise/predwi_*_denoised.mif files"
      echo "lack of anatomy in noise_res is a marker of accuracy"
      mrview $PRD/connectivity/predwi.mif \
             $PRD/connectivity/predwi_denoised.mif \
             $PRD/connectivity/noise.mif \
             $PRD/connectivity/noise_res.mif  
    fi
  fi
fi

# topup/eddy corrections
if [ ! -f $PRD/connectivity/predwi_denoised_preproc.mif ]
then
  if [ "$TOPUP" = "reversed" ] || [ "$TOPUP" = "eddy_correct" ]; then
    # eddy maybe topup corrections depending of the encoding scheme
    echo "apply eddy and maybe topup"
    dwipreproc $PRD/connectivity/predwi_denoised.mif \
               $PRD/connectivity/predwi_denoised_preproc.mif \
               -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final \
               -rpe_header -eddy_options ' --repol' -cuda -force \
               -nthreads "$NB_THREADS"    
    if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
        echo "check topup/eddy preprocessed mif file"
        mrview $PRD/connectivity/predwi.mif \
               $PRD/connectivity/predwi_denoised.mif \
               $PRD/connectivity/predwi_denoised_preproc.mif 
    fi
  else # no topup/eddy
    echo "no topup/eddy applied"
    mrconvert $PRD/connectivity/predwi_denoised.mif \
              $PRD/connectivity/predwi_denoised_preproc.mif \
              -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final \
              -force -nthreads "$NB_THREADS"
  fi
  # check preproc files
  if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
    echo "check preprocessed mif file (no topup/no eddy)"
    mrview $PRD/connectivity/predwi.mif \
           $PRD/connectivity/predwi_denoised.mif \
           $PRD/connectivity/predwi_denoised_preproc.mif 
  fi
fi

# TOCHECK: Masking step before or after biascorrect?
# Native-resolution mask creation
if [ ! -f $PRD/connectivity/mask_native.mif ]; then
  echo "create dwi mask"
  dwi2mask $PRD/connectivity/predwi_denoised_preproc.mif \
           $PRD/connectivity/mask_native.mif -nthreads "$NB_THREADS"
  # check mask file
  if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
    echo "check native mask mif file"
    mrview $PRD/connectivity/predwi_denoised_preproc.mif \
           -overlay.load $PRD/connectivity/mask_native.mif \
           -overlay.opacity 0.5
  fi
fi

# Bias field correction
if [ ! -f $PRD/connectivity/predwi_denoised_preproc_bias.mif ]; then
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
  # check bias field correction
  if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
    echo "check native mask mif file"
    mrview $PRD/connectivity/predwi.mif \
           $PRD/connectivity/predwi_denoised_preproc.mif \
           $PRD/connectivity/predwi_denoised_preproc_bias.mif 
  fi
fi

# TOCHECK: why not upsampling to vox=1.25 as recommended in mrtrix?
# upsampling and reorienting a la fsl
# reorienting from DiCOM to FSL, RAS to LAS, means -stride -1,+2,+3,+4
# see: http://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html
# upsampling (Dyrby TB. Neuroimage. 2014 Dec;103:202-13.) can help registration
# with structural and is common with mrtrix3 fixel analysis pipeline
# see: http://community.mrtrix.org/t/upsampling-dwi-vs-tckgen-defaults/998/2
if [ ! -f $PRD/connectivity/dwi.mif ]; then
  echo "upsample dwi"
  mrresize $PRD/connectivity/predwi_denoised_preproc_bias.mif - -scale 2 -force | \
  mrconvert - -datatype float32 -stride -1,+2,+3,+4 $PRD/connectivity/dwi.mif -force 
fi

if [ ! -f $PRD/connectivity/mask.mif ]; then
  # for dwi2fod step, a permissive, dilated mask can be used to minimize
  # streamline premature termination, see BIDS protocol: 
  # https://github.com/BIDS-Apps/MRtrix3_connectome/blob/master/run.py
  echo "upsample mask"
  mrresize $PRD/connectivity/mask_native.mif - -scale 2 -force | \
  mrconvert - $PRD/connectivity/mask.mif -datatype bit -stride -1,+2,+3 \
            -force -nthreads "$NB_THREADS"
  maskfilter $PRD/connectivity/mask.mif dilate \
             $PRD/connectivity/mask_dilated.mif -npass 2 -force \
             -nthreads "$NB_THREADS" 
  # check upsampled files
  if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
    echo "check upsampled mif files"
    mrview $PRD/connectivity/dwi.mif \
           -overlay.load $PRD/connectivity/mask.mif \
           -overlay.load $PRD/connectivity/mask_dilated.mif \
           -overlay.opacity 0.5 -norealign 
  fi
fi

## FLIRT registration
# a comparison of registration methods is available in:
# Ou Y, et al. IEEE Trans Med Imaging. 2014 Oct;33(10):2039-65
# other potentials methods for registration TOCHECK include
# FLIRT -bbr (Kerstin's protocol): http://community.mrtrix.org/t/registration-of-structural-and-diffusion-weighted-data/203/8
# bbregister (FS)

# low b extraction to FSL
if [ ! -f $PRD/connectivity/lowb.nii.gz ]; then
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
fi

# generating FSl brain.mgz
if [ ! -f $PRD/connectivity/brain.nii.gz ]; then
  # brain.mgz seems to be superior to diff to T1
  # as the main problem for registration is the wmgm interface that we want to
  # remove and BET stripping is unfortunate in many situations, 
  # and FS pial eddited volumes already present
  # stride from FS to FSL: RAS to LAS
  # see: http://www.grahamwideman.com/gw/brain/fs/coords/fscoords.htm
  echo "generating FSL orientation for masked brain"
  mrconvert $FS/$SUBJ_ID/mri/brain.mgz $PRD/connectivity/brain.nii.gz \
            -datatype float32 -stride -1,+2,+3,+4 -force -nthreads "$NB_THREADS" 
fi

# TOCHECK: test and compare to lowb method, or leave as an option
# # Generate transform image (dwi) for alternative registration method: 
# # replace lowb.nii.gz with output lowb_pseudobrain.nii.gz in the subsequent 
# # registration steps
# # see: Bhushan C, et al. Neuroimage. 2015 Jul 15;115:269-8
# # used in: https://github.com/BIDS-Apps/MRtrix3_connectome/blob/master/run.py
# if [ ! -f $PRD/connectivity/lowb_pseudobrain.nii.gz ]; then
#   echo "extracting b0 vols for registration: pseudostructural"
#   dwiextract $PRD/connectivity/dwi.mif -bzero - \
#   | mrmath - mean - -axis 3 \
#   | mrcalc 1 - -divide $PRD/connectivity/mask_upsampled.mif -multiply - \
#   | mrconvert - - -stride -1,+2,+3 \
#   | mrhistmatch - $PRD/connectivity/brain.nii.gz \
#                 $PRD/connectivity/lowb_pseudobrain.nii.gz
#   if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
#     echo "check pseudo lowb files"
#     mrview $PRD/connectivity/lowb_pseudobrain.nii.gz \
#            $PRD/connectivity/lowb.nii.gz -overlay.load \
#            $PRD/connectivity/lowb_pseudobrain.nii.gz \
#            -overlay.opacity 0.5 -norealign
#   fi
# fi

# aparc+aseg to FSL
if [ ! -f $PRD/connectivity/aparc+aseg.nii.gz ]; then
  echo "generating FSL orientation for aparc+aseg"
  # stride from FS to FSL: RAS to LAS
  mrconvert $FS/$SUBJ_ID/mri/aparc+aseg.mgz \
            $PRD/connectivity/aparc+aseg.nii.gz -stride -1,+2,+3,+4 -force \
            -nthreads "$NB_THREADS" 
fi

# check orientations
if [ ! -f $PRD/connectivity/aparc+aseg_reorient.nii.gz ]; then
  echo "reorienting the region parcellation"
  "$FSL"fslreorient2std $PRD/connectivity/aparc+aseg.nii.gz \
                  $PRD/connectivity/aparc+aseg_reorient.nii.gz
  # check parcellation to brain.mgz
  if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
    # TODO: mrview discrete colour scheme?
    echo "check parcellation"
    echo "if it's correct, just close the window." 
    echo "Otherwise... well, it should be correct anyway"
    mrview $PRD/connectivity/brain.nii.gz \
           -overlay.load $PRD/connectivity/aparc+aseg_reorient.nii.gz \
           -overlay.opacity 0.5 -norealign
  fi
fi

# aparcaseg to diff by inverser transform
if [ ! -f $PRD/connectivity/aparcaseg_2_diff.nii.gz ]; then
  # 6 dof; see:
  # http://web.mit.edu/fsl_v5.0.8/fsl/doc/wiki/FLIRT(2f)FAQ.html#What_cost_function_and.2BAC8-or_degrees_of_freedom_.28DOF.29_should_I_use_in_FLIRT.3F
  echo "register aparc+aseg to diff"
  "$FSL"flirt -in $PRD/connectivity/lowb.nii.gz \
              -ref $PRD/connectivity/brain.nii.gz \
              -omat $PRD/connectivity/diffusion_2_struct.mat \
              -out $PRD/connectivity/lowb_2_struct.nii.gz -dof 6 \
              -searchrx -180 180 -searchry -180 180 -searchrz -180 180 \
              -cost mutualinfo
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
if [ ! -f $PRD/connectivity/brain_2_diff.nii.gz ]; then
  echo "register brain to diff"
  mrtransform $PRD/connectivity/brain.nii.gz \
              $PRD/connectivity/brain_2_diff.nii.gz \
              -linear $PRD/connectivity/diffusion_2_struct_mrtrix.txt \
              -inverse -force 
  # check parcellation to diff
  if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]; then
    echo "check parcellation registration to diffusion space"
    echo "if it's correct, just close the window."
    echo "Otherwise you will have to do the registration by hand"
    mrview $PRD/connectivity/brain_2_diff.nii.gz \
           $PRD/connectivity/lowb.nii.gz \
           -overlay.load $PRD/connectivity/aparcaseg_2_diff.nii.gz \
           -overlay.opacity 0.5 -norealign
  fi
fi


# prepare file for act
if [ "$ACT" = "yes" ] && [ ! -f $PRD/connectivity/act.mif ]; then
  echo "prepare files for act"
  5ttgen fsl $PRD/connectivity/brain_2_diff.nii.gz $PRD/connectivity/act.mif \
         -premasked -force  -nthreads "$NB_THREADS"
  if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
  then
      echo "check tissue segmented image"
      5tt2vis $PRD/connectivity/act.mif $PRD/connectivity/act_vis.mif -force \
              -nthreads "$NB_THREADS"
      mrview $PRD/connectivity/act_vis.mif -colourmap 4
  fi
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
  if [ ! -f $PRD/connectivity/response_wm.txt ]; then
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
      if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]; then
          echo "check ODF image iwth the selected voxels"
          mrview $PRD/connectivity/meanlowb.mif \
                 -overlay.load $PRD/connectivity/RF_voxels.mif \
                 -overlay.opacity 0.5
      fi
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
  if [ ! -f $PRD/connectivity/response_wm.txt ]; then
    echo "estimating response using dhollander algorithm"
    # dwi2response tournier $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt -force -voxels $PRD/connectivity/RF_voxels.mif -mask $PRD/connectivity/mask.mif
    dwi2response dhollander $PRD/connectivity/dwi.mif \
                 $PRD/connectivity/response_wm.txt \
                 $PRD/connectivity/response_gm.txt \
                 $PRD/connectivity/response_csf.txt \
                 -voxels $PRD/connectivity/RF_voxels.mif \
                 -mask $PRD/connectivity/mask.mif -force \
                 -nthreads "$NB_THREADS"
    if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]; then
      echo "check ODF image"
      mrview $PRD/connectivity/meanlowb.mif \
             -overlay.load $PRD/connectivity/RF_voxels.mif \
             -overlay.opacity 0.5
    fi
  fi
fi

# Fibre orientation distribution estimation
if [ ! -f $PRD/connectivity/wm_CSD$lmax.mif ]; then
  # Both for multishell and single shell since we use dhollander in the 
  # single shell case
  # see: http://community.mrtrix.org/t/wm-odf-and-response-function-with-dhollander-option---single-shell-versus-multi-shell/572/4
  echo "calculating fod on multishell or single shell data"
  dwi2fod msmt_csd $PRD/connectivity/dwi.mif \
          $PRD/connectivity/response_wm.txt \
          $PRD/connectivity/wm_CSD$lmax.mif \
          $PRD/connectivity/response_gm.txt \
          $PRD/connectivity/gm_CSD$lmax.mif \
          $PRD/connectivity/response_csf.txt \
          $PRD/connectivity/csf_CSD$lmax.mif \
          -mask $PRD/connectivity/mask_dilated.mif -force \
          -nthreads "$NB_THREADS"
  if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]; then
    echo "check ODF image"
    mrconvert $PRD/connectivity/wm_CSD$lmax.mif - -coord 3 0 \
    -nthreads "$NB_THREADS" \
    | mrcat $PRD/connectivity/csf_CSD$lmax.mif \
            $PRD/connectivity/gm_CSD$lmax.mif - \
            $PRD/connectivity/tissueRGB.mif -axis 3 -nthreads "$NB_THREADS"
    mrview $PRD/connectivity/tissueRGB.mif \
           -odf.load_sh $PRD/connectivity/wm_CSD$lmax.mif 
  fi
fi

# tractography
if [ ! -f $PRD/connectivity/whole_brain.tck ]; then
  if [ "$SIFT" = "sift"]; then
    # temporarily change number of tracks for sift
    number_tracks=$(($NUMBER_TRACKS*$SIFT_MULTIPLIER))
  fi
  native_voxelsize=$(mrinfo $PRD/connectivity/mask_native.mif -vox \
                   | cut -f 1 -d " " | xargs printf "%1.f")
  stepsize=$( bc -l <<< "scale=2; "$native_voxelsize"/2" )
  angle=$( bc -l <<< "scale=2; 90*"$stepsize"/"$native_voxelsize"" )
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
      tckgen $PRD/connectivity/wm_CSD"$lmax".mif \
             $PRD/connectivity/whole_brain.tck \
             -seed_gmwmi $PRD/connectivity/gmwmi_mask.mif 
             -act $PRD/connectivity/act.mif -select "$NUMBER_TRACKS" \
             -seed_unidirectional -crop_at_gmwmi -backtrack \
             -minlength 4 -maxlength 250 -step "$stepsize" -angle "$angle" \
             -cutoff 0.06 -force -nthreads "$NB_THREADS"
    elif [ "$SEED" = "dynamic" ]; then
       # -dynamic seeding may work slightly better than gmwmi, 
       # see Smith RE Neuroimage. 2015 Oct 1;119:338-51.
      echo "seeding dynamically"   
      tckgen $PRD/connectivity/wm_CSD"$lmax".mif \
             $PRD/connectivity/whole_brain.tck \
             -seed_dynamic $PRD/connectivity/wm_CSD$lmax.mif \
             -act $PRD/connectivity/act.mif -select "$NUMBER_TRACKS" \
             -crop_at_gmwmi -backtrack -minlength 4 -maxlength 250 \
             -step "$stepsize" -angle "$angle" -cutoff 0.06 -force \
             -nthreads "$NB_THREADS"
    fi  
  else
    echo "generating tracks without using act" 
    echo "seeding dynamically" 
    tckgen $PRD/connectivity/wm_CSD"$lmax".mif \
           $PRD/connectivity/whole_brain.tck \
           -seed_dynamic $PRD/connectivity/wm_CSD"$lmax".mif \
           -mask $PRD/connectivity/mask.mif -select "$NUMBER_TRACKS" \
           -maxlength 250 -step "$stepsize" -angle "$angle" -cutoff 1  \
           -force -nthreads "$NB_THREADS"
  fi
fi

# postprocessing
if [ ! -f $PRD/connectivity/whole_brain_post.tck ]; then
  if [ "$SIFT" = "sift" ]; then
    echo "using sift"
    number_tracks=$(($NUMBER_TRACKS/$SIFT_MULTIPLIER))
    if [ "$ACT" = "yes" ]; then
        echo "trimming tracks using sift/act" 
        tcksift $PRD/connectivity/whole_brain.tck \
                $PRD/connectivity/wm_CSD"$lmax".mif \
                $PRD/connectivity/whole_brain_post.tck \
                -act $PRD/connectivity/act.mif \
                -out_mu $PRD/connectivity/mu.txt \
                -term_number $NUMBER_TRACKS -fd_scale_gm -force \
                -nthreads "$NB_THREADS"
    else
        echo "trimming tracks using sift/without act" 
        tcksift $PRD/connectivity/whole_brain.tck \
                $PRD/connectivity/wm_CSD"$lmax".mif \
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
               $PRD/connectivity/wm_CSD"$lmax".mif \
               $PRD/connectivity/streamline_weights.csv\
               -act $PRD/connectivity/act.mif \
               -out_mu $PRD/connectivity/mu.txt \
               -out_coeffs $PRD/connectivity/streamline_coeffs.csv \
               -fd_scale_gm -force -nthreads "$NB_THREADS"
    else
      tcksift2 $PRD/connectivity/whole_brain.tck \
               $PRD/connectivity/wm_CSD"$lmax".mif \
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
if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif ]; then
  echo " compute FS labels"
  labelconvert $PRD/connectivity/aparcaseg_2_diff.nii.gz \
               $FREESURFER_HOME/FreeSurferColorLUT.txt \
               fs_region.txt $PRD/connectivity/aparcaseg_2_diff_fs.mif \
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
                $PRD/connectivity/brain_2_diff.nii.gz fs_region.txt \
                $PRD/connectivity/aparcaseg_2_diff_fsl.mif -premasked \
                -force -nthreads "$NB_THREADS"   
  fi 
fi

if [ ! -f $PRD/connectivity/weights.csv ]; then
  echo "compute connectivity matrix weights"
  if [ "$SIFT" = "sift2" ]; then
    # -tck_weights_in flag only needed for sift2 but not for sift/no processing
    # TOCHECK:  mrtrix3 currently generates upper_triangular weights 
    # matrices, need to add -symmetric flag if needed, also -zero_diagonal 
    # if needed (did not see that in the original code)
    # I think I do the symmetric in the compute_connectivity files.py
    # diagonal we want to keep it
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

if [ ! -f $PRD/connectivity/tract_lengths.csv ]; then
  echo "compute connectivity matrix edge lengths"
  # TODO: I don't think it makes sense for the length to use the SIFT2 weighting
  #if [ "$SIFT" = "sift2" ]; then
  #  echo "df"
  #  # mean length result: weight by the length, then average
  #  # see: http://community.mrtrix.org/t/tck2connectome-edge-statistic-sift2-questions/1059/2 
  #  # TOCHECK: be careful when applying sift2, as here the mean is 
  #  # sum(streamline length * streamline weight)/no streamlines, a bit more
  #  # fuzzy to interpret than with sift, however left it as option
  #  tck2connectome $PRD/connectivity/whole_brain_post.tck \
  #                 $PRD/connectivity/aparcaseg_2_diff_$ASEG.mif \
  #                 $PRD/connectivity/tract_lengths.csv \
  #                 -tck_weights_in $PRD/connectivity/streamline_weights.csv \
  #                 -assignment_radial_search 2 -zero_diagonal -scale_length \
  #                 -stat_edge mean -force -nthreads "$NB_THREADS"
  #else
  tck2connectome $PRD/connectivity/whole_brain_post.tck \
                 $PRD/connectivity/aparcaseg_2_diff_"$ASEG".mif \
                 $PRD/connectivity/tract_lengths.csv \
                 -assignment_radial_search 2 -zero_diagonal -scale_length \
                 -stat_edge mean -force -nthreads "$NB_THREADS"
  #fi
fi

# view connectome
if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
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
  # TOCHECK: in mrview, load the lut table (fs_region.txt) for node correspondence, 
  # and exemplars.tck if wanting to see edges as streamlines 
  mrview $PRD/connectivity/aparcaseg_2_diff_$ASEG.mif \
         -connectome.init $PRD/connectivity/aparcaseg_2_diff_$ASEG.mif \
         -connectome.load $PRD/connectivity/weights.csv 
fi

# view tractogram and tdi
if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
  echo "view tractogram and tdi image"
  if [ ! -f $PRD/connectivity/whole_brain_post_decimated.tck ]; then
    # $(( number_tracks/100)) this follows some recommendations by 
    # JD Tournier to avoid mrview to be to slow
    # (for visualization no more than 100-200K streamlines)
    # => min(100k, number_tracks/100)
    if [ "$SIFT" = "sift2" ]; then
        tckedit $PRD/connectivity/whole_brain_post.tck \
                $PRD/connectivity/whole_brain_post_decimated.tck \
                -tck_weights_in $PRD/connectivity/streamline_weights.csv \
                -number $(($NUMBER_TRACKS<100000?$NUMBER_TRACKS:100000))
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
cp cortical.txt $PRD/$SUBJ_ID/connectivity/cortical.txt

# compute centers, areas and orientations
if [ ! -f $PRD/$SUBJ_ID/connectivity/weights.txt ]; then
  echo "generate useful files for TVB"
  native_voxelsize=$(mrinfo $PRD/connectivity/mask_native.mif -vox \
                 | cut -f 1 -d " " | xargs printf "%1.f")
  export stepsize=$( bc -l <<< "scale=2; "$native_voxelsize"/2" )
  python compute_connectivity_files.py
fi

# zip to put in final format
pushd . > /dev/null
cd $PRD/$SUBJ_ID/connectivity > /dev/null
zip $PRD/$SUBJ_ID/connectivity.zip areas.txt average_orientations.txt \
  weights.txt tract_lengths.txt cortical.txt centres.txt -q
popd > /dev/null 



# Done 
read -p "Press [Enter] key to continue..." 

# TODO : update sub parcellations for mrtrix3
###################################################
# compute sub parcellations connectivity if asked
if [ -n "$K_LIST" ]; then
  for K in $K_LIST; do
    export curr_K=$(( 2**K ))
    echo $curr_K
    mkdir -p $PRD/$SUBJ_ID/connectivity_"$curr_K"
    if [ -n "$MATLAB" ]; then
      if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii.gz ]; then
      $MATLAB -r "run subparcel.m; quit;" -nodesktop -nodisplay 
      gzip $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii
      fi
    else
      if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii.gz ]
      then
      sh subparcel/distrib/run_subparcel.sh $MCR  
      gzip $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii
      fi
    fi
    if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif ]
    then
      labelconfig $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii.gz $PRD/connectivity/corr_mat_"$curr_K".txt $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif  -lut_basic $PRD/connectivity/corr_mat_"$curr_K".txt
    fi
    if [ ! -f $PRD/connectivity/weights_$curr_K.csv ]
    then
      echo "compute connectivity sub matrix using act"
      tck2connectome $PRD/connectivity/whole_brain_post.tck $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif $PRD/connectivity/weights_"$curr_K".csv -assignment_radial_search 2
      tck2connectome  $PRD/connectivity/whole_brain_post.tck $PRD/connectivity/aparcaseg_2_diff_"$curr_K".mif $PRD/connectivity/tract_lengths_"$curr_K".csv -metric meanlength -assignment_radial_search 2 -zero_diagonal 
    fi
    if [ ! -f $PRD/$SUBJ_ID/connectivity_"$curr_K"/weights.txt ]
    then
      echo "generate files for TVB subparcellations"
      python compute_connectivity_sub.py $PRD/connectivity/weights_"$curr_K".csv $PRD/connectivity/tract_lengths_"$curr_K".csv $PRD/$SUBJ_ID/connectivity_"$curr_K"/weights.txt $PRD/$SUBJ_ID/connectivity_"$curr_K"/tract_lengths.txt
    fi
    pushd . > /dev/null
    cd $PRD/$SUBJ_ID/connectivity_"$curr_K" > /dev/null
    zip $PRD/$SUBJ_ID/connectivity_"$curr_K".zip weights.txt tract_lengths.txt centres.txt average_orientations.txt -q 
    popd > /dev/null
  done
fi

######################## compute MEG and EEG forward projection matrices
# make BEM surfaces
if [ ! -h ${FS}/${SUBJ_ID}/bem/inner_skull.surf ]; then
  echo "generating bem surfaces"
  mne_watershed_bem --subject ${SUBJ_ID}
  ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_inner_skull_surface ${FS}/${SUBJ_ID}/bem/inner_skull.surf
  ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_outer_skin_surface  ${FS}/${SUBJ_ID}/bem/outer_skin.surf
  ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_outer_skull_surface ${FS}/${SUBJ_ID}/bem/outer_skull.surf
fi

# export to ascii
if [ ! -f ${FS}/${SUBJ_ID}/bem/inner_skull.asc ]; then
  echo "importing bem surface from freesurfer"
  mris_convert $FS/$SUBJ_ID/bem/inner_skull.surf $FS/$SUBJ_ID/bem/inner_skull.asc
  mris_convert $FS/$SUBJ_ID/bem/outer_skull.surf $FS/$SUBJ_ID/bem/outer_skull.asc
  mris_convert $FS/$SUBJ_ID/bem/outer_skin.surf $FS/$SUBJ_ID/bem/outer_skin.asc
fi

# triangles and vertices bem
if [ ! -f $PRD/$SUBJ_ID/surface/inner_skull_vertices.txt ]; then
  echo "extracting bem vertices and triangles"
  python extract_bem.py inner_skull 
  python extract_bem.py outer_skull 
  python extract_bem.py outer_skin 
fi

if [ ! -f ${FS}/${SUBJ_ID}/bem/${SUBJ_ID}-head.fif ]; then
  echo "generating head bem"
  mkheadsurf -s $SUBJ_ID
  mne_surf2bem --surf ${FS}/${SUBJ_ID}/surf/lh.seghead --id 4 --check --fif ${FS}/${SUBJ_ID}/bem/${SUBJ_ID}-head.fif 
fi

if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]; then
  echo "check bem surfaces"
  freeview -v ${FS}/${SUBJ_ID}/mri/T1.mgz -f ${FS}/${SUBJ_ID}/bem/inner_skull.surf:color=yellow:edgecolor=yellow ${FS}/${SUBJ_ID}/bem/outer_skull.surf:color=blue:edgecolor=blue ${FS}/${SUBJ_ID}/bem/outer_skin.surf:color=red:edgecolor=red
fi

# Setup BEM
if [ ! -f ${FS}/${SUBJ_ID}/bem/*-bem.fif ]; then
  worked=0
  outershift=0
  while [ "$worked" == 0 ]; do
    echo "try generate forward model with 0 shift"
    worked=1
    mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 --outershift $outershift || worked=0 
    if [ "$worked" == 0 ]; then
      echo "try generate foward model with 1 shift"
      worked=1
      mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 --outershift 1 || worked=0 
    fi
    if [ "$worked" == 0 ] && [ "$CHECK" = "yes" ]; then
      echo 'you can try using a different shifting value for outer skull, please enter a value in mm'
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

