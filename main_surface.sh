#!/usr/bin/env bash

# TODO: add explicits echo for what to check in the figures
######### import config
while getopts ":c:" opt; do
     case $opt in
     c)
         export CONFIG=$OPTARG
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

if [ ! -n "$CONFIG" ]
then
    echo "you must provide a config file"
    exit 1
fi

if [ ! -n "$number_tracks" ]
then
    echo "config file not correct"
    exit 1
fi

if [ ! -n "$SUBJECTS_DIR" ]
then
    echo "you have to set the SUBJECTS_DIR environnement
    variable for FreeSurfer"
    exit 1
else
    export FS=$SUBJECTS_DIR
fi

######### build cortical surface and region mapping
if [ ! -f $PRD/data/T1/T1.nii.gz ]
then
    echo "generating T1 from DICOM"
    mrconvert $PRD/data/T1/ $PRD/data/T1/T1.nii.gz
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
if [ ! -f $PRD/surface/lh.pial.asc ]
then
    echo "importing left pial surface from freesurfer"
    mris_convert $FS/$SUBJ_ID/surf/lh.pial $PRD/surface/lh.pial.asc
    # take care of the c_(ras) shift which is not done by FS (thks FS!)
    mris_info $FS/$SUBJ_ID/surf/lh.pial >& $PRD/surface/lhinfo.txt
fi

# triangles and vertices high
if [ ! -f $PRD/surface/lh_vertices_high.txt ]
then
    echo "extracting left vertices and triangles"
    python extract_high.py lh
fi

# decimation using brainvisa
if [ ! -f $PRD/surface/lh_vertices_low.txt ]
then
    echo "left decimation using remesher"
    # -> to mesh
    python txt2off.py $PRD/surface/lh_vertices_high.txt $PRD/surface/lh_triangles_high.txt $PRD/surface/lh_high.off
    #  decimation
    ./remesher/cmdremesher/cmdremesher $PRD/surface/lh_high.off $PRD/surface/lh_low.off
    # export to list vertices triangles
    python off2txt.py $PRD/surface/lh_low.off $PRD/surface/lh_vertices_low.txt $PRD/surface/lh_triangles_low.txt
fi

# create left the region mapping
if [ ! -f $PRD/surface/lh_region_mapping_low_not_corrected.txt ]
then
    echo "generating the left region mapping on the decimated surface"
    if [ -n "$matlab" ]
    then
        $matlab -r "rl='lh';run region_mapping.m; quit;" -nodesktop -nodisplay
    else
        sh region_mapping/distrib/run_region_mapping.sh $MCR
    fi
fi

# correct
if [ ! -f $PRD/surface/lh_region_mapping_low.txt ]
then
    echo "correct the left region mapping"
    python correct_region_mapping.py lh
    echo "check left region mapping"
    python check_region_mapping.py lh
fi

###################################### right hemisphere
# export pial into text file
if [ ! -f $PRD/surface/rh.pial.asc ]
then
    echo "importing right pial surface from freesurfer"
    mris_convert $FS/$SUBJ_ID/surf/rh.pial $PRD/surface/rh.pial.asc
    # take care of the c_(ras) shift which is not done by FS (thks FS!)
    mris_info $FS/$SUBJ_ID/surf/rh.pial >& $PRD/surface/rhinfo.txt
fi

# triangles and vertices high
if [ ! -f $PRD/surface/rh_vertices_high.txt ]
then
    echo "extracting right vertices and triangles"
    python extract_high.py rh
fi

# decimation using brainvisa
if [ ! -f $PRD/surface/rh_vertices_low.txt ]
then
    echo "right decimation using remesher"
    # -> to mesh
    python txt2off.py $PRD/surface/rh_vertices_high.txt $PRD/surface/rh_triangles_high.txt $PRD/surface/rh_high.off
    #  decimation
    ./remesher/cmdremesher/cmdremesher $PRD/surface/rh_high.off $PRD/surface/rh_low.off
    # export to list vertices triangles
    python off2txt.py $PRD/surface/rh_low.off $PRD/surface/rh_vertices_low.txt $PRD/surface/rh_triangles_low.txt
fi

if [ ! -f $PRD/surface/rh_region_mapping_low_not_corrected.txt ]
then
    echo "generating the right region mapping on the decimated surface"
    # create left the region mapping
    if [ -n "$matlab" ]
    then
        $matlab -r "rl='rh'; run region_mapping.m; quit;" -nodesktop -nodisplay
    else
        sh region_mapping/distrib/run_region_mapping.sh $MCR
    fi
fi

# correct
if [ ! -f $PRD/surface/rh_region_mapping_low.txt ]
then
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
if [ ! -f $PRD/$SUBJ_ID/surface/region_mapping.txt ]
then
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
if [ ! -f $PRD/surface/subcortical/aseg_058_vert.txt ]
then
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

# if single acquisition  with reversed directions
function mrchoose () {
    choice=$1
    shift
    $@ << EOF
$choice
EOF
}

# handle encoding scheme
if [ "$topup" = "reversed" ]
then
    if [ ! -f $PRD/connectivity/predwi_1.mif ] || [ ! -f $PRD/connectivity/predwi_2.mif ]
    then # TODO: find a way for HCP dataset
         mrchoose 0 mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi_1.mif \
                              -datatype float32 -stride 0,0,0,1 -force
         mrchoose 1 mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi_2.mif \
                              -force -datatype float32 -stride 0,0,0,1 
        # check mif files
        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
        then
            echo "check predwi_*.mif files"
            mrview $PRD/connectivity/predwi_1.mif $PRD/connectivity/predwi_2.mif
        fi
    fi
    # recombining PE dir files 
    if [ ! -f $PRD/connectivity/predwi.mif ]
    then
        mrcat $PRD/connectivity/predwi_1.mif $PRD/connectivity/predwi_2.mif \
              $PRD/connectivity/predwi.mif -axis 3
        mrinfo $PRD/connectivity/predwi.mif \
               -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init \
               -export_pe_table $PRD/connectivity/pe_table -force 
    fi
else
    echo "generate dwi mif file for use without topup (fsl)"
    if [ ! -f $PRD/connectivity/predwi.mif ]
    then # TODO: check the strides
        mrconvert $PRD/data/DWI/ $PRD/connectivity/predwi.mif\
                  -export_pe_table $PRD/connectivity/pe_table 
                  -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_init
                  -datatype float32 -stride 0,0,0,1 -force
        # check mif file
        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
        then
            echo "check predwi_*.mif file"
            mrview $PRD/connectivity/predwi.mif 
        fi
    fi
fi

# denoising the volumes
if [ ! -f $PRD/connectivity/predwi_denoised.mif ]
then # denoising the combined-directions file is preferable to denoising \
     # predwi1 and 2 separately because of a higher no of volumes
    dwidenoise $PRD/connectivity/predwi.mif \
               $PRD/connectivity/predwi_denoised.mif 
               -noise $PRD/connectivity/noise.mif -force
    if [ ! -f $PRD/connectivity/noise_res.mif ]
    then # calculate residuals noise
        mrcalc $PRD/connectivity/predwi.mif \
               $PRD/connectivity/predwi_denoised.mif \
               -subtract $PRD/connectivity/noise_res.mif
        # check noise file: lack of anatomy is a marker of accuracy
        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
        then # noise.mif can also be used for SNR calculation
            echo "check noise/predwi_*_denoised.mif files"
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
    if [ "$topup" = "reversed" ]
    then # topup and eddy corrections
         # TOCHECK: use this hacked version (repol)necesarily if having an equal no. of forward
         # and reverse volumes, otherwise use dwipreproc command
         # TODO: compare with 
         # http://mrtrix.readthedocs.io/en/latest/dwi_preprocessing/dwipreproc.html#dwipreproc-page
        echo "apply topup/eddy"
        dwipreproc $PRD/connectivity/predwi_denoised.mif \
                   $PRD/connectivity/predwi_denoised_preproc.mif \
                   -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final \
                   -rpe_header -eddy_options ' --repol' -cuda -force    
        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
        then
            echo "check topup/eddy preprocessed mif file"
            mrview $PRD/connectivity/predwi.mif \
                   $PRD/connectivity/predwi_denoised.mif \
                   $PRD/connectivity/predwi_denoised_preproc.mif 
        fi
    elif [ "$topup" = "eddy_correct" ]
    then # eddy only 
         # TODO: check that the header contains non for rpe encoding
        echo "use eddy (fsl); no topup applied"
        dwipreproc $PRD/connectivity/predwi_denoised.mif \
                   $PRD/connectivity/predwi_denoised_preproc.mif 
                   -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final \
                   -rpe_header -eddy_options ' --repol' -cuda -force
        # check preproc files
        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
        then
            echo "check eddy (no topup) preprocessed mif file"
            mrview $PRD/connectivity/predwi.mif \
                   $PRD/connectivity/predwi_denoised.mif \
                   $PRD/connectivity/predwi_denoised_preproc.mif 
        fi
    else # no topup/eddy
        echo "no topup/eddy applied"
        mrconvert $PRD/connectivity/predwi_denoised.mif \
                  $PRD/connectivity/predwi_denoised_preproc.mif \
                  -export_grad_mrtrix $PRD/connectivity/bvecs_bvals_final -force
        # check preproc files
        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
        then
            echo "check preprocessed mif file (no topup/no eddy)"
            mrview $PRD/connectivity/predwi.mif \
                   $PRD/connectivity/predwi_denoised.mif \
                   $PRD/connectivity/predwi_denoised_preproc.mif 
        fi
    fi
fi

# Native-resolution mask creation
if [ ! -f $PRD/connectivity/mask_native.mif ]
then
    echo "create dwi mask"
    dwi2mask $PRD/connectivity/predwi_denoised_preproc.mif \
             $PRD/connectivity/mask_native.mif
###    check mask file
    if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
    then
        echo "check native mask mif file"
        mrview $PRD/connectivity/predwi_denoised_preproc.mif \
               -overlay.load $PRD/connectivity/mask_native.mif \
               -overlay.opacity 0.5
    fi
fi

# Bias field correction
if [ ! -f $PRD/connectivity/predwi_denoised_preproc_bias.mif ]
then
    echo "bias correct"
    # ANTS seems better than FSL
    # see http://mrtrix.readthedocs.io/en/0.3.16/workflows/DWI_preprocessing_for_quantitative_analysis.html
    # TODO: check if ANTS is installed, otherwise use FSL
    dwibiascorrect $PRD/connectivity/predwi_denoised_preproc.mif \
                   $PRD/connectivity/predwi_denoised_preproc_bias.mif \
                   -mask $PRD/connectivity/mask_native.mif \
                   -bias $PRD/connectivity/B1_bias.mif -ants -force
fi

## FLIRT registration
# low b extraction to FSL
if [ ! -f $PRD/connectivity/lowb.nii.gz ]
then
    echo "extracting b0 vols for registration"
    dwiextract $PRD/connectivity/dwi.mif $PRD/connectivity/lowb.mif \
               -bzero -force 
    # note: stride from mrtrix to FSL, RAS to LAS
    # see: http://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html
    mrconvert $PRD/connectivity/lowb.mif $PRD/connectivity/lowb.nii.gz \
              -stride -1,+2,+3,+4
    # for visualization 
    mrmath $PRD/connectivity/lowb.mif mean $PRD/connectivity/meanlowb.mif \
           -axis 3 -force
fi

# generating FSl brain.mgz
if [ ! -f $PRD/connectivity/T1.nii.gz ]
then # note: brain.mgz seems to be superior to diff to T1
     # as BET stripping is unfortunate in many situations, and FS pial eddited volumes already present
     # TODO: ref? T1 option?
     # note: stride from FS to FSL: RAS to LAS
     # see: http://www.grahamwideman.com/gw/brain/fs/coords/fscoords.htm
    echo "generating FSL orientation for masked brain"
    mrconvert $FS/$SUBJ_ID/mri/brain.mgz $PRD/connectivity/brain.nii.gz
              -datatype float32 -stride -1,+2,+3,+4 -force 
fi

# TODO
## Generate transform image (dwi) for alternative registration method: replace lowb.nii.gz with output lowb_pseudobrain.nii.gz in the subsequent registration steps
##    if [ ! -f $PRD/connectivity/lowb_pseudobrain.nii.gz ]
##    then
##        echo "extracting b0 vols for registration: pseudostructural"
##        dwiextract $PRD/connectivity/dwi.mif -bzero - | mrmath - mean - -axis 3 | mrcalc 1 - -divide $PRD/connectivity/mask_upsampled.mif -multiply - | mrconvert - - -stride -1,+2,+3 | mrhistmatch - $PRD/connectivity/brain.nii.gz $PRD/connectivity/lowb_pseudobrain.nii.gz
##        if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
##        then
##            echo "check pseudo lowb files"
##            mrview $PRD/connectivity/lowb_pseudobrain.nii.gz $PRD/connectivity/lowb.nii.gz -overlay.load $PRD/connectivity/lowb_pseudobrain.nii.gz -overlay.opacity 0.5 -norealign
##        fi
##    fi

# aparc+aseg to FSL
if [ ! -f $PRD/connectivity/aparc+aseg.nii.gz ]
then
    echo "generating FSL orientation for aparc+aseg"
    # note: stride from FS to FSL: RAS to LAS
    mrconvert $FS/$SUBJ_ID/mri/aparc+aseg.mgz \
              $PRD/connectivity/aparc+aseg.nii.gz -stride -1,+2,+3,+4 
fi

# check orientations
if [ ! -f $PRD/connectivity/aparc+aseg_reorient.nii.gz ]
then
    echo "reorienting the region parcellation"
    fslreorient2std $PRD/connectivity/aparc+aseg.nii.gz \
                    $PRD/connectivity/aparc+aseg_reorient.nii.gz
    # check parcellation to brain.mgz
    if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
    then
        echo "check parcellation"
        echo "if it's correct, just close the window." 
        echo "Otherwise... well, it should be correct anyway"
        mrview $PRD/connectivity/brain.nii.gz \
               -overlay.load $PRD/connectivity/aparc+aseg_reorient.nii.gz \
               -overlay.opacity 0.5 -norealign
    fi
fi

# aparcaseg to diff by inverser transform
if [ ! -f $PRD/connectivity/aparcaseg_2_diff.nii.gz ]
then # TOCHECK:6 dof vs 12 dof
    echo " register aparc+aseg to diff"
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
    mrtransform $PRD/connectivity/aparcaseg_2_diff.nii.gz \
                -linear $PRD/connectivity/diffusion_2_struct_mrtrix.txt \
                -inverse $PRD/connectivity/aparc+aseg_reorient.nii.gz \
                -datatype uint32 -force 
fi

# brain to diff by inverse transform
if [ ! -f $PRD/connectivity/brain_2_diff.nii.gz ]
then
    echo "register brain to diff"
    mrtransform $PRD/connectivity/brain_2_diff.nii.gz \
                -linear $PRD/connectivity/diffusion_2_struct_mrtrix.txt \
                -inverse $PRD/connectivity/brain.nii.gz -force 

    # check parcellation to diff
    if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
    then
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
### TODO: test
if [ "$act" = "yes" ] && [ ! -f $PRD/connectivity/act.mif ]
then
    echo "prepare files for act"
    5ttgen fsl $PRD/connectivity/brain_2_diff.nii.gz $PRD/connectivity/act.mif \
           -premasked -force       
    if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
    then
        echo "check tissue segmented image"
        5tt2vis -force $PRD/connectivity/act.mif $PRD/connectivity/act_vis.mif
        mrview $PRD/connectivity/act_vis.mif
    fi
fi

# Response function estimation
# Check if multi or single shell
shells=$(mrinfo -shells $PRD/connectivity/dwi.mif)
echo "shell b values are $shells"
nshells=($shells)
no_shells=${#nshells[@]}
echo "no of shells are $no_shells"

if [ "$no_shells" -gt 2 ] 
then
# Multishell
    if [ ! -f $PRD/connectivity/response_wm.txt ] 
    then
        if [ "$act" = "yes" ]
        then 
            echo "estimating response using msmt algorithm"
            dwi2response msmt_5tt $PRD/connectivity/dwi.mif \
                         $PRD/connectivity/act.mif \
                         $PRD/connectivity/response_wm.txt \
                         $PRD/connectivity/response_gm.txt \
                         $PRD/connectivity/response_csf.txt \
                         -voxels $PRD/connectivity/RF_voxels.mif \
                         -mask $PRD/connectivity/mask.mif -force
            if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
            then
                echo "check ODF image"
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
                         -mask $PRD/connectivity/mask.mif -force 
        fi
    fi
else
# Single shell only
    if [ ! -f $PRD/connectivity/response.txt ]
    then
        echo "estimating response using dhollander algorithm"
        # dwi2response tournier $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt -force -voxels $PRD/connectivity/RF_voxels.mif -mask $PRD/connectivity/mask.mif
        dwi2response dhollander $PRD/connectivity/dwi.mif \
                     $PRD/connectivity/response_wm.txt \
                     $PRD/connectivity/response_gm.txt \
                     $PRD/connectivity/response_csf.txt \
                     -voxels $PRD/connectivity/RF_voxels.mif \
                     -mask $PRD/connectivity/mask.mif -force 
        if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
        then
            echo "check ODF image"
            mrview $PRD/connectivity/meanlowb.mif \
                   -overlay.load $PRD/connectivity/RF_voxels.mif \
                   -overlay.opacity 0.5
        fi
    fi
fi


# Fibre orientation distribution estimation
if [ ! -f $PRD/connectivity/wm_CSD$lmax.mif ]
then
# Multishell
    if [ "$no_shells" -gt 2 ] 
    then # TOCHECK: lmax?
        echo "calculating fod on multishell data"
        dwi2fod msmt_csd $PRD/connectivity/dwi.mif \
                $PRD/connectivity/response_wm.txt \
                $PRD/connectivity/wm_CSD$lmax.mif \
                $PRD/connectivity/response_gm.txt \
                $PRD/connectivity/gm_CSD$lmax.mif \
                $PRD/connectivity/response_csf.txt \
                $PRD/connectivity/csf_CSD$lmax.mif \
                -mask $PRD/connectivity/mask_dilated.mif -force 
    else
        # Single shell only
        echo "calculating fod on single-shell data"
        # performing msmt_csd on single shell data
        # see: http://community.mrtrix.org/t/msmt-csd-for-single-shell-data/1052
        dwiextract $PRD/connectivity/dwi.mif - 
        | dwi2fod msmt_csd - $PRD/connectivity/response.txt \
                  $PRD/connectivity/wm_CSD$lmax.mif \
                  -mask $PRD/connectivity/mask_dilated.mif -lmax $lmax
    fi
    if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
    then
        if [ "$no_shells" -gt 2 ] 
        then
            echo "check ODF image"
            mrconvert $PRD/connectivity/wm_CSD$lmax.mif - -coord 3 0 \
            | mrcat $PRD/connectivity/csf_CSD$lmax.mif \
                    $PRD/connectivity/gm_CSD$lmax.mif - \
                    $PRD/connectivity/tissueRGB.mif -axis 3
            mrview $PRD/connectivity/tissueRGB.mif \
                   -odf.load_sh $PRD/connectivity/wm_CSD$lmax.mif 
        else
            echo "check ODF image"
            mrview $PRD/connectivity/dwi.mif \
            -odf.load_sh $PRD/connectivity/wm_CSD$lmax.mif
        fi
    fi
fi

# tractography
if [ ! -f $PRD/connectivity/whole_brain.tck ]
then
    if [ "$act" = "yes" ]
    then
        echo "generating tracks using act" 
        if [ "$seed" = "gmwmi" ]
        then
            echo "seeding from gmwmi" 
            5tt2gmwmi $PRD/connectivity/act.mif \
                      $PRD/connectivity/gmwmi_mask.mif -force 
            # TODO: cutoff add not msmt csd back to default?; min length check andreas paper; angle
            tckgen $PRD/connectivity/wm_CSD"$lmax".mif \
                   $PRD/connectivity/whole_brain.tck \
                   -seed_gmwmi $PRD/connectivity/gmwmi_mask.mif 
                   -act $PRD/connectivity/act.mif -select $number_tracks \
                   -seed_unidirectional -crop_at_gmwmi -backtrack \
                   -minlength 4 -maxlength 250 -step 1 -angle 45 -cutoff 0.06 \
                   -force
        else # [ "$seed" = "dynamic" ] default. TODO: check if good without SIFT ; see_unidirectional?
             # -dynamic seeding may work slightly better than gmwmi, 
             # see Smith RE Neuroimage. 2015 Oct 1;119:338-51.
            echo "seeding dynamically"   
            tckgen $PRD/connectivity/wm_CSD"$lmax".mif \
                   $PRD/connectivity/whole_brain.tck \
                   -seed_dynamic $PRD/connectivity/wm_CSD$lmax.mif \
                   -act $PRD/connectivity/act.mif -select $number_tracks \
                   -crop_at_gmwmi -backtrack -minlength 4 -maxlength 250 \
                   -step 1 -angle 45 -cutoff 0.06 -force 
            fi
        fi  
    else
        echo "generating tracks without using act" 
        echo "seeding dynamically" 
        tckgen $PRD/connectivity/wm_CSD"$lmax".mif \
               $PRD/connectivity/whole_brain.tck \
               -seed_dynamic $PRD/connectivity/wm_CSD"$lmax".mif \
               -mask $PRD/connectivity/mask.mif -select $number_tracks \
               -maxlength 250 -step 1 -angle 45 -cutoff 0.06  -force 
        fi
    fi
fi

# postprocessing
if [ ! -f $PRD/connectivity/whole_brain_post.tck ]
then
    if [ "$sift" = "sift2" ] 
    then 
        echo "using sift2"
        ln -s $PRD/connectivity/whole_brain.tck $PRD/connectivity/whole_brain_post.tck
        if [ "$act" = "yes" ]
        then
            echo "using act" 
            tcksift2 $PRD/connectivity/whole_brain.tck \
                     $PRD/connectivity/wm_CSD"$lmax".mif \
                     $PRD/connectivity/streamline_weights.csv\
                     -act $PRD/connectivity/act.mif \
                     -out_mu $PRD/connectivity/mu.txt \
                     -out_coeffs $PRD/connectivity/streamline_coeffs.csv \
                     -fd_scale_gm -force
        else
            tcksift2 $PRD/connectivity/whole_brain.tck \
                     $PRD/connectivity/wm_CSD"$lmax".mif \
                     $PRD/connectivity/streamline_weights.csv \
                     -out_mu $PRD/connectivity/mu.txt \
                     -out_coeffs $PRD/connectivity/streamline_coeffs.csv -force
        fi
    else 
        echo "not using sift2"
        ln -s $PRD/connectivity/whole_brain.tck \
              $PRD/connectivity/whole_brain_post.tck
    fi
fi

####TODO # now compute connectivity and length matrix
if [ ! -f $PRD/connectivity/aparcaseg_2_diff.mif ]
then
    echo " compute labels"
    labelconfig $PRD/connectivity/aparcaseg_2_diff.nii.gz fs_region.txt $PRD/connectivity/aparcaseg_2_diff.mif -lut_freesurfer $FREESURFER_HOME/FreeSurferColorLUT.txt
fi

if [ ! -f $PRD/connectivity/weights.csv ]
then
    echo "compute connectivity matrix"
    tck2connectome $PRD/connectivity/whole_brain_post.tck $PRD/connectivity/aparcaseg_2_diff.mif $PRD/connectivity/weights.csv -assignment_radial_search 2
    tck2connectome $PRD/connectivity/whole_brain_post.tck $PRD/connectivity/aparcaseg_2_diff.mif $PRD/connectivity/tract_lengths.csv -metric meanlength -zero_diagonal -assignment_radial_search 2 
fi

# Compute other files
# we do not compute hemisphere
# subcortical is already done
cp cortical.txt $PRD/$SUBJ_ID/connectivity/cortical.txt

# compute centers, areas and orientations
if [ ! -f $PRD/$SUBJ_ID/connectivity/weights.txt ]
then
    echo " generate useful files for TVB"
    python compute_connectivity_files.py
fi

# zip to put in final format
pushd . > /dev/null
cd $PRD/$SUBJ_ID/connectivity > /dev/null
zip $PRD/$SUBJ_ID/connectivity.zip areas.txt average_orientations.txt weights.txt tract_lengths.txt cortical.txt centres.txt -q
popd > /dev/null 

###################################################
# compute sub parcellations connectivity if asked
if [ -n "$K_list" ]
then
    for K in $K_list
    do
        export curr_K=$(( 2**K ))
        mkdir -p $PRD/$SUBJ_ID/connectivity_"$curr_K"

        if [ -n "$matlab" ]  
        then
            if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii.gz ]
            then
            $matlab -r "run subparcel.m; quit;" -nodesktop -nodisplay 
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
if [ ! -h ${FS}/${SUBJ_ID}/bem/inner_skull.surf ]
then
    echo "generating bem surfaces"
    mne_watershed_bem --subject ${SUBJ_ID}
    ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_inner_skull_surface ${FS}/${SUBJ_ID}/bem/inner_skull.surf
    ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_outer_skin_surface  ${FS}/${SUBJ_ID}/bem/outer_skin.surf
    ln -s ${FS}/${SUBJ_ID}/bem/watershed/${SUBJ_ID}_outer_skull_surface ${FS}/${SUBJ_ID}/bem/outer_skull.surf
fi

# export to ascii
if [ ! -f ${FS}/${SUBJ_ID}/bem/inner_skull.asc ]
then
    echo "importing bem surface from freesurfer"
    mris_convert $FS/$SUBJ_ID/bem/inner_skull.surf $FS/$SUBJ_ID/bem/inner_skull.asc
    mris_convert $FS/$SUBJ_ID/bem/outer_skull.surf $FS/$SUBJ_ID/bem/outer_skull.asc
    mris_convert $FS/$SUBJ_ID/bem/outer_skin.surf $FS/$SUBJ_ID/bem/outer_skin.asc
fi

# triangles and vertices bem
if [ ! -f $PRD/$SUBJ_ID/surface/inner_skull_vertices.txt ]
then
    echo "extracting bem vertices and triangles"
    python extract_bem.py inner_skull 
    python extract_bem.py outer_skull 
    python extract_bem.py outer_skin 
fi

if [ ! -f ${FS}/${SUBJ_ID}/bem/${SUBJ_ID}-head.fif ]
then
    echo "generating head bem"
    mkheadsurf -s $SUBJ_ID
    mne_surf2bem --surf ${FS}/${SUBJ_ID}/surf/lh.seghead --id 4 --check --fif ${FS}/${SUBJ_ID}/bem/${SUBJ_ID}-head.fif 
fi

if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
then
    echo "check bem surfaces"
    freeview -v ${FS}/${SUBJ_ID}/mri/T1.mgz -f ${FS}/${SUBJ_ID}/bem/inner_skull.surf:color=yellow:edgecolor=yellow ${FS}/${SUBJ_ID}/bem/outer_skull.surf:color=blue:edgecolor=blue ${FS}/${SUBJ_ID}/bem/outer_skin.surf:color=red:edgecolor=red
fi

# Setup BEM
if [ ! -f ${FS}/${SUBJ_ID}/bem/*-bem.fif ]
then
    worked=0
    outershift=0
    while [ "$worked" == 0 ]
    do
        echo "try generate forward model with 0 shift"
        worked=1
        mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 --outershift $outershift || worked=0 
        if [ "$worked" == 0 ]
        then
            echo "try generate foward model with 1 shift"
            worked=1
            mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 --outershift 1 || worked=0 
        fi
        if [ "$worked" == 0 ] && [ "$CHECK" = "yes" ]
        then
            echo 'you can try using a different shifting value for outer skull, please enter a value in mm'
            read outershift;
            echo $outershift
        elif [ "$worked" == 0 ]
        then
            echo "bem did not worked"
            worked=1
        elif [ "$worked" == 1 ]
        then
            echo "success!"
        fi
    done
fi

