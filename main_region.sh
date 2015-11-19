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

######### build cortical surface and region mapping
if [ ! -f $PRD/data/T1/T1.nii.gz ]
then
echo "generating T1 from DICOM"
mrconvert $PRD/data/T1/ $PRD/data/T1/T1.nii.gz
fi


# mrtrix
mkdir -p $PRD/connectivity_regions
mkdir -p $PRD/"$SUBJ_ID"_regions/connectivity
# mrconvert
if [ ! -f $PRD/connectivity_regions/dwi.mif ]
then
if [ -f $PRD/data/DWI/*.nii.gz ]
then
ls $PRD/data/DWI/ | grep .nii.gz | xargs -I {} mrconvert $PRD/data/DWI/{} $PRD/connectivity_regions/dwi.mif 
else
mrconvert $PRD/data/DWI/ $PRD/connectivity_regions/dwi.mif
fi
fi
# brainmask 
if [ ! -f $PRD/connectivity_regions/lowb.nii.gz  ]
then
mrmath $PRD/connectivity_regions/dwi.mif -axis 3 mean $PRD/connectivity_regions/lowb.nii.gz
fi
if [ ! -f $PRD/connectivity_regions/mask_not_checked.mif ]
then
mrthreshold -percent $percent_value_mask $PRD/connectivity_regions/lowb.nii.gz - | mrfilter - median  - | mrfilter - median $PRD/connectivity_regions/mask_not_checked.mif
fi

# check the mask
if [ ! -f $PRD/connectivity_regions/mask_checked.mif ]
then
cp $PRD/connectivity_regions/mask_not_checked.mif $PRD/connectivity_regions/mask_checked.mif
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
echo "check the mask and then close the window"
while true; do
mrview $PRD/connectivity_regions/mask_checked.mif
read -p "was the mask good?" yn
case $yn in
[Yy]* ) break;;
[Nn]* ) read -p "enter new threshold value" percent_value_mask; echo $percent_value_mask; rm $PRD/connectivity_regions/mask_checked.mif; 
	mrthreshold -percent $percent_value_mask $PRD/connectivity_regions/lowb.nii.gz - | mrfilter - median - | mrfilter - median $PRD/connectivity_regions/mask_checked.mif;;
 * ) echo "Please answer y or n.";;
esac
done
fi
fi

if [ -f $PRD/connectivity_regions/mask_checked.mif ]
then
cp $PRD/connectivity_regions/mask_checked.mif $PRD/connectivity_regions/mask.mif
elif [ -f $PRD/connectivity_regions/mask_not_checked.mif ]
then
cp $PRD/connectivity_regions/mask_not_checked.mif $PRD/connectivity_regions/mask.mif
fi

# tensor imaging
if [ ! -f $PRD/connectivity_regions/dt.mif ]
then
if [ -f $PRD/data/DWI/*.nii.gz ]
then
dwi2tensor $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/dt.mif -fslgrad $PRD/data/DWI/bvecs $PRD/data/DWI/bvals
else
dwi2tensor $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/dt.mif
fi
fi
if [ ! -f $PRD/connectivity_regions/fa.mif ]
then
tensor2metric $PRD/connectivity_regions/dt.mif -fa - | mrcalc - $PRD/connectivity_regions/mask.mif -mult $PRD/connectivity_regions/fa.mif
fi
if [ ! -f $PRD/connectivity_regions/ev.mif ]
then
tensor2metric $PRD/connectivity_regions/dt.mif -vector - | mrcalc - $PRD/connectivity_regions/fa.mif -mult $PRD/connectivity_regions/ev.mif
fi

#constrained spherical decconvolution
if [ ! -f $PRD/connectivity_regions/sf.mif ]
then
maskfilter $PRD/connectivity_regions/mask.mif -npass 3 erode - | mrcalc $PRD/connectivity_regions/fa.mif - -mult - | mrthreshold - -abs 0.7 $PRD/connectivity_regions/sf.mif
fi
if [ ! -f $PRD/connectivity_regions/response.txt ]
then
if [ -f $PRD/data/DWI/*.nii.gz ]
then
dwi2response $PRD/connectivity_regions/dwi.mif -mask $PRD/connectivity_regions/sf.mif -lmax $lmax $PRD/connectivity_regions/response.txt -fslgrad $PRD/data/DWI/bvecs $PRD/data/DWI/bvals
else
dwi2response $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/sf.mif -lmax $lmax $PRD/connectivity_regions/response.txt
fi
fi
if  [ -n "$DISPLAY" ]  &&  [ "$CHECK" = "yes" ]
then
echo "check the response function"
echo "it should be broadest in the axial plane, and have low amplitude along the z-axis."
shview -response $PRD/connectivity_regions/response.txt
fi

if [ ! -f $PRD/connectivity_regions/CSD$lmax.mif ]
then
if [ -f $PRD/data/DWI/*.nii.gz ]
then
dwi2fod $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/response.txt -lmax $lmax -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/CSD$lmax.mif -fslgrad $PRD/data/DWI/bvecs $PRD/data/DWI/bvals
else
dwi2fod $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/response.txt -lmax $lmax -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/CSD$lmax.mif
fi
fi


# FLIRT registration
# parcellation MNI to T1 using FNIRT
echo "first preregistration T1 to MNI using FLIRT"
#bet $PRD/data/T1/T1.nii $PRD/data/T1/T1_bet.nii
flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in $PRD/data/T1/T1.nii.gz -omat $PRD/connectivity_regions/t1_2_mni_transf.mat
echo "then register T1 to MNI using FNIRT"
fnirt --in=$PRD/data/T1/T1.nii.gz --aff=$PRD/connectivity_regions/t1_2_mni_transf.mat --cout=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf --config=T1_2_MNI152_2mm
echo "apply the warp to check the registration"
applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --in=$PRD/data/T1/T1.nii.gz --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf --out=$PRD/connectivity_regions/t1_2_mni_warped_structural_2mm
echo "inverse the linear registration"
convert_xfm -omat $PRD/connectivity_regions/t1_2_mni_transf_inverse.mat -inverse $PRD/connectivity_regions/t1_2_mni_transf.mat
echo "inverse the warp"
invwarp --ref=$PRD/data/T1/T1.nii.gz --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf.nii.gz --out=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf_inverse
echo "Bring the chosen parcellation to T1-Space"
applywarp --ref=$PRD/data/T1/T1.nii.gz --in=parcellations/$parcel --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf_inverse --out=$PRD/connectivity_regions/region_parcellation

#reorient parcellation to standard orientation to match to T1
if [ ! -f $PRD/connectivity_regions/region_parcellation_reorient.nii.gz ]
then
echo "reorienting the region parcellation"
fslreorient2std $PRD/connectivity_regions/region_parcellation $PRD/connectivity_regions/region_parcellation_reorient
fi

#reorient T1 to standard orientation
if [ ! -f $PRD/connectivity_regions/T1.nii.gz ]
then
echo "reorienting the T1"
fslreorient2std $PRD/data/T1/T1.nii.gz $PRD/connectivity_regions/T1.nii.gz
fi

# check parcellation to T1
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
echo "check parcellation"
echo "if it's correct, just close the window. Otherwise you will have to do the registration by hand"
"$FSL"fslview $PRD/connectivity_regions/T1.nii.gz $PRD/connectivity_regions/region_parcellation_reorient.nii.gz -l "Cool"
fi

if [ ! -f $PRD/connectivity_regions/region_parcellation_2_diff.nii.gz ]
then
	echo "register parcellation to diff"
	"$FSL"flirt -in $PRD/connectivity_regions/lowb.nii.gz -ref $PRD/connectivity_regions/T1.nii.gz -omat $PRD/connectivity_regions/diffusion_2_struct.mat -out $PRD/connectivity_regions/lowb_2_struct.nii.gz -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -cost mutualinfo
	"$FSL"convert_xfm -omat $PRD/connectivity_regions/diffusion_2_struct_inverse.mat -inverse $PRD/connectivity_regions/diffusion_2_struct.mat
	"$FSL"flirt -applyxfm -in $PRD/connectivity_regions/region_parcellation_reorient.nii.gz -ref $PRD/connectivity_regions/lowb.nii.gz -out $PRD/connectivity_regions/region_parcellation_2_diff.nii.gz -init $PRD/connectivity_regions/diffusion_2_struct_inverse.mat -interp nearestneighbour

	if [ ! -f $PRD/connectivity_regions/T1_2_diff.nii.gz ]
    then
        echo " register parcellation to diff"
        "$FSL"flirt -applyxfm -in $PRD/connectivity/T1.nii.gz -ref $PRD/connectivity/lowb.nii.gz -out $PRD/connectivity_regions/T1_2_diff.nii.gz -init $PRD/connectivity_regions/diffusion_2_struct_inverse.mat -interp nearestneighbour
    fi

	# check parcellation to diff
	if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
	then
		echo "check parcellation"
		echo "if it's correct, just close the window. Otherwise you will have to do the registration by hand"
		"$FSL"fslview $PRD/connectivity_regions/lowb.nii.gz $PRD/connectivity_regions/region_parcellation_2_diff -l "Cool"
	fi
fi



#prepare file for act
if [ "$act" = "yes" ] && [ ! -f $PRD/connectivity_regions/act.mif ]
then
echo "prepare files for act"
act_anat_prepare_fsl $PRD/connectivity/T1_2_diff.nii.gz $PRD/connectivity_regions/act.mif
fi

# tractography
if [ ! -f $PRD/connectivity_regions/whole_brain.tck ]
then
	if [ "$act" = "yes" ]
	then
		echo "generating tracks using act"
		5tt2gmwmi $PRD/connectivity_regions/act.mif $PRD/connectivity_regions/gmwmi_mask.mif
		tckgen $PRD/connectivity_regions/CSD$lmax.mif $PRD/connectivity_regions/whole_brain.tck -unidirectional -seed_gmwmi $PRD/connectivity_regions/gmwmi_mask.mif -num $number_tracks -act $PRD/connectivity_regions/act.mif -maxlength 250 -step 0.5
	else
		echo "generating tracks without using act"
		tckgen $PRD/connectivity_regions/CSD$lmax.mif $PRD/connectivity_regions/whole_brain.tck -unidirectional -algorithm iFOD2 -seed_image $PRD/connectivity/mask.mif -mask $PRD/connectivity/mask.mif -maxlength 250 -step 0.5 -num $number_tracks
	fi
fi

if [ "$sift" = "yes" ] && [ ! -f $PRD/connectivity_regions/whole_brain_post.tck ]
then
	echo "using sift"
	if [ "$act" = "yes" ]
	then
		tcksift $PRD/connectivity_regions/whole_brain.tck $PRD/connectivity_regions/CSD$lmax.mif $PRD/connectivity_regions/whole_brain_post.tck -act $PRD/connectivity_regions/act.mif -term_number $(( number_tracks/2 ))
	else
		tcksift $PRD/connectivity_regions/whole_brain.tck $PRD/connectivity_regions/CSD$lmax.mif $PRD/connectivity_regions/whole_brain_post.tck  -term_number $(( number_tracks/2 ))
	fi
	elif [ ! -f $PRD/connectivity_regions/whole_brain_post.tck ]
	then
		echo "not using SIFT"
		ln -s $PRD/connectivity_regions/whole_brain.tck $PRD/connectivity_regions/whole_brain_post.tck
fi

#####now compute connectivity and length matrix
if [ ! -f $PRD/connectivity_regions/region_parcellation_2_diff.mif ]
then
    echo " compute labels"
    labelconfig $PRD/connectivity_regions/region_parcellation_2_diff.nii.gz fs_region.txt $PRD/connectivity_regions/region_parcellation_2_diff.mif -lut_freesurfer $FREESURFER_HOME/FreeSurferColorLUT.txt
fi

if [ ! -f $PRD/"$SUBJ_ID"_regions/connectivity/weights.csv ]
then
    echo "compute connectivity matrix"
    tck2connectome $PRD/connectivity/whole_brain_post.tck $PRD/connectivity_regions/region_parcellation_2_diff.mif $PRD/"$SUBJ_ID"_regions/connectivity/weights.csv -assignment_radial_search 2
    tck2connectome $PRD/connectivity/whole_brain_post.tck $PRD/connectivity_regions/region_parcellation_2_diff.mif $PRD/"$SUBJ_ID"_regions/connectivity/tract_lengths.csv -metric meanlength -zero_diagonal -assignment_radial_search 2 
fi

#compute centres.txt
if [ ! -f $PRD/"$SUBJ_ID"_regions/connectivity/centres.txt ]
then
echo "compute region centres"
#needed to unzip region_parcellation.nii.gz if the parameter of load_nii matlab function (from compute_region_centres(.m)) not modified
#gunzip $PRD/connectivity_regions/region_parcellation.nii.gz
if [ -n "$matlab" ]
then
$matlab -r "run compute_region_centres.m; quit;" -nodesktop -nodisplay
else
sh compute_region_centres/distrib/run_compute_region_centres.sh $MCR
fi
fi

#compute weights.txt and tract_lengths.txt
if [ ! -f $PRD/"$SUBJ_ID"_regions/connectivity/weights.txt ]
then
    echo " generate useful files for TVB"
    python compute_connectivity_files.py
fi

# zip to put in final format
pushd .
cd $PRD/"$SUBJ_ID"_regions/connectivity
zip $PRD/"$SUBJ_ID"_regions/connectivity.zip weights.txt tract_lengths.txt centres.txt
popd
