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
if [ ! -f $PRD/data/T1/T1.nii ]
then
echo "generating T1 from DICOM"
mrconvert $PRD/data/T1/ $PRD/data/T1/T1.nii
fi

########################## build connectivity
# mrtrix
mkdir -p $PRD/connectivity_regions
mkdir -p $PRD/"$SUBJ_ID"_regions/connectivity
# mrconvert
if [ ! -f $PRD/connectivity_regions/dwi.mif ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .nii | xargs -I {} mrconvert $PRD/data/DWI/{} $PRD/connectivity_regions/dwi.mif 
else
mrconvert $PRD/data/DWI/ $PRD/connectivity_regions/dwi.mif
fi
fi
# brainmask 
if [ ! -f $PRD/connectivity_regions/lowb.nii  ]
then
average $PRD/connectivity_regions/dwi.mif -axis 3 $PRD/connectivity_regions/lowb.nii
fi
if [ ! -f $PRD/connectivity_regions/mask_not_checked.mif ]
then
threshold -percent $percent_value_mask $PRD/connectivity_regions/lowb.nii - | median3D - - | median3D - $PRD/connectivity_regions/mask_not_checked.mif
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
	threshold -percent $percent_value_mask $PRD/connectivity_regions/lowb.nii - | median3D - - | median3D - $PRD/connectivity_regions/mask_checked.mif;;
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
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .b | xargs -I {} dwi2tensor $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/dt.mif -grad $PRD/data/DWI/{}
else
dwi2tensor $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/dt.mif
fi
fi
if [ ! -f $PRD/connectivity_regions/fa.mif ]
then
tensor2FA $PRD/connectivity_regions/dt.mif - | mrmult - $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/fa.mif
fi
if [ ! -f $PRD/connectivity_regions/ev.mif ]
then
tensor2vector $PRD/connectivity_regions/dt.mif - | mrmult - $PRD/connectivity_regions/fa.mif $PRD/connectivity_regions/ev.mif
fi
# constrained spherical decconvolution
if [ ! -f $PRD/connectivity_regions/sf.mif ]
then
erode $PRD/connectivity_regions/mask.mif -npass 3 - | mrmult $PRD/connectivity_regions/fa.mif - - | threshold - -abs 0.7 $PRD/connectivity_regions/sf.mif
fi
if [ ! -f $PRD/connectivity_regions/response.txt ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .b | xargs -I {} estimate_response $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/sf.mif -lmax $lmax $PRD/connectivity_regions/response.txt -grad $PRD/data/DWI/{}
else
estimate_response $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/sf.mif -lmax $lmax $PRD/connectivity_regions/response.txt
fi
fi
if  [ -n "$DISPLAY" ]  &&  [ "$CHECK" = "yes" ]
then
echo "check the response function"
echo "it should be broadest in the axial plane, and have low amplitude along the z-axis."
disp_profile -response $PRD/connectivity_regions/response.txt
fi
if [ ! -f $PRD/connectivity_regions/CSD6.mif ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .b | xargs -I {} csdeconv $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/response.txt -lmax $lmax -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/CSD6.mif -grad $PRD/data/DWI/{}
else
csdeconv $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/response.txt -lmax $lmax -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/CSD6.mif
fi
fi

# tractography
if [ ! -f $PRD/connectivity_regions/whole_brain_1.tck ]
then
echo "generating tracks"
for I in $(seq 1 $number_tracks)
do
streamtrack SD_PROB $PRD/connectivity_regions/CSD6.mif -seed $PRD/connectivity_regions/mask.mif -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/whole_brain_$I.tck -num 100000
done
fi

# FLIRT registration
# parcellation MNI to T1 using FNIRT
# first preregistration T1 to MNI using FLIRT
# flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in $PRD/data/T1/T1.nii -omat $PRD/connectivity_regions/t1_2_mni_transf.mat
# then register T1 to MNI using FNIRT
# fnirt --in=$PRD/data/T1/T1.nii --aff=$PRD/connectivity_regions/t1_2_mni_transf.mat --cout=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf --config=T1_2_MNI152_2mm
# apply the warp to check the registration 
# applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --in=$PRD/data/T1/T1.nii --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf --out=$PRD/connectivity_regions/t1_2_mni_warped_structural_2mm
# inverse the linear registration
# convert_xfm -omat $PRD/connectivity_regions/t1_2_mni_transf_inverse.mat -inverse $PRD/connectivity_regions/t1_2_mni_transf.mat
# inverse the warp
# invwarp --ref=$PRD/data/T1/T1.nii --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf.nii.gz --out=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf_inverse
#Bring the chosen parcellation to T1-Space
# applywarp --ref=$PRD/data/T1/T1.nii --in=parcellations/$parcel --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf_inverse --out=$PRD/connectivity_regions/region_parcellation

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
fslreorient2std $PRD/data/T1/T1.nii $PRD/connectivity_regions/T1.nii.gz
fi

# check parcellation to T1
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
echo "check parcellation"
echo "if it's correct, just close the window. Otherwise you will have to
do the registration by hand"
fslview $PRD/connectivity_regions/T1.nii.gz $PRD/connectivity_regions/region_parcellation_reorient -l "Cool"
fi

if [ ! -f $PRD/connectivity_regions/region_parcellation_2_diff.nii.gz ]
then
echo "register parcellation to diff"
flirt -in $PRD/connectivity_regions/lowb.nii -ref $PRD/connectivity_regions/T1.nii.gz -omat $PRD/connectivity_regions/diffusion_2_struct.mat -out $PRD/connectivity_regions/lowb_2_struct.nii -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -cost mutualinfo
convert_xfm -omat $PRD/connectivity_regions/diffusion_2_struct_inverse.mat -inverse $PRD/connectivity_regions/diffusion_2_struct.mat
flirt -applyxfm -in $PRD/connectivity_regions/region_parcellation_reorient.nii -ref $PRD/connectivity_regions/lowb.nii -out $PRD/connectivity_regions/region_parcellation_2_diff.nii -init $PRD/connectivity_regions/diffusion_2_struct_inverse.mat -interp nearestneighbour
fi

# check parcellation to diff
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
echo "check parcellation"
echo "if it's correct, just close the window. Otherwise you will have to
do the registration by hand"
fslview $PRD/connectivity_regions/lowb.nii $PRD/connectivity_regions/region_parcellation_2_diff -l "Cool"
fi

# now compute connectivity and length matrix
if [ ! -f $PRD/"$SUBJ_ID"_regions/connectivity/weights.txt ]
then
echo "compute connectivity region"
if [ -n "$matlab" ]
then
$matlab -r "run compute_connectivity_region.m; quit;" -nodesktop -nodisplay
else
sh compute_connectivity_region/distrib/run_compute_connectivity_region.sh $MCR
fi
fi

# now compute connectivity and length matrix
if [ ! -f $PRD/"$SUBJ_ID"_regions/connectivity/centres.txt ]
then
echo "compute region centres"
if [ -n "$matlab" ]
then
$matlab -r "run compute_region_centres.m; quit;" -nodesktop -nodisplay
else
sh compute_region_centres/distrib/run_compute_region_centres.sh $MCR
fi
fi


# zip to put in final format
pushd .
cd $PRD/"$SUBJ_ID"_regions/connectivity
zip $PRD/"$SUBJ_ID"_regions/connectivity.zip weights.txt tract_lengths.txt centres.txt
bzip2 *
popd
