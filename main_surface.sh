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

###################### freesurfer
if [ ! -d $FS/$SUBJ_ID ] 
then
echo "running recon-all of freesurfer"
recon-all -i $PRD/data/T1/T1.nii -s $SUBJ_ID -all
fi

###################################### left hemisphere
# export pial into text file
mkdir -p $PRD/surface
if [ ! -f $PRD/surface/lh.pial.asc ]
then
echo "importing left pial surface from freesurfer"
mris_convert $FS/$SUBJ_ID/surf/lh.pial $PRD/surface/lh.pial.asc
fi

# triangles and vertices high
if [ ! -f $PRD/surface/lh_vertices_high.txt ]
then
echo "extracting left vertices and triangles"
python left_extract_high.py
fi

# decimation using brainvisa
if [ ! -f $PRD/surface/lh_vertices_low.txt ]
then
echo "left decimation using brainvisa"
# -> to mesh
$BV/bin/python left_transform_mesh_high.py
#  decimation
$BV/bin/AimsMeshDecimation $PRD/surface/lh_mesh_high.mesh $PRD/surface/lh_mesh_low.mesh
# export to list vertices triangles
$BV/bin/python left_export_to_vertices.py
fi

# create left the region mapping
if [ ! -f $PRD/surface/lh_region_mapping_low_not_corrected.txt ]
then
echo "generating the left region mapping on the decimated surface"
if [ -n "$matlab" ]
then
$matlab -r "run left_region_mapping.m; quit;" -nodesktop -nodisplay
else
sh left_region_mapping/distrib/run_left_region_mapping.sh $MCR
fi
fi

# check
if [ ! -f $PRD/surface/lh_region_mapping_low.txt ]
then
if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ] 
then
echo "check left region mapping"
python check_left_region_mapping.py
else
echo "correct the left region mapping"
# correct
python correct_left_region_mapping.py
fi
fi

###################################### right hemisphere
# export pial into text file
if [ ! -f $PRD/surface/rh.pial.asc ]
then
echo "importing right pial surface from freesurfer"
mris_convert $FS/$SUBJ_ID/surf/rh.pial $PRD/surface/rh.pial.asc
fi

# triangles and vertices high
if [ ! -f $PRD/surface/rh_vertices_high.txt ]
then
echo "extracting right vertices and triangles"
python right_extract_high.py
fi

# decimation using brainvisa
if [ ! -f $PRD/surface/rh_vertices_low.txt ]
then
echo "right decimation using brainvisa"
# -> to mesh
$BV/bin/python right_transform_mesh_high.py
#  decimation
$BV/bin/AimsMeshDecimation $PRD/surface/rh_mesh_high.mesh $PRD/surface/rh_mesh_low.mesh
# export to list vertices triangles
$BV/bin/python right_export_to_vertices.py
fi

if [ ! -f $PRD/surface/rh_region_mapping_low_not_corrected.txt ]
then
echo "generating the right region mapping on the decimated surface"
# create left the region mapping
if [ -n "$matlab" ]
then
$matlab -r "run right_region_mapping.m; quit;" -nodesktop -nodisplay
else
sh right_region_mapping/distrib/run_right_region_mapping.sh $MCR
fi
fi

# check
if [ ! -f $PRD/surface/rh_region_mapping_low.txt ]
then
if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
then
echo "check right region mapping"
python check_right_region_mapping.py
else
echo " correct the right region mapping"
# correct
python correct_right_region_mapping.py
fi
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
pushd .
cd $PRD/$SUBJ_ID/surface
zip $PRD/$SUBJ_ID/surface.zip vertices.txt triangles.txt
cp region_mapping.txt ..
popd

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

########################## build connectivity
# mrtrix
mkdir -p $PRD/connectivity
mkdir -p $PRD/$SUBJ_ID/connectivity
# mrconvert
if [ ! -f $PRD/connectivity/dwi.mif ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .nii | xargs -I {} mrconvert $PRD/data/DWI/{} $PRD/connectivity/dwi.mif 
else
mrconvert $PRD/data/DWI/ $PRD/connectivity/dwi.mif
fi
fi

# brainmask 
if [ ! -f $PRD/connectivity/lowb.nii  ]
then
average $PRD/connectivity/dwi.mif -axis 3 $PRD/connectivity/lowb.nii
fi
if [ ! -f $PRD/connectivity/mask_not_checked.mif ]
then
threshold -percent $percent_value_mask $PRD/connectivity/lowb.nii - | median3D - - | median3D - $PRD/connectivity/mask_not_checked.mif
fi

# check the mask
if [ ! -f $PRD/connectivity/mask_checked.mif ]
then
cp $PRD/connectivity/mask_not_checked.mif $PRD/connectivity/mask_checked.mif
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
while true; do
mrview $PRD/connectivity/mask_checked.mif
read -p "was the mask good?" yn
case $yn in
[Yy]* ) break;;
[Nn]* ) read -p "enter new threshold value" percent_value_mask; echo $percent_value_mask; rm $PRD/connectivity/mask_checked.mif; 
	threshold -percent $percent_value_mask $PRD/connectivity/lowb.nii - | median3D - - | median3D - $PRD/connectivity/mask.mif;;
 * ) echo "Please answer y or n.";;
esac
done
fi
fi

if [ -f $PRD/connectivity/mask_checked.mif ]
then
cp $PRD/connectivity/mask_checked.mif $PRD/connectivity/mask.mif
elif [ -f $PRD/connectivity/mask_not_checked.mif ]
then
cp $PRD/connectivity/mask_not_checked.mif $PRD/connectivity/mask.mif
fi

# tensor imaging
if [ ! -f $PRD/connectivity/dt.mif ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .b | xargs -I {} dwi2tensor $PRD/connectivity/dwi.mif $PRD/connectivity/dt.mif -grad $PRD/data/DWI/{}
else
dwi2tensor $PRD/connectivity/dwi.mif $PRD/connectivity/dt.mif
fi
fi
if [ ! -f $PRD/connectivity/fa.mif ]
then
tensor2FA $PRD/connectivity/dt.mif - | mrmult - $PRD/connectivity/mask.mif $PRD/connectivity/fa.mif
fi
if [ ! -f $PRD/connectivity/ev.mif ]
then
tensor2vector $PRD/connectivity/dt.mif - | mrmult - $PRD/connectivity/fa.mif $PRD/connectivity/ev.mif
fi
# constrained spherical deconvolution
if [ ! -f $PRD/connectivity/sf.mif ]
then
erode $PRD/connectivity/mask.mif -npass 3 - | mrmult $PRD/connectivity/fa.mif - - | threshold - -abs 0.7 $PRD/connectivity/sf.mif
fi
if [ ! -f $PRD/connectivity/response.txt ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .b | xargs -I {} estimate_response $PRD/connectivity/dwi.mif $PRD/connectivity/sf.mif -lmax $lmax $PRD/connectivity/response.txt -grad $PRD/data/DWI/{}
else
estimate_response $PRD/connectivity/dwi.mif $PRD/connectivity/sf.mif -lmax $lmax $PRD/connectivity/response.txt
fi

if  [ -n "$DISPLAY" ]  &&  [ "$CHECK" = "yes" ]
then
disp_profile -response $PRD/connectivity/response.txt
fi
fi
if [ ! -f $PRD/connectivity/CSD6.mif ]
then
if [ -f $PRD/data/DWI/*.nii ]
then
ls $PRD/data/DWI/ | grep .b | xargs -I {} csdeconv $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt -lmax $lmax -mask $PRD/connectivity/mask.mif $PRD/connectivity/CSD6.mif -grad $PRD/data/DWI/{}
else
csdeconv $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt -lmax $lmax -mask $PRD/connectivity/mask.mif $PRD/connectivity/CSD6.mif
fi
fi

# tractography
if [ ! -f $PRD/connectivity/whole_brain_1.tck ]
then
echo "generating tracks"
for I in $(seq 1 $number_tracks)
do
streamtrack SD_PROB $PRD/connectivity/CSD6.mif -seed $PRD/connectivity/mask.mif -mask $PRD/connectivity/mask.mif $PRD/connectivity/whole_brain_$I.tck -num 100000
done
fi

# FLIRT registration
#Diff to T1
if [ ! -f $PRD/connectivity/T1.nii ]
then
echo "generating good orientation for T1"
mri_convert --in_type mgz --out_type nii --out_orientation RAS $FS/$SUBJ_ID/mri/T1.mgz $PRD/connectivity/T1.nii
fi
if [ ! -f $PRD/connectivity/aparc+aseg.nii ]
then
echo " getting aparc+aseg"
mri_convert --in_type mgz --out_type nii --out_orientation RAS $FS/$SUBJ_ID/mri/aparc+aseg.mgz $PRD/connectivity/aparc+aseg.nii
fi


if [ ! -f $PRD/connectivity/aparc+aseg_reorient.nii ]
then
echo "reorienting the region parcellation"
fslreorient2std $PRD/connectivity/aparc+aseg.nii $PRD/connectivity/aparc+aseg_reorient.nii
# check parcellation to T1
if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
then
echo "check parcellation"
echo " if it's correct, just close the window. Otherwise... well, it should be correct anyway"
fslview $PRD/connectivity/T1.nii $PRD/connectivity/aparc+aseg_reorient -l "Cool"
fi
fi

if [ ! -f $PRD/connectivity/aparcaseg_2_diff.nii.gz ]
then
echo " register aparc+aseg to diff"
flirt -in $PRD/connectivity/lowb.nii -ref $PRD/connectivity/T1.nii -omat $PRD/connectivity/diffusion_2_struct.mat -out $PRD/connectivity/lowb_2_struct.nii -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -cost mutualinfo
convert_xfm -omat $PRD/connectivity/diffusion_2_struct_inverse.mat -inverse $PRD/connectivity/diffusion_2_struct.mat
flirt -applyxfm -in $PRD/connectivity/aparc+aseg_reorient.nii -ref $PRD/connectivity/lowb.nii -out $PRD/connectivity/aparcaseg_2_diff.nii -init $PRD/connectivity/diffusion_2_struct_inverse.mat -interp nearestneighbour
#flirt -in $PRD/connectivity/aparc+aseg.nii -ref $PRD/connectivity/lowb.nii -out $PRD/connectivity/aparcaseg_2_diff.nii -interp nearestneighbour 
# check parcellation to diff
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
echo "check parcellation registration to diffusion space"
echo "if it's correct, just close the window. Otherwise you will have to
do the registration by hand"
fslview $PRD/connectivity/lowb.nii $PRD/connectivity/aparcaseg_2_diff -l "Cool"
fi
fi

# compute sub parcellations connectivity if asked
if [ -n "$K" ]
then
export curr_K=$(( 2**K ))
mkdir -p $PRD/$SUBJ_ID/connectivity_"$curr_K"
if [ -n "$matlab" ]  
then
if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii ]
then
$matlab -r "run subparcel.m; quit;" -nodesktop -nodisplay 
fi
if [ ! -f $PRD/$SUBJ_ID/connectivity_"$curr_K"/weights.txt ]
then
$matlab -r "run compute_connectivity.m; quit;" -nodesktop -nodisplay
fi
else
if [ ! -f $PRD/connectivity/aparcaseg_2_diff_"$curr_K".nii ]
then
sh subparcel/distrib/run_subparcel.sh $MCR  
fi
if [ ! -f $PRD/$SUBJ_ID/connectivity_"$curr_K"/weights.txt ]
then
sh compute_connectivity/distrib/run_compute_connectivity.sh $MCR
fi
fi
pushd .
cd $PRD/$SUBJ_ID/connectivity_"$curr_K"
zip $PRD/$SUBJ_ID/connectivity_"$curr_K".zip weights.txt tract_lengths.txt centres.txt
popd
fi


# now compute connectivity and length matrix
if [ ! -f $PRD/$SUBJ_ID/connectivity/weights.txt ]
then
echo "compute connectivity matrix"
export curr_K=""
if [ -n "$matlab" ]
then
$matlab -r "run compute_connectivity.m; quit;" -nodesktop -nodisplay
else
sh compute_connectivity/distrib/run_compute_connectivity.sh $MCR
fi
fi

########
# we do not compute hemisphere
# subcortical is already done
cp cortical.txt $PRD/$SUBJ_ID/connectivity/cortical.txt

# # compute centers, areas and orientations
if [ ! -f $PRD/$SUBJ_ID/connectivity/centres.txt ]
then
echo " generate useful files for TVB"
python compute_other_files.py
fi

# zip to put in final format
pushd .
cd $PRD/$SUBJ_ID/connectivity
zip $PRD/$SUBJ_ID/connectivity.zip areas.txt average_orientations.txt weights.txt tract_lengths.txt cortical.txt centres.txt
bzip2 *
popd
