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

# mrconvert
if [ ! -f $PRD/connectivity/dwi.mif ]
then
if [ -f $PRD/data/DWI/*.nii.gz ]
then
ls $PRD/data/DWI/ | grep '.nii.gz$' | xargs -I {} mrconvert $PRD/data/DWI/{} $PRD/connectivity/dwi.mif -fslgrad $PRD/data/DWI/bvecs $PRD/data/DWI/bvals
else
mrconvert $PRD/data/DWI/ $PRD/connectivity/dwi.mif
fi
fi

# topup and eddy correction
if [ "$topup" = "yes" ]
then
echo "topup and eddy not implemented yes"
fi

if [ ! -f $PRD/connectivity/mask.mif ]
then
dwi2mask $PRD/connectivity/dwi.mif $PRD/connectivity/mask.mif
fi

if [ ! -f $PRD/connectivity/lowb.nii.gz ]
then
dwiextract $PRD/connectivity/dwi.mif $PRD/connectivity/lowb.mif -bzero
mrconvert $PRD/connectivity/lowb.mif $PRD/connectivity/lowb.nii.gz 
fi

# FLIRT registration
#Diff to T1
if [ ! -f $PRD/connectivity/T1.nii.gz ]
then
echo "generating good orientation for T1"
mri_convert --in_type mgz --out_type nii --out_orientation RAS $FS/$SUBJ_ID/mri/T1.mgz $PRD/connectivity/T1.nii.gz
fi
if [ ! -f $PRD/connectivity/aparc+aseg.nii.gz ]
then
echo " getting aparc+aseg"
mri_convert --in_type mgz --out_type nii --out_orientation RAS $FS/$SUBJ_ID/mri/aparc+aseg.mgz $PRD/connectivity/aparc+aseg.nii.gz
fi


if [ ! -f $PRD/connectivity/aparc+aseg_reorient.nii.gz ]
then
echo "reorienting the region parcellation"
fslreorient2std $PRD/connectivity/aparc+aseg.nii.gz $PRD/connectivity/aparc+aseg_reorient.nii.gz
# check parcellation to T1
if [ -n "$DISPLAY" ] && [ "$CHECK" = "yes" ]
then
echo "check parcellation"
echo " if it's correct, just close the window. Otherwise... well, it should be correct anyway"
fslview $PRD/connectivity/T1.nii.gz $PRD/connectivity/aparc+aseg_reorient.nii.gz -l "Cool"
fi
fi

if [ ! -f $PRD/connectivity/aparcaseg_2_diff.nii.gz ]
then
echo " register aparc+aseg to diff"
flirt -in $PRD/connectivity/lowb.nii.gz -ref $PRD/connectivity/T1.nii.gz -omat $PRD/connectivity/diffusion_2_struct.mat -out $PRD/connectivity/lowb_2_struct.nii.gz -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -cost mutualinfo
convert_xfm -omat $PRD/connectivity/diffusion_2_struct_inverse.mat -inverse $PRD/connectivity/diffusion_2_struct.mat
flirt -applyxfm -in $PRD/connectivity/aparc+aseg_reorient.nii.gz -ref $PRD/connectivity/lowb.nii.gz -out $PRD/connectivity/aparcaseg_2_diff.nii.gz -init $PRD/connectivity/diffusion_2_struct_inverse.mat -interp nearestneighbour

if [ ! -f $PRD/connectivity/T1_2_diff.nii.gz ]
then
echo " register aparc+aseg to diff"
flirt -applyxfm -in $PRD/connectivity/T1.nii.gz -ref $PRD/connectivity/lowb.nii.gz -out $PRD/connectivity/T1_2_diff.nii.gz -init $PRD/connectivity/diffusion_2_struct_inverse.mat -interp nearestneighbour
fi

# check parcellation to diff
if [ -n "$DISPLAY" ]  && [ "$CHECK" = "yes" ]
then
echo "check parcellation registration to diffusion space"
echo "if it's correct, just close the window. Otherwise you will have to
do the registration by hand"
fslview $PRD/connectivity/lowb.nii.gz $PRD/connectivity/aparcaseg_2_diff -l "Cool"
fi
fi

# response function estimation
if [ ! -f $PRD/connectivity/response.txt ]
then
echo "estimating response"
dwi2response $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt -mask $PRD/connectivity/mask.mif 
fi


# fibre orientation distribution estimation
if [ ! -f $PRD/connectivity/CSD$lmax.mif ]
then
echo "calculating fod"
dwi2fod $PRD/connectivity/dwi.mif $PRD/connectivity/response.txt $PRD/connectivity/CSD$lmax.mif -lmax $lmax -mask $PRD/connectivity/mask.mif
fi

# prepare file for act
if [ "$act" = "yes" ] && [ ! -f $PRD/connectivity/act.mif ]
then
echo "prepare files for act"
act_anat_prepare_fsl $PRD/connectivity/T1_2_diff.nii.gz $PRD/connectivity/act.mif
fi

# tractography
if [ ! -f $PRD/connectivity/whole_brain.tck ]
then
if [ "$act" = "yes" ]
then
echo "generating tracks using act" 
5tt2gmwmi $PRD/connectivity/act.mif $PRD/connectivity/gmwmi_mask.mif
tckgen $PRD/connectivity/CSD$lmax.mif $PRD/connectivity/whole_brain.tck -unidirectional -seed_gmwmi $PRD/connectivity/gmwmi_mask.mif -num $number_tracks -act $PRD/connectivity/act.mif -maxlength 150
else
echo "generating tracks without using act" 
tckgen $PRD/connectivity/CSD$lmax.mif $PRD/connectivity/whole_brain.tck -unidirectional -algorithm iFOD2 -seed_image $PRD/connectivity/aparcaseg_2_diff.nii.gz -mask $PRD/connectivity/mask.mif -maxlength 150 -num $number_tracks
fi
fi

if [ "$sift" = "yes" ] && [ ! -f $PRD/connectivity/whole_brain_post.tck ]
then
echo "using sift"
if [ "$act" = "yes" ]
then
tcksift $PRD/connectivity/whole_brain.tck $PRD/connectivity/CSD"$lmax".mif  $PRD/connectivity/whole_brain_post.tck -act $PRD/connectivity/act.mif -term_number $(( number_tracks/2 ))
else
tcksift $PRD/connectivity/whole_brain.tck $PRD/connectivity/CSD"$lmax".mif  $PRD/connectivity/whole_brain_post.tck -term_number $(( number_tracks/2 ))
fi
elif [ ! -f $PRD/connectivity/whole_brain_post.tck ]
then
echo "not using SIFT"
    ln -s $PRD/connectivity/whole_brain.tck  $PRD/connectivity/whole_brain_post.tck
fi


# now compute connectivity and length matrix
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
if [ -n "$K" ]
then
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
python compute_connectivity_sub.py
fi

pushd . > /dev/null
cd $PRD/$SUBJ_ID/connectivity_"$curr_K" > /dev/null
zip $PRD/$SUBJ_ID/connectivity_"$curr_K".zip weights.txt tract_lengths.txt centres.txt average_orientations.txt -q 
popd > /dev/null
fi

######################## compute MEG and EEG forward projection matrices
# make BEM surfaces
if [ ! -f ${FS}/${SUBJ_ID}/bem/inner_skull.surf ]
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
if [ ! -f ${FS}/${SUBJ_ID}/bem/*-bem-sol.fif ]
then
worked=0
outershift=0
while [ "$worked" == 0 ]
do
worked=1
mne_setup_forward_model --subject ${SUBJ_ID} --surf --ico 4 --outershift $outershift || worked=0 
if [ "$worked" == 0 ]
then
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
else
    echo "success!"
fi
done
fi

