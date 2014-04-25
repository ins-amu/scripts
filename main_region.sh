# Script for regions simulations 
################ 
#the directory with all files
export PRD=/disk2/Work/Processed_data/jonathan_subjects/ab/
# subject name
export SUBJ_ID=ab
# matlab path
alias matlab='/home/tim/Matlab/bin/matlab'
# error handling
set -e

########## Important parameters
# This parameter modify the mask for diffusion processing
# check for different values with mrview
export percent_value_mask=10

# This parameter is the maximum harmonic order for spherical deconvolution
# It depends of the number of directions used during acquisition
# Please refer to the following table 
# Maximum harmonic order (lmax)	Number of parameters required
#             2				     6
#             4				     15
#             6				     28
#             8				     45
#             10			     66
#             12			     91
#             n				½ (n+1)(n+2)
export lmax = 6 

######### build cortical surface and region mapping
# cd $PRD/scripts
# mrconvert $PRD/data/T1/ $PRD/data/T1.nii

########################## build connectivity
# mrtrix
mkdir -p $PRD/connectivity_regions
mkdir -p $PRD/"$SUBJ_ID"_regions/connectivity
# mrconvert
mrconvert $PRD/data/DWI/ $PRD/connectivity_regions/dwi.mif
# brainmask 
average $PRD/connectivity_regions/dwi.mif -axis 3 $PRD/connectivity_regions/lowb.nii
threshold -percent $percent_value_mask $PRD/connectivity_regions/lowb.nii - | median3D - - | median3D - $PRD/connectivity_regions/mask.mif
# check the mask
if [ -n "$DISPLAY" ]; then mrview $PRD/connectivity_regions/mask.mif; fi
# tensor imaging
dwi2tensor $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/dt.mif
tensor2FA $PRD/connectivity_regions/dt.mif - | mrmult - $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/fa.mif
tensor2vector $PRD/connectivity_regions/dt.mif - | mrmult - $PRD/connectivity_regions/fa.mif $PRD/connectivity_regions/ev.mif
# constrained spherical decconvolution
erode $PRD/connectivity_regions/mask.mif -npass 3 - | mrmult $PRD/connectivity_regions/fa.mif - - | threshold - -abs 0.7 $PRD/connectivity_regions/sf.mif
estimate_response $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/sf.mif -lmax $lmax $PRD/connectivity_regions/response.txt
if [ -n "$DISPLAY" ]; then disp_profile -response $PRD/connectivity_regions/response.txt; fi
csdeconv $PRD/connectivity_regions/dwi.mif $PRD/connectivity_regions/response.txt -lmax $lmax -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/CSD6.mif
# tractography
for I in 1 2 3 4 5 6 7 8 9 10
do
streamtrack SD_PROB $PRD/connectivity_regions/CSD6.mif -seed $PRD/connectivity_regions/mask.mif -mask $PRD/connectivity_regions/mask.mif $PRD/connectivity_regions/whole_brain_$I.tck -num 100000
done

# FLIRT registration
# parcellation (T1 space) to diff
flirt -in $PRD/connectivity_regions/region_parcellation.nii -ref $PRD/connectivity_regions/lowb.nii -out $PRD/connectivity_regions/region_parcellation_2_diff.nii -interp nearestneighbour 

# now compute connectivity and length matrix
matlab -r "run compute_connectivity_region.m; quit;" -nodesktop -nodisplay

########
# we do not compute hemisphere
# subcortical is already done
cp cortical.txt $PRD/$SUBJ_ID/connectivity/cortical.txt

# # compute centers, areas and orientations
matlab -r "compute_regions_centres.m; quit;" -nodesktop -nodisplay

# zip to put in final format
cd $PRD/"$SUBJ_ID"_regions/connectivity
zip $PRD/"$SUBJ_ID"_regions/connectivity.zip weights.txt tracts.txt centres.txt
cd $PRD/scripts
