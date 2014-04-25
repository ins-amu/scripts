# Script for regions simulations 
################ 
#the directory with all files
export PRD=/disk2/Work/Processed_data/jonathan_subjects/ab/
# subject name
export SUBJ_ID=ab
# FSL_DIR
export FSL_DIR=/usr/share/fsl/5.0
# matlab path
alias matlab='/home/tim/Matlab/bin/matlab'
# error handling
#set -e

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
export lmax=6 

# chosen parcellation
export parcel=AAL

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
# parcellation MNI to T1 using FNIRT
# first preregistration T1 to MNI using FLIRT
#flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in $PRD/data/T1/T1.nii -omat $PRD/connectivity_regions/t1_2_mni_transf.mat
## then register T1 to MNI using FNIRT
#fnirt --in=$PRD/data/T1/T1.nii --aff=$PRD/connectivity_regions/t1_2_mni_transf.mat --cout=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf --config=T1_2_MNI152_2mm
## apply the warp to check the registration 
#applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --in=$PRD/data/T1/T1.nii --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf --out=$PRD/connectivity_regions/t1_2_mni_warped_structural_2mm
## inverse the linear registration
#convert_xfm -omat $PRD/connectivity_regions/t1_2_mni_transf_inverse.mat -inverse $PRD/connectivity_regions/t1_2_mni_transf.mat
## inverse the warp
#invwarp --ref=$PRD/data/T/T1.nii --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf.nii.gz --out=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf_inverse
##Bring the chosen parcellation to T1-Space
#applywarp --ref=$PRD/data/T1/T1.nii --in=parcellations/$parcel --warp=$PRD/connectivity_regions/t1_2_mni_nonlinear_transf_inverse --out=$PRD/connectivity_regions/region_parcellation
# parcellation (T1 space) to diff
flirt -in $PRD/connectivity_regions/region_parcellation.nii -ref $PRD/connectivity_regions/lowb.nii -out $PRD/connectivity_regions/region_parcellation_2_diff.nii -interp nearestneighbour 

# now compute connectivity and length matrix
matlab -r "run compute_connectivity_region.m; quit;" -nodesktop -nodisplay

########
# we do not compute hemisphere
# subcortical is already done
cp cortical.txt $PRD/$SUBJ_ID/connectivity/cortical.txt

# # compute centers
matlab -r "run compute_region_centres.m; quit;" -nodesktop -nodisplay

# zip to put in final format
cd $PRD/"$SUBJ_ID"_regions/connectivity
zip $PRD/"$SUBJ_ID"_regions/connectivity.zip weights.txt tracts.txt centres.txt
cd $PRD/scripts
