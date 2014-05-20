#######################
# Example config file #
#######################
####################################
# Copy this file wherever you want #
# (in the root directory for files #
# would be a good idea)            #
# and change the relevant options  #
####################################

#### Running the pipeline ####
# to run the pipeline, do:
# for region pipeline:
# sh path_to_scripts/main_region.sh -c path_to_config/config.sh
# for surface pipeline: 
# sh path_to_scripts/main_surface.sh -c path_to_config/config.sh

#the root directory for files 
# this is where all data and processed data are
# we advice to also put this config file in this directory
export PRD=/path_to_root_dir/

# path for subjects for freesurfer 
# $SUBJECTS_DIR should be in your .bashrc as you put
# it during freesurfer installation
export FS=$SUBJECTS_DIR

# subject name
# this will determine the name of your subject
# in brainvisa and in the final directory
export SUBJ_ID=name_subj

# brainvisa directory
export BV=/pat_to_brainvisa/

# Matlab Runtime Compiler path
# if you have matlab, do mcrinstaller at the matlab prompt
# to find the location of the MCR
# alternatively you can uncomment the following line and comment the MCR line
# export matlab=/path_to_matlab/
# if you don't have matlab, download the MCR (linux 64 bits) here: 
# http://www.mathworks.com/products/compiler/mcr/index.html
export MCR=/path_to_matlab_runtime_compiler/

# error handling: in case of error, the pipeline 
# stops immediately
#set -e

# check the processed data when the pipeline is running
# (you need a display) (no/yes)
export CHECK=no

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

# This parameter is important for the correction of the region mapping. Between 0 and 1. The bigger it is, the bigger is the correction. (only import for the surface pipeline: main_surface.sh)
export region_mapping_corr="0.42"

# number of tracks (divided by 100.000) used in the tractography step.
# note that you will have less tracks at the end.
export number_tracks=10

# if you want subdivided parcellations, you can set the folowing value
# according to the following table
# K:                7    8    9    10    11   12
# Number of Nodes:  128  256  512  1024  2048 4096
# export K=9

# uniquely for region pipeline (main_region.sh)
# choice of the parcellation
export parcel=AAL

