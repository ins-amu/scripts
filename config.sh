#######################
# Example config file #
#######################
######################################
# Copy this file wherever you want
# (in the root directory for processed
# data would be a good idea)
# and change the relevant options
######################################
# to run the pipeline, do:
# for region pipeline:
# sh path_to_scripts/main_region.sh -config path_to_config/config.sh
# for surface pipeline: 
# sh path_to_scripts/main_surface.sh -config path_to_config/config.sh
#####################################

#the root directory for files 
export PRD=/home/tim/Work/Processed_data/tim_pipeline/TREC/
# freesurfer 
export FS=$SUBJECTS_DIR
# subject name
export SUBJ_ID=TREC
# brainvisa directory
export BV=/home/tim/Work/Soft/brainvisa-4.3.0/
# matlab path
alias matlab='/usr/local/MATLAB/R2013a/bin/matlab'
# error handling
#set -e
# checking
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

