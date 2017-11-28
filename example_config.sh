#######################
# Example config file #
#######################

####################################
# Copy this file wherever you want #
# (in the root directory for your  #
# data would be a good idea)       #
# and change the relevant options  #
####################################

#### Path configurations

# to run the pipeline, do:
# for region pipeline:
# bash path_to_scripts/main_region.sh -c path_to_config/config.sh
# for surface pipeline: 
# bash path_to_scripts/main_surface.sh -c path_to_config/config.sh

# The defaults indicates the options chosen by the pipeline if the option is 
# not set

#### Mandatory path parameters

# the root directory for files 
# this is where all data and processed data are
# we advice to also put this config file in this directory
export PRD=/path_to_root_dir/

# subject name
# this will determine the name of your subject
# in brainvisa and in the final directory
export SUBJ_ID=name_subj

# Matlab path
# export MATLAB=/path_to_matlab/
export MATLAB=$(which matlab)

#### Standard options

# error handling: in case of error, the pipeline stops immediately
# Can also be set with option -e to main_surface.sh
# default: not set; 
# set -e

# Force: all choices are made automatically without asking the user
# WARNING: can lead to issues
# options: ["no", "yes"]; default: "no"
# Can also be set with option -f to main_surface.sh
export FORCE="no"

# Quiet: run the pipeline without any output
# options: ["no", "yes"]; default: "no"
# Can also be set with option -q to main_surface.sh
export QUIET="no"


#### Optional parameters
# FSL prefix in case of use of fsl5.0 and fsl 4 is present
# for instance FSL="fsl5.0' or FSL="" otherwise
# if only fsl5.0 is installed, leave empty
# default: empty
export FSL=""

# HCP option
# Set to yes if your data are coming from the Human Connectome Project
# options: ["no", "yes"]; default: "no"
export HCP="no"

#### Pipeline parameters
# The defaults indicates the options chosen by the pipeline if the option is 
# not set

# check the processed data when the pipeline is running
# (you need a display and mrview installed) 
# options: ["yes"/"no"/"force"]; default: "no"
export CHECK="no"

# Methods for the FLIRT registration
# options: ["regular"/"boundary"/"pseudo"]; default: "regular"
export REGISTRATION="regular"

# This parameter is important for the correction of the region mapping. 
# Between 0 and 1. The bigger it is, the bigger is the correction. 
# (only import for the surface pipeline: main_surface.sh)
# default: 0.42
export REGION_MAPPING_CORR="0.42"

# for computing subconnectivity
# if you want subdivided parcellations, you can set the folowing value
# according to the following table
# K:                 0   1   2    3    4    5
# Number of Nodes:  70  140 280  560  1120 2240
# default: ""
# Needs to be a list of integers
export K_LIST="0 1 2 3 4 5"

# number of tracks used in the tractography step.
# default: 10.000.000
# Needs to be an integer
export NUMBER_TRACKS=10000000

# choice of the parcellation
# options ["desikan", "destrieux", "HCP-MMP"]; default: "desikan" 
export parcel="desikan"

# use topup and eddy distortion correction
# this depends of you images
# options: ["no", "eddy_correct"], default: "eddy_correct"
export TOPUP="eddy_correct"

# use Anatomically Constrained Tractography (yes/no)
# options ["yes", "no"]; default: "yes"
export ACT="yes"

# using Spherical-deconvolution informed filtering of tractograms 
# options: ["sift", "sift2", "no"]; default: SIFT2
export SIFT="sift2"

# if using SIFT, you can set the sift_multiplier variable:
# the number of tracks generated will be NUMBER_TRACKS*SIFT_MULTIPLIER
# default: 10
# Needs to be an integer
# export SIFT_MULTIPLIER=10

# seeding mechanism for tckgen if using act, otherwise default to dynamic
# options: ["gmwmi", "dynamic"]; default: "gmwmi"
export SEED="gmwmi"

# subcortical segmentation correction
# options: ["fs", "fsl"]; default: "fsl"
export ASEG="fsl" 

# 5ttgen
# options: ["fs", "fsl"]; default: "fsl"
export 5TTGEN="fsl" 

# compute forward model for MEG and EEG
# options: ["yes", "no"]; default: "yes"
export FORWARD_MODEL="yes"

# number of threads
# default: value in ~/.mrtrix.conf file if present, or 1 if not present
# Needs to be an integer
export NB_THREADS=1
