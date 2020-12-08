#!/usr/bin/env bash
# This is the separate  config file for computing connectome harmonics
# It assumes the main SCRIPTS config file has been set up and the main_surface.sh ran without errors.

# By default, the routine compute subject specific harmonics based on all the tracks computed by the tractography


#######################
## Setup environment ##
#######################

# Path to SCRIPTS root directory
export SCRIPTS_DIR=`pwd`

# Path to Matlab Utils (for dependencies)
export PMU=$SCRIPTS_DIR/matlab_utils/

# number of parallel cores to compute high-resolution connectome (different than NTHREADS because more memory consuming)
export NWORKERS=4


