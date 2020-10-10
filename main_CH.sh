#!/usr/bin/env bash

source config_CH.sh

#####################################
## Convert tracks into ascii files ##
#####################################
if [ ! -d $PRD/connectivity/tmp_ascii_tck ]; then
  mkdir -p $PRD/connectivity/tmp_ascii_tck
  echo "exporting tracts in ascii... this can take a while"
  tckconvert $PRD/connectivity/whole_brain_post.tck \
             $PRD/connectivity/tmp_ascii_tck/output-[].txt \
             -scanner2voxel $PRD/connectivity/brain_2_diff.nii.gz \
             -nthreads "$NB_THREADS"
fi

#######################################
## Run Connectome Harmonics pipeline ##
#######################################
echo "Run Pipeline matlab routine..."
$MATLAB -r "run CH_pipeline_integrated.m; quit;" -nodesktop
