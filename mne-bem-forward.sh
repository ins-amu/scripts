#!/bin/bash

# this script should do most of the anatomical pipeline, and the only
# user intervention is the create the correct mri-meg trans file.

# args: input mri, subject id

# could set the subjects_dir according to project, etc. 

export INPUT=$1
export SUBJECTS_DIR=$2
export SUBJECT=$3

# setup freesurfer & mne environments
source $FREESURFER_HOME/FreeSurferEnv.sh
source $MNE_ROOT/bin/mne_setup_sh 

# Run freesurfer reconstruction
recon-all -all -i ${INPUT} -s $SUBJECT

# remove fsaverage symlinks
#find ${SUBJECTS_DIR} -type l -exec rm -rf {} \;

# use watershed algorithm to make BEM surfaces
if [ ! -f ${SUBJECTS_DIR}/${SUBJECT}/bem/${SUBJECT}-inner-skull.surf ];
then
    cd ${SUBJECTS_DIR}/${SUBJECT}/bem
    mne_watershed_bem
    ln -s watershed/${SUBJECT}_inner_skull_surface ${SUBJECT}-inner_skull.surf
    ln -s watershed/${SUBJECT}_outer_skin_surface ${SUBJECT}-outer_skin.surf
    ln -s watershed/${SUBJECT}_outer_skull_surface ${SUBJECT}-outer_skull.surf
fi

# Make high resolution head surface
HEAD_FNAME=${SUBJECTS_DIR}/${SUBJECT}/bem/${SUBJECT}-head.fif
if [ ! -f ${HEAD_FNAME} ];
then
    mkheadsurf -s ${SUBJECT}
    mne_surf2bem --surf ${SUBJECTS_DIR}/${SUBJECT}/surf/lh.seghead --id 4 --check --fif ${HEAD_FNAME}
fi

# Setup BEM
mne_setup_forward_model --surf --ico 4

