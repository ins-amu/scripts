#!/usr/bin/env bash

setUp() {
  source /disk2/Work/Processed_data/brown/scripts/cf-p08/config_cf-p08.sh
  set +e
  name_dir=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 10)
  mkdir -p "$PRD"/test_"$name_dir"
  cp -r "$PRD"/data "$PRD"/test_"$name_dir"/data
  PRD="$PRD"/test_"$name_dir"
}

Teardown() {
	rm -r "$PRD"
}

test_fsl_5() {
  export FSL="fsl5.0-"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_registration_boundary() {
  export REGISTRATION="boundary"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_registration_pseudo() {
  export REGISTRATION="pseudo"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_region_mapping_corr() {
  export REGION_MATPPING_COOR="0.5"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_k_list() {
  export K_LIST="0 2 5"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_no_k_list() {
  export K_LIST=""
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_number_tracks() {
  export NUMBER_TRACKS=10000
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_parcel_destrieux() {
  export PARCEL="destrieux"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_parcel_HCP() {
  export PARCEL="HCP"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_no_topup() {
  export TOPUP="no"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_no_act() {
  export ACT="no"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_no_sift() {
  export SIFT="no"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_sift() {
  export SIFT="yes"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_seed_dynamic() {
  export SEED="dynamic"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_aseg_fs() {
  export ASEG="fs"
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

test_nb_threads_2() {
  export NB_THREADS=2
  bash ../main_surface.sh -c "test" -e > /dev/null && out=0 || out=1
  echo ">>> Test output is : "$out""
}

setUp; test_registration_boundary; Teardown &
setUp; test_registration_pseudo; Teardown &
setUp; test_registration_regular; Teardown &
setUp; test_number_tracks; Teardown &
setUp; test_no_topup; Teardown &
setUp; test_topup; Teardown &
setUp; test_topup; Teardown &
setUp; test_topup; Teardown &

exit
