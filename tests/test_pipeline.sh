#!/usr/bin/env bash

# import and check config
while getopts "c:e" opt; do
  case $opt in
  c)
    CONFIG=$OPTARG
    if [ ! -f "$CONFIG" -a "$CONFIG" != "test" ];then
      echo "Config file "$CONFIG" does not exist." >&2
      exit 1
    elif [ $CONFIG = "test" ]; then
      echo "test mode"
    else
      echo "Using config file $CONFIG." >&2
      source "$CONFIG"
    fi
    ;;
  e) 
    set -e 
    ;;
  \?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

if [ -z "$CONFIG" ]; then
  echo "You must provide a config file."
  exit 1
fi

setUp() {
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
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_fsl_5 output is : "$out" <<< \n"
}

test_registration_boundary() {
  export REGISTRATION="boundary"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_registration_boundary output is : "$out" <<< \n"
}

test_registration_pseudo() {
  export REGISTRATION="pseudo"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_registration_pseudo output is : "$out" <<< \n"
}

test_region_mapping_corr() {
  export REGION_MATPPING_COOR="0.5"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_region_mapping_corr output is : "$out" <<< \n"
}

test_k_list() {
  export K_LIST="0 2 5"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_k_list output is : "$out" <<< \n"
}

test_no_k_list() {
  export K_LIST=""
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_no_k_list output is : "$out" <<< \n"
}

test_number_tracks() {
  export NUMBER_TRACKS=10000
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_number_tracks output is : "$out" <<< \n"
}

test_parcel_destrieux() {
  export PARCEL="destrieux"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_parcel_destrieux output is : "$out" <<< \n"
}

test_parcel_HCP() {
  export PARCEL="HCP-MMP"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_parcel_HCP output is : "$out" <<< \n"
}

test_parcel_Yeo7() {
  export PARCEL="HCP-MMP"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_parcel_HCP output is : "$out" <<< \n"
}

test_parcel_Yeo17() {
  export PARCEL="HCP-MMP"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_parcel_HCP output is : "$out" <<< \n"
}

test_no_topup() {
  export TOPUP="no"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_no_topup output is : "$out" <<< \n"
}

test_no_act() {
  export ACT="no"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_no_act output is : "$out" <<< \n"
}

test_no_sift() {
  export SIFT="no"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_no_sift output is : "$out" <<< \n"
}

test_sift() {
  export SIFT="sift"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_no_sift output is : "$out" <<< \n"
}

test_sift_multiplier() {
  export SIFT="sift"
  export SIFT_MULTIPLIER=2
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_no_sift output is : "$out" <<< \n"
}

test_seed_dynamic() {
  export SEED="dynamic"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_seed_dynamic output is : "$out" <<< \n"
}

test_aseg_fs() {
  export ASEG="fs"
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_aseg_fs output is : "$out" <<< \n"
}

test_nb_threads_2() {
  export NB_THREADS=2
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_nb_threads_2 output is : "$out" <<< \n"
}

#( setUp; test_fsl_5 ) &
#( setUp; test_registration_boundary; Teardown ) &
#( setUp; test_registration_pseudo; Teardown ) &
#( setUp; test_region_mapping_corr; Teardown ) &
#( setUp; test_k_list; Teardown  ) &
#( setUp; test_no_k_list; Teardown ) &
#( setUp; test_number_tracks; Teardown ) &
( setUp; test_parcel_destrieux ) &
( setUp; test_parcel_HCP ) &
( setUp; test_parcel_Yeo7 ) &
( setUp; test_parcel_Yeo17 ) &
#( setUp; test_no_topup; Teardown ) &
#( setUp; test_no_act; Teardown ) &
#( setUp; test_no_sift; Teardown ) &
#( setUp; test_sift; Teardown ) &
#( setUp; test_sift_multiplier; Teardown ) &
#( setUp; test_seed_dynamic; Teardown ) &
#( setUp; test_aseg_fs; Teardown ) &
#( setUp; test_nb_threads_2; Teardown ) &

exit