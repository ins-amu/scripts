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
  mkdir -p "$PRD"/test_"$NAME_TEST"
  cp -r "$PRD"/data "$PRD"/test_"$NAME_TEST"/data
  PRD="$PRD"/test_"$NAME_TEST"
}

Teardown() {
	rm -r "$PRD"
}


test_function() {
  bash ./main_surface.sh -c "test" -e -q -f > /dev/null && out="success" || out="fail"
  printf "\n >>> Test "$NAME_TEST" output is : "$out" <<< \n"
}

#( export NAME_TEST="fsl_5"; export FSL="fsl5.0-"; setUp; test_function ) &
#( export NAME_TEST="registration_boundary"; export REGISTRATION="boundary"; setUp; test_function; Teardown ) &
#( export NAME_TEST="registration_pseudo"; export REGISTRATION="pseudo"; setUp; test_function; Teardown ) &
#( export NAME_TEST="region_mapping_corr"; export REGION_MAPPING_COOR="0.5"; setUp; test_function; Teardown ) &
#( export NAME_TEST="k_list"; export K_LIST="0 2 5"; setUp; test_function; Teardown  ) &
#( export NAME_TEST="no_k_list"; export K_LIST=""; setUp; test_function; Teardown ) &
#( export NAME_TEST="number_tracks"; export NUMBER_TRACKS=10000"; setUp; test_function; Teardown ) &
( export NAME_TEST="parcel_destrieux"; PARCEL="destrieux"; setUp; test_function ) &
#( export NAME_TEST="parcel_HCP"; export PARCEL="HCP-MMP"; setUp; test_function ) &
#( export NAME_TEST="parcel_Yeo7"; export PARCEL="Yeo-7nets"; setUp; test_function) &
#( export NAME_TEST="parcel_Yeo17"; export PARCEL="Yeo-17nets"; setUp; test_function; Teardown ) &
#( export NAME_TEST="parcel_no_topup"; export TOPUP="no"; setUp; test_function ) &
#( export NAME_TEST="no_act"; export ACT="no"; setUp; test_function; Teardown ) &
#( export NAME_TEST="no_sift"; export SIFT="no""; setUp; test_function; Teardown ) &
#( export NAME_TEST="sift"; export SIFT="sift""; setUp; test_function; Teardown ) &
#( export NAME_TEST="sift_multiplier"; export SIFT="sift"; export SIFT_MULTIPLIER=2; setUp; test_function; Teardown ) &
#( export NAME_TEST="seed_dynamic"; export SEED="dynamic"; setUp; test_function; Teardown ) &
#( export NAME_TEST="aseg_fs"; export PARCEL="destrieux"; setUp; test_function; Teardown ) &
#( export NAME_TEST="nb_threads_2"; export ASEG="fs"; setUp; test_function; Teardown ) &

exit