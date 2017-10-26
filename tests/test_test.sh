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

test_test_false() {
  false > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_test_false output is : "$out" <<< \n"
}

test_test_true() {
  true > /dev/null && out="success" || out="fail"
  printf "\n >>> Test test_test_true output is : "$out" <<< \n"
}

( setUp; test_test_false; Teardown ) &
( setUp; test_test_true ; Teardown ) &

exit
