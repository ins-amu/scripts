#!/usr/bin/env bash

setUp() {
  source /disk2/Work/Processed_data/brown/scripts/cf-p08/config_cf-p08.sh
  set +e
  name_dir=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 10)
  mkdir -p "$PRD"/test_"$name_dir"
  cp -r "$PRD"/data "$PRD"/test_"$name_dir"/data
  PRD="$PRD"/test_"$name_dir"
}

test_1() {
  bash ../main_surface.sh -c "test" -e && out=0 || out=1
  assertEquals "$out" 0
}

test_2() {
  bash ../main_surface.sh -c "test" -e  && out=0 || out=1
  assertEquals "$out" 0
}

. ../shunit2/shunit2