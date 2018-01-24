#!/bin/bash - 
# set -x

export FC=$(which gfortran)
export CC=$(which gcc)
export CXX=$(which g++)

export CORES=$(grep -c ^processor /proc/cpuinfo)
echo "Number cores: ${CORES}"

build_types=(Debug Release)
static_types=(TRUE FALSE)

pwd=$PWD

failed=false
results=()
for bt in "${build_types[@]}"; do
  for st in "${static_types[@]}"; do
    dir_name="gnu_${bt}_${st}"

    mkdir ${dir_name}
    ./.AWS/build.sh ${dir_name} ${bt} ${st}
    make_result=$?
    if [ ${make_result} -ne 0 ]; then
      failed=true
      results+=("(${dir_name}:\t${make_result})")
      continue
    fi

    ./.AWS/test.sh ${dir_name} ${bt} ${st}
    test_result=$?
    if [ ${test_result} -ne 0 ]; then
      failed=true
      results+=("(${dir_name}:\t${test_result})")
    fi
    cd ${pwd}
  done
done

echo
echo "Results:"
for res in ${results[@]}; do
  echo -e $res
done
echo

if [ ${failed} == true ]; then
  exit 1
fi
