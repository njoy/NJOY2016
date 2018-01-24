#!/bin/bash - 
set -x

export FC=$(which gfortran)
export CC=$(which gcc)
export CXX=$(which g++)

export CORES=$(grep -c ^processor /proc/cpuinfo)
echo "Number cores: ${CORES}"

build_types=(Debug Release)
static_types=(TRUE False)

pwd=$PWD

failed=false
results=()
for bt in "${build_types[@]}"; do
  for st in "${static_types[@]}"; do
    dir_name="gnu_${bt}_${st}"
    mkdir ${dir_name}
    ./.AWS/build.sh ${dir_name} ${bt} ${st}
    ./.AWS/test.sh ${dir_name} ${bt} ${st}
    res=$?
    echo "res: ${res}"
    if [ ${res} -ne 0 ]; then
      failed=true
    fi
    results+=("(${dir_name}:\t${res})")
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
