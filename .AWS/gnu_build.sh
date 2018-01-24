#!/bin/bash - 

export FC=$(which gfortran)
export CC=$(which gcc)
export CXX=$(which g++)

export CORES=$(grep -c ^processor /proc/cpuinfo)
echo "Number cores: ${CORES}"

# build_types=(Debug Release)
# static_types=(TRUE FALSE)
build_types=(Debug)
static_types=(TRUE)

pwd=$PWD

results=()
for bt in "${build_types[@]}"; do
  for st in "${static_types[@]}"; do
    dir_name="gnu_${bt}_${st}"
    mkdir ${dir_name}
    ./.AWS/build.sh ${dir_name} ${bt} ${st}
    ./.AWS/test.sh ${dir_name} ${bt} ${st}
    results+=(${dir_name}, $?)
    cd ${pwd}
  done
done

echo ${results}
