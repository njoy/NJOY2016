#!/bin/bash - 
# set -x

compiler_family=$1

echo -e
echo "-----------------------------------------------------"
echo "Compiler Family: ${compiler_family}"
echo -e "FC:\t${FC}"
echo -e "CC:\t${CC}"
echo -e "CXX:\t${CXX}"
echo
echo "Compiler Versions:"
echo -e "FC version:\t" `${FC} --version`
echo -e "CC: version\t" `${CC} --version`
echo -e "CXX: version\t" `${CXX} --version`
echo "-----------------------------------------------------"
echo

export CORES=$(grep -c ^processor /proc/cpuinfo)
echo "Number cores: ${CORES}"

build_types=(Debug Release)
static_types=(TRUE FALSE)

pwd=$PWD

failed=false
results=()
for bt in "${build_types[@]}"; do
  for st in "${static_types[@]}"; do
    dir_name="${compiler_family}_${bt}_${st}"

    mkdir ${dir_name}
    ./.AWS/build.sh ${dir_name} ${bt} ${st}
    make_result=$?
    if [ ${make_result} -ne 0 ]; then
      failed=true
      results+=("${dir_name}:\t${make_result}")
      continue
    fi

    ./.AWS/test.sh ${dir_name} ${bt} ${st}
    test_result=$?
    if [ ${test_result} -ne 0 ]; then
      failed=true
      results+=("${dir_name}:\t${test_result}")
    fi
    results+=("${dir_name}:\t 0")
    cd ${pwd}
  done
done

echo
echo "-----------------------------"
echo "Results:"
for res in ${results[@]}; do
  echo -e $res
done
echo "-----------------------------"
echo

if [ ${failed} == true ]; then
  exit 1
fi
