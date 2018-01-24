#!/bin/bash
set -x

dir=$1
build_type=$2
static_libs=$3

echo
echo "--------------------------------------"
echo "build_type: $build_type"
echo "static_libs: $static_libs"
echo "--------------------------------------"
echo

echo "fortran compiler version:" 
$FC --version

cd ${dir}
echo `pwd`

cmake -D CMAKE_BUILD_TYPE=$build_type \
      -D static_libraries=$static_libs \
      ../

make -j${CORES}
