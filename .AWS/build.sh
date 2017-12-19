#! /bin/bash

gcc --version

CORES=$(grep -c ^processor /proc/cpuinfo)

echo "Number cores: ${CORES}"

echo `pwd`
mkdir bin
cd bin
cmake -D link_time_optimization=ON ../

make VERBOSE=1 -j${CORES}

