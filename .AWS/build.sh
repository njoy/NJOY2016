#! /bin/bash

gcc --version

CORES=$(grep -c ^processor /proc/cpuinfo)

echo "Number cores: ${CORES}"

mkdir bin
cd bin
cmake -D link_time_optimization=ON ../

make VERBOSE=1 -j${CORES}

ctest --output-on-failure -j${CORES}

export TEST_FAILURE=$?
if [ $TEST_FAILURE -ne 0 ];
then
    exit 1
fi
