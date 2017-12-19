#!/bin/bash
CORES=$(grep -c ^processor /proc/cpuinfo)

ctest --output-on-failure -j${CORES}

export TEST_FAILURE=$?
if [ $TEST_FAILURE -ne 0 ];
then
    exit 1
fi
