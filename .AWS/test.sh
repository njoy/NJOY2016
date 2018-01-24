#!/bin/bash - 

cd ${dir}
echo `pwd`

ctest --output-on-failure -j${CORES}

export TEST_FAILURE=$?
if [ $TEST_FAILURE -ne 0 ];
then
    exit 1
fi
