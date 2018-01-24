#!/bin/bash - 

cd ${dir}
echo `pwd`

ctest --output-on-failure -j${CORES}

export TEST_FAILURE=$?
return $?
