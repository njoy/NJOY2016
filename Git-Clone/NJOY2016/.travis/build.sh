#!/bin/bash
set -x

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 90 \
       --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-6 \
       --slave /usr/bin/gcov gcov /usr/bin/gcov-6 \
       --slave /usr/bin/gcc-ar ar /usr/bin/gcc-ar-6 \
       --slave /usr/bin/gcc-nm nm /usr/bin/gcc-nm-6 \
       --slave /usr/bin/gcc-ranlib ranlib /usr/bin/gcc-ranlib-6
  sudo update-alternatives --config gcc
  export CUSTOM=('-D CMAKE_AR=/usr/bin/gcc-ar' '-D CMAKE_NM=/usr/bin/gcc-nm' '-D CMAKE_RANLIB=/usr/bin/gcc-ranlib')
else
  export CUSTOM=("-D no_link_time_optimization=TRUE")
fi

export FC=$(which gfortran)

mkdir build
cd build
cmake ${CUSTOM[@]}\
      -D CMAKE_BUILD_TYPE=$build_type \
      -D static_libraries=$static_libraries \
      -D appended_flags="$appended_flags" ..
make -j2
export COMPILATION_FAILURE=$?
if [ $COMPILATION_FAILURE -ne 0 ];
then
  exit 1
fi

ctest --output-on-failure -j2
export TEST_FAILURE=$?
if [ $TEST_FAILURE -ne 0 ];
then
    exit 1
fi
if [ "$build_type" = "coverage" ]
then
  pip install --upgrade pip
  pip install --user cpp-coveralls
  echo "loading coverage information"
  coveralls  --exclude-pattern "/usr/include/.*|.*/CMakeFiles/.*|.*subprojects.*|.*dependencies.*|.*test\.cpp" --root ".." --build-root "." --gcov-options '\-lp'
fi
