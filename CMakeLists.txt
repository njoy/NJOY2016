
cmake_minimum_required( VERSION 3.15 )

find_package( Python3 3.5 REQUIRED )

set( CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported configuration types" FORCE )

project( njoy LANGUAGES Fortran )

message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )
message( STATUS "" )
message( STATUS "njoy" )
message( STATUS "Git current branch: ${GIT_BRANCH}" )
message( STATUS "Git commit hash: ${GIT_HASH}" )
message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )

add_executable( njoy_executable
             "${CMAKE_CURRENT_SOURCE_DIR}/src/main.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/acecm.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/acedo.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/acefc.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/acepa.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/acepn.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/acer.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/aceth.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/broadr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/ccccr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/covr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/dtfr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/endf.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/errorr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/gaminr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/gaspr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/graph.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/groupr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/heatr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/leapr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/locale.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/main.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/mainio.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/mathm.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/matxsr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/mixr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/moder.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/phys.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/plotr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/powr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/purr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/reconr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resxsr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/samm.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/thermr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/unresr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/util.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/vers.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/viewr.f90"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/wimsr.f90" )
set_target_properties( njoy_executable PROPERTIES OUTPUT_NAME njoy )

enable_testing()
add_subdirectory( tests )
