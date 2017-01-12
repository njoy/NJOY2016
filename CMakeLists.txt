cmake_minimum_required( VERSION 3.2 ) 

if( NOT DEFINED build_type )
    if( VERBOSE )
        message( STATUS "build_type variable not specified")
        message( STATUS "build_type defaulted to debug")
    endif()
    set ( build_type "debug" )
endif()
if( NOT DEFINED njoy_build_type )
    if( VERBOSE )
        message( STATUS "njoy_build_type not specified")
        message( STATUS "njoy_build_type defaulted to value of build_type variable")
    endif()
    set( njoy_build_type "${build_type}" )
 endif()

if( NOT DEFINED static_libraries )
    if( VERBOSE )
        message( STATUS "static_libraries variable not specified")
        message( STATUS "static_libraries defaulted to FALSE")
    endif()
    set( build_static_libraries FALSE )
else()
    set( build_static_libraries ${static_libraries} )
endif()

if( NOT DEFINED static_njoy )
    if( VERBOSE )
        message( STATUS "static_njoy variable not specified")
        message( STATUS "static_njoy defaulted to value of static_libraries variable")
    endif()
    set( static_njoy ${build_static_libraries} )
endif()

if( ${static_njoy} )
    set( njoy_policy STATIC )
else()
    set( njoy_policy SHARED )
endif()

project( njoy VERSION 2016 LANGUAGES Fortran )

if( DEFINED Fortran_compiler_flags )
    if( ( NOT DEFINED njoy_compiler_flags ) AND ( NOT DEFINED njoy_use_default_compiler_flags ) )
        set( njoy_compiler_flags "${Fortran_compiler_flags}" )
    endif()
endif()
if( NOT DEFINED njoy_compiler_flags )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
        if( njoy_build_type STREQUAL "debug" )
            set( njoy_compiler_flags " -Wall -Wextra -Wpedantic -g -gdwarf-3 -fdefault-real-8" )
        elseif( njoy_build_type STREQUAL "coverage" )
            set( njoy_compiler_flags " -Wall -Wextra -Wpedantic -g -gdwarf-3 -fprofile-arcs -ftest-coverage -fno-inline" )
        elseif( njoy_build_type STREQUAL "release" )
            set( njoy_compiler_flags " -Wall -Wextra -Wpedantic -O3 -DNDEBUG" )
            if( NOT no_link_time_optimization )
                set( njoy_compiler_flags "${njoy_compiler_flags} -flto" )
            endif()
        elseif( njoy_build_type STREQUAL "native" )
            set( njoy_compiler_flags " -Wall -Wextra -Wpedantic -O3 -DNDEBUG -march=native" )
            if( NOT no_link_time_optimization )
                set( njoy_compiler_flags "${njoy_compiler_flags} -flto" )
            endif()
        endif()
    else()
        message( WARNING "Compiler vendor not recognized. No compilation flags set" )
    endif()
endif()
if( DEFINED appended_flags AND NOT njoy_no_appended_flags )
    set( njoy_compiler_flags "${njoy_compiler_flags} ${appended_flags}" )
endif()
if( DEFINED njoy_appended_flags )
    set( njoy_compiler_flags "${njoy_compiler_flags} ${njoy_appended_flags}" )
endif()

if ( NOT GIT_EXECUTABLE )
    find_package( Git )
    if ( NOT GIT_FOUND )
        message( FATAL_ERROR "git installation was not found." )
    endif()
endif()
execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )
message( STATUS "" )
message( STATUS "njoy Version: 2016" )
message( STATUS "Git current branch: ${GIT_BRANCH}" )
message( STATUS "Git commit hash: ${GIT_HASH}" )
message( STATUS "" )
message( STATUS "njoy flags: ${njoy_compiler_flags}" ) 
message( STATUS "" ) 
message( STATUS "-----------------------------------------------------------" ) 
message( STATUS "" ) 

add_library( njoy ${njoy_policy}
             src/acecm.f90
             src/acedo.f90
             src/acefc.f90
             src/acepa.f90
             src/acepn.f90
             src/acer.f90
             src/aceth.f90
             src/broadr.f90
             src/ccccr.f90
             src/covr.f90
             src/dtfr.f90
             src/endf.f90
             src/errorr.f90
             src/gaminr.f90
             src/gaspr.f90
             src/graph.f90
             src/groupr.f90
             src/heatr.f90
             src/leapr.f90
             src/locale.f90
             src/mainio.f90
             src/mathm.f90
             src/matxsr.f90
             src/mixr.f90
             src/moder.f90
             src/phys.f90
             src/plotr.f90
             src/powr.f90
             src/purr.f90
             src/reconr.f90
             src/resxsr.f90
             src/samm.f90
             src/thermr.f90
             src/unresr.f90
             src/util.f90
             src/vers.f90
             src/viewr.f90
             src/wimsr.f90 )
separate_arguments( njoy_compiler_flags_list UNIX_COMMAND "${njoy_compiler_flags}" )
foreach( flag IN LISTS njoy_compiler_flags_list )
    target_compile_options( njoy PUBLIC ${flag} )
endforeach( flag )           
set_target_properties( njoy PROPERTIES LINK_FLAGS "${njoy_compiler_flags}" )
target_include_directories( njoy PUBLIC ${PROJECT_BINARY_DIRECTORY} ) 

add_executable( njoy_executable src/main.f90 )

target_link_libraries( njoy_executable PUBLIC njoy )
foreach( flag IN LISTS $njoy_compiler_flags_list )
    target_compile_options( njoy_executable PUBLIC ${flag} )
endforeach( flag )           
set_target_properties( njoy_executable PROPERTIES LINK_FLAGS "${njoy_compiler_flags}" )
set_target_properties( njoy_executable PROPERTIES OUTPUT_NAME njoy )
target_include_directories( njoy_executable PUBLIC ${PROJECT_BINARY_DIRECTORY} )

enable_testing()
add_subdirectory(test)
