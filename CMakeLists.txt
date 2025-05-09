cmake_minimum_required(VERSION 3.15)

project(NJOY LANGUAGES Fortran)
string(TOLOWER ${CMAKE_PROJECT_NAME} CMAKE_PROJECT_NAME_LOWER) # For install paths like "njoy"

# Include standard modules
include(GNUInstallDirs)
include(FeatureSummary)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported configuration types" FORCE)

# Determine if this is being built as a subproject or the main project
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    set(NJOY_IS_SUBPROJECT FALSE) # This is the main/top-level project
else()
    set(NJOY_IS_SUBPROJECT TRUE)  # This is being included as a subproject
endif()
add_feature_info("Building as Subproject" NJOY_IS_SUBPROJECT "NJOY is being built as part of a larger CMake project.")

# Compiler Version Check
set(NJOY_GNU_MINIMUM_VERSION 5.1)

if(CMAKE_Fortran_COMPILER_ID AND NJOY_${CMAKE_Fortran_COMPILER_ID}_MINIMUM_VERSION)
    if(CMAKE_Fortran_COMPILER_VERSION AND
       CMAKE_Fortran_COMPILER_VERSION VERSION_LESS NJOY_${CMAKE_Fortran_COMPILER_ID}_MINIMUM_VERSION)
        message(FATAL_ERROR
            "${CMAKE_Fortran_COMPILER_ID} version ${CMAKE_Fortran_COMPILER_VERSION} "
            "is less than required minimum ${NJOY_${CMAKE_Fortran_COMPILER_ID}_MINIMUM_VERSION}")
    endif()
endif()

# Build Options
option(NJOY_STRICT "Compile time warnings are converted to errors" OFF)
option(NJOY_COVERAGE "Enable binary instrumentation for test coverage (DEBUG configuration)" OFF)
option(NJOY_PROFILE_GENERATE "Enable binary instrumentation for execution profiling (RELEASE configuration)" OFF)
option(NJOY_PROFILE_USE "Leverage previously generated execution profile for optimization (RELEASE configuration)" OFF)
option(NJOY_ENABLE_LTO "Enable Link Time Optimization (IPO)" ON)
option(NJOY_NONPORTABLE_OPTIMIZATION "Enable optimizations compromising binary portability (RELEASE configuration)" OFF)
option(NJOY_LINK_EXECUTABLE_STATICALLY "Attempt to link the NJOY executable statically" OFF)

# Standard CMake option to control library type (SHARED vs STATIC)
# Default to ON (shared libraries) if not set, to match previous default behavior
# where NJOY_STATIC_SELF_LIBRARY was OFF by default.
if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")
endif()

if(NJOY_PROFILE_GENERATE AND NJOY_PROFILE_USE)
    message(FATAL_ERROR "Cannot both generate (NJOY_PROFILE_GENERATE) and use (NJOY_PROFILE_USE) execution profiles in the same configuration.")
endif()

if(NJOY_PROFILE_GENERATE OR NJOY_PROFILE_USE)
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/profiling")
endif()

# Interprocedural Optimization (LTO)
if(NJOY_ENABLE_LTO)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

# Use PROJECT_BINARY_DIR for the module directory path.
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME_LOWER}_fortran_modules" CACHE PATH "Directory for NJOY Fortran modules")
file(MAKE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")

set(NJOY_APPENDED_COMPILE_FLAGS "" CACHE STRING "User-specified compile flags appended to NJOY targets")
set(NJOY_APPENDED_LINK_FLAGS "" CACHE STRING "User-specified link flags appended to NJOY targets")

# Compiler Specific Flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(NJOY_GNU_COMMON_FLAGS "-Wall;-Wextra")
    set(NJOY_GNU_DEBUG_FLAGS "-O0;-g;-gdwarf-3;-frounding-math;-fsignaling-nans;-fcheck=all;-ffpe-trap=invalid,zero,overflow")
    set(NJOY_GNU_RELEASE_FLAGS "-O3;-DNDEBUG") # CMake's IPO will add -flto if NJOY_ENABLE_LTO is ON
    set(NJOY_GNU_STRICT_FLAGS "-Werror")
    set(NJOY_GNU_COVERAGE_FLAGS "--coverage")
    set(NJOY_GNU_PROFILE_GENERATE_FLAGS "-fprofile-generate=${PROJECT_BINARY_DIR}/profiling")
    set(NJOY_GNU_PROFILE_USE_FLAGS "-fprofile-use=${PROJECT_BINARY_DIR}/profiling")
    set(NJOY_GNU_NONPORTABLE_OPTIMIZATION_FLAGS "-march=native")
    set(NJOY_GNU_STATIC_LINK_FLAGS "-static") # For njoy_executable if NJOY_LINK_EXECUTABLE_STATICALLY is ON
    set(NJOY_GNU_SUBPROJECT_COMPILE_FLAGS "-Wno-maybe-uninitialized;-Wno-conversion;-Wno-unused-parameter;-Wno-compare-reals;-Wno-unused-variable;-Wno-intrinsic-shadow;-Wno-unused-dummy-argument;-Wno-uninitialized;-Wno-unused-function;-Wno-unused-label;-Wno-character-truncation")
    set(NJOY_GNU_BASE_PROJECT_COMPILE_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(WARNING "Intel Fortran compiler flags are placeholders. Please verify and complete them for ${CMAKE_Fortran_COMPILER_ID}.")
    set(NJOY_INTEL_COMMON_FLAGS "-Wall")
    set(NJOY_INTEL_DEBUG_FLAGS "-O0;-g;-debug;-check all;-fpe0;-traceback")
    set(NJOY_INTEL_RELEASE_FLAGS "-O3;-DNDEBUG") # CMake's IPO will add -ipo if NJOY_ENABLE_LTO is ON
    set(NJOY_INTEL_STRICT_FLAGS "-warn errors")
    set(NJOY_INTEL_COVERAGE_FLAGS "-prof-gen=srcpos")
    set(NJOY_INTEL_PROFILE_GENERATE_FLAGS "-prof-gen=${PROJECT_BINARY_DIR}/profiling")
    set(NJOY_INTEL_PROFILE_USE_FLAGS "-prof-use=${PROJECT_BINARY_DIR}/profiling")
    set(NJOY_INTEL_NONPORTABLE_OPTIMIZATION_FLAGS "-xHost")
    set(NJOY_INTEL_STATIC_LINK_FLAGS "-static-intel") # For njoy_executable if NJOY_LINK_EXECUTABLE_STATICALLY is ON
    set(NJOY_INTEL_SUBPROJECT_COMPILE_FLAGS "-diag-disable remark")
    set(NJOY_INTEL_BASE_PROJECT_COMPILE_FLAGS "")
else()
    message(WARNING "Compiler ${CMAKE_Fortran_COMPILER_ID} is not explicitly supported with optimized flags. Using default CMake flags.")
endif()

function(njoy_apply_build_options TARGET_NAME)
    set(COMPILE_OPTIONS_LIST "")
    set(LINK_OPTIONS_LIST "")
    set(COMPILER_FLAGS_PREFIX "NJOY_${CMAKE_Fortran_COMPILER_ID}")

    if(DEFINED ${COMPILER_FLAGS_PREFIX}_COMMON_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST ${${COMPILER_FLAGS_PREFIX}_COMMON_FLAGS})
    endif()
    if(NJOY_STRICT AND DEFINED ${COMPILER_FLAGS_PREFIX}_STRICT_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST ${${COMPILER_FLAGS_PREFIX}_STRICT_FLAGS})
    endif()

    if(NJOY_IS_SUBPROJECT)
        if(DEFINED ${COMPILER_FLAGS_PREFIX}_SUBPROJECT_COMPILE_FLAGS)
            list(APPEND COMPILE_OPTIONS_LIST ${${COMPILER_FLAGS_PREFIX}_SUBPROJECT_COMPILE_FLAGS})
        endif()
    else() # Base project
        if(DEFINED ${COMPILER_FLAGS_PREFIX}_BASE_PROJECT_COMPILE_FLAGS)
            list(APPEND COMPILE_OPTIONS_LIST ${${COMPILER_FLAGS_PREFIX}_BASE_PROJECT_COMPILE_FLAGS})
        endif()
    endif()

    # Debug specific flags
    if(DEFINED ${COMPILER_FLAGS_PREFIX}_DEBUG_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST $<$<CONFIG:Debug>:${${COMPILER_FLAGS_PREFIX}_DEBUG_FLAGS}>)
    endif()
    if(NJOY_COVERAGE AND DEFINED ${COMPILER_FLAGS_PREFIX}_COVERAGE_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST $<$<CONFIG:Debug>:${${COMPILER_FLAGS_PREFIX}_COVERAGE_FLAGS}>)
        list(APPEND LINK_OPTIONS_LIST $<$<CONFIG:Debug>:${${COMPILER_FLAGS_PREFIX}_COVERAGE_FLAGS}>)
    endif()

    # Release specific flags
    if(DEFINED ${COMPILER_FLAGS_PREFIX}_RELEASE_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_RELEASE_FLAGS}>)
    endif()
    if(NJOY_PROFILE_GENERATE AND DEFINED ${COMPILER_FLAGS_PREFIX}_PROFILE_GENERATE_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_PROFILE_GENERATE_FLAGS}>)
        list(APPEND LINK_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_PROFILE_GENERATE_FLAGS}>)
    endif()
    if(NJOY_PROFILE_USE AND DEFINED ${COMPILER_FLAGS_PREFIX}_PROFILE_USE_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_PROFILE_USE_FLAGS}>)
        list(APPEND LINK_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_PROFILE_USE_FLAGS}>)
    endif()

    # Non-portable optimization flags
    if(NJOY_NONPORTABLE_OPTIMIZATION AND DEFINED ${COMPILER_FLAGS_PREFIX}_NONPORTABLE_OPTIMIZATION_FLAGS)
        list(APPEND COMPILE_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_NONPORTABLE_OPTIMIZATION_FLAGS}>)
        if(NJOY_ENABLE_LTO) # If LTO is enabled, -march=native might also be needed at final link time
            list(APPEND LINK_OPTIONS_LIST $<$<CONFIG:Release>:${${COMPILER_FLAGS_PREFIX}_NONPORTABLE_OPTIMIZATION_FLAGS}>)
        endif()
    endif()

    if(NJOY_APPENDED_COMPILE_FLAGS)
        string(REPLACE " " ";" _temp_flags "${NJOY_APPENDED_COMPILE_FLAGS}")
        list(APPEND COMPILE_OPTIONS_LIST ${_temp_flags})
    endif()
    if(NJOY_APPENDED_LINK_FLAGS)
        string(REPLACE " " ";" _temp_flags "${NJOY_APPENDED_LINK_FLAGS}")
        list(APPEND LINK_OPTIONS_LIST ${_temp_flags})
    endif()

    target_compile_options(${TARGET_NAME} PRIVATE ${COMPILE_OPTIONS_LIST})
    if(LINK_OPTIONS_LIST) # Only call if list is not empty
        target_link_options(${TARGET_NAME} PRIVATE ${LINK_OPTIONS_LIST})
    endif()
endfunction()

# RPATH settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) # Important: RPATH is added at install time, not build time
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(NJOY_RPATH_PREFIX "@loader_path")
else()
    set(NJOY_RPATH_PREFIX "$ORIGIN")
endif()
file(RELATIVE_PATH RPATH_TO_LIB_DIR "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_LIBDIR}")
set(CMAKE_INSTALL_RPATH "${NJOY_RPATH_PREFIX}/${RPATH_TO_LIB_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) # Ensures a target's link interface RPATHs are used

# Feature Summary Setup for NJOY Options
add_feature_info("Strict Warnings As Errors" NJOY_STRICT "Compile time warnings are converted to errors.")
add_feature_info("Test Coverage (Debug)" NJOY_COVERAGE "Enable binary instrumentation for test coverage.")
add_feature_info("Profiling Generation (Release)" NJOY_PROFILE_GENERATE "Enable binary instrumentation for execution profiling.")
add_feature_info("Profiling Use (Release)" NJOY_PROFILE_USE "Leverage previously generated execution profile for optimization.")
add_feature_info("CMake LTO/IPO (Release)" NJOY_ENABLE_LTO "Enable Link Time Optimization via CMake's INTERPROCEDURAL_OPTIMIZATION property (Release config).")
add_feature_info("Non-Portable Optimizations (Release)" NJOY_NONPORTABLE_OPTIMIZATION "Enable optimizations compromising binary portability (e.g., -march=native).")
add_feature_info("Build NJOY as Shared Library" BUILD_SHARED_LIBS "Build NJOY as a shared library.")
add_feature_info("Link Executable Statically" NJOY_LINK_EXECUTABLE_STATICALLY "Attempt to link the NJOY executable statically (e.g., using -static).")

# Determine NJOY Library type for messaging
if(BUILD_SHARED_LIBS)
    set(NJOY_LIBRARY_TYPE_MESSAGE "SHARED (BUILD_SHARED_LIBS=ON)")
else()
    set(NJOY_LIBRARY_TYPE_MESSAGE "STATIC (BUILD_SHARED_LIBS=OFF)")
endif()

# Configuration Summary
message(STATUS "")
message(STATUS "-----------------------------------------------------------")
message(STATUS " ${CMAKE_PROJECT_NAME} Project Configuration")
message(STATUS "-----------------------------------------------------------")
message(STATUS "  Project Source Directory: ${PROJECT_SOURCE_DIR}")
message(STATUS "  Project Binary Directory: ${PROJECT_BINARY_DIR}")
message(STATUS "  Installation Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "    Install Binaries to: ${CMAKE_INSTALL_FULL_BINDIR}")
message(STATUS "    Install Libraries to: ${CMAKE_INSTALL_FULL_LIBDIR}")
message(STATUS "    Install Modules/Includes to: ${CMAKE_INSTALL_FULL_INCLUDEDIR}/${CMAKE_PROJECT_NAME_LOWER}")
message(STATUS "  Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  Build types: ${CMAKE_CONFIGURATION_TYPES}")
message(STATUS "  NJOY Library Type: ${NJOY_LIBRARY_TYPE_MESSAGE}")
message(STATUS "  NJOY LTO/IPO Setting: ${NJOY_ENABLE_LTO}")
message(STATUS "-----------------------------------------------------------")
feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
message(STATUS "-----------------------------------------------------------")

set(NJOY_LIBRARY_SOURCES
    "src/acecm.f90" "src/acedo.f90" "src/acefc.f90" "src/acepa.f90" "src/acepn.f90"
    "src/acer.f90" "src/aceth.f90" "src/broadr.f90" "src/ccccr.f90" "src/covr.f90"
    "src/dtfr.f90" "src/endf.f90" "src/errorr.f90" "src/gaminr.f90" "src/gaspr.f90"
    "src/graph.f90" "src/groupr.f90" "src/heatr.f90" "src/leapr.f90" "src/locale.f90"
    "src/mainio.f90" "src/mathm.f90" "src/matxsr.f90" "src/mixr.f90" "src/moder.f90"
    "src/phys.f90" "src/plotr.f90" "src/powr.f90" "src/purr.f90" "src/reconr.f90"
    "src/resxsr.f90" "src/samm.f90" "src/thermr.f90" "src/unresr.f90" "src/util.f90"
    "src/vers.f90" "src/viewr.f90" "src/wimsr.f90"
)

# Add NJOY library. Type (SHARED/STATIC) is controlled by BUILD_SHARED_LIBS.
add_library(${CMAKE_PROJECT_NAME} ${NJOY_LIBRARY_SOURCES})
set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}> # For build tree consumption
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME_LOWER}> # For installed consumption
)
njoy_apply_build_options(${CMAKE_PROJECT_NAME})

if(NOT NJOY_IS_SUBPROJECT)
    add_executable(${CMAKE_PROJECT_NAME}_executable src/main.f90)
    set_target_properties(${CMAKE_PROJECT_NAME}_executable PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME_LOWER})
    target_link_libraries(${CMAKE_PROJECT_NAME}_executable PRIVATE ${CMAKE_PROJECT_NAME})
    njoy_apply_build_options(${CMAKE_PROJECT_NAME}_executable)

    if(NJOY_LINK_EXECUTABLE_STATICALLY)
        if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
            message(WARNING "NJOY_LINK_EXECUTABLE_STATICALLY=ON: Fully static executables are not well-supported on macOS and may lead to issues or failures. Prefer dynamic linking or linking a static NJOY library to a dynamically linked executable.")
        endif()

        set(_compiler_static_link_flags_var "NJOY_${CMAKE_Fortran_COMPILER_ID}_STATIC_LINK_FLAGS")
        if(DEFINED ${_compiler_static_link_flags_var} AND NOT "${${_compiler_static_link_flags_var}}" STREQUAL "")
            message(STATUS "Applying static link flags to executable ${CMAKE_PROJECT_NAME}_executable: ${${_compiler_static_link_flags_var}}")
            target_link_options(${CMAKE_PROJECT_NAME}_executable PRIVATE ${${_compiler_static_link_flags_var}})
        else()
            message(WARNING "NJOY_LINK_EXECUTABLE_STATICALLY is ON, but no specific static link flags (NJOY_${CMAKE_Fortran_COMPILER_ID}_STATIC_LINK_FLAGS) are defined for compiler ${CMAKE_Fortran_COMPILER_ID}. The executable might not be fully static.")
        endif()
    endif()

    enable_testing()
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests")
        add_subdirectory(tests)
    else()
        message(WARNING "Tests directory not found at ${CMAKE_CURRENT_SOURCE_DIR}/tests. Skipping tests setup.")
    endif()
endif()

# Installation
set(NJOY_INSTALLATION_TARGETS ${CMAKE_PROJECT_NAME})
if(NOT NJOY_IS_SUBPROJECT)
    list(APPEND NJOY_INSTALLATION_TARGETS ${CMAKE_PROJECT_NAME}_executable)
endif()

install(TARGETS ${NJOY_INSTALLATION_TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
)
install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME_LOWER}
    FILES_MATCHING PATTERN "*.mod"
    PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)