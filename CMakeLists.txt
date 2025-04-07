#################
# Project setup #
#################

# Camke version
cmake_minimum_required(VERSION 3.0)
# Project name
project(CafeMol VERSION 3.2.1)
# Language support
enable_language(Fortran)
# Testing support
enable_testing()
# Define executable name
set(CAFEMOL_EXE "cafemol-bin")

###############
# Build types #
###############
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "Check")
endif()

#########################
# Paths and directories #
#########################

# Define directories
set(SRC "${CMAKE_SOURCE_DIR}/src")
set(BIN "${CMAKE_SOURCE_DIR}/bin")
set(TEST "${CMAKE_SOURCE_DIR}/test")
set(CMAKE "${CMAKE_SOURCE_DIR}/cmake")

# Store .mod files in the mods folder
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mods/cafemol")

# Include cmake helper files
include("${CMAKE}/fortran.cmake")

# Add directories
add_subdirectory("${SRC}")
add_subdirectory("${TEST}")
