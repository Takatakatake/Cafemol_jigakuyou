project(pdcp)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

find_package(Eigen3 REQUIRED)
include_directories("${EIGEN3_INCLUDE_DIR}")

find_library(LIBBOOST_PROGRAM_OPTIONS boost_program_options)

add_subdirectory(src)
