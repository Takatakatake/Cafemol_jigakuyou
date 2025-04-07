project(surface)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

find_library(LIBBOOST_PROGRAM_OPTIONS boost_program_options)

add_subdirectory(src)
