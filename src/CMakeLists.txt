################
# Binary setup #
################

# Configure the header file
configure_file("config.h.in" "config.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Define directories
set(CAFE_STREAM "${CMAKE_CURRENT_SOURCE_DIR}/mt_stream_cafe")
set(MAIN_SRC "${CMAKE_CURRENT_SOURCE_DIR}/main.F90")

# Define files
file(GLOB FILE_SRC "*.F90" "${CAFE_STREAM}/*.F90")
list(REMOVE_ITEM FILE_SRC
    "${MAIN_SRC}"
    "${CAFE_STREAM}/check_stream.F90"
    "${CAFE_STREAM}/f_jump_coeff.F90"
    "${CAFE_STREAM}/mt_cafe.F90")

add_library(cafemol STATIC ${FILE_SRC})

# Define project executable
add_executable("${CAFEMOL_EXE}" "${MAIN_SRC}")
target_link_libraries("${CAFEMOL_EXE}" cafemol)

# Add flag to specific source files
set_source_files_properties("${CAFE_STREAM}/mt_stream.F90" PROPERTIES COMPILE_FLAGS "-fno-range-check")
set_source_files_properties("${CAFE_STREAM}/gf2xe.F90" PROPERTIES COMPILE_FLAGS "-fno-range-check")

# Add properties to the target binary
set_target_properties("${CAFEMOL_EXE}"
    PROPERTIES
    OUTPUT_NAME cafemol
    RUNTIME_OUTPUT_DIRECTORY "${BIN}"
    COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
)

# Link MPI libraries
if(MPI_PAR)
    target_link_libraries("${CAFEMOL_EXE}" "${MPI_Fortran_LIBRARIES}")
endif()

add_subdirectory(utils)
