################
# Binary setup #
################

# Define directories
set(CAFE_MT_STREAM "${SRC}/mt_stream_cafe")

# Define files
file(GLOB CAFE_MT_STREAM_SRC "${CAFE_MT_STREAM}/*.F90")
list(REMOVE_ITEM CAFE_MT_STREAM_SRC
    "${CAFE_MT_STREAM}/check_stream.F90"
    "${CAFE_MT_STREAM}/f_jump_coeff.F90"
    "${CAFE_MT_STREAM}/mt_cafe.F90")

# Add flag to specific source files
set_source_files_properties("${CAFE_MT_STREAM}/mt_stream.F90" PROPERTIES 
    COMPILE_FLAGS "-fno-range-check")
set_source_files_properties("${CAFE_MT_STREAM}/gf2xe.F90" PROPERTIES 
    COMPILE_FLAGS "-fno-range-check")

# ------------------------------------------------------------------------------
# Test for mod_unit
# ------------------------------------------------------------------------------
# Define test executable
add_executable(test_mod_unit
    "${SRC}/mod_assertion.F90"
    "${SRC}/mod_error.F90"
    "${SRC}/mod_unit.F90"
    test_mod_unit.F90
)

# Store .mod files in the mods folder
set_target_properties(test_mod_unit PROPERTIES Fortran_MODULE_DIRECTORY 
    "${CMAKE_BINARY_DIR}/mods/test_mod_unit")

# Add test
add_test("Test_mod_unit" test_mod_unit)

# ------------------------------------------------------------------------------
# Test for mod_file
# ------------------------------------------------------------------------------
# Define test executable
add_executable(test_mod_file
    "${SRC}/mod_assertion.F90"
    "${SRC}/mod_error.F90"
    "${SRC}/mod_unit.F90"
    "${SRC}/mod_file.F90"
    test_mod_file.F90
)

# Store .mod files in the mods folder
set_target_properties(test_mod_file PROPERTIES Fortran_MODULE_DIRECTORY 
    "${CMAKE_BINARY_DIR}/mods/test_mod_file")

# Add test
add_test("Test_mod_file" test_mod_file)

# ------------------------------------------------------------------------------
# Test for mod_logger
# ------------------------------------------------------------------------------
# Define test executable
add_executable(test_mod_logger
    "${SRC}/mod_assertion.F90"
    "${SRC}/mod_error.F90"
    "${SRC}/mod_unit.F90"
    "${SRC}/mod_file.F90"
    "${SRC}/mod_logger.F90"
    test_mod_logger.F90
)

# Store .mod files in the mods folder
set_target_properties(test_mod_logger PROPERTIES Fortran_MODULE_DIRECTORY
    "${CMAKE_BINARY_DIR}/mods/test_mod_logger")

# Add test
add_test("Test_mod_logger" test_mod_logger)

