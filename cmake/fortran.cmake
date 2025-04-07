#####################
# Compiling options #
#####################

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    # Debug flags
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fpe0")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -debug all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn all")
    # Check flags
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -fpe0 ")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -traceback")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -check all")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -DREAL_128")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -prec-div")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -prec-sqrt")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -gen-interfaces nosource")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -warn interfaces")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -fp-model strict")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -fp-speculation=off")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    # General flags
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0")
    # Debug flags
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=invalid,zero,overflow")
    # Check flags
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -g")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -Wall")
    set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} -fbounds-check")
endif()

# General flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Build Types
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DLOG_DEBUG)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

IF(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
ENDIF(CMAKE_BUILD_TYPE MATCHES Release)

# Architecture dependent flags
if(APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000")
endif()

# MPI option
option(MPI_PAR "Compile with MPI" OFF)
if(MPI_PAR)
    # MPI setup
    find_package(MPI REQUIRED)
    #include_directories("${MPI_Fortran_INCLUDE_PATH}")
    add_definitions(-DMPI_PAR)
endif()

# OpenMP
option(OMP_PAR "Compile with OpenMP" OFF)
if(OMP_PAR)
    # Add OpenMP flags
    if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
    elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()

# TIME option
option(TIME "Measure processing time" ON)
if(TIME)
    add_definitions(-DTIME)
endif()

# Precision
set(PREC AUTO CACHE STRING "Set precision for REAL")
set_property(CACHE PREC PROPERTY STRINGS AUTO REAL32 REAL64 REAL128)
if(PREC STREQUAL "AUTO")
    add_definitions(-DREAL_64)
    message("-- Compiling CafeMol with REAL64 precision")
elseif(PREC STREQUAL "REAL32")
    add_definitions(-DREAL_32)
    message("-- Compiling CafeMol with REAL32 precision")
elseif(PREC STREQUAL "REAL64")
    add_definitions(-DREAL_64)
    message("-- Compiling CafeMol with REAL64 precision")
elseif(PREC STREQUAL "REAL128")
    add_definitions(-DREAL_128)
    message("-- Compiling CafeMol with REAL128 precision")
endif()
