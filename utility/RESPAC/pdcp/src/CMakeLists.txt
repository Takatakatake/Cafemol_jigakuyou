set(CPP_FILES
    Charge.cpp Dx.cpp Input.cpp Parameter.cpp Pdc.cpp Pqr.cpp util.cpp)

add_definitions("-O3 -std=c++11 -Wall")

add_executable(pdcp
    pdcp.cpp ${CPP_FILES})

target_link_libraries(pdcp
    ${LIBBOOST_PROGRAM_OPTIONS})
