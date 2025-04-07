set(CPP_FILES
    Grid.cpp Input.cpp Pqr.cpp util.cpp)

add_definitions("-O3 -std=c++11 -Wall")

add_executable(surface
    surface.cpp ${CPP_FILES})

target_link_libraries(surface
    ${LIBBOOST_PROGRAM_OPTIONS})
