file(GLOB FILE_SRC "*.F90")

foreach(SRC IN LISTS FILE_SRC)
    get_filename_component(EXE ${SRC} NAME_WE)
    add_executable(${EXE} ${SRC})
    target_link_libraries(${EXE} cafemol)
    set_target_properties(${EXE}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${BIN}"
        COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
    )
endforeach()

