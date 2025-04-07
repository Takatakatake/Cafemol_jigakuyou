program test_mod_file

    use, intrinsic :: iso_fortran_env
    use mod_assertion
    use mod_error
    use mod_file
    implicit none

    ! Internal variables
    type(file_t) :: foo, bar
    type(error_t) :: error
    integer :: exit_status
    character(len=LINE_BUFFER_LEN) :: line

    !//////////////////////////////////////////////////////////////////////////
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test mod_file start"
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test initial status"
    call assert_false(foo%is_open())
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test open file"
    call foo%open_file(file_name="../../test/files/sample1.txt", mode="r")
    call assert_true(foo%is_open())
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test read line"
    write (output_unit, "(A)") "Read first line"
    call foo%read_line(line, error)
    call assert_eq(error%id, NO_ERROR)
    call assert_eq(line, "abcde")

    write (output_unit, "(A)") "Read second line"
    call foo%read_line(line, error)
    call assert_eq(error%id, NO_ERROR)
    call assert_eq(line, "12345")

    write (output_unit, "(A)") "End of file"
    call foo%read_line(line, error)
    call assert_eq(error%id, ERROR_END_OF_FILE)

    write (output_unit, "(A)") "Keep reading eof"
    call foo%read_line(line, error)
    call assert_eq(error%id, ERROR_END_OF_FILE)

    write (output_unit, "(A)") "Write an only read file"
    call foo%write_line(line, error)
    call assert_eq(error%id, ERROR_WRITE)

    write (output_unit, "(A)") "Rewind file"
    call foo%rewind_file(error)
    call assert_eq(error%id, NO_ERROR)

    write (output_unit, "(A)") "Read first line"
    call foo%read_line(line, error)
    call assert_eq(error%id, NO_ERROR)
    call assert_eq(line, "abcde")
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test close file"
    call foo%close_file()
    call assert_eq(foo%is_open(), .false.)
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test write file"
    write (output_unit, "(A)") "Open file"
    call foo%open_file(file_name="../../test/files/foo.txt", mode="w", &
                       new_file=.true.)
    call assert_eq(foo%is_open(), .true.)

    write (output_unit, "(A)") "Write line"
    call foo%write_line("Hello, world!", error)
    call assert_eq(error%id, NO_ERROR)

    write (output_unit, "(A)") "Read an only write file"
    call foo%read_line(line, error)
    call assert_eq(error%id, ERROR_READ)

    write (output_unit, "(A)") "Close file"
    call foo%close_file()
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test read/write file"
    write (output_unit, "(A)") "Open file"
    call bar%open_file(file_name="../../test/files/bar.txt", mode="rw", &
                       new_file=.true.)
    call assert_true(bar%is_open())

    write (output_unit, "(A)") "Write line"
    call bar%write_line("Hello, world!", error)
    call assert_eq(error%id, NO_ERROR)

    write (output_unit, "(A)") "Rewind file"
    call bar%rewind_file(error)
    call assert_eq(error%id, NO_ERROR)

    write (output_unit, "(A)") "Read line"
    call bar%read_line(line, error)
    call assert_eq(error%id, NO_ERROR)
    call assert_eq(line, "Hello, world!")

    write (output_unit, "(A)") "Close file"
    call bar%close_file()
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test mod_file end"
    !--------------------------------------------------------------------------
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

end program test_mod_file
