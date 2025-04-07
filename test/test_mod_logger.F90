program test_mod_logger

    use, intrinsic :: iso_fortran_env
    use mod_assertion
    use mod_logger
    implicit none

    ! Internal variables
    type(logger_t) :: foo, bar
    integer :: exit_status
    logical :: file_exists

    !//////////////////////////////////////////////////////////////////////////
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test mod_logger start"
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test get logger level"
    call foo%open_log_file()
    call assert_eq(foo%get_logger_level(), LOG_LEVEL%WARNING)
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test set logger level"
    call foo%set_logger_level(LOG_LEVEL%DEBUG)
    call assert_eq(foo%get_logger_level(), LOG_LEVEL%DEBUG)
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test automatic naming"
    file_exists = .false.
    inquire (file="unit10.log", exist=file_exists)
    call assert_true(file_exists)

    call bar%open_log_file()
    file_exists = .false.
    inquire (file="unit11.log", exist=file_exists)
    call assert_true(file_exists)
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Clean up"
    call foo%close_log_file()
    call bar%close_log_file()
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test mod_logger end"
    !--------------------------------------------------------------------------
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

end program test_mod_logger
