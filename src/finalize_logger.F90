subroutine finalize_logger()

    use, intrinsic :: iso_fortran_env
    use var_logger

    implicit none

    ! Close logger
    call logger%close_log_file()

#ifdef LOG_DEBUG
    ! Close logger
    call debug_logger%close_log_file()
#endif

end subroutine finalize_logger
