subroutine initialize_logger()

    use, intrinsic :: iso_fortran_env
    use var_logger

    implicit none

    ! Initialize logger
    call logger%set_logger_level(LOG_LEVEL%INFO)
    ! Open log file
    call logger%open_log_file(log_unit=998, path="./", file_name="info.log")

#ifdef LOG_DEBUG
    ! Initialize debug logger
    call debug_logger%set_logger_level(LOG_LEVEL%DEBUG)
    ! Open log file
    call debug_logger%open_log_file(log_unit=999, path="./", file_name="debug.log")
#endif

end subroutine initialize_logger
