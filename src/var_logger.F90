module var_logger

    use, intrinsic :: iso_fortran_env
    use mod_logger

    implicit none

    ! Public variables
    character(len=100) :: log_msg_buffer
    type(logger_t) logger
#ifdef LOG_DEBUG
    type(logger_t) debug_logger
#endif

end module var_logger
