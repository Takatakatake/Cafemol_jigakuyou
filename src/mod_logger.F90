module mod_logger

    use, intrinsic :: iso_fortran_env
    use mod_unit, only: new_unit
    use mod_error
    use mod_file
    implicit none

    private

    ! Private parameters
    integer, parameter :: BAD_UNIT = -1

    ! Private types
    type :: log_level_t
        integer :: DEBUG
        integer :: INFO
        integer :: WARNING
        integer :: ERROR
        integer :: CRITICAL
    end type log_level_t

    ! Public variables
    type(log_level_t), parameter, public :: LOG_LEVEL = &
        log_level_t(1, 2, 3, 4, 5)

    ! Public types
    type, public :: logger_t
        private
        integer :: logger_level = LOG_LEVEL%WARNING
        type(file_t) :: log_file
    contains
        ! I/O control
        procedure :: open_log_file
        procedure :: close_log_file
        ! Level property setter/getter procedures
        procedure :: set_logger_level
        procedure :: get_logger_level
        ! Logging procedures
        procedure :: log_debug
        procedure :: log_info
        procedure :: log_warning
        procedure :: log_error
        procedure :: log_critical
    end type logger_t

contains

    subroutine open_log_file(this, log_unit, path, file_name)

        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        integer, optional, intent(in) :: log_unit
        character(len=*), optional, intent(in) :: path
        character(len=*), optional, intent(in) :: file_name

        ! Internal variables
        integer :: open_status, log_unit_
        character(len=:), allocatable :: path_
        character(len=:), allocatable :: file_name_
        character(len=100) :: buffer
        type(error_t) :: error

        ! Check optional variables
        if (present(log_unit)) then
            log_unit_ = log_unit
        else
            log_unit_ = new_unit()
        end if

        if (present(path)) then
            path_ = path
        else
            path_ = "./"
        end if

        if (present(file_name)) then
            file_name_ = file_name
        else
            write (buffer, "(A,I0,A)") "unit", log_unit_, ".log"
            file_name_ = trim(adjustl(buffer))
        end if

        ! Open file
        call this%log_file%open_file(file_name=file_name_, mode="w", &
            new_file=.true., file_unit=log_unit_, error=error)
        call check_error(error)

        return

    end subroutine open_log_file

    subroutine close_log_file(this)

        ! Subroutine arguments
        class(logger_t), intent(inout) :: this

        ! Internal variables
        integer :: close_status
        type(error_t) :: error

        ! Close file
        call this%log_file%close_file(error)
        call check_error(error)

    end subroutine close_log_file

    subroutine set_logger_level(this, logger_level)

        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        integer, intent(in) :: logger_level
        type(error_t) :: error

        ! Check logger level value
        if (logger_level < 1 .or. logger_level > 5) then
            error%id = ERROR_BAD_VALUE
            call check_error(error)
        end if

        ! Set logger level
        this%logger_level = logger_level

        return

    end subroutine set_logger_level

    function get_logger_level(this) result(logger_level)

        ! Function arguments
        class(logger_t), intent(in) :: this
        integer :: logger_level

        ! Get level
        logger_level = this%logger_level

    end function get_logger_level

    subroutine log_message(message, message_level, logger)

        ! Subroutine arguments
        character(len=*), intent(in) :: message
        integer, intent(in) :: message_level
        type(logger_t), intent(inout) :: logger

        ! Internal variables
        type(error_t) :: error
        integer :: date_time_val(8)
        character(len=20) :: buffer
        character(len=60) :: name
        character(len=:), allocatable :: header
        character(len=:), allocatable :: line

        ! Check level
        if (logger%get_logger_level() > message_level) return 

        ! Get date and time
        call date_and_time(values=date_time_val)

        ! Write date and time info
        write (buffer, "(I4,A1,5(I0.2,A1))") &
            date_time_val(1), "-", & ! Year
            date_time_val(2), "-", & ! Month
            date_time_val(3), " ", & ! Day
            date_time_val(5), ":", & ! Hour
            date_time_val(6), ":", & ! Minutes
            date_time_val(7), " "    ! Seconds

        ! Get host name
        call hostnm(name)

        ! Prepare header
        select case(message_level)
            case(LOG_LEVEL%DEBUG)
                header = buffer // trim(adjustl(name)) // " DEBUG "
            case(LOG_LEVEL%INFO)
                header = buffer // trim(adjustl(name)) // " INFO "
            case(LOG_LEVEL%WARNING)
                header = buffer // trim(adjustl(name)) // " WARNING "
            case(LOG_LEVEL%ERROR)
                header = buffer // trim(adjustl(name)) // " ERROR "
            case(LOG_LEVEL%CRITICAL)
                header = buffer // trim(adjustl(name)) // " CRITICAL "
            case default
                error%id = ERROR_BAD_VALUE
                call check_error(error)
        end select

        ! Write log
        line = header // trim(adjustl(message))
        call logger%log_file%write_line(line, error)
        call check_error(error)

        return

    end subroutine log_message

    subroutine log_debug(this, message)
        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        character(len=*), intent(in) :: message

        call log_message(message, LOG_LEVEL%DEBUG, this)

    end subroutine log_debug

    subroutine log_info(this, message)
        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        character(len=*), intent(in) :: message

        call log_message(message, LOG_LEVEL%INFO, this)

    end subroutine log_info

    subroutine log_warning(this, message)
        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        character(len=*), intent(in) :: message

        call log_message(message, LOG_LEVEL%WARNING, this)

    end subroutine log_warning

    subroutine log_error(this, message)
        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        character(len=*), intent(in) :: message

        call log_message(message, LOG_LEVEL%ERROR, this)

    end subroutine log_error

    subroutine log_critical(this, message)
        ! Subroutine arguments
        class(logger_t), intent(inout) :: this
        character(len=*), intent(in) :: message

        call log_message(message, LOG_LEVEL%CRITICAL, this)

    end subroutine log_critical

end module mod_logger
