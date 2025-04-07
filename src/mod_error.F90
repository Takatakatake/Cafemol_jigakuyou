module mod_error

    use, intrinsic :: iso_fortran_env
    implicit none

    private

    ! Public variables
    integer, parameter, public :: NO_ERROR = 0
    integer, parameter, public :: ERROR_ALLOCATION = 1
    integer, parameter, public :: ERROR_DEALLOCATION = 2
    integer, parameter, public :: ERROR_OPEN_FILE = 3
    integer, parameter, public :: ERROR_CLOSE_FILE = 4
    integer, parameter, public :: ERROR_READ = 5
    integer, parameter, public :: ERROR_WRITE = 6
    integer, parameter, public :: ERROR_REWIND = 7
    integer, parameter, public :: ERROR_NOT_FOUND = 8
    integer, parameter, public :: ERROR_MAX_UNIT = 9
    integer, parameter, public :: ERROR_BAD_UNIT = 10
    integer, parameter, public :: ERROR_BAD_CALL = 11
    integer, parameter, public :: ERROR_BAD_VALUE = 12
    integer, parameter, public :: ERROR_BAD_FORMAT = 13
    integer, parameter, public :: ERROR_END_OF_FILE = 14
    integer, parameter, public :: ERROR_UNKNOWN = 15

    ! Public types
    type, public :: error_t
        integer :: id
        character(len=:), allocatable :: message
    end type error_t

    ! Public methods
    public :: check_error

contains

    subroutine check_error(error, line_number, file_name)

        ! Subroutine arguments
        type(error_t), intent(in) :: error
        integer, optional, intent(in) :: line_number
        character(len=*), optional, intent(in) :: file_name

        if (error%id == NO_ERROR) return

        if (present(line_number) .and. present(file_name)) then
            write (output_unit, "(A,I0,A,A)") "An error ocurred at line ", &
                line_number, " in file ", file_name
        else
            write (output_unit, "(A)") "An error ocurred."
        end if

        write (output_unit, "(A,I0)") "Error id: ", error%id

        if (len_trim(error%message) == 0) then
            write (output_unit, "(A)") explain_me_this(error)
        else
            write (output_unit, "(A)") error%message
        end if

        stop 1

    end subroutine check_error

    function explain_me_this(error) result(error_message)

        ! Subroutine arguments
        type(error_t), intent(in) :: error
        character(len=:), allocatable :: error_message

        select case (error%id)
        case (ERROR_ALLOCATION)
            error_message = "Error - Couldn't allocate array"
        case (ERROR_DEALLOCATION)
            error_message = "Error - Couldn't deallocate array"
        case (ERROR_OPEN_FILE)
            error_message = "Error - Couldn't open file"
        case (ERROR_CLOSE_FILE)
            error_message = "Error - Couldn't close file"
        case (ERROR_READ)
            error_message = "Error - Couldn't read value"
        case (ERROR_WRITE)
            error_message = "Error - Couldn't write value"
        case (ERROR_REWIND)
            error_message = "Error - Couldn't rewind unit"
        case (ERROR_NOT_FOUND)
            error_message = "Error - Couldn't find value"
        case (ERROR_MAX_UNIT)
            error_message = "Error - No more units are available"
        case (ERROR_BAD_UNIT)
            error_message = "Error - Bad unit number"
        case (ERROR_BAD_CALL)
            error_message = "Error - Bad subroutine call"
        case (ERROR_BAD_VALUE)
            error_message = "Error - Bad value encountered"
        case (ERROR_BAD_FORMAT)
            error_message = "Error - Bad format encountered"
        case (ERROR_UNKNOWN)
            error_message = "Error - Uknown error ocurred"
        case default
            error_message = "Error - Bad error id"
        end select

    end function explain_me_this

end module mod_error
