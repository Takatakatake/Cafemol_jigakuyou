module mod_file

    use, intrinsic :: iso_fortran_env
    use mod_unit, only: new_unit
    use mod_error
    implicit none

    private

    ! Public parameters
    integer, parameter, public :: LINE_BUFFER_LEN = 120

    ! Public types
    type, public :: file_t
        private
        logical :: is_opened = .false.
        logical :: is_binary = .false.
        logical :: reach_eof = .false.
        integer :: file_unit = -1
        character(len=:), allocatable :: file_mode
        character(len=:), allocatable :: file_name
    contains
        ! Public procedures
        procedure :: is_open
        procedure :: open_file
        procedure :: close_file
        procedure :: read_line
        procedure :: write_line
        procedure :: rewind_file
        procedure :: get_file_unit
        procedure :: get_file_name
        procedure :: write_str
        procedure :: read_str
        generic :: write_int => write_int8, write_int16, write_int32, write_int64
        generic :: read_int => read_int8, read_int16, read_int32, read_int64
        generic :: write_real => write_real32, write_real64, write_real128
        generic :: read_real => read_real32, read_real64, read_real128
        ! Private procedures
        procedure, private :: write_int8
        procedure, private :: write_int16
        procedure, private :: write_int32
        procedure, private :: write_int64
        procedure, private :: read_int8
        procedure, private :: read_int16
        procedure, private :: read_int32
        procedure, private :: read_int64
        procedure, private :: write_real32
        procedure, private :: write_real64
        procedure, private :: write_real128
        procedure, private :: read_real32
        procedure, private :: read_real64
        procedure, private :: read_real128
    end type file_t

contains

    function is_open(this) result(rslt)

        ! Function arguments
        class(file_t), intent(in) :: this
        logical :: rslt

        if (this%is_opened) then
            rslt = .true.
        else
            rslt = .false.
        end if

    end function is_open

    subroutine open_file(this, file_name, mode, new_file, binary, file_unit, &
                         error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        character(len=*), intent(in) :: file_name
        character(len=*), intent(in) :: mode
        logical, optional, intent(in) :: new_file
        logical, optional, intent(in) :: binary
        integer, optional, intent(in) :: file_unit
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        logical :: binary_, in_use, read_mode, write_mode
        integer :: open_status, file_unit_
        character(len=:), allocatable :: action_val
        character(len=:), allocatable :: status_val
        character(len=:), allocatable :: access_val
        character(len=:), allocatable :: file_name_
        type(error_t) :: error_

        ! Check if file is already opened
        if (this%is_open()) return

        if (present(new_file)) then
            if (new_file) then
                status_val = "replace"
            else
                status_val = "unknown"
            end if
        else
            status_val = "unknown"
        end if

        if (present(binary)) then
            binary_ = binary
            if (binary) then
                access_val = "stream"
            else
                access_val = "sequential"
            end if
        else
            binary_ = .false.
            access_val = "sequential"
        end if

        if (present(file_unit)) then
            file_unit_ = file_unit
        else
            file_unit_ = new_unit()
        end if

        ! Check unit value
        inquire (unit=file_unit_, opened=in_use)
        if (file_unit_ < 0 .or. in_use) then
            if (present(error)) then
                error%id = ERROR_BAD_UNIT
                return
            else
                error_%id = ERROR_BAD_UNIT
                call check_error(error_)
            end if
        end if

        ! Parse mode
        read_mode = .false.
        if (index(mode, 'r') /= 0) read_mode = .true.

        write_mode = .false.
        if (index(mode, 'w') /= 0) write_mode = .true.

        if (read_mode .and. .not. write_mode) then
            action_val = "read"
        else if (.not. read_mode .and. write_mode) then
            action_val = "write"
        else if (read_mode .and. write_mode) then
            action_val = "readwrite"
        end if

        ! Remove blanks from file_name
        file_name_ = trim(adjustl(file_name))

        ! Open file
        open (unit=file_unit_, file=file_name_, access=access_val, &
            action=action_val, status=status_val, iostat=open_status)

        ! Check open_status
        if (open_status > 0) then
            if (present(error)) then
                error%id = ERROR_OPEN_FILE
                error%message = "Couldn't open " // file_name_
                return
            else
                error_%id = ERROR_OPEN_FILE
                error_%message = "Couldn't open " // file_name_
                call check_error(error_)
            end if
        end if

        ! Update file property
        this%file_unit = file_unit_
        this%is_opened = .true.
        this%file_name = file_name_
        this%file_mode = action_val
        if (binary_) this%is_binary = .true.

        if (present(error)) error%id = NO_ERROR

    end subroutine open_file

    subroutine close_file(this, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: close_status
        type(error_t) :: error_

        ! If it's already closed, return
        if (.not. this%is_open()) return

        ! Close file
        close (unit=this%file_unit, iostat=close_status)

        ! Check close_status
        if (close_status > 0) then
            if (present(error)) then
                error%id = ERROR_CLOSE_FILE
                error%message = "Couldn't close "//this%file_name
                return
            else
                error_%id = ERROR_CLOSE_FILE
                error_%message = "Couldn't close "//this%file_name
                call check_error(error_)
            end if
        end if

        ! Update file properties
        this%is_opened = .false.
        this%file_unit = -1
        this%reach_eof = .false.

        if (present(error)) error%id = NO_ERROR

    end subroutine close_file

    subroutine read_line(this, line, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        character(len=LINE_BUFFER_LEN), intent(out) :: line
        type(error_t), optional, intent(out) :: error

        ! Internal variables
        integer :: read_status
        logical :: read_mode
        type(error_t) :: error_

        ! Check if binary
        if (this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_[int,real,...] method with binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_[int,real,...] method with binary files."
                call check_error(error_)
            end if
        end if

        ! Check if eof
        if (this%reach_eof) then
            if (present(error)) then
                error%id = ERROR_END_OF_FILE
                return
            else
                error_%id = ERROR_END_OF_FILE
                call check_error(error_)
            end if
        end if

        ! Check mode
        read_mode = .true.
        if (index(this%file_mode, 'r') /= 0) read_mode = .true.

        if (.not. read_mode) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = "Can't read line if mode does not contain 'r'."
                return
            else
                error_%id = ERROR_READ
                error_%message = "Can't read line if mode does not contain 'r'."
                call check_error(error_)
            end if
        end if

        ! Read line
        read (this%file_unit, "(A)", iostat=read_status) line

        ! Check status
        if (read_status < 0) then
            this%reach_eof = .true.
            if (present(error)) then
                error%id = ERROR_END_OF_FILE
            else
                error_%id = ERROR_END_OF_FILE
                call check_error(error_)
            end if
        else if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_line

    subroutine write_line(this, line, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        character(len=*), intent(in) :: line
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        logical :: write_mode
        type(error_t) :: error_

        ! Check if binary
        if (this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_[int,real,...] method with binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_[int,real,...] method with binary files."
                call check_error(error_)
            end if
        end if

        ! Check mode
        write_mode = .false.
        if (index(this%file_mode, 'w') /= 0) write_mode = .true.

        if (.not. write_mode) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = "Can't write line if mode does not contain 'w'."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = "Can't write line if mode does not contain 'w'."
                call check_error(error_)
            end if
        end if

        ! Write line
        write (this%file_unit, "(A)", iostat=write_status) line

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_line

    subroutine rewind_file(this, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        type(error_t) :: error_
        integer :: rewind_status

        ! Check if file is opened
        if (.not. this%is_open()) then
            if (present(error)) then
                error%id = ERROR_REWIND
                error%message = "Can't rewind if file is not opened."
                return
            else
                error_%id = ERROR_REWIND
                error_%message = "Can't rewind if file is not opened."
                call check_error(error_)
            end if
        end if

        ! Rewind file
        rewind (this%file_unit, iostat=rewind_status)

        ! Check status
        if (rewind_status > 0) then
            if (present(error)) then
                error%id = ERROR_REWIND
                return
            else
                error_%id = ERROR_REWIND
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

        ! Update file properties
        this%reach_eof = .false.

    end subroutine rewind_file

    function get_file_unit(this) result(rslt)

        ! Function arguments
        class(file_t), intent(in) :: this
        integer :: rslt

        rslt = this%file_unit

    end function get_file_unit

    function get_file_name(this) result(rslt)

        ! Function arguments
        class(file_t), intent(in) :: this
        character(len=:), allocatable :: rslt

        rslt = trim(adjustl(this%file_name))

    end function get_file_name

    subroutine write_str(this, str, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        character(len=*), intent(in) :: str
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) str

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_str

    subroutine read_str(this, str, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        character(len=*), intent(out) :: str
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) str

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_str

    subroutine write_int8(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int8), intent(in) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) int

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_int8

    subroutine write_int16(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int16), intent(in) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) int

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_int16

    subroutine write_int32(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int32), intent(in) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) int

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_int32

    subroutine write_int64(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int64), intent(in) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) int

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_int64

    subroutine read_int8(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int8), intent(out) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) int

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_int8

    subroutine read_int16(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int16), intent(out) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) int

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_int16

    subroutine read_int32(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int32), intent(out) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) int

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_int32

    subroutine read_int64(this, int, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        integer(int64), intent(out) :: int
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) int

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_int64

    subroutine write_real32(this, float, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        real(real32), intent(in) :: float
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) float

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_real32

    subroutine write_real64(this, float, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        real(real64), intent(in) :: float
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) float

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_real64

    subroutine write_real128(this, float, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        real(real128), intent(in) :: float
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: write_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_WRITE
                error%message = &
                    "Use write_line method with non-binary files."
                return
            else
                error_%id = ERROR_WRITE
                error_%message = &
                    "Use write_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Write value
        write (unit=this%file_unit, iostat=write_status) float

        if (write_status > 0) then
            if (present(error)) then
                error%id = ERROR_WRITE
            else
                error_%id = ERROR_WRITE
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine write_real128

    subroutine read_real32(this, float, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        real(real32), intent(out) :: float
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) float

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_real32

    subroutine read_real64(this, float, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        real(real64), intent(out) :: float
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) float

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_real64

    subroutine read_real128(this, float, error)

        ! Subroutine arguments
        class(file_t), intent(inout) :: this
        real(real128), intent(out) :: float
        type(error_t), optional, intent(inout) :: error

        ! Internal variables
        integer :: read_status
        type(error_t) :: error_

        ! Check if binary
        if (.not. this%is_binary) then
            if (present(error)) then
                error%id = ERROR_READ
                error%message = &
                    "Use read_line method with non-binary files."
                return
            else
                error_%id = ERROR_READ
                error_%message = &
                    "Use read_line method with non-binary files."
                call check_error(error_)
            end if
        end if

        ! Read value
        read (unit=this%file_unit, iostat=read_status) float

        if (read_status > 0) then
            if (present(error)) then
                error%id = ERROR_READ
            else
                error_%id = ERROR_READ
                call check_error(error_)
            end if
        else
            if (present(error)) error%id = NO_ERROR
        end if

    end subroutine read_real128

end module mod_file
