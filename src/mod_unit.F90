module mod_unit

    use, intrinsic :: iso_fortran_env
    implicit none

    private
    integer(int32), parameter :: MIN_UNIT = 10
    integer(int32), parameter :: MAX_UNIT = 999

    ! Public procedures
    public :: new_unit

contains

    function new_unit(unit_id) result(unit_number)

        ! Function arguments
        integer(int32), intent(out), optional :: unit_id
        integer(int32) :: unit_number

        ! Internal variables
        integer(int32) :: iunit
        logical :: in_use

        do iunit = MIN_UNIT, MAX_UNIT
            inquire (unit=iunit, opened=in_use)
            if (.not. in_use) then
                unit_number = iunit
                if (present(unit_id)) unit_id = iunit
                return
            end if
        end do

        write (output_unit, "(A)") "Error: Unit number not available"
        stop 1

    end function

end module mod_unit
