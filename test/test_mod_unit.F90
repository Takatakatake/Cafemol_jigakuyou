program test_mod_unit

    use, intrinsic :: iso_fortran_env
    use mod_assertion
    use mod_error
    use mod_unit
    implicit none

    ! Internal variables
    integer(int32) :: unit_id, exit_status

    !//////////////////////////////////////////////////////////////////////////
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test mod_unit start"
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test units start at 10"
    unit_id = new_unit()
    call assert_eq(unit_id, 10)
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test units are not repeated"
    open (unit=unit_id, file="foo.txt")
    call assert_eq(new_unit(), 11)
    close (unit=unit_id)
    !--------------------------------------------------------------------------
    write (output_unit, "(A)") "Test mod_unit end"
    !--------------------------------------------------------------------------
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

end program test_mod_unit
