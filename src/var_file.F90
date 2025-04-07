module var_file

    use, intrinsic :: iso_fortran_env
    use mod_file

    implicit none

    ! Public variables
    type(file_t) :: input_file
    type(file_t) :: restart_file

end module var_file
