! setp_pro_dna_pwm
!> @brief This subroutine is to read and set the parameters for protein-DNA
!> sequence specific interacions.

subroutine setp_pro_dna_pwm()

    use const_maxsize
    use const_index
    use var_inp, only: infile, outfile
    use var_setp, only: inpdpwm

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! --------------------------------------------------------------------
    ! local variables
    integer :: luninp, lunout
    integer :: iline, nlines, iequa, nequat

    character(4) :: kfind
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MSG_ERROR) :: error_message

    ! --------------------------------------------------------------------
    luninp = infile%inp
    lunout = outfile%data
    inpdpwm%pro_dna_pwm_en0 = 0.0e0_PREC
    inpdpwm%pro_dna_pwm_factor = 1.0e0_PREC

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (luninp)
        call ukoto_uiread2(luninp, lunout, 'in_pdpwm         ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)
        if (kfind /= 'FIND') then
            ! error_message = 'Error: cannot find "in_pdpwm" field in the input file'
            ! call util_error(ERROR%STOP_ALL, error_message)
        else

            do iline = 1, nlines
                call ukoto_uiequa2(cwkinp(iline), nequat, csides)

                do iequa = 1, nequat
                    cvalue = 'pro_dna_pwm_en0'
                    call ukoto_rvalue2(lunout, csides(1, iequa), inpdpwm%pro_dna_pwm_en0, cvalue)

                    cvalue = 'pro_dna_pwm_factor'
                    call ukoto_rvalue2(lunout, csides(1, iequa), inpdpwm%pro_dna_pwm_factor, cvalue)
                end do
            end do
        end if

        ! -----------------------------------------------------------------
        ! checking input variables
        if (inpdpwm%pro_dna_pwm_en0 >= 1.0e5_PREC) then
            error_message = 'Error: invalid value for pro_dna_pwm_en0 (should be <= 100000).'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (inpdpwm%pro_dna_pwm_factor >= 1.0e5_PREC) then
            error_message = 'Error: invalid value for pro_dna_pwm_factor (should be <= 100000).'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef MPI_PAR
    end if
    call MPI_Bcast(inpdpwm, inpdpwm%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_pro_dna_pwm
