! setp_mapara_pro_dna_pwm

! @brief Read parameters from pro_dna_pwm.para file. The parameters are &
!      used for protein-DNA sequence specific integrations

subroutine setp_mapara_pro_dna_pwm()
    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: infile, outfile
    use var_setp, only: inpdpwm

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! -------------------------------------------------------------------
    ! Local variables
    integer :: lunout, lunpara
    integer :: iline, nlines, iequa, nequat

    character(4) :: kfind
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: cline
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MSG_ERROR) :: error_message

    ! -------------------------------------------------------------------
    ! Initialize parameters

    inpdpwm%pro_dna_pwm_sigma = 1.0
    inpdpwm%pro_dna_pwm_phi = 10.0

    ! -------------------------------------------------------------------
    ! Initialize file handle
    lunpara = infile%para_pdpwm
    lunout = outfile%data

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (lunpara)

        call ukoto_uiread2(lunpara, lunout, 'pro_dna_pwm_para ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)

        if (kfind /= 'FIND') then
            error_message = 'Error: cannot find "pro_dna_pwm_para" in the pro_dna_pwm.para file'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
            cline = cwkinp(iline)

            call ukoto_uiequa2(cline, nequat, csides)
            do iequa = 1, nequat

                cvalue = 'pro_dna_pwm_sigma'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdpwm%pro_dna_pwm_sigma, cvalue)

                cvalue = 'pro_dna_pwm_phi'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdpwm%pro_dna_pwm_phi, cvalue)

                cvalue = 'energy_unit_pro_dna_pwm'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdpwm%energy_unit_pro_dna_pwm, cvalue)
            end do
        end do
        inpdpwm%pro_dna_pwm_phi = inpdpwm%pro_dna_pwm_phi*F_PI/1.800e2_PREC

#ifdef MPI_PAR
    end if
    call MPI_Bcast(inpdpwm, inpdpwm%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_mapara_pro_dna_pwm
