! setp_mapara_pro_dna_nonspec

! @brief Read parameters from pro_dna_nonspec.para file. The parameters are &
!      used for protein-DNA sequence specific integrations

subroutine setp_mapara_pro_dna_nonspec()
    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: infile, outfile
    use var_setp, only: inpdns

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

    inpdns%pro_dna_nonspec_sigma = 1.0
    inpdns%pro_dna_nonspec_phi = 10.0

    ! -------------------------------------------------------------------
    ! Initialize file handle
    lunpara = infile%para_pdns
    lunout = outfile%data

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (lunpara)

        call ukoto_uiread2(lunpara, lunout, 'pro_dna_nonspec_para ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)

        if (kfind /= 'FIND') then
            error_message = 'Error: cannot find "pro_dna_nonspec_para" in the pro_dna_nonspec.para file'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
            cline = cwkinp(iline)

            call ukoto_uiequa2(cline, nequat, csides)
            do iequa = 1, nequat

                cvalue = 'pro_dna_nonspec_sigma'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdns%pro_dna_nonspec_sigma, cvalue)

                cvalue = 'pro_dna_nonspec_phi'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdns%pro_dna_nonspec_phi, cvalue)

                cvalue = 'energy_unit_pro_dna_nonspec'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdns%energy_unit_pro_dna_nonspec, cvalue)
            end do
        end do
        inpdns%pro_dna_nonspec_phi = inpdns%pro_dna_nonspec_phi*F_PI/1.800e2_PREC

#ifdef MPI_PAR
    end if
    call MPI_Bcast(inpdns, inpdns%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_mapara_pro_dna_nonspec
