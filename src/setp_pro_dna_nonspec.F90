! setp_pro_dna_nonspec
!> @brief This subroutine is to read and set the parameters for protein-DNA
!> sequence specific interacions.

subroutine setp_pro_dna_nonspec()

    use const_maxsize
    use const_index
    use var_inp, only: infile, outfile
    use var_setp, only: inpdns

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
    inpdns%pro_dna_nonspec_factor = 1.0e30_PREC

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (luninp)
        call ukoto_uiread2(luninp, lunout, 'in_pdns         ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)
        if (kfind /= 'FIND') then
            error_message = 'Error: cannot find "in_pdns" field in the input file'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
            call ukoto_uiequa2(cwkinp(iline), nequat, csides)

            do iequa = 1, nequat
                cvalue = 'pro_dna_nonspec_factor'
                call ukoto_rvalue2(lunout, csides(1, iequa), inpdns%pro_dna_nonspec_factor, cvalue)
            end do
        end do

        ! -----------------------------------------------------------------
        ! checking input variables
        if (inpdns%pro_dna_nonspec_factor >= 1.0e5_PREC) then
            error_message = 'Error: invalid value for pro_dna_nonspec_factor (should be <= 100000).'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef MPI_PAR
    end if
    call MPI_Bcast(inpdns, inpdns%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_pro_dna_nonspec

