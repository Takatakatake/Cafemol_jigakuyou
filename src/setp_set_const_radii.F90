! setp_set_const_radii
!> @brief Read the "<<<< set_const_radii" field from the .inp file &
!>        and set exv radii for some groups

!#define _DEBUG
subroutine setp_set_const_radii()

    use const_maxsize
    use const_index
    use var_inp, only: infile, outfile
    use var_setp, only: inmisc
    use var_struct, only: exv_radius_mp

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! --------------------------------------------------------------------
    ! local variables
    integer :: icol, icol_ini
    integer :: luninp, lunout
    integer :: iline, nlines, iequa, nequat
    integer :: inunit(2), instate
    integer :: imp, jmp, iexv

    character(4) :: kfind
    character(12) :: char12
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MXCOLM) :: ctmp00
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MSG_ERROR) :: error_message

    ! --------------------------------------------------------------------
    luninp = infile%inp
    lunout = outfile%data

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (luninp)
        call ukoto_uiread2(luninp, lunout, 'set_const_radii ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)
        if (kfind /= 'FIND') then
            error_message = 'Error: cannot find "set_const_radii" field in the input file'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
            call ukoto_uiequa2(cwkinp(iline), nequat, csides)
            do iequa = 1, nequat
                cvalue = 'const_exv_radii'
                call ukoto_rvalue2(lunout, csides(1, iequa), inmisc%const_exv_radii, cvalue)
            end do
        end do

        do iline = 1, nlines
            ctmp00 = cwkinp(iline)

            if (ctmp00(1:13) == 'SET_EXV_RADII') then
                icol_ini = 15
                do icol = icol_ini, CARRAY_MXCOLM
                    if (ctmp00(icol:icol) == ')') exit
                end do

                read (ctmp00(icol_ini:icol - 1), *) char12
                call util_unitstate(char12, inunit, instate)
                write (lunout, '(3a)') '---reading residue for setting constant exv radius: ', ctmp00(1:13), ' ', trim(char12)

                imp = inunit(1)
                jmp = inunit(2)

                do iexv = imp, jmp
                    exv_radius_mp(iexv) = inmisc%const_exv_radii*0.5e0_PREC
                end do
            end if
        end do

        ! ----- CHECK =====
        ! write (*,*) " EXV_RADII check: "
        ! do imp = 1, nmp_all
        !    write (*,*) imp, "    ", exv_radius_mp(imp)
        ! end do

#ifdef MPI_PAR
    end if

    call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(exv_radius_mp, MXMP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

#endif

end subroutine setp_set_const_radii
