!setp_modified_muca
!> @brief Reads "<<<< modified_muca" field. Data are stored in  &
!>        "inmmc" struct.

subroutine setp_modified_muca()

    use const_maxsize
    use const_index
    use var_inp, only: infile, outfile
    use var_setp, only: inmmc

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! --------------------------------------------------------------------
    ! intent(out) :: xbox, ybox, zbox

    ! --------------------------------------------------------------------
    ! local variables
    integer :: luninp, lunout
    integer :: iline, nlines, iequa, nequat
    character(4) :: kfind
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MSG_ERROR) :: error_message

    character(CARRAY_MXCOLM) :: ctmp00

    ! --------------------------------------------------------------------
    luninp = infile%inp
    lunout = outfile%data
    inmmc%em_depth = -1.0
    inmmc%em_mid = -1.0
    inmmc%em_sigma = -1.0

    ! for Wang-Landau MuCa
    inmmc%muca_cv_min = -1200.0
    inmmc%muca_cv_max = 800.0
    inmmc%muca_increment = 0.00001
    inmmc%muca_stop = 2000000000
    inmmc%muca_bins = 100

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (luninp)
        call ukoto_uiread2(luninp, lunout, 'modified_muca     ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)
        if (kfind /= 'FIND') then
            error_message = 'Error: cannot find "modified_muca" field in the input file'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
            call ukoto_uiequa2(cwkinp(iline), nequat, csides)
            ctmp00 = cwkinp(iline)

            do iequa = 1, nequat
                if (inmmc%i_modified_muca == 1) then
                    cvalue = 'em_depth'
                    call ukoto_rvalue2(lunout, csides(1, iequa), &
                                       inmmc%em_depth, cvalue)

                    cvalue = 'em_mid'
                    call ukoto_rvalue2(lunout, csides(1, iequa), &
                                       inmmc%em_mid, cvalue)

                    cvalue = 'em_sigma'
                    call ukoto_rvalue2(lunout, csides(1, iequa), &
                                       inmmc%em_sigma, cvalue)

                end if

                if (inmmc%i_modified_muca == 2) then

                    ! Wang-Landau MuCa
                    cvalue = 'muca_cv_max'
                    call ukoto_rvalue2(lunout, csides(1, iequa), &
                                       inmmc%muca_cv_max, cvalue)
                    cvalue = 'muca_cv_min'
                    call ukoto_rvalue2(lunout, csides(1, iequa), &
                                       inmmc%muca_cv_min, cvalue)
                    cvalue = 'muca_increment'
                    call ukoto_rvalue2(lunout, csides(1, iequa), &
                                       inmmc%muca_increment, cvalue)
                    cvalue = 'muca_stop'
                    call ukoto_ivalue2(lunout, csides(1, iequa), &
                                       inmmc%muca_stop, cvalue)
                    cvalue = 'muca_bins'
                    call ukoto_ivalue2(lunout, csides(1, iequa), &
                                       inmmc%muca_bins, cvalue)
                    cvalue = 'i_muca_init'
                    call ukoto_ivalue2(lunout, csides(1, iequa), &
                                       inmmc%i_muca_init, cvalue)

                    if (csides(1, iequa) == 'muca_init') then
                        inmmc%muca_init = csides(2, iequa)
                    endif

                end if
            end do
        end do

        ! -----------------------------------------------------------------
        ! checking input variables
        if (inmmc%i_modified_muca == 1) then
            if (inmmc%em_depth <= 0.0) then
                error_message = 'Error: invalid value for em_depth'
                call util_error(ERROR%STOP_ALL, error_message)
            else if (inmmc%em_sigma <= 0.0) then
                error_message = 'Error: invalid value for em_sigma'
                call util_error(ERROR%STOP_ALL, error_message)
            end if
        end if

        if (inmmc%i_modified_muca == 2) then
            if (inmmc%muca_bins <= 0) then
                error_message = 'Error: muca_bins must be larger than 0'
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            if (inmmc%muca_cv_max < inmmc%muca_cv_min) then
                error_message = 'Error: muca_cv_max must be larger than muca_cv_min'
                call util_error(ERROR%STOP_ALL, error_message)
            end if
        end if

#ifdef MPI_PAR
    end if

    call MPI_Bcast(inmmc, inmmc%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)

#endif

end subroutine setp_modified_muca

