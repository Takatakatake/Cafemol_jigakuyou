! setp_stage
!> @brief This subroutine is to read and set the parameters for stage.
! vim: tabstop=3 shiftwidth=3

subroutine setp_stage()

    use const_maxsize
    use const_index
    use var_inp, only: infile, outfile
    use var_stage

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    integer      :: inp_id, out_id ! file id
    character(4) :: block_found    ! ukoto status
    integer      :: iline, nlines  ! # of lines and counter
    character(CARRAY_MXCOLM) :: input_lines(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: bufline

    ! to parse `x = y`
    integer                  :: i_equation, n_equation
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MXCOLM) :: cvalue

    ! a buffer for error message
    character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
    write (*, *) "DEBUG: subroutine setp_stage start"
#endif

    inp_id = infile%inp
    out_id = outfile%data

    ! set default values
    stage%group_id = -1
    stage%eps = 1.2
    stage%height = 0.0
    stage%cutoff_ratio = 2.5
    stage%distance_offset = 2.0
    
#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        !--------------------------------------------------------------------------
        ! find <<<< afm_fitting block
        rewind (inp_id)
        call ukoto_uiread2(inp_id, out_id, 'stage           ', block_found, &
                           CARRAY_MXLINE, nlines, input_lines)

        if (block_found /= 'FIND') then
            error_message = 'Error: "stage" field not found while i_stage = 1'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        !--------------------------------------------------------------------------
        ! read n_x, n_y first to allocate
        do iline = 1, nlines
            bufline = input_lines(iline)
            call ukoto_uiequa2(input_lines(iline), n_equation, csides)

            do i_equation = 1, n_equation
                cvalue = 'eps'
                call ukoto_rvalue2(out_id, csides(1, i_equation), stage%eps, cvalue)
                cvalue = 'height'
                call ukoto_rvalue2(out_id, csides(1, i_equation), stage%height, cvalue)
                cvalue = 'cutoff_ratio'
                call ukoto_rvalue2(out_id, csides(1, i_equation), stage%cutoff_ratio, cvalue)
                cvalue = 'distance_offset'
                call ukoto_rvalue2(out_id, csides(1, i_equation), stage%distance_offset, cvalue)

                cvalue = 'group_id'
                call ukoto_ivalue2(out_id, csides(1, i_equation), stage%group_id, cvalue)
            end do
        end do

        !--------------------------------------------------------------------------
        ! check input and allocate pixels
        if (stage%group_id == -1) then
            error_message = 'Error: "group_id" field not found in "stage" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef MPI_PAR
    end if
#endif

#ifdef _DEBUG
    write (*, *) "DEBUG: subroutine setp_stage end"
#endif

end subroutine setp_stage
