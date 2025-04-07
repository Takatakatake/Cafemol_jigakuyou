! setp_afm_fitting
!> @brief This subroutine is to read and set the parameters for afm fitting.
! vim: tabstop=3 shiftwidth=3

subroutine setp_afm_fitting()

    use const_maxsize
    use const_index
    use var_inp, only: infile, outfile
    use var_struct, only: grp
    use var_afmfit ! write into this
!$  use omp_lib

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    integer      :: inp_id, out_id ! file id
    character(4) :: block_found    ! ukoto status
    integer      :: iline, nlines  ! # of lines and counter
    character(CARRAY_MXCOLM) :: input_lines(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: bufline

    ! buffer to read AFM_IMAGE lines
    character(9) :: ident
    integer      :: i_pixel

    ! to parse `x = y`
    integer                  :: i_equation, n_equation
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MXCOLM) :: cvalue

    ! a buffer for error message
    character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
    write (*, *) "DEBUG: subroutine setp_afm_fitting start"
#endif

    inp_id = infile%inp
    out_id = outfile%data

    ! reset if already set
    if (allocated(afm%x_at))                deallocate (afm%x_at)     !&
    if (allocated(afm%y_at))                deallocate (afm%y_at)     !&
    if (allocated(afm%z_sim_at))            deallocate (afm%z_sim_at) !&
    if (allocated(afm%z_ref_at))            deallocate (afm%z_ref_at) !&
    if (allocated(afm%imp2pxl))             deallocate (afm%imp2pxl)  !&
    if (allocated(afm%sumexp))              deallocate (afm%sumexp)   !&
    if (allocated(afm%r_sumexp))            deallocate (afm%r_sumexp) !&
    if (allocated(afm%sumexp_thread_local)) deallocate (afm%sumexp_thread_local) !&
    if (allocated(afm%z_ref_sim_sum_thloc)) deallocate (afm%z_ref_sim_sum_thloc) !&
    if (allocated(afm%z_sim_sim_sum_thloc)) deallocate (afm%z_sim_sim_sum_thloc) !&

    ! set default values
    afm%n_x          = -1        !&
    afm%n_y          = -1        !&
    afm%n_pixel      =  0        !&
    afm%pixel_size_x =  0.0_PREC !&
    afm%pixel_size_y =  0.0_PREC !&
    afm%k            =  0.0_PREC !&
    afm%beta         =  0.0_PREC !&
    afm%sigma_x      =  0.0_PREC !&
    afm%sigma_y      =  0.0_PREC !&
    afm%n_max_pixel  = -1        !&
    afm%cutoff_sigma =  5.0_PREC !&
    afm%cutoff_exp   = 20.0_PREC !&
    afm%z0           =  0.0_PREC !&
    afm%num_max_threads = 1      !&

!$  afm%num_max_threads = omp_get_max_threads()

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        !--------------------------------------------------------------------------
        ! find <<<< afm_fitting block
        rewind (inp_id)
        call ukoto_uiread2(inp_id, out_id, 'afm_fitting     ', block_found, &
                           CARRAY_MXLINE, nlines, input_lines)

        if (block_found /= 'FIND') then
            error_message = 'Error: "afm_fitting" field not found while i_afm_fitting = 1'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        !--------------------------------------------------------------------------
        ! read n_x, n_y first to allocate
        do iline = 1, nlines
            bufline = input_lines(iline)
            if (bufline(1:9) /= 'AFM_IMAGE') then
                call ukoto_uiequa2(input_lines(iline), n_equation, csides)

                do i_equation = 1, n_equation
                    cvalue = 'k'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%k, cvalue)
                    cvalue = 'beta'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%beta, cvalue)
                    cvalue = 'sigma_x'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%sigma_x, cvalue)
                    cvalue = 'sigma_y'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%sigma_y, cvalue)
                    cvalue = 'z0'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%z0, cvalue)
                    cvalue = 'cutoff_sigma'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%cutoff_sigma, cvalue)
                    cvalue = 'cutoff_exp'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%cutoff_exp, cvalue)

                    cvalue = 'n_x'
                    call ukoto_ivalue2(out_id, csides(1, i_equation), afm%n_x, cvalue)
                    cvalue = 'n_y'
                    call ukoto_ivalue2(out_id, csides(1, i_equation), afm%n_y, cvalue)
                    cvalue = 'pixel_size_x'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%pixel_size_x, cvalue)
                    cvalue = 'pixel_size_y'
                    call ukoto_rvalue2(out_id, csides(1, i_equation), afm%pixel_size_y, cvalue)

                    cvalue = 'group_id'
                    call ukoto_ivalue2(out_id, csides(1, i_equation), afm%group_id, cvalue)
                end do
            end if
        end do

        !--------------------------------------------------------------------------
        ! check input and allocate pixels
        if (afm%n_x == -1) then
            error_message = 'Error: "n_x" field not found in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (afm%n_y == -1) then
            error_message = 'Error: "n_y" field not found in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        afm%n_pixel = afm%n_x*afm%n_y
        write (out_id, '(1a, 1i10, 1a, 1i10, 1a, 1i10)') '---afm image size: ', &
            afm%n_x, ' * ', afm%n_y, ' = ', afm%n_pixel

        allocate (afm%x_at(afm%n_pixel))
        allocate (afm%y_at(afm%n_pixel))
        allocate (afm%z_ref_at(afm%n_pixel))
        allocate (afm%z_sim_at(afm%n_pixel))
        allocate (afm%sumexp(afm%n_pixel))
        allocate (afm%r_sumexp(afm%n_pixel))
        allocate (afm%sumexp_thread_local(afm%n_pixel, afm%num_max_threads))
        allocate (afm%z_sim_sim_sum_thloc(afm%num_max_threads))
        allocate (afm%z_ref_sim_sum_thloc(afm%num_max_threads))

        !--------------------------------------------------------------------------
        ! check input and setup parameters
        if (afm%beta == 0.0_PREC) then
            error_message = 'Error: "beta" not found (or = 0.0) in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (afm%sigma_x == 0.0_PREC) then
            error_message = 'Error: "sigma_x" not found (or = 0.0) in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (afm%sigma_y == 0.0_PREC) then
            error_message = 'Error: "sigma_y" not found (or = 0.0) in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (afm%pixel_size_x == 0.0_PREC) then
            error_message = 'Error: "pixel_size_x" not found (or = 0.0) in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (afm%pixel_size_y == 0.0_PREC) then
            error_message = 'Error: "pixel_size_y" not found (or = 0.0) in "afm_fitting" block'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        afm%rbeta = 1.0_PREC/afm%beta
        afm%rsigma_x = 1.0_PREC/afm%sigma_x
        afm%rsigma_y = 1.0_PREC/afm%sigma_y
        afm%rsigma_x_sq = 1.0_PREC/(afm%sigma_x*afm%sigma_x)
        afm%rsigma_y_sq = 1.0_PREC/(afm%sigma_y*afm%sigma_y)
        afm%n_max_pixel = & ! for cutoff
            (1 + 2*int(ceiling(afm%sigma_x*afm%cutoff_sigma/afm%pixel_size_x)))* &
            (1 + 2*int(ceiling(afm%sigma_y*afm%cutoff_sigma/afm%pixel_size_y)))
        allocate (afm%imp2pxl(grp%nmp(afm%group_id)*afm%n_max_pixel))
        afm%imp2pxl = -1

        !--------------------------------------------------------------------------
        ! write debug informations
#ifdef _DEBUG
        write (*, *) "DEBUG: omp_get_max_threads            = ", afm%num_max_threads
        write (*, *) "DEBUG: afm%sigma_x * afm%cutoff_sigma = ", afm%sigma_x*afm%cutoff_sigma
        write (*, *) "DEBUG: afm%sigma_y * afm%cutoff_sigma = ", afm%sigma_y*afm%cutoff_sigma
        write (*, *) "DEBUG: pixel_size_x / r_c             = ", &
            int(ceiling(afm%sigma_x*afm%cutoff_sigma/afm%pixel_size_x))
        write (*, *) "DEBUG: pixel_size_y / r_c             = ", &
            int(ceiling(afm%sigma_y*afm%cutoff_sigma/afm%pixel_size_y))
        write (*, *) "DEBUG: n_max_pixel                    = ", afm%n_max_pixel
        write (*, *) "DEBUG: grp%nmp(afm%group_id)          = ", grp%nmp(afm%group_id)
#endif

        !--------------------------------------------------------------------------
        ! read AFM_IMAGE section. the arrays are already allocated!
        afm%z_ref_sq_sum = 0.0_PREC
        afm%r_z_ref_sq_sum = 0.0_PREC
        i_pixel = 0
        do iline = 1, nlines
            bufline = input_lines(iline)
            if (bufline(1:9) == 'AFM_IMAGE') then
                i_pixel = i_pixel + 1
                if (i_pixel > afm%n_pixel) then
                    error_message = 'Error: too many AFM_IMAGE lines in "afm_fitting" block'
                    call util_error(ERROR%STOP_ALL, error_message)
                end if

                read (bufline, *) ident, afm%x_at(i_pixel), afm%y_at(i_pixel), &
                    afm%z_ref_at(i_pixel)
                if (ident /= 'AFM_IMAGE') then
                    error_message = 'Error: internal error while reading "afm_fitting" block'
                    call util_error(ERROR%STOP_ALL, error_message)
                end if

                afm%z_ref_sq_sum = afm%z_ref_sq_sum + afm%z_ref_at(i_pixel)**2

                write (out_id, '(1a, 3g12.4)') '---reading afm pixels: ', &
                    afm%x_at(i_pixel), afm%y_at(i_pixel), afm%z_ref_at(i_pixel)
            end if
        end do
        write (out_id, '(1a, 1i10, 1a)') '---there were ', i_pixel, ' AFM_IMAGE lines'

        if (afm%z_ref_sq_sum == 0.0_PREC) then
            error_message = 'Error: input AFM image is 0.0 everywhere. the correlation become undefined.'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
        afm%r_z_ref_sq_sum = 1.0_PREC/afm%z_ref_sq_sum

        if (i_pixel /= afm%n_pixel) then
            error_message = 'Error: number of "AFM_IMAGE" lines differs from n_x * n_y'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef MPI_PAR
    end if
#endif

#ifdef _DEBUG
    write (*, *) "DEBUG: subroutine setp_afm_fitting end"
#endif

end subroutine setp_afm_fitting

subroutine deallocate_afmfit()
    use var_afmfit
    implicit none

    integer :: i_pixel

    if (allocated(afm%x_at))                deallocate (afm%x_at)     !&
    if (allocated(afm%y_at))                deallocate (afm%y_at)     !&
    if (allocated(afm%z_sim_at))            deallocate (afm%z_sim_at) !&
    if (allocated(afm%z_ref_at))            deallocate (afm%z_ref_at) !&
    if (allocated(afm%imp2pxl))             deallocate (afm%imp2pxl)  !&
    if (allocated(afm%sumexp))              deallocate (afm%sumexp)   !&
    if (allocated(afm%r_sumexp))            deallocate (afm%r_sumexp) !&
    if (allocated(afm%sumexp_thread_local)) deallocate (afm%sumexp_thread_local) !&
    if (allocated(afm%z_ref_sim_sum_thloc)) deallocate (afm%z_ref_sim_sum_thloc) !&
    if (allocated(afm%z_sim_sim_sum_thloc)) deallocate (afm%z_sim_sim_sum_thloc) !&

end subroutine deallocate_afmfit
