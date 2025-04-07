! simu_muca_wl
!> @brief Wang-Landau Multicanonical

subroutine simu_muca_wl(irep)
    use const_maxsize, only: PREC
    use const_physical
    use const_index
    use mpiconst
    use if_mloop ! without this segmentation fault will be occured when calling simu_energy
    use var_inp, only: outfile
    use var_setp, only: inmmc, insimu
    use var_simu, only: muca_cv, pnlet_muca, pnle_unit_muca, velo_mp, force_mp, &
        istep, muca_bias, muca_bias_last, muca_bias_buf, muca_increment, muca_coef

    implicit none
    integer, intent(in) :: irep

    integer :: i_bins
    real(PREC), save :: cv_max, cv_min, mu, sigma, normalizer, threshold, w_bin
    real(PREC) :: muca_bias_sum(inmmc%muca_bins)

    ! initialize

    if (istep == 1) then
        ! for all MPI process
        cv_max = inmmc%muca_cv_max
        cv_min = inmmc%muca_cv_min
        sigma = 3.0*(cv_max - cv_min)/inmmc%muca_bins
        normalizer = 1.0/(sqrt(2*F_PI)*sigma)
        muca_increment = inmmc%muca_increment
        threshold = 0.99 ! timing of increment decay
        w_bin = (cv_max - cv_min)/inmmc%muca_bins

        ! load initial histogram if exist (TODO)
        ! write(*,*) "i_muca_init", inmmc%i_muca_init

        if (inmmc%i_muca_init == 1) then
            write (*, *) "open muca file (32):", inmmc%muca_init
            open (32, file=inmmc%muca_init, status='old')
            read (32, *) muca_bias(1:inmmc%muca_bins)
            close (32)

        else
            muca_bias = 0
        endif
        muca_bias_last = muca_bias

        ! once only for one task
        if ((local_rank_mpi == 0) .and. (local_rank_rep == 0)) then
            ! write histogram metadata
            write (outfile%muca, *) "#cv_max", cv_max
            write (outfile%muca, *) "#cv_min", cv_min
            write (outfile%muca, *) "#bins", inmmc%muca_bins

            ! write initial histogram
            do i_bins = 1, inmmc%muca_bins
                write (outfile%muca, '(20f15.6)', advance='no') muca_bias(i_bins)
            enddo
        endif
    endif

    !--------------------------------------------------
    ! calculate CV
    !--------------------------------------------------
    call simu_energy(irep, velo_mp(:, :, irep), pnlet_muca, pnle_unit_muca)
    muca_cv = pnlet_muca(E_TYPE%TOTAL)

    !--------------------------------------------------
    ! increment histogram
    !--------------------------------------------------
    if (istep < inmmc%muca_stop) then
        i_bins = cv2histo(muca_cv)
#ifdef MPI_PAR
        if ((i_bins >= 1) .and. (i_bins <= inmmc%muca_bins)) then
            muca_bias_buf(i_bins) = muca_bias_buf(i_bins) + muca_increment
        endif
        call MPI_ALLREDUCE(muca_bias_buf, muca_bias_sum, inmmc%muca_bins, PREC_MPI, MPI_SUM, MPI_COMM_REP, ierr)
        muca_bias(1:inmmc%muca_bins) = muca_bias(1:inmmc%muca_bins) + muca_bias_sum
        muca_bias_buf = 0
        muca_bias_sum = 0
#else
        if ((i_bins >= 1) .and. (i_bins <= inmmc%muca_bins)) then
            muca_bias(i_bins) = muca_bias(i_bins) + muca_increment
        endif
#endif
    endif

    !--------------------------------------------------
    ! calculate coef to force
    ! (gaussian deposit.)
    !--------------------------------------------------
    muca_coef = 1

    if (muca_cv > cv_min) then
        do i_bins = 1, inmmc%muca_bins
            mu = i2mu(i_bins)
            muca_coef = muca_coef - muca_bias(i_bins)*normalizer*(muca_cv - mu)/sigma**2 &
                        *exp(-(muca_cv - mu)**2/(2*sigma**2))*w_bin
        enddo

        ! to flatten outer histogram
        do i_bins = inmmc%muca_bins, 2*inmmc%muca_bins
            mu = i2mu(i_bins)
            muca_coef = muca_coef - muca_bias(inmmc%muca_bins)*normalizer*(muca_cv - mu)/sigma**2 &
                        *exp(-(muca_cv - mu)**2/(2*sigma**2))*w_bin
        enddo
    endif

    !--------------------------------------------------
    ! multiple force
    !--------------------------------------------------
    force_mp = force_mp*muca_coef

    !--------------------------------------------------
    ! write histogram and other info to muca file
    !--------------------------------------------------
    if ((local_rank_mpi + local_rank_rep == 0) .and. (mod(istep, insimu%n_step_save) == 0)) then
        write (outfile%muca, *) "#istep ", istep
        write (outfile%muca, *) "#E_total", muca_cv
        write (outfile%muca, *) "#muca_coef", muca_coef
        if (istep < inmmc%muca_stop) then
            do i_bins = 1, inmmc%muca_bins
                write (outfile%muca, '(20f15.6)', advance='no') muca_bias(i_bins)
            enddo
        endif
    endif

contains
    integer function cv2histo(cv)
        real(PREC) cv
        cv2histo = 0
        if ((cv < cv_max) .and. (cv > cv_min)) then
            cv2histo = int((cv - cv_min)/(cv_max - cv_min)*inmmc%muca_bins) + 1
        endif
        if (cv >= cv_max) then
            cv2histo = inmmc%muca_bins + 1
        endif
    end function
    real(PREC) function i2mu(i)
        integer i
        i2mu = (cv_max - cv_min)/inmmc%muca_bins*(i - 0.5) + cv_min
    end function
end subroutine
