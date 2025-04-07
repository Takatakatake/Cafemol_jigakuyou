! simu_force_afm_fitting
!> @brief Calculates and adds the force related to afmfit interaction
! vim: tabstop=3 shiftwidth=3

subroutine simu_force_afm_fitting(irep, force_mp)

    use const_maxsize
    use const_index
    use var_struct, only: xyz_mp_rep, nmp_all, grp, exv_radius_mp
    use var_afmfit
!$  use omp_lib
    implicit none

    ! -------------------------------------------------------------------------
    integer,    intent(in)    :: irep                 !&
    real(PREC), intent(inout) :: force_mp(3, nmp_all) !&

    ! -------------------------------------------------------------------------
    integer    :: idx, imp, i_pxl, i_list, thread_id
    real(PREC) ::   z_ref_sq_sum,   z_sim_sq_sum,   z_ref_sim_sum !&
    real(PREC) :: r_z_ref_sq_sum, r_z_sim_sq_sum, r_z_ref_sim_sum !&
    real(PREC) :: stage_term, cc, coef, term
    real(PREC) :: dx, dy, dz, tmp
    real(PREC) :: pos(3), force(3)

    character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
    write (*, *) "DEBUG: simu_force_afm_fitting start"
    write (*, *) "DEBUG: n_pixel = ", afm%n_pixel
#endif

    ! -------------------------------------------------------------------------
    ! clear/initialize values

    z_ref_sq_sum    = afm%z_ref_sq_sum                   !&
    z_sim_sq_sum    = 0.0_PREC                           !&
    z_ref_sim_sum   = 0.0_PREC                           !&
    r_z_ref_sq_sum  = afm%r_z_ref_sq_sum                 !&
    r_z_sim_sq_sum  = 0.0_PREC                           !&
    r_z_ref_sim_sum = 0.0_PREC                           !&
    stage_term      = exp((0.0_PREC - afm%z0)*afm%rbeta) !&
    thread_id       = 1                                  !&

    afm%imp2pxl             = -1       !&
    afm%sumexp              = 0.0_PREC !&
    afm%r_sumexp            = 0.0_PREC !&
    afm%sumexp_thread_local = 0.0_PREC !&
    afm%z_sim_sim_sum_thloc = 0.0_PREC !&
    afm%z_ref_sim_sum_thloc = 0.0_PREC !&

    ! -------------------------------------------------------------------------
    ! search pixels that affects each particle and calculate contributions

!$omp parallel
!$omp do private(thread_id, idx, imp, i_list, i_pxl, pos, dx, dy, dz, tmp)
    do idx = 1, grp%nmp(afm%group_id)
!$      thread_id = omp_get_thread_num() + 1
        imp = grp%implist(idx, afm%group_id)
        pos = xyz_mp_rep(1:3, imp, irep)

        i_list = 0
        do i_pxl = 1, afm%n_pixel
            dx = pos(1) - afm%x_at(i_pxl)
            dy = pos(2) - afm%y_at(i_pxl)
            if (abs(dx) > afm%cutoff_sigma * afm%sigma_x) cycle !&
            if (abs(dy) > afm%cutoff_sigma * afm%sigma_y) cycle !&

            dz = pos(3) - afm%z0
            if (dz * afm%rbeta < -afm%cutoff_exp) cycle !&

            tmp = exp(-dx * dx * afm%rsigma_x_sq * 0.5_PREC - & !&
                       dy * dy * afm%rsigma_y_sq * 0.5_PREC + & !&
                      (dz + exv_radius_mp(imp)) * afm%rbeta)    !&

            afm%sumexp_thread_local(i_pxl, thread_id) = &
                afm%sumexp_thread_local(i_pxl, thread_id) + tmp

            i_list = i_list + 1
            afm%imp2pxl((idx - 1) * afm%n_max_pixel + i_list) = i_pxl
        end do

        if (i_list > afm%n_max_pixel) then
            error_message = 'Error: internal error: i_list exceeds n_max_pixel'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef _DEBUG
        write (*, *) "DEBUG: i_list for ", imp, "-th particle = ", i_list
#endif
    end do
!$omp end do
! never add `nowait` to any of the parallel do in this subroutine.

    ! -------------------------------------------------------------------------
    ! merge thread local pixel contributions to sumexp

!$omp do private(i_pxl)
    do i_pxl = 1, afm%n_pixel
        afm%sumexp(i_pxl) = stage_term + sum(afm%sumexp_thread_local(i_pxl, :))
    end do
!$omp end do

    ! -------------------------------------------------------------------------
    ! calculate height in simulation and correlation coefficient

    ! CafeMol requires OpenMP implementation should be deterministic to run
    ! its standard test script. But the order of omp reduction is not specified
    ! (OpenMP 4.5 2.15.3.6, p.205 l.23). This means that it is not guaranteed
    ! that the bitwise identical results will be obtained when comparing one
    ! parallel run to another because it uses floating point operation.
    ! Therefore, here it does not use `omp reduction(z_sim_sq_sum:+)`.

!$omp do private(i_pxl, thread_id), schedule(static)
    do i_pxl = 1, afm%n_pixel
!$      thread_id = omp_get_thread_num() + 1
        afm%z_sim_at(i_pxl) = afm%beta * log(stage_term + afm%sumexp(i_pxl)) !&
        afm%r_sumexp(i_pxl) = 1.0_PREC / afm%sumexp(i_pxl)                   !&

        afm%z_sim_sim_sum_thloc(thread_id) =          &
           afm%z_sim_sim_sum_thloc(thread_id) +       &
           afm%z_sim_at(i_pxl) * afm%z_sim_at(i_pxl) !&
        afm%z_ref_sim_sum_thloc(thread_id) =          &
           afm%z_ref_sim_sum_thloc(thread_id) +       &
           afm%z_ref_at(i_pxl) * afm%z_sim_at(i_pxl) !&
    end do
!$omp end do

!$omp single
    z_sim_sq_sum  = sum(afm%z_sim_sim_sum_thloc(:)) !&
    z_ref_sim_sum = sum(afm%z_ref_sim_sum_thloc(:)) !&

    r_z_sim_sq_sum  = 1.0_PREC / z_sim_sq_sum            !&
    r_z_ref_sim_sum = 1.0_PREC / z_ref_sim_sum           !&
    cc = z_ref_sim_sum / sqrt(z_sim_sq_sum*z_ref_sq_sum) !&

    coef = afm%k * cc !&
!$omp end single

    ! -------------------------------------------------------------------------
    ! calculate force

!$omp do private(idx, imp, i_list, i_pxl, pos, force, dx, dy, dz, term)
    do idx = 1, grp%nmp(afm%group_id)
        imp = grp%implist(idx, afm%group_id)
        pos = xyz_mp_rep(1:3, imp, irep)

        force(1:3) = 0.0_PREC
        do i_list = 1, afm%n_max_pixel
            i_pxl = afm%imp2pxl((idx - 1) * afm%n_max_pixel + i_list) !&
            if (i_pxl == -1) exit

            dx = pos(1) - afm%x_at(i_pxl)
            dy = pos(2) - afm%y_at(i_pxl)
            dz = pos(3) - afm%z0

            term = coef * (afm%z_ref_at(i_pxl) * r_z_ref_sim_sum - & !&
                           afm%z_sim_at(i_pxl) * r_z_sim_sq_sum) * & !&
                   afm%beta * afm%r_sumexp(i_pxl) *                & !&
                   exp(-0.5_PREC * dx * dx * afm%rsigma_x_sq -     & !&
                        0.5_PREC * dy * dy * afm%rsigma_y_sq +     & !&
                       dz * afm%rbeta)

            force(1) = force(1) - term * dx * afm%rsigma_x_sq !&
            force(2) = force(2) - term * dy * afm%rsigma_x_sq !&
            force(3) = force(3) + term * afm%rbeta            !&
        end do
        force_mp(1:3, imp) = force_mp(1:3, imp) + force(1:3)
    end do
!$omp end do
!$omp end parallel

#ifdef _DEBUG
    write (*, *) "DEBUG: simu_force_afm_fitting end"
#endif

end subroutine simu_force_afm_fitting
