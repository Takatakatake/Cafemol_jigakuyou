! simu_energy_afm_fitting
!> @brief Calculate energy of afm fitting interaction
! vim: tabstop=3 shiftwidth=3

subroutine simu_energy_afm_fitting(irep, pnlet)

    use const_index
    use const_maxsize
    use var_struct, only: xyz_mp_rep, grp, exv_radius_mp
    use var_afmfit
!$  use omp_lib
    implicit none

    ! -------------------------------------------------------------------------
    integer,    intent(in)    :: irep     !&
    real(PREC), intent(inout) :: pnlet(:) !&

    ! -------------------------------------------------------------------------
    integer    :: idx, imp, i_pxl, thread_id
    real(PREC) :: z_ref_sq_sum, z_sim_sq_sum, z_ref_sim_sum
    real(PREC) :: stage_term, cc, dx, dy, dz, tmp
    real(PREC) :: pos(3)

#ifdef _DEBUG
    write (*, *) "DEBUG: simu_energy_afm_fitting start"
#endif

    ! -------------------------------------------------------------------------
    ! clear values
    z_ref_sq_sum  = afm%z_ref_sq_sum                     !&
    z_sim_sq_sum  = 0.0_PREC                             !&
    z_ref_sim_sum = 0.0_PREC                             !&
    stage_term    = exp((0.0_PREC - afm%z0) * afm%rbeta) !&
    thread_id     = 1                                    !&

    afm%sumexp              = 0.0_PREC !&
    afm%r_sumexp            = 0.0_PREC !&
    afm%sumexp_thread_local = 0.0_PREC !&
    afm%z_sim_sim_sum_thloc = 0.0_PREC !&
    afm%z_ref_sim_sum_thloc = 0.0_PREC !&

    ! -------------------------------------------------------------------------
    ! search pixels that affects to each particle and calculate contributions

!$omp parallel
!$omp do private(thread_id, idx, imp, pos, i_pxl, dx, dy, dz, tmp)
    do idx = 1, grp%nmp(afm%group_id)
!$      thread_id = omp_get_thread_num() + 1
        imp = grp%implist(idx, afm%group_id)
        pos = xyz_mp_rep(1:3, imp, irep)

        do i_pxl = 1, afm%n_pixel
            dx = pos(1) - afm%x_at(i_pxl)
            dy = pos(2) - afm%y_at(i_pxl)
            if (abs(dx) > afm%cutoff_sigma * afm%sigma_x .or. & !&
                abs(dy) > afm%cutoff_sigma * afm%sigma_y) cycle !&
            dz = pos(3) - afm%z0
            if (dz*afm%rbeta < -afm%cutoff_exp) cycle

            tmp = exp(-dx * dx * afm%rsigma_x_sq * 0.5_PREC - & !&
                       dy * dy * afm%rsigma_y_sq * 0.5_PREC + & !&
                      (dz + exv_radius_mp(imp)) * afm%rbeta)    !&

            afm%sumexp_thread_local(i_pxl, thread_id) = &
                afm%sumexp_thread_local(i_pxl, thread_id) + tmp
        end do
    end do
!$omp end do

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

        afm%z_sim_sim_sum_thloc(thread_id) =    &
           afm%z_sim_sim_sum_thloc(thread_id) + &
           afm%z_sim_at(i_pxl) * afm%z_sim_at(i_pxl)  !&
        afm%z_ref_sim_sum_thloc(thread_id) =    &
           afm%z_ref_sim_sum_thloc(thread_id) + &
           afm%z_ref_at(i_pxl) * afm%z_sim_at(i_pxl)  !&
    end do
!$omp end do
!$omp end parallel

    z_sim_sq_sum  = sum(afm%z_sim_sim_sum_thloc(:)) !&
    z_ref_sim_sum = sum(afm%z_ref_sim_sum_thloc(:)) !&

#ifdef _DEBUG
    write (*, *) "DEBUG: z_sim_sq_sum  = ", z_sim_sq_sum
    write (*, *) "DEBUG: z_ref_sq_sum  = ", z_ref_sq_sum
    write (*, *) "DEBUG: z_ref_sim_sum = ", z_ref_sim_sum
#endif

    ! -------------------------------------------------------------------------
    ! submit energy into pnlet

    cc = z_ref_sim_sum / sqrt(z_sim_sq_sum * z_ref_sq_sum) !&

#ifdef _DEBUG
    write (*, *) "DEBUG: correlation coefficient = ", cc                      !&
    write (*, *) "DEBUG: energy for AFM fitting  = ", afm%k * (1.0_PREC - cc) !&
#endif

    pnlet(E_TYPE%AFM_CORRELATION) = cc                      !&
    pnlet(E_TYPE%AFM_FITTING)     = afm%k * (1.0_PREC - cc) !&

#ifdef _DEBUG
    write (*, *) "DEBUG: simu_energy_afm_fitting end"
#endif

end subroutine simu_energy_afm_fitting
