! simu_energy_stage
!> @brief Lennard-Jones form interaction between a group of particle and a plane
!>        that is perpendicular to the Z axis. all the radiis are the same value
!>        as EXV. epsilon values are uniform. The main objective is afm_fitting.
! vim: tabstop=3 shiftwidth=3

subroutine simu_energy_stage(irep, pnlet)

    use const_index
    use var_struct, only: xyz_mp_rep, grp, exv_radius_mp
    use var_stage
    implicit none

    ! ------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: pnlet(:)

    ! ------------------------------------------------------------
    integer    :: idx, imp
    real(PREC) :: dz, sgm, sz, sz3, sz6

#ifdef _DEBUG
    write (*, *) "DEBUG: simu_energy_stage start"
#endif

    ! ------------------------------------------------------------
!$omp do private(imp, dz, sgm, sz, sz3, sz6)
    do idx = 1, grp%nmp(stage%group_id)
        imp = grp%implist(idx, stage%group_id)
        dz = xyz_mp_rep(3, imp, irep) - stage%height
        sgm = exv_radius_mp(imp) + stage%distance_offset
        if (sgm*stage%cutoff_ratio < dz) cycle
        sz = sgm/dz
        sz3 = sz*sz*sz
        sz6 = sz3*sz3

        pnlet(E_TYPE%Z_STAGE) = pnlet(E_TYPE%Z_STAGE) + &
                                4.0_PREC*stage%eps*sz6*(sz6 - 1.0_PREC)

    end do
!$omp end do nowait

#ifdef _DEBUG
    write (*, *) "DEBUG: simu_energy_stage end"
#endif

    return
end subroutine simu_energy_stage
