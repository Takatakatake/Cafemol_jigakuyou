! simu_energy_hp
!> @brief Calculate hydorophobic energy

! **************************************************************************
! Calculation of hydrophobic energy

subroutine simu_energy_hp(irep, pnlet, pnle_unit)

    use const_maxsize
    use const_physical
    use const_index
    use var_neighbor_list, only: hp_list
    use var_setp, only: inhp
    use var_struct, only: imp2unit, xyz_mp_rep, &
        nhp, ihp2mp, ncoor_hp, ncoor_max_hp, coef_aa_hp
#ifdef MPI_PAR3
    use mpiconst
#endif

    implicit none

    ! ------------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(out)   :: pnlet(:)         ! (E_TYPE%MAX)
    real(PREC), intent(out)   :: pnle_unit(:, :, :) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

    ! ------------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: iunit
    integer :: ineigh, ihp, jhp
    real(PREC) :: dist
    real(PREC) :: v21(SPACE_DIM)
    real(PREC) :: rho(nhp)
    real(PREC) :: pnl
#ifdef MPI_PAR3
    integer :: klen
#endif

    ! ------------------------------------------------------------------------
    ! hydrophobic

    rho(:) = 0.0e0_PREC

#ifdef MPI_PAR3
    klen = (hp_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, hp_list(irep)%num_pairs)
#else
    ksta = 1
    kend = hp_list(irep)%num_pairs
#endif

!$omp do private(ihp,jhp,v21,dist)
    do ineigh = ksta, kend
        ihp = hp_list(irep)%pairs(1, ineigh)
        jhp = hp_list(irep)%pairs(2, ineigh)

        v21(1:3) = xyz_mp_rep(1:3, ihp2mp(jhp), irep) - xyz_mp_rep(1:3, ihp2mp(ihp), irep)
        dist = sqrt(dot_product(v21, v21))

        if (dist <= hp_list(irep)%coefs(1, ineigh)) then
            rho(ihp) = rho(ihp) + ncoor_hp(jhp)
        else if (dist < hp_list(irep)%coefs(2, ineigh)) then
            rho(ihp) = rho(ihp) + ncoor_hp(jhp)* &
                       calc_uhp(dist, hp_list(irep)%coefs(1, ineigh), hp_list(irep)%coefs(2, ineigh))
        end if
    end do
!$omp end do

    rho(1:nhp) = rho(1:nhp)/ncoor_max_hp(1:nhp)

#ifdef MPI_PAR3
    klen = (nhp - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, nhp)
#else
    ksta = 1
    kend = nhp
#endif

!$omp do private(pnl, iunit)
    do ihp = ksta, kend
        pnl = -inhp%coef_hp*coef_aa_hp(ihp)*calc_shp(inhp%coef_rho_hp, rho(ihp), inhp%rho_min_hp)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%HPENE) = pnlet(E_TYPE%HPENE) + pnl
        iunit = imp2unit(ihp2mp(ihp))
        pnle_unit(iunit, iunit, E_TYPE%HPENE) = pnle_unit(iunit, iunit, E_TYPE%HPENE) + pnl
    end do
!$omp end do nowait

contains

    pure function calc_uhp(r, r_min, r_max) result(uhp)
        implicit none
        real(PREC), intent(in) :: r, r_min, r_max
        real(PREC) :: uhp

        uhp = 0.5e0_PREC*(1.0e0_PREC + cos(F_PI*(r - r_min)/(r_max - r_min)))
    end function

    pure function calc_shp(coef, rho, rho_min) result(shp)
        implicit none
        real(PREC), intent(in) :: coef, rho, rho_min
        real(PREC) :: shp

        if (rho >= 1.0e0_PREC) then
            shp = 1.0e0_PREC
        else if (rho <= rho_min) then
            shp = coef*rho
        else
            shp = coef*rho + 0.5e0_PREC*(1.0e0_PREC - coef)* &
                  (1.0e0_PREC + cos(F_PI*(1.0e0_PREC - rho)/(1.0e0_PREC - rho_min)))
        end if
    end function

end subroutine simu_energy_hp
