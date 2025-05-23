! simu_energy_velo
!> @brief Calculate the kinetic energy

subroutine simu_energy_velo(velo_mp, pnle_unit, pnlet)

    use const_maxsize
    use const_index
    use var_setp, only: ifix_mp
    use var_struct, only: nmp_real, cmass_mp, imp2unit
#ifdef MPI_PAR3
    use mpiconst
#endif

    implicit none

    ! ---------------------------------------------------------------------
    real(PREC), intent(in)    :: velo_mp(:, :)     ! (SPACE_DIM, nmp_real)
    real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX)
    real(PREC), intent(inout) :: pnle_unit(:, :, :) ! (nunit_all, nunit_all, E_TYPE%MAX)

    ! ---------------------------------------------------------------------
    ! local variables
    integer :: imp, iunit
    integer :: ksta, kend
    real(PREC) :: ev
#ifdef MPI_PAR3
    integer :: klen
#endif

    ! ---------------------------------------------------------------------
#ifdef MPI_PAR3
    klen = (nmp_real - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, nmp_real)
#else
    ksta = 1
    kend = nmp_real
#endif
!$omp do private(ev,iunit)
    do imp = ksta, kend

        if (ifix_mp(imp) == 1) cycle

        ev = 0.5e0_PREC*cmass_mp(imp) &
             *(velo_mp(1, imp)**2 + velo_mp(2, imp)**2 + velo_mp(3, imp)**2)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%VELO) = pnlet(E_TYPE%VELO) + ev

        iunit = imp2unit(imp)
        pnle_unit(iunit, iunit, E_TYPE%VELO) = pnle_unit(iunit, iunit, E_TYPE%VELO) + ev
    end do
!$omp end do nowait

end subroutine simu_energy_velo
