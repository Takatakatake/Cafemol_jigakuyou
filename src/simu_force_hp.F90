!>simu_force_hp
!> @brief Calculates and adds the hydrophobic force.

subroutine simu_force_hp(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_neighbor_list, only: hp_list
    use var_setp, only: inhp
    use var_struct, only: xyz_mp_rep, &
        nhp, ihp2mp, ncoor_hp, ncoor_max_hp, coef_aa_hp, nmp_all
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! ------------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! ------------------------------------------------------------------------
    ! local variables
    integer :: imp, jmp
    integer :: ineigh, ihp, jhp
    integer :: ksta, kend
    real(PREC) :: dist
    real(PREC) :: v21(SPACE_DIM), force(SPACE_DIM)
    real(PREC) :: cduhp_dr
    real(PREC) :: drho_dr(1:SPACE_DIM, 1:hp_list(irep)%num_pairs)
    real(PREC) :: uhp, rho(nhp), coefs(nhp)
#ifdef MPI_PAR
    integer :: klen
#endif

#ifdef _DEBUG
    write (*, *) '#### start simu_force_hp'
#endif

    ! ------------------------------------------------------------------------
    ! hydrophobic

#ifdef MPI_PAR
    klen = (nhp - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, nhp)
#else
    ksta = 1
    kend = nhp
#endif

#ifdef MPI_PAR
    klen = (hp_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, hp_list(irep)%num_pairs)
#else
    ksta = 1
    kend = hp_list(irep)%num_pairs
#endif

    rho(:) = 0.0
    drho_dr(:, :) = 0.0
!$omp do private(ihp,jhp,v21,dist,uhp,cduhp_dr)
    do ineigh = ksta, kend
        ihp = hp_list(irep)%pairs(1, ineigh)
        jhp = hp_list(irep)%pairs(2, ineigh)

        v21(1:3) = xyz_mp_rep(1:3, ihp2mp(jhp), irep) - xyz_mp_rep(1:3, ihp2mp(ihp), irep)
        dist = sqrt(dot_product(v21, v21))

        if (dist <= hp_list(irep)%coefs(1, ineigh)) then
            rho(ihp) = rho(ihp) + ncoor_hp(jhp)
        else if (dist < hp_list(irep)%coefs(2, ineigh)) then
            call calc_uhp_and_diff(dist, hp_list(irep)%coefs(1, ineigh), &
                                   hp_list(irep)%coefs(2, ineigh), uhp, cduhp_dr)
            rho(ihp) = rho(ihp) + ncoor_hp(jhp)*uhp
            drho_dr(1:3, ineigh) = ncoor_hp(jhp)*cduhp_dr*v21(1:3)
        end if
    end do
!$omp end do

    coefs(:) = 0.0
    rho(1:nhp) = rho(1:nhp)/ncoor_max_hp(1:nhp)
    coefs(1:nhp) = inhp%coef_hp*coef_aa_hp(1:nhp)* &
                   calc_dshp_drho(inhp%coef_rho_hp, rho(1:nhp), inhp%rho_min_hp)/ncoor_max_hp(1:nhp)

!$omp do private(ihp,jhp,imp,jmp,force)
    do ineigh = ksta, kend
        ihp = hp_list(irep)%pairs(1, ineigh)
        jhp = hp_list(irep)%pairs(2, ineigh)
        imp = ihp2mp(ihp)
        jmp = ihp2mp(jhp)
        force(1:3) = coefs(ihp)*drho_dr(1:3, ineigh)
        force_mp(1:3, jmp) = force_mp(1:3, jmp) + force(1:3)
        force_mp(1:3, imp) = force_mp(1:3, imp) - force(1:3)
    end do
!$omp end do nowait

contains

    pure subroutine calc_uhp_and_diff(r, r_min, r_max, uhp, diff)
        implicit none
        real(PREC), intent(in) :: r, r_min, r_max
        real(PREC), intent(out) :: uhp, diff
        real(PREC) :: r_interval, r_offset, theta

        r_interval = r_max - r_min
        r_offset = r - r_min
        theta = F_PI*r_offset/r_interval

        uhp = 0.5*(1.0 + cos(theta))
        diff = -0.5*sin(theta)*F_PI/r_interval/r
    end subroutine

    elemental function calc_dshp_drho(coef, rho, rho_min) result(dshp_drho)
        implicit none
        real(PREC), intent(in) :: coef, rho, rho_min
        real(PREC) :: dshp_drho
        real(PREC) :: rho_interval, rho_offset

        if (rho >= 1.0) then
            dshp_drho = 0.0
        else if (rho <= rho_min) then
            dshp_drho = coef
        else
            rho_interval = 1.0 - rho_min
            rho_offset = 1.0 - rho

            dshp_drho = coef + 0.5*(1.0 - coef)* &
                        sin(F_PI*rho_offset/rho_interval)*F_PI/rho_interval
        end if
    end function

end subroutine simu_force_hp
