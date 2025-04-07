!simu_force_exv_wca
!> @brief Calculates the force related to excluded volume

subroutine simu_force_exv_wca(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: exv_wca_list
    use var_setp, only: indtrna
    use var_struct, only: nmp_all, xyz_mp_rep, pxyz_mp_rep
    use mpiconst

    implicit none

    ! --------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: klen, ksta, kend
    integer :: imp1, imp2
    integer :: ipnl, imirror
    real(PREC) :: dist2
    real(PREC) :: coef
    real(PREC) :: cdist2
    real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist14
    real(PREC) :: dvdw_dr
    real(PREC) :: v21(SPACE_DIM), for(SPACE_DIM)

    ! --------------------------------------------------------------------
    !! Currently this potential is available noly for RNA.
    cdist2 = indtrna%cdist_exvwca**2
    coef = 12.0e0_PREC*indtrna%coef_exvwca/cdist2

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
    klen = (exv_wca_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, exv_wca_list(irep)%num_pairs)
#else
    ksta = 1
    kend = exv_wca_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
    print *, "exv_wca      = ", kend - ksta + 1
#endif
#else
    ksta = 1
    kend = exv_wca_list(irep)%num_pairs
#endif
!$omp do private(imp1,imp2,v21,dist2,&
!$omp&           roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist14,dvdw_dr,for,imirror)
    do ipnl = ksta, kend

        imp1 = exv_wca_list(irep)%pairs(1, ipnl)
        imp2 = exv_wca_list(irep)%pairs(2, ipnl)

!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        if (inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
            imirror = exv_wca_list(irep)%pairs(3, ipnl)
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

        dist2 = dot_product(v21, v21)

        if (dist2 > cdist2) cycle

        ! -----------------------------------------------------------------
        roverdist2 = cdist2/dist2
        roverdist4 = roverdist2*roverdist2
        roverdist8 = roverdist4*roverdist4
        roverdist14 = roverdist2*roverdist4*roverdist8
        dvdw_dr = coef*(roverdist14 - roverdist8)
        if (dvdw_dr > DE_MAX) then
            !write (*, *) "exvol protein", imp1, imp2, dvdw_dr
            dvdw_dr = DE_MAX
        end if
!     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

        for(1:3) = dvdw_dr*v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
    end do
!$omp end do nowait

end subroutine simu_force_exv_wca
