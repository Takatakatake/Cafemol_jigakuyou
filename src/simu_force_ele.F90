! simu_force_ele
!> @brief This subroutine calculates the force of electrostatic interaction.

subroutine simu_force_ele(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: ele_list
    use var_setp, only: inmisc, inele, inion
    use var_struct, only: xyz_mp_rep, pxyz_mp_rep, &
        nmp_all, iontype_mp, imp2type
    use var_replica, only: irep2grep
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! --------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: imp1, imp2
    integer :: ksta, kend
    integer :: grep, iele1, itype1, itype2
    integer :: imptype1, imptype2
    integer :: imirror
    real(PREC) :: dist1, dist2, rdist1, xtanh, rsig, rek_corr
    real(PREC) :: dvdw_dr, rcdist, cutoff2
    real(PREC) :: ek, ek_simu
    real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
    integer :: klen
#endif

    ! --------------------------------------------------------------------
#ifdef _DEBUG
    write (*, *) '#### start simu_force_ele'
#endif

    ! for speed up
    grep = irep2grep(irep)
    cutoff2 = (inele%cutoff_ele*inele%cdist(grep))**2
    rcdist = 1.0e0_PREC/inele%cdist(grep)

    ek_simu = inele%diele
    ek = inele%diele_water

#ifdef _DEBUG
    write (*, *) 'lele(irep), ', ele_list(irep)%num_pairs
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH
    klen = (ele_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, ele_list(irep)%num_pairs)
#else
    ksta = 1
    kend = ele_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
    print *, "pnl2_6       = ", kend - ksta + 1
#endif

#else
    ksta = 1
    kend = ele_list(irep)%num_pairs
#endif

!$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,itype1,itype2, imptype1, imptype2,&
!$omp&           rsig,xtanh,rek_corr,dvdw_dr,for,imirror)
    do iele1 = ksta, kend
        imp1 = ele_list(irep)%pairs(1, iele1)
        imp2 = ele_list(irep)%pairs(2, iele1)

        if (inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
            imirror = ele_list(irep)%pairs(3, iele1)
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

        ! v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        if (dist2 > cutoff2) cycle

        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        rdist1 = 1.0/dist1

        itype1 = iontype_mp(imp1)
        itype2 = iontype_mp(imp2)

        if ((.not. inmisc%class_flag(CLASS%ION)) .or. &
            itype1 <= 0 .or. itype1 > IONTYPE%MAX_ALL .or. &
            itype2 <= 0 .or. itype2 > IONTYPE%MAX_ALL) then
            dvdw_dr = ele_list(irep)%coefs(1, iele1)*rdist1*rdist1* &
                      (1.0e0_PREC*rdist1 + rcdist)* &
                      exp(-dist1*rcdist)
        else
            rsig = 1.0/inion%csigmame(itype1, itype2)
            xtanh = tanh((dist1 - inion%cdistme(itype1, itype2))*rsig)
            rek_corr = 1.0/(0.5*(ek + 5.2) + 0.5*(ek - 5.2)*xtanh)

            dvdw_dr = ek_simu*rek_corr*ele_list(irep)%coefs(1, iele1)*rdist1*rdist1 &
                      *(1.0e0_PREC*rdist1 + rcdist &
                        + rek_corr*rsig*0.5*(ek - 5.2)*(1.0 - xtanh**2))* &
                      exp(-dist1*rcdist)
        end if

        if (dvdw_dr > DE_MAX) then
            dvdw_dr = DE_MAX
        end if
        ! if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

        imptype1 = imp2type(imp1)
        imptype2 = imp2type(imp2)
        if (imptype1 == MPTYPE%DNA2_PHOS .AND. imptype2 == MPTYPE%PRO) then
            dvdw_dr = dvdw_dr*(-inele%dna2_phos_pro_charge/0.6)
        else if (imptype2 == MPTYPE%DNA2_PHOS .AND. imptype1 == MPTYPE%PRO) then
            dvdw_dr = dvdw_dr*(-inele%dna2_phos_pro_charge/0.6)
        end if

        for(1:3) = dvdw_dr*v21(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
    end do
!$omp end do nowait

end subroutine simu_force_ele
