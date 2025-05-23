! simu_energy_ele
!> @brief Calculate the energy of electrostatic interaction

subroutine simu_energy_ele(irep, pnlet, pnle_unit)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: ele_list
    use var_setp, only: inmisc, inele, inion
    use var_struct, only: imp2unit, xyz_mp_rep, pxyz_mp_rep, &
        iontype_mp, imp2type
    use var_replica, only: irep2grep
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
    integer :: imp1, imp2, iunit, junit, grep, iele1, imirror
    integer :: itype1, itype2
    integer :: imptype1, imptype2
    real(PREC) :: dist1, dist2
    real(PREC) :: pnl, rcdist, cutoff2, xtanh, rsig, rek_corr
    real(PREC) :: ew, ek, rk
    real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
    integer :: klen
#endif

    grep = irep2grep(irep)
    cutoff2 = (inele%cutoff_ele*inele%cdist(grep))**2
    rcdist = 1.0/inele%cdist(grep)

    ew = inele%diele_water
    ek = inele%diele

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH
    klen = (ele_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, ele_list(irep)%num_pairs)
#else
    ksta = 1
    kend = ele_list(irep)%num_pairs
#endif
#else
    ksta = 1
    kend = ele_list(irep)%num_pairs
#endif
!$omp do private(imp1,imp2,v21,dist2,dist1,itype1,itype2,imptype1,imptype2, &
!$omp&           rsig,xtanh,rek_corr,pnl,iunit,junit,imirror)
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

        ! --------------------------------------------------------------------
        dist1 = sqrt(dist2)

        itype1 = iontype_mp(imp1)
        itype2 = iontype_mp(imp2)

     if((.not. inmisc%class_flag(CLASS%ION)) .or. itype1 <= 0 .or. itype1 > IONTYPE%MAX_ALL .or. itype2 <= 0 .or. itype2 > IONTYPE%MAX_ALL) then

            if (inmisc%i_temp_independent == 0) then
                pnl = ele_list(irep)%coefs(1, iele1)/dist1*exp(-dist1*rcdist)

            else if (inmisc%i_temp_independent == 1) then
                pnl = ele_list(irep)%coefs(1, iele1)*(1.0e0_PREC/dist1 - 0.5e0_PREC*rcdist) &
                      *exp(-dist1*rcdist)

            else if (inmisc%i_temp_independent == 2) then
                rk = dist1*rcdist
                pnl = ele_list(irep)%coefs(1, iele1)/dist1*exp(-rk) &
                      *(-(1.0e0_PREC + 0.5e0_PREC*rk)*inele%diele_dTcoef)

            endif

        else
            rsig = 1.0/inion%csigmame(itype1, itype2)
            xtanh = tanh((dist1 - inion%cdistme(itype1, itype2))*rsig)
            rek_corr = 1.0/(0.5*(ew + 5.2) + 0.5*(ew - 5.2)*xtanh)

            pnl = ek*rek_corr*ele_list(irep)%coefs(1, iele1)/dist1*exp(-dist1*rcdist)
        end if

        ! ------ reset charge for phosphate in 3SPN2 ------
        imptype1 = imp2type(imp1)
        imptype2 = imp2type(imp2)
        if (imptype1 == MPTYPE%DNA2_PHOS .AND. imptype2 == MPTYPE%PRO) then
            pnl = pnl*(-inele%dna2_phos_pro_charge/0.6)
        else if (imptype2 == MPTYPE%DNA2_PHOS .AND. imptype1 == MPTYPE%PRO) then
            pnl = pnl*(-inele%dna2_phos_pro_charge/0.6)
        end if
        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%ELE) = pnlet(E_TYPE%ELE) + pnl

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%ELE) = pnle_unit(iunit, junit, E_TYPE%ELE) + pnl
    end do
!$omp end do nowait

end subroutine simu_energy_ele

