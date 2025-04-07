!simu_energy_exv_wca
!> @brief Calculates the energy related to excluded volume.   &
!>        The values are added into "pnlet(ENERGY%EXV_WCA)" and  &
!>        "pnle_unit(ENERGY%EXV_WCA)".
!
! Weeks−Chandler−Andersen (WCA) potential
!
! Reference:
!   Equation (2) in
!   N.A. Denesyuk and D. Thirumalai, J Phys. Chem. B (2013)
!   The original paper is
!   D. Chandler, J.D. Weeks and H.C. Andersen, Science (1983)
!
! Potential function:
!    U_ev = epsilon0 [ (D0/r)^12 - 2(D0/r)^6 + 1 ]   (r <= D0)
!         = 0    (r > D0)
!
! Coefficients:
!     epsilon0 = coef_exvwca
!        D0    = cdist_exvwca

subroutine simu_energy_exv_wca(irep, pnle_unit, pnlet)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: exv_wca_list
    use var_setp, only: indtrna
    use var_struct, only: imp2unit, xyz_mp_rep, pxyz_mp_rep
#ifdef MPI_PAR3
    use mpiconst
#endif

    implicit none

    ! ------------------------------------------------------------------------
    integer, intent(in)  :: irep
    real(PREC), intent(out) :: pnlet(:)         ! (E_TYPE%MAX)
    real(PREC), intent(out) :: pnle_unit(:, :, :) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

    ! ------------------------------------------------------------------------
    ! local variables
    integer :: klen, ksta, kend
    integer :: imp1, imp2, iunit, junit
    integer :: ipnl, imirror
    real(PREC) :: dist2
    real(PREC) :: coef, cdist2
    !real(PREC) :: coef, coef_pro, coef_rna, coef_rna_pro, coef_lig
    !real(PREC) :: cdist2, cdist2_pro, cdist2_rna, cdist2_rna_pro, &
    !              cdist2_llig, cdist2_lpro
    real(PREC) :: roverdist2, roverdist4, roverdist6, roverdist12
    real(PREC) :: ene
    real(PREC) :: v21(SPACE_DIM)

    ! ------------------------------------------------------------------------
    !! Currently this potential is available noly for RNA.
    cdist2 = indtrna%cdist_exvwca**2
    coef = indtrna%coef_exvwca
    !cdist2_rna = inrna%cdist_exvwca**2
    !cdist2_rna_pro = (inrna%cdist_exvwca + inpro%cdist_exvwca) ** 2 / 4
    !cdist2_pro = inpro%cdist_exvwca * inpro%cdist_exvwca
    !cdist2_llig = inligand%cdist_exvwca_llig * inligand%cdist_exvwca_llig
    !cdist2_lpro = inligand%cdist_exvwca_lpro * inligand%cdist_exvwca_lpro
    !coef_rna = inrna%coef_exvwca
    !coef_pro = inpro%coef_exvwca
    !coef_rna_pro = inrna%coef_exvwca
    !coef_lig = inligand%coef_exvwca

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
    klen = (exv_wca_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, exv_wca_list(irep)%num_pairs)
#else
    ksta = 1
    kend = exv_wca_list(irep)%num_pairs
#endif
#else
    ksta = 1
    kend = exv_wca_list(irep)%num_pairs
#endif
!$omp do private(imp1,imp2,v21,dist2, &
!$omp&           roverdist2,roverdist4,roverdist6,roverdist12, &
!$omp&           ene,iunit,junit,imirror)
    do ipnl = ksta, kend

        imp1 = exv_wca_list(irep)%pairs(1, ipnl)
        imp2 = exv_wca_list(irep)%pairs(2, ipnl)

        if (inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
            imirror = exv_wca_list(irep)%pairs(3, ipnl)
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

        dist2 = dot_product(v21, v21)

        !
        !! Currently this potential is available noly for RNA.
        !
        !if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        !   cdist2  = cdist2_rna
        !   coef    = coef_rna
        !else if ((iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%PRO) .OR. &
        !         (iclass_mp(imp1) == CLASS%PRO .AND. iclass_mp(imp2) == CLASS%RNA)) then
        !   cdist2  = cdist2_rna_pro
        !   coef    = coef_rna_pro
        !else if(iclass_mp(imp1) == CLASS%LIG .OR. iclass_mp(imp2) == CLASS%LIG) then
        !   cdist2  = cdist2_llig
        !   coef    = coef_lig
        !   if(iclass_mp(imp1) == CLASS%PRO .OR. iclass_mp(imp2) == CLASS%PRO) then
        !      cdist2  = cdist2_lpro
        !   end if
        !else
        !   cdist2  = cdist2_pro
        !   coef    = coef_pro
        !endif

        if (dist2 > cdist2) cycle

        ! --------------------------------------------------------------------
        roverdist2 = cdist2/dist2
        roverdist4 = roverdist2*roverdist2
        roverdist6 = roverdist2*roverdist4
        roverdist12 = roverdist6*roverdist6
        ene = coef*(roverdist12 - 2*roverdist6 + 1.0e0_PREC)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%EXV_WCA) = pnlet(E_TYPE%EXV_WCA) + ene

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%EXV_WCA) = pnle_unit(iunit, junit, E_TYPE%EXV_WCA) + ene
    end do
!$omp end do nowait

end subroutine simu_energy_exv_wca
