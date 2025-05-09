!simu_force_pnl
!> @brief Calculates the force related to excluded volume    &
!>        excepting DNA particle.

subroutine simu_force_pnl(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: exv_list
    use var_setp, only: inpro, inrna, inligand, inmisc
    use var_struct, only: nmp_all, xyz_mp_rep, pxyz_mp_rep, iclass_mp
    use mpiconst

    implicit none

    ! --------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: imp1, imp2
    integer :: ipnl, imirror
    real(PREC) :: dist2
    real(PREC) :: coef, coef_pro, coef_rna, coef_rna_pro, coef_llig, coef_lpro
    real(PREC) :: cdist2, cdist2_pro, cdist2_rna, cdist2_rna_pro, &
                  cdist2_llig, cdist2_lpro
    real(PREC) :: cutoff2, cutoff2_pro, cutoff2_rna, cutoff2_rna_pro, &
                  cutoff2_llig, cutoff2_lpro
    real(PREC) :: roverdist2, roverdist4, roverdist8
    real(PREC) :: roverdist14
    real(PREC) :: dvdw_dr
    real(PREC) :: v21(SPACE_DIM), for(SPACE_DIM)
#ifdef MPI_PAR
    integer :: klen
#endif

    ! --------------------------------------------------------------------
    ! exvol protein
    ! for speed up
    cutoff2_pro = (inpro%cutoff_exvol*inpro%cdist_rep12)**2
    cdist2_pro = inpro%cdist_rep12**2
    coef_pro = 12.0e0_PREC*inpro%crep12/cdist2_pro
    if (inmisc%class_flag(CLASS%RNA)) then
        cutoff2_rna = (inrna%cutoff_exvol*inrna%cdist_rep12)**2
        cdist2_rna = inrna%cdist_rep12**2
        coef_rna = 12.0e0_PREC*inrna%crep12/cdist2_rna
        cutoff2_rna_pro = (inrna%cutoff_exvol*inrna%cdist_rep12 + inpro%cutoff_exvol*inpro%cdist_rep12)**2/4
        cdist2_rna_pro = (inrna%cdist_rep12 + inpro%cdist_rep12)**2/4
        coef_rna_pro = 12.0e0_PREC*inrna%crep12/cdist2_rna_pro
    endif
    if (inmisc%class_flag(CLASS%LIG)) then
        cutoff2_llig = (inligand%cutoff_exvol*inligand%cdist_rep12_llig)**2
        cdist2_llig = inligand%cdist_rep12_llig**2
        coef_llig = 12.0e0_PREC*inligand%crep12/cdist2_llig
        cutoff2_lpro = (inligand%cutoff_exvol*inligand%cdist_rep12_lpro)**2
        cdist2_lpro = inligand%cdist_rep12_lpro**2
        coef_lpro = 12.0e0_PREC*inligand%crep12/cdist2_lpro
    endif

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
    klen = (exv_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, exv_list(irep)%num_pairs)
#else
    ksta = 1
    kend = exv_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
    print *, "pnl          = ", kend - ksta + 1
#endif
#else
    ksta = 1
    kend = exv_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,cutoff2,cdist2,coef,&
!$OMP&           roverdist2,roverdist4, &
!$OMP&           roverdist8,roverdist14,dvdw_dr,for,imirror)
    do ipnl = ksta, kend

        imp1 = exv_list(irep)%pairs(1, ipnl)
        imp2 = exv_list(irep)%pairs(2, ipnl)

        if (inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
            imirror = exv_list(irep)%pairs(3, ipnl)
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
            cutoff2 = cutoff2_rna
            cdist2 = cdist2_rna
            coef = coef_rna
        else if ((iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%PRO) .OR. &
                 (iclass_mp(imp1) == CLASS%PRO .AND. iclass_mp(imp2) == CLASS%RNA)) then
            cutoff2 = cutoff2_rna_pro
            cdist2 = cdist2_rna_pro
            coef = coef_rna_pro
        else if (iclass_mp(imp1) == CLASS%LIG .OR. iclass_mp(imp2) == CLASS%LIG) then
            cutoff2 = cutoff2_llig
            cdist2 = cdist2_llig
            coef = coef_llig
            if (iclass_mp(imp1) == CLASS%PRO .OR. iclass_mp(imp2) == CLASS%PRO) then
                cutoff2 = cutoff2_lpro
                cdist2 = cdist2_lpro
                coef = coef_lpro
            end if
        else
            cutoff2 = cutoff2_pro
            cdist2 = cdist2_pro
            coef = coef_pro
        endif

        if (dist2 > cutoff2) cycle

        ! -----------------------------------------------------------------
        roverdist2 = cdist2/dist2
        roverdist4 = roverdist2*roverdist2
        roverdist8 = roverdist4*roverdist4
        roverdist14 = roverdist2*roverdist4*roverdist8
        dvdw_dr = coef*roverdist14
        if (dvdw_dr > DE_MAX) then
            !write (*, *) "exvol protein", imp1, imp2, dvdw_dr
            dvdw_dr = DE_MAX
        end if
!     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

        for(1:3) = dvdw_dr*v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
    end do
!$OMP end do nowait

end subroutine simu_force_pnl
