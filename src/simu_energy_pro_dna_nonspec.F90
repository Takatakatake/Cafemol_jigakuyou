!simu_energy_pro_dna_nonspec
!>@brief Calculates protein-DNA sequence specific interacions.
!> The values are added into "pnlet(ENERGY%PRO_DNA_NONSPEC)" and "pnle_unit(ENERGY%PRO_DNA_NONSPEC)".

subroutine simu_energy_pro_dna_nonspec(irep, pnle_unit, pnlet)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: pro_dna_nonspec_list
    use var_setp, only: inpdns
    use var_struct, only: xyz_mp_rep, pxyz_mp_rep, &
        pro_mp2pdns_para, pro_mp2pdns_r0_cutoff, &
        pdns_mp_N_ca, pdns_mp_C_ca, &
        pdns_para_r0, pdns_para_aNC0, pdns_para_aSB0, &
        pdns_coef, imp2unit
#ifdef MPI_PAR3
    use mpiconst
#endif

    implicit none

    ! ------------------------------------------------------------------------
    integer, intent(in)  :: irep
    real(PREC), intent(out) :: pnlet(:)         ! (E_TYPE%MAX)
    real(PREC), intent(out) :: pnle_unit(:, :, :) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: ipdns, imirror
    integer :: ica_0, ica_N, ica_C
    integer :: iphos_0
    integer :: isugar
    integer :: itmp1, itmp2, itmp3
    integer :: iunit, junit, kunit
    real(PREC) :: sigma, phi, phi2, ene_factor
    character(CARRAY_MSG_ERROR) :: error_message

    ! force calculation variables ------------
    ! native values:
    real(PREC) :: r0_tmp, t10_tmp, t30_tmp
    ! vectors:
    real(PREC) :: vbcx, vbcy, vbcz
    real(PREC) :: vsbx, vsby, vsbz
    real(PREC) :: vncx, vncy, vncz
    ! inner-products:
    real(PREC) :: ip_bc_sb
    real(PREC) :: ip_bc_nc
    real(PREC) :: ip_bc2
    real(PREC) :: ip_sb2
    real(PREC) :: ip_nc2
    ! distances:
    real(PREC) :: nm_bc, nm_nc, nm_sb
    ! angles:
    real(PREC) :: cost1, cost3, t1, t3
    ! Energy modulating factors:
    real(PREC) :: dr, sig_sqr, sig_sqr_2, f
    real(PREC) :: dt1, dt3, cosdt1, cosdt3
    real(PREC) :: ktheta_2, ktheta, g1, g3
    ! energy:
    real(PREC) :: energy_local, e_tmp
#ifdef MPI_PAR3
    integer :: klen
#endif

    ! write (*,*) " simu_energy_pro_dna_nonspec begins here: "
    ! write (*,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

    sigma = inpdns%pro_dna_nonspec_sigma
    phi = inpdns%pro_dna_nonspec_phi
    phi2 = phi*2.0e0_PREC
    ene_factor = inpdns%pro_dna_nonspec_factor

    sig_sqr = sigma*sigma
    sig_sqr_2 = 2.0e0_PREC*sig_sqr
    ktheta_2 = F_PI/phi
    ktheta = ktheta_2/2.0e0_PREC

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
    klen = (pro_dna_nonspec_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, pro_dna_nonspec_list(irep)%num_pairs)
#else
    ksta = 1
    kend = pro_dna_nonspec_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
    print *, "pnl          = ", kend - ksta + 1
#endif
#else
    ksta = 1
    kend = pro_dna_nonspec_list(irep)%num_pairs
#endif

!$omp   do private(ica_0,ica_N,ica_C,iphos_0,isugar,imirror, &
!$omp  & itmp1,itmp2,itmp3,                                                              &
!$omp  & vbcx, vbcy, vbcz, vsbx, vsby, vsbz,                                             &
!$omp  & vncx, vncy, vncz, ip_bc_sb, ip_bc_nc, ip_bc2, ip_sb2,                           &
!$omp  & ip_nc2, nm_bc, nm_nc, nm_sb, cost1, cost3, t1, t3, dr,                          &
!$omp  & f, dt1, dt3, cosdt1, cosdt3,                                                    &
!$omp  & g1, g3, energy_local, e_tmp,                                                    &
!$omp  & r0_tmp, t10_tmp, t30_tmp,                                                       &
!$omp  & error_message)

    do ipdns = ksta, kend
        ica_0 = pro_dna_nonspec_list(irep)%pairs(1, ipdns)
        iphos_0 = pro_dna_nonspec_list(irep)%pairs(2, ipdns)
        ! write (*,*) ica_0, iphos_0

        ! =========================== compute distance from phos to Ca ===================
        if (inperi%i_periodic == 0) then
            vbcx = xyz_mp_rep(1, ica_0, irep) - xyz_mp_rep(1, iphos_0, irep)
            vbcy = xyz_mp_rep(2, ica_0, irep) - xyz_mp_rep(2, iphos_0, irep)
            vbcz = xyz_mp_rep(3, ica_0, irep) - xyz_mp_rep(3, iphos_0, irep)
        else
            imirror = pro_dna_nonspec_list(irep)%pairs(3, ipdns)
            vbcx = pxyz_mp_rep(1, ica_0, irep) - pxyz_mp_rep(1, iphos_0, irep) + inperi%d_mirror(1, imirror)
            vbcy = pxyz_mp_rep(2, ica_0, irep) - pxyz_mp_rep(2, iphos_0, irep) + inperi%d_mirror(2, imirror)
            vbcz = pxyz_mp_rep(3, ica_0, irep) - pxyz_mp_rep(3, iphos_0, irep) + inperi%d_mirror(3, imirror)
        end if
        ip_bc2 = vbcx*vbcx + vbcy*vbcy + vbcz*vbcz
        ! ================================================================================

        ! ============================= cycle if dist > cutoff ===========================
        if (ip_bc2 >= pro_mp2pdns_r0_cutoff(ica_0)) then
            cycle
        end if
        ! ================================================================================

        ! ======================== assign index for other particles ======================
        ica_N = pdns_mp_N_ca(ica_0)
        ica_C = pdns_mp_C_ca(ica_0)
        isugar = iphos_0 + 1
        ! ================================================================================

        ! ======================== compute vector coordinates ======================
        vsbx = xyz_mp_rep(1, isugar, irep) - xyz_mp_rep(1, iphos_0, irep)
        vsby = xyz_mp_rep(2, isugar, irep) - xyz_mp_rep(2, iphos_0, irep)
        vsbz = xyz_mp_rep(3, isugar, irep) - xyz_mp_rep(3, iphos_0, irep)

        vncx = xyz_mp_rep(1, ica_N, irep) - xyz_mp_rep(1, ica_C, irep)
        vncy = xyz_mp_rep(2, ica_N, irep) - xyz_mp_rep(2, ica_C, irep)
        vncz = xyz_mp_rep(3, ica_N, irep) - xyz_mp_rep(3, ica_C, irep)

        ! ==================== inner-products ====================
        ip_bc_sb = vbcx*vsbx + vbcy*vsby + vbcz*vsbz
        ip_bc_nc = vbcx*vncx + vbcy*vncy + vbcz*vncz
        ip_sb2 = vsbx*vsbx + vsby*vsby + vsbz*vsbz
        ip_nc2 = vncx*vncx + vncy*vncy + vncz*vncz

        ! ==================== norm of vectors ====================
        nm_bc = sqrt(ip_bc2)
        nm_sb = sqrt(ip_sb2)
        nm_nc = sqrt(ip_nc2)

        ! ============================== angles ==============================
        cost1 = ip_bc_sb/(nm_bc*nm_sb)
        cost3 = ip_bc_nc/(nm_bc*nm_nc)
        if (cost1 > 1.0e0_PREC) then
            cost1 = 1.0e0_PREC
        else if (cost1 < -1.0e0_PREC) then
            cost1 = -1.0e0_PREC
        end if
        if (cost3 > 1.0e0_PREC) then
            cost3 = 1.0e0_PREC
        else if (cost3 < -1.0e0_PREC) then
            cost3 = -1.0e0_PREC
        end if
        t1 = acos(cost1)
        t3 = acos(cost3)

        ! ========== guess interaction_pair ==========
        energy_local = 1.0e30_PREC
        do itmp1 = 1, MXMPPDNS
            itmp2 = pro_mp2pdns_para(itmp1, ica_0)
            if (itmp2 == 0) then
                exit
            end if
            r0_tmp = pdns_para_r0(itmp2)
            t10_tmp = pdns_para_aSB0(itmp2)
            t30_tmp = pdns_para_aNC0(itmp2)
            dt1 = t1 - t10_tmp; dt3 = t3 - t30_tmp
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ range of angles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (dt1 >= phi2 .or. dt1 <= -phi2 .or. &
                dt3 >= phi2 .or. dt3 <= -phi2) then
                cycle
            end if
            if (dt1 < phi .and. dt1 > -phi) then
                g1 = 1.0e0_PREC
            else
                cosdt1 = cos(ktheta*dt1)
                g1 = 1.0e0_PREC - cosdt1*cosdt1
            end if
            if (dt3 < phi .and. dt3 > -phi) then
                g3 = 1.0e0_PREC
            else
                cosdt3 = cos(ktheta*dt3)
                g3 = 1.0e0_PREC - cosdt3*cosdt3
            end if
            ! =============================================================================
            dr = nm_bc - r0_tmp
            f = exp(-(dr*dr)/sig_sqr_2)
            e_tmp = ene_factor*pdns_coef(itmp2)*f*g1*g3
            ! write (*,*) e_tmp
            if (e_tmp < energy_local) then
                energy_local = e_tmp
            end if
        end do
        if (energy_local > 1.0e5_PREC) cycle
        ! write (*,*) " ene_local from simu_ene_pdns: ", energy_local

        pnlet(E_TYPE%PRO_DNA_NONSPEC) = pnlet(E_TYPE%PRO_DNA_NONSPEC) + energy_local

        iunit = imp2unit(ica_0)
        junit = imp2unit(iphos_0)
        if (iunit > junit) then
            kunit = junit; junit = iunit; iunit = kunit
        end if
        pnle_unit(iunit, junit, E_TYPE%PRO_DNA_NONSPEC) = pnle_unit(iunit, junit, E_TYPE%PRO_DNA_NONSPEC) + energy_local
    end do
!$omp   end do nowait

end subroutine simu_energy_pro_dna_nonspec
