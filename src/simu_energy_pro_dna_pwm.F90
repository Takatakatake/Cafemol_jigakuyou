!simu_energy_pro_dna_pwm
!>@brief Calculates protein-DNA sequence specific interacions.
!> The values are added into "pnlet(ENERGY%PRO_DNA_PWM)" and "pnle_unit(ENERGY%PRO_DNA_PWM)".

subroutine simu_energy_pro_dna_pwm(irep, pnle_unit, pnlet)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: pro_dna_pwm_list
    use var_setp, only: inpdpwm
    use var_struct, only: xyz_mp_rep, pxyz_mp_rep, &
        pro_mp2pdpwm_para, pro_mp2pdpwm_r0_cutoff, &
        pdpwm_mp_N_ca, pdpwm_mp_C_ca, pdpwm_mp_5_base, pdpwm_mp_3_base, &
        pdpwm_para_r0, pdpwm_para_aNC0, pdpwm_para_aSB0, pdpwm_para_a530, &
        pdpwm_local_coef, pdpwm_local_shift, pdpwm_para_pwm_A, pdpwm_para_pwm_C, &
        pdpwm_para_pwm_G, pdpwm_para_pwm_T, cmp2seq, imp2unit
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
    integer :: ipdpwm, imirror
    integer :: ica_0, ica_N, ica_C
    integer :: ibase_0, ibase_5, ibase_3
    integer :: isugar
    integer :: itmp1, itmp2, itmp3
    integer :: iunit, junit, kunit
    real(PREC) :: sigma, phi, phi2, ene_shift, ene_factor

    ! force calculation variables ------------
    ! native values:
    real(PREC) :: r0, t10, t20, t30
    real(PREC) :: ene_pwm
    real(PREC) :: local_coef
    real(PREC) :: local_shift
    ! vectors:
    real(PREC) :: vbcx, vbcy, vbcz
    real(PREC) :: vsbx, vsby, vsbz
    real(PREC) :: v53x, v53y, v53z
    real(PREC) :: vncx, vncy, vncz
    ! inner-products:
    real(PREC) :: ip_bc_sb
    real(PREC) :: ip_bc_53
    real(PREC) :: ip_bc_nc
    real(PREC) :: ip_bc2
    real(PREC) :: ip_sb2
    real(PREC) :: ip_532
    real(PREC) :: ip_nc2
    ! distances:
    real(PREC) :: nm_bc, nm_nc, nm_53, nm_sb
    ! angles:
    real(PREC) :: cost1, cost2, cost3, t1, t2, t3
    ! Energy modulating factors:
    real(PREC) :: dr, sig_sqr, sig_sqr_2, f
    real(PREC) :: dt1, dt2, dt3, cosdt1, cosdt2, cosdt3
    real(PREC) :: ktheta_2, ktheta, g1, g2, g3
    ! energy:
    real(PREC) :: energy_local
    ! ----------------------------------------
    character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR3
    integer :: klen
#endif

    ! write (*,*) " simu_energy_pro_dna_pwm begins here: "
    ! write (*,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

    sigma = inpdpwm%pro_dna_pwm_sigma
    phi = inpdpwm%pro_dna_pwm_phi
    phi2 = phi*2.0e0_PREC
    ene_shift = inpdpwm%pro_dna_pwm_en0
    ene_factor = inpdpwm%pro_dna_pwm_factor

    sig_sqr = sigma*sigma
    sig_sqr_2 = 2.0e0_PREC*sig_sqr
    ktheta_2 = F_PI/phi
    ktheta = ktheta_2/2.0e0_PREC

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
    klen = (pro_dna_pwm_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, pro_dna_pwm_list(irep)%num_pairs)
#else
    ksta = 1
    kend = pro_dna_pwm_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
    print *, "pnl          = ", kend - ksta + 1
#endif
#else
    ksta = 1
    kend = pro_dna_pwm_list(irep)%num_pairs
#endif

!$omp   do private(ica_0,ica_N,ica_C,ibase_0,ibase_5,ibase_3,isugar,imirror, &
!$omp  & itmp1,itmp2,itmp3,                                                              &
!$omp  & vbcx, vbcy, vbcz, vsbx, vsby, vsbz, v53x, v53y,                                 &
!$omp  & v53z, vncx, vncy, vncz, ip_bc_sb, ip_bc_53, ip_bc_nc, ip_bc2, ip_sb2, ip_532,   &
!$omp  & ip_nc2, nm_bc, nm_nc, nm_53, nm_sb, cost1, cost2, cost3, t1, t2, t3, dr,        &
!$omp  & f, dt1, dt2, dt3, cosdt1, cosdt2, cosdt3,                                       &
!$omp  & g1, g2, g3, energy_local,                                                       &
!$omp  & r0, t10, t20, t30,                                                              &
!$omp  & ene_pwm,                                                                        &
!$omp  & local_shift,                                                                    &
!$omp  & local_coef,                                                                     &
!$omp  & error_message)

    do ipdpwm = ksta, kend
        ica_0 = pro_dna_pwm_list(irep)%pairs(1, ipdpwm)
        ibase_0 = pro_dna_pwm_list(irep)%pairs(2, ipdpwm)
        ! write (*,*) ica_0, ibase_0

        ! =========================== compute distance from base to Ca ===================
        if (inperi%i_periodic == 0) then
            vbcx = xyz_mp_rep(1, ica_0, irep) - xyz_mp_rep(1, ibase_0, irep)
            vbcy = xyz_mp_rep(2, ica_0, irep) - xyz_mp_rep(2, ibase_0, irep)
            vbcz = xyz_mp_rep(3, ica_0, irep) - xyz_mp_rep(3, ibase_0, irep)
        else
            imirror = pro_dna_pwm_list(irep)%pairs(3, ipdpwm)
            vbcx = pxyz_mp_rep(1, ica_0, irep) - pxyz_mp_rep(1, ibase_0, irep) + inperi%d_mirror(1, imirror)
            vbcy = pxyz_mp_rep(2, ica_0, irep) - pxyz_mp_rep(2, ibase_0, irep) + inperi%d_mirror(2, imirror)
            vbcz = pxyz_mp_rep(3, ica_0, irep) - pxyz_mp_rep(3, ibase_0, irep) + inperi%d_mirror(3, imirror)
        end if
        ip_bc2 = vbcx*vbcx + vbcy*vbcy + vbcz*vbcz
        ! ================================================================================

        ! ============================= cycle if dist > cutoff ===========================
        if (ip_bc2 >= pro_mp2pdpwm_r0_cutoff(ica_0)) then
            cycle
        end if
        ! ================================================================================

        ! ======================== assign index for other particles ======================
        ica_N = pdpwm_mp_N_ca(ica_0)
        ica_C = pdpwm_mp_C_ca(ica_0)
        ibase_5 = pdpwm_mp_5_base(ibase_0)
        ibase_3 = pdpwm_mp_3_base(ibase_0)
        isugar = ibase_0 - 1
        ! ================================================================================
        ! =========================== cycle if at the end of DNA =========================
        if (ibase_5 == ibase_0 .or. ibase_3 == ibase_0) then
            cycle
        end if
        ! ================================================================================

        ! ======================== compute vector coordinates ======================
        vsbx = xyz_mp_rep(1, isugar, irep) - xyz_mp_rep(1, ibase_0, irep)
        vsby = xyz_mp_rep(2, isugar, irep) - xyz_mp_rep(2, ibase_0, irep)
        vsbz = xyz_mp_rep(3, isugar, irep) - xyz_mp_rep(3, ibase_0, irep)

        v53x = xyz_mp_rep(1, ibase_3, irep) - xyz_mp_rep(1, ibase_5, irep)
        v53y = xyz_mp_rep(2, ibase_3, irep) - xyz_mp_rep(2, ibase_5, irep)
        v53z = xyz_mp_rep(3, ibase_3, irep) - xyz_mp_rep(3, ibase_5, irep)

        vncx = xyz_mp_rep(1, ica_N, irep) - xyz_mp_rep(1, ica_C, irep)
        vncy = xyz_mp_rep(2, ica_N, irep) - xyz_mp_rep(2, ica_C, irep)
        vncz = xyz_mp_rep(3, ica_N, irep) - xyz_mp_rep(3, ica_C, irep)

        ! ==================== inner-products ====================
        ip_bc_sb = vbcx*vsbx + vbcy*vsby + vbcz*vsbz
        ip_bc_53 = vbcx*v53x + vbcy*v53y + vbcz*v53z
        ip_bc_nc = vbcx*vncx + vbcy*vncy + vbcz*vncz
        ip_sb2 = vsbx*vsbx + vsby*vsby + vsbz*vsbz
        ip_532 = v53x*v53x + v53y*v53y + v53z*v53z
        ip_nc2 = vncx*vncx + vncy*vncy + vncz*vncz

        ! ==================== norm of vectors ====================
        nm_bc = sqrt(ip_bc2)
        nm_sb = sqrt(ip_sb2)
        nm_53 = sqrt(ip_532)
        nm_nc = sqrt(ip_nc2)

        ! ============================== angles ==============================
        cost1 = ip_bc_sb/(nm_bc*nm_sb)
        cost2 = ip_bc_53/(nm_bc*nm_53)
        cost3 = ip_bc_nc/(nm_bc*nm_nc)
        if (cost1 > 1.0e0_PREC) then
            cost1 = 1.0e0_PREC
        else if (cost1 < -1.0e0_PREC) then
            cost1 = -1.0e0_PREC
        end if
        if (cost2 > 1.0e0_PREC) then
            cost2 = 1.0e0_PREC
        else if (cost2 < -1.0e0_PREC) then
            cost2 = -1.0e0_PREC
        end if
        if (cost3 > 1.0e0_PREC) then
            cost3 = 1.0e0_PREC
        else if (cost3 < -1.0e0_PREC) then
            cost3 = -1.0e0_PREC
        end if
        t1 = acos(cost1)
        t2 = acos(cost2)
        t3 = acos(cost3)

        ! ========== looping interaction_pair ==========
        ! energy_local = 1.0e30_PREC
        do itmp1 = 1, MXMPPDPWM
            itmp2 = pro_mp2pdpwm_para(itmp1, ica_0)
            if (itmp2 == 0) then
                exit
            end if
            r0 = pdpwm_para_r0(itmp2)
            t10 = pdpwm_para_aSB0(itmp2)
            t20 = pdpwm_para_a530(itmp2)
            t30 = pdpwm_para_aNC0(itmp2)
            dt1 = t1 - t10; dt2 = t2 - t20; dt3 = t3 - t30
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ range of angles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (dt1 >= phi2 .or. dt1 <= -phi2 .or. &
                dt2 >= phi2 .or. dt2 <= -phi2 .or. &
                dt3 >= phi2 .or. dt3 <= -phi2) then
                cycle
            end if
            if (dt1 < phi .and. dt1 > -phi) then
                g1 = 1.0e0_PREC
            else
                cosdt1 = cos(ktheta*dt1)
                g1 = 1.0e0_PREC - cosdt1*cosdt1
            end if
            if (dt2 < phi .and. dt2 > -phi) then
                g2 = 1.0e0_PREC
            else
                cosdt2 = cos(ktheta*dt2)
                g2 = 1.0e0_PREC - cosdt2*cosdt2
            end if
            if (dt3 < phi .and. dt3 > -phi) then
                g3 = 1.0e0_PREC
            else
                cosdt3 = cos(ktheta*dt3)
                g3 = 1.0e0_PREC - cosdt3*cosdt3
            end if
            ! =============================================================================
            dr = nm_bc - r0
            f = exp(-(dr*dr)/sig_sqr_2)
            ! -------------------- read pwm energy --------------------
            if (cmp2seq(ibase_0) == 'DA ') then
                ene_pwm = pdpwm_para_pwm_A(itmp2)
            else if (cmp2seq(ibase_0) == 'DG ') then
                ene_pwm = pdpwm_para_pwm_G(itmp2)
            else if (cmp2seq(ibase_0) == 'DT ') then
                ene_pwm = pdpwm_para_pwm_T(itmp2)
            else if (cmp2seq(ibase_0) == 'DC ') then
                ene_pwm = pdpwm_para_pwm_C(itmp2)
            else
                write (error_message, *) ' Name of DNA particle in PDPWM is wrong! WTF! => ', ibase_0
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            local_coef = pdpwm_local_coef(itmp2)
            local_shift = pdpwm_local_shift(itmp2)
            ! energy_local = ene_factor * (ene_pwm + ene_shift) * f * g1 * g2 * g3
            energy_local = ene_factor*local_coef*(ene_pwm + ene_shift + local_shift)*f*g1*g2*g3
            ! write (*,*) " ene_local from simu_ene_pdpwm: ", energy_local

            pnlet(E_TYPE%PRO_DNA_PWM) = pnlet(E_TYPE%PRO_DNA_PWM) + energy_local

            iunit = imp2unit(ica_0)
            junit = imp2unit(ibase_0)
            if (iunit > junit) then
                kunit = junit; junit = iunit; iunit = kunit
            end if
            pnle_unit(iunit, junit, E_TYPE%PRO_DNA_PWM) = pnle_unit(iunit, junit, E_TYPE%PRO_DNA_PWM) + energy_local
        end do
    end do
!$omp   end do nowait

end subroutine simu_energy_pro_dna_pwm
