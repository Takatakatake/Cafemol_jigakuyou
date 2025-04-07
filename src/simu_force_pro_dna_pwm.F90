!simu_force_pnl_restype
!> @brief Calculates the force related to protein-DNA sequence-specific intaractions

subroutine simu_force_pro_dna_pwm(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: pro_dna_pwm_list
    use var_setp, only: inpdpwm

    use var_struct, only: nmp_all, xyz_mp_rep, pxyz_mp_rep, &
        pro_mp2pdpwm_para, pro_mp2pdpwm_r0_cutoff, &
        pdpwm_mp_N_ca, pdpwm_mp_C_ca, pdpwm_mp_5_base, pdpwm_mp_3_base, &
        pdpwm_para_r0, pdpwm_para_aNC0, pdpwm_para_aSB0, pdpwm_para_a530, &
        pdpwm_para_pwm_A, pdpwm_para_pwm_C, pdpwm_para_pwm_G, pdpwm_para_pwm_T, &
        pdpwm_local_coef, pdpwm_local_shift, cmp2seq

    use mpiconst

    implicit none

    ! --------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: ipdpwm, imirror
    integer :: ica_0, ica_N, ica_C
    integer :: ibase_0, ibase_5, ibase_3
    integer :: isugar
    integer :: itmp1, itmp2
    real(PREC) :: sigma, phi, phi2, ene_shift, ene_factor

    ! force calculation variables ------------
    ! native values:
    real(PREC) :: r0, t10, t20, t30
    real(PREC) :: ene_pwm
    real(PREC) :: local_shift
    real(PREC) :: local_coef
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
    real(PREC) :: ktheta_2, ktheta, g1, g2, g3, gg1, gg2, gg3
    ! energy:
    real(PREC) :: energy_local
    ! force coefficients:
    real(PREC) :: k0, k1, k2, k3
    real(PREC) :: q0, q1bc, q1sb, q2bc, q253, q3bc, q3nc
    real(PREC) :: sqrt_small_value
    ! forces blocks:
    real(PREC) :: f_BC_0_x, f_BC_0_y, f_BC_0_z
    real(PREC) :: f_BC_1_x, f_BC_1_y, f_BC_1_z
    real(PREC) :: f_BC_2_x, f_BC_2_y, f_BC_2_z
    real(PREC) :: f_BC_3_x, f_BC_3_y, f_BC_3_z
    real(PREC) :: f_SB_x, f_SB_y, f_SB_z
    real(PREC) :: f_53_x, f_53_y, f_53_z
    real(PREC) :: f_NC_x, f_NC_y, f_NC_z
    ! forces:
    real(PREC) :: f_C0_x, f_C0_y, f_C0_z
    real(PREC) :: f_B0_x, f_B0_y, f_B0_z
    ! ----------------------------------------
    character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR
    integer :: klen
#endif

    ! write (*,*) " simu_force_pro_dna_pwm begins here: "
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

#ifdef MPI_PAR
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
!$omp  & itmp1,itmp2,                                                                    &
!$omp  & r0, t10, t20, t30, ene_pwm, vbcx, vbcy, vbcz, vsbx, vsby, vsbz, v53x, v53y,     &
!$omp  & v53z, vncx, vncy, vncz, ip_bc_sb, ip_bc_53, ip_bc_nc, ip_bc2, ip_sb2, ip_532,   &
!$omp  & ip_nc2, nm_bc, nm_nc, nm_53, nm_sb, cost1, cost2, cost3, t1, t2, t3, dr,        &
!$omp  & f, dt1, dt2, dt3, cosdt1, cosdt2, cosdt3,                                       &
!$omp  & g1, g2, g3, energy_local, k0, k1, k2, k3, q0, q1bc, q1sb, q2bc, q253,           &
!$omp  & q3bc, q3nc, f_BC_0_x, f_BC_0_y, f_BC_0_z, f_BC_1_x, f_BC_1_y, f_BC_1_z,         &
!$omp  & f_BC_2_x, f_BC_2_y, f_BC_2_z, f_BC_3_x, f_BC_3_y, f_BC_3_z, f_SB_x, f_SB_y,     &
!$omp  & f_SB_z, f_53_x, f_53_y, f_53_z, f_NC_x, f_NC_y, f_NC_z, f_C0_x, f_C0_y, f_C0_z, &
!$omp  & f_B0_x, f_B0_y, f_B0_z, gg1, gg2, gg3, sqrt_small_value,                        &
!$omp  & local_shift, local_coef,                                                        &
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

        ! ========== guess interaction_pair ==========
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
            if (dt1 <= phi .and. dt1 >= -phi) then
                g1 = 1.0e0_PREC
            else
                cosdt1 = cos(ktheta*dt1)
                g1 = 1.0e0_PREC - cosdt1*cosdt1
            end if
            if (dt2 <= phi .and. dt2 >= -phi) then
                g2 = 1.0e0_PREC
            else
                cosdt2 = cos(ktheta*dt2)
                g2 = 1.0e0_PREC - cosdt2*cosdt2
            end if
            if (dt3 <= phi .and. dt3 >= -phi) then
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
            ! write (*,*) " ene_local from simu_force_pdpwm: ", energy_local

            ! ================================================================================
            ! ============================== forces calculation ==============================
            ! ---------- k0, k1, k2, k3 ----------
            k0 = energy_local*dr/sig_sqr
            if (g1 <= ZERO_JUDGE) then
                gg1 = ene_factor*local_coef*(ene_pwm + ene_shift + local_shift)*f*g2*g3
            else
                gg1 = energy_local/g1
            end if
            if (g2 <= ZERO_JUDGE) then
                gg2 = ene_factor*local_coef*(ene_pwm + ene_shift + local_shift)*f*g1*g3
            else
                gg2 = energy_local/g2
            end if
            if (g3 <= ZERO_JUDGE) then
                gg3 = ene_factor*local_coef*(ene_pwm + ene_shift + local_shift)*f*g2*g1
            else
                gg3 = energy_local/g3
            end if
            if (dt1 <= phi .and. dt1 >= -phi) then
                f_BC_1_x = 0.0e0_PREC; f_BC_1_y = 0.0e0_PREC; f_BC_1_z = 0.0e0_PREC
                f_SB_x = 0.0e0_PREC; f_SB_y = 0.0e0_PREC; f_SB_z = 0.0e0_PREC
            else
                sqrt_small_value = sqrt(ip_bc2*ip_sb2 - ip_bc_sb*ip_bc_sb)
                if (sqrt_small_value <= ZERO_JUDGE) then
                    sqrt_small_value = ZERO_JUDGE
                end if
                k1 = (gg1*ktheta*sin(ktheta_2*dt1))/sqrt_small_value
                q1bc = k1*ip_bc_sb/ip_bc2; q1sb = k1*ip_bc_sb/ip_sb2
                f_BC_1_x = k1*vsbx - q1bc*vbcx; f_BC_1_y = k1*vsby - q1bc*vbcy; f_BC_1_z = k1*vsbz - q1bc*vbcz
                f_SB_x = k1*vbcx - q1sb*vsbx; f_SB_y = k1*vbcy - q1sb*vsby; f_SB_z = k1*vbcz - q1sb*vsbz
            end if
            if (dt2 <= phi .and. dt2 >= -phi) then
                f_BC_2_x = 0.0e0_PREC; f_BC_2_y = 0.0e0_PREC; f_BC_2_z = 0.0e0_PREC
                f_53_x = 0.0e0_PREC; f_53_y = 0.0e0_PREC; f_53_z = 0.0e0_PREC
            else
                sqrt_small_value = sqrt(ip_bc2*ip_532 - ip_bc_53*ip_bc_53)
                if (sqrt_small_value <= ZERO_JUDGE) then
                    sqrt_small_value = ZERO_JUDGE
                end if
                k2 = (gg2*ktheta*sin(ktheta_2*dt2))/sqrt_small_value
                q2bc = k2*ip_bc_53/ip_bc2; q253 = k2*ip_bc_53/ip_532
                f_BC_2_x = k2*v53x - q2bc*vbcx; f_BC_2_y = k2*v53y - q2bc*vbcy; f_BC_2_z = k2*v53z - q2bc*vbcz
                f_53_x = k2*vbcx - q253*v53x; f_53_y = k2*vbcy - q253*v53y; f_53_z = k2*vbcz - q253*v53z
            end if
            if (dt3 <= phi .and. dt3 >= -phi) then
                f_BC_3_x = 0.0e0_PREC; f_BC_3_y = 0.0e0_PREC; f_BC_3_z = 0.0e0_PREC
                f_NC_x = 0.0e0_PREC; f_NC_y = 0.0e0_PREC; f_NC_z = 0.0e0_PREC
            else
                sqrt_small_value = sqrt(ip_bc2*ip_nc2 - ip_bc_nc*ip_bc_nc)
                if (sqrt_small_value <= ZERO_JUDGE) then
                    sqrt_small_value = ZERO_JUDGE
                end if
                k3 = (gg3*ktheta*sin(ktheta_2*dt3))/sqrt_small_value
                q3bc = k3*ip_bc_nc/ip_bc2; q3nc = k3*ip_bc_nc/ip_nc2
                f_BC_3_x = k3*vncx - q3bc*vbcx; f_BC_3_y = k3*vncy - q3bc*vbcy; f_BC_3_z = k3*vncz - q3bc*vbcz
                f_NC_x = k3*vbcx - q3nc*vncx; f_NC_y = k3*vbcy - q3nc*vncy; f_NC_z = k3*vbcz - q3nc*vncz
            end if
            q0 = k0/nm_bc
            f_BC_0_x = q0*vbcx; f_BC_0_y = q0*vbcy; f_BC_0_z = q0*vbcz
            ! ~~~~~~~~~~~~~~~~~~~~ oh, at last! ~~~~~~~~~~~~~~~~~~~~
            f_C0_x = f_BC_0_x + f_BC_1_x + f_BC_2_x + f_BC_3_x
            f_C0_y = f_BC_0_y + f_BC_1_y + f_BC_2_y + f_BC_3_y
            f_C0_z = f_BC_0_z + f_BC_1_z + f_BC_2_z + f_BC_3_z
            f_B0_x = -f_C0_x - f_SB_x
            f_B0_y = -f_C0_y - f_SB_y
            f_B0_z = -f_C0_z - f_SB_z

            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ update global force list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! write (*,*) " ~~~~~~~~~~ force check ~~~~~~~~~~"
            ! write (*,*) " pair: ", ica_0, ibase_0
            ! write (*,*) ibase_0, " | => ", f_B0_x, f_B0_y, f_B0_z
            ! write (*,*) ica_0, " | => ", f_C0_x, f_C0_y, f_C0_z
            ! write (*,*) isugar, " | => ", f_SB_x, f_SB_y, f_SB_z
            ! write (*,*) ibase_5, " | => ", - f_53_x, - f_53_y, - f_53_z
            ! write (*,*) ibase_3, " | => ", f_53_x, f_53_y, f_53_z
            ! write (*,*) ica_N, " | => ", f_NC_x, f_NC_y, f_NC_z
            ! write (*,*) ica_C, " | => ", -f_NC_x, -f_NC_y, -f_NC_z
            ! write (*,*) " ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~"

            force_mp(1, ibase_0) = force_mp(1, ibase_0) + f_B0_x
            force_mp(2, ibase_0) = force_mp(2, ibase_0) + f_B0_y
            force_mp(3, ibase_0) = force_mp(3, ibase_0) + f_B0_z
            force_mp(1, ica_0) = force_mp(1, ica_0) + f_C0_x
            force_mp(2, ica_0) = force_mp(2, ica_0) + f_C0_y
            force_mp(3, ica_0) = force_mp(3, ica_0) + f_C0_z
            force_mp(1, isugar) = force_mp(1, isugar) + f_SB_x
            force_mp(2, isugar) = force_mp(2, isugar) + f_SB_y
            force_mp(3, isugar) = force_mp(3, isugar) + f_SB_z
            force_mp(1, ibase_3) = force_mp(1, ibase_3) + f_53_x
            force_mp(2, ibase_3) = force_mp(2, ibase_3) + f_53_y
            force_mp(3, ibase_3) = force_mp(3, ibase_3) + f_53_z
            force_mp(1, ibase_5) = force_mp(1, ibase_5) - f_53_x
            force_mp(2, ibase_5) = force_mp(2, ibase_5) - f_53_y
            force_mp(3, ibase_5) = force_mp(3, ibase_5) - f_53_z
            force_mp(1, ica_N) = force_mp(1, ica_N) + f_NC_x
            force_mp(2, ica_N) = force_mp(2, ica_N) + f_NC_y
            force_mp(3, ica_N) = force_mp(3, ica_N) + f_NC_z
            force_mp(1, ica_C) = force_mp(1, ica_C) - f_NC_x
            force_mp(2, ica_C) = force_mp(2, ica_C) - f_NC_y
            force_mp(3, ica_C) = force_mp(3, ica_C) - f_NC_z

            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DONE! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        end do
    end do
!$omp   end do nowait

end subroutine simu_force_pro_dna_pwm
