!mloop_nativeinfo
!> @brief Reads the filename of .ninfo by calling "mloop_ninfo_inp"    &
!>        in each simulation stage, and control native-info reading procedure.

subroutine mloop_nativeinfo(istep_sim)

    use const_maxsize
    use const_index
    use var_inp, only: iopen_lunnum
    use var_setp, only: inpro, indna, indna2, inrna, inligand, inmisc
    use var_struct, only: nmp_all, nunit_all, imp2unit, iclass_unit, &
        nbd, ibd2mp, bd_nat, factor_bd, coef_bd, nba, iba2mp, ba_nat, factor_ba, &
        coef_ba, ndih, idih2mp, dih_nat, dih_sin_nat, dih_cos_nat, factor_dih, &
        coef_dih, ncon, icon2mp, icon2unit, go_nat, go_nat2, factor_go, &
        icon_dummy_mgo, coef_go, get_icon_type, ncon_unit, imp2type, &
        iallcon2unit, nmorse, imorse2mp, imorse2unit, imorse_dummy_mgo, &
        morse_nat, morse_nat2, factor_morse, nrna_bp, irna_bp2mp, irna_bp2unit, &
        rna_bp_nat, rna_bp_nat2, factor_rna_bp, irna_bp_dummy_mgo, nhb_bp, &
        coef_rna_bp, nrna_st, irna_st2mp, irna_st2unit, rna_st_nat, rna_st_nat2, &
        factor_rna_st, irna_st_dummy_mgo, coef_rna_st, rna_base_type, ndtrna_st, &
        idtrna_st2mp, dtrna_st_nat, coef_dtrna_st, idtrna_st2nn, ndtrna_hb, &
        idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, ibd2type, iba2type, idih2type, &
        icon2type, imorse2type, pro_mp2pdpwm_r0_cutoff, &
        ipdpwm2pro_mp, pro_mp2pdpwm_para, pdpwm_para_r0, pdpwm_para_a530, & ! PDPWM
        pdpwm_mp_N_ca, pdpwm_mp_C_ca, pdpwm_mp_5_base, pdpwm_mp_3_base, & ! PDPWM
        pdpwm_para_aNC0, pdpwm_para_aSB0, pdpwm_para_pwm_A, pdpwm_para_pwm_C, & ! PDPWM
        pdpwm_local_coef, pdpwm_local_shift, & ! PDPWM
        pdpwm_para_pwm_G, pdpwm_para_pwm_T, npdpwm_pro, npdpwm_para, & ! PDPWM
        pro_mp2pdns_r0_cutoff, & ! PDNNNNNNNNNS
        ipdns2pro_mp, pro_mp2pdns_para, pdns_para_r0, & ! PDNNNNNNNNNS
        pdns_mp_N_ca, pdns_mp_C_ca, & ! PDNNNNNNNNNS
        pdns_para_aNC0, pdns_para_aSB0, pdns_coef, & ! PDNNNNNNNNNS
        npdns_pro, npdns_para, & ! PDNNNNNNNNNS
        coef_aicg13_gauss, wid_aicg13_gauss, aicg13_nat, factor_aicg13, & ! AICG2
        coef_aicg14_gauss, wid_aicg14_gauss, aicg14_nat, factor_aicg14, & ! AICG2
        coef_dih_gauss, wid_dih_gauss                                               ! AICG2
    use var_enm, only: inenm
    use var_replica, only: n_replica_all
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! ---------------------------------------------------------------------
    integer, intent(in) :: istep_sim

    ! ---------------------------------------------------------------------
    ! local variables
    integer :: ier
    integer :: iunit, junit, icon, ibd, iba, idih
    integer :: ibp, ist
    integer :: imp1, imp2, imp3, imp4
    integer :: itype1, itype2, itype3, itype4, iConType
    integer :: lun
    integer :: iopen_status
    integer :: nld, inat, i_ninfo_type
    integer :: inat_unit(MXUNIT, MXUNIT)
    integer :: ipost2pre_con(nmp_all*MXMPCON)
    integer :: ipost2pre_rna_bp(nmp_all*MXMPRNABP)
    integer :: ipost2pre_rna_st(nmp_all*MXMPRNABP)
    real(PREC) :: cgo_pro
    character(CARRAY_MXFILE) :: cnat_fname(MXUNIT*MXUNIT)
    character(CARRAY_MSG_ERROR) :: error_message
    logical :: flg_enm

    ! ---------------------------------------------------------------------
    ! read native information from nativeinfo in input file
    ! ---------------------------------------------------------------------

    ! ---------------------------------------------------------------------
    ! i_ninfo_type : shows whether native_info input type used
    call mloop_ninfo_inp(istep_sim, i_ninfo_type, inat_unit, cnat_fname)

    nbd = 0
    nba = 0
    ndih = 0
    ncon = 0
    nmorse = 0
    nrna_bp = 0
    nrna_st = 0
    ndtrna_st = 0
    ndtrna_hb = 0
    npdpwm_pro = 0
    npdpwm_para = 0
    npdns_pro = 0
    npdns_para = 0

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        if (i_ninfo_type == 1) then

            inat = 1
            lun = iopen_lunnum + 1
            nld = index(cnat_fname(inat), ' ')
            open (lun, file=cnat_fname(inat) (1:nld - 1), &
                  status='old', action='read', iostat=iopen_status)
            write (*, '(a16,i3,a3,a)') "open ninfo file(", lun, "): ", trim(cnat_fname(inat))
            if (iopen_status > 0) then
                error_message = 'Error: cannot open '//trim(cnat_fname(inat))
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            call read_nativeinfo(lun, i_ninfo_type, 1, 1)
            close (lun)

        else

            do iunit = 1, nunit_all
                do junit = iunit, nunit_all

                    inat = inat_unit(iunit, junit)
                    if (inat == 0) cycle

                    lun = iopen_lunnum + 1
                    nld = index(cnat_fname(inat), ' ')
                    open (lun, file=cnat_fname(inat) (1:nld - 1), &
                          status='old', action='read', iostat=iopen_status)
                    write (*, '(a16,i3,a3,a)') "open ninfo file(", lun, "): ", trim(cnat_fname(inat))
                    if (iopen_status > 0) then
                        error_message = 'Error: cannot open '//trim(cnat_fname(inat))
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if

                    call read_nativeinfo(lun, i_ninfo_type, iunit, junit)
                    close (lun)

                end do
            end do
        endif

        ! ------------------------------------------------------------
#ifdef MPI_PAR
    end if

    call MPI_Bcast(imp2type, MXMP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! bond
    call MPI_Bcast(nbd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ibd2mp, 2*MXMPBD*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ibd2type, MXMPBD*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(bd_nat, MXMPBD*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_bd, MXMPBD*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_bd, 2*MXMPBD*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! angle
    call MPI_Bcast(nba, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(iba2mp, 3*MXMPBA*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(iba2type, MXMPBA*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ba_nat, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_ba, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_ba, 2*MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! aicg13(L_AICG2)
    call MPI_Bcast(aicg13_nat, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_aicg13, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(wid_aicg13_gauss, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_aicg13_gauss, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! dihedral
    call MPI_Bcast(ndih, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(idih2mp, 4*MXMPDIH*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(idih2type, MXMPDIH*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dih_nat, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dih_sin_nat, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dih_cos_nat, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_dih, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_dih, 2*MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! aicg14(L_AICG2)
    call MPI_Bcast(aicg14_nat, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_aicg14, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(wid_aicg14_gauss, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_aicg14_gauss, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! aicgdih(L_AICG2_PLUS)
    call MPI_Bcast(dih_nat, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_aicg14, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(wid_dih_gauss, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_dih_gauss, MXMPDIH*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! go (LJ1210)
    call MPI_Bcast(ncon, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(icon2mp, 2*MXMPCON*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(icon2type, MXMPCON*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(icon2unit, 2*MXMPCON*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(icon_dummy_mgo, MXMPCON*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(go_nat, MXMPCON*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(go_nat2, MXMPCON*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(factor_go, MXMPCON*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(coef_go, MXMPCON*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

    ! go (morse)
    call MPI_Bcast(nmorse, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (nmorse > 0) then
        call MPI_Bcast(imorse2mp, 2*MXMPMORSE*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(imorse2type, MXMPMORSE*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(imorse2unit, 2*MXMPMORSE*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(imorse_dummy_mgo, MXMPMORSE*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(morse_nat, MXMPMORSE*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(morse_nat2, MXMPMORSE*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(factor_morse, MXMPMORSE*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif

    ! RNA basepair
    call MPI_Bcast(nrna_bp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (nrna_bp > 0) then
        call MPI_Bcast(irna_bp2mp, 2*MXRNABP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(irna_bp2unit, 2*MXRNABP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(nhb_bp, MXRNABP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(irna_bp_dummy_mgo, MXRNABP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(rna_bp_nat, MXRNABP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(rna_bp_nat2, MXRNABP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(factor_rna_bp, MXRNABP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(coef_rna_bp, MXRNABP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif
    ! RNA stack
    call MPI_Bcast(nrna_st, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (nrna_st > 0) then
        call MPI_Bcast(irna_st2mp, 2*MXRNAST, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(irna_st2unit, 2*MXRNAST, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(irna_st_dummy_mgo, MXRNAST, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(rna_st_nat, MXRNAST, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(rna_st_nat2, MXRNAST, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(factor_rna_st, MXRNAST, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(coef_rna_st, MXRNAST, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif

    ! DT_RNA stack
    call MPI_Bcast(ndtrna_st, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ndtrna_st > 0) then
        call MPI_Bcast(idtrna_st2mp, 7*MXDTRNAST, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(idtrna_st2nn, MXDTRNAST, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(dtrna_st_nat, 3*MXDTRNAST, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(coef_dtrna_st, 4*MXDTRNAST*n_replica_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif

    ! DT_RNA hydrogen bond (hbond)
    call MPI_Bcast(ndtrna_hb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ndtrna_hb > 0) then
        call MPI_Bcast(idtrna_hb2mp, 6*MXDTRNAHB, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(dtrna_hb_nat, 6*MXDTRNAHB, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(coef_dtrna_hb, 7*MXDTRNAHB, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif

    call MPI_Bcast(npdpwm_pro, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(npdpwm_para, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (inmisc%force_flag(INTERACT%PRO_DNA_PWM)) then
        call MPI_Bcast(ipdpwm2pro_mp, nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pro_mp2pdpwm_para, MXMPPDPWM*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pro_mp2pdpwm_r0_cutoff, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_mp_N_ca, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_mp_C_ca, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_mp_5_base, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_mp_3_base, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_r0, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_a530, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_aNC0, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_aSB0, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_local_coef, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_local_shift, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_pwm_A, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_pwm_C, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_pwm_G, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdpwm_para_pwm_T, MXMPPDPWM*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif

    call MPI_Bcast(npdns_pro, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(npdns_para, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (inmisc%force_flag(INTERACT%PRO_DNA_NONSPEC)) then
        call MPI_Bcast(ipdns2pro_mp, nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pro_mp2pdns_para, MXMPPDNS*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pro_mp2pdns_r0_cutoff, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdns_mp_N_ca, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdns_mp_C_ca, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdns_para_r0, MXMPPDNS*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdns_para_aNC0, MXMPPDNS*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdns_para_aSB0, MXMPPDNS*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(pdns_coef, MXMPPDNS*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
    endif

#endif

    ! ------------------------------------------------------------
    call util_sort_contact(ipost2pre_con)
    if (inmisc%class_flag(CLASS%RNA)) then
        call util_sort_rna_bp(ipost2pre_rna_bp)
        call util_sort_rna_st(ipost2pre_rna_st)
    endif

    ! ------------------------------------------------------------
    ! calc coef_go and ncon_unit
    if (inmisc%force_flag(INTERACT%ENM)) then
        cgo_pro = inenm%cenm
        flg_enm = .true.
    else
        cgo_pro = inpro%cgo1210
        flg_enm = .false.
    endif

    ncon_unit(1:nunit_all, 1:nunit_all) = 0
    do icon = 1, ncon
        iunit = icon2unit(1, icon)
        junit = icon2unit(2, icon)

        if (.not. inmisc%flg_coef_from_ninfo) then
            iConType = get_icon_type(icon2mp(1, icon), icon2mp(2, icon))

            selectcase (iConType)
            case (CONTYPE%RP_RP)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_P_P
            case (CONTYPE%RP_RB)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_P_S
            case (CONTYPE%RP_RS)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_P_B
            case (CONTYPE%RS_RS)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_S_S
            case (CONTYPE%RS_RB)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_S_B
            case (CONTYPE%RB_RB)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_B_B
            case (CONTYPE%PRO_RP)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_pro_P
            case (CONTYPE%PRO_RS)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_pro_B
            case (CONTYPE%PRO_RB)
                coef_go(icon) = factor_go(icon)*inrna%cgo1210_pro_S
            case (CONTYPE%PRO_PRO)
                coef_go(icon) = factor_go(icon)*inpro%cgo1210
            case (CONTYPE%DNA_DNA)
                write (error_message, *) 'Error: logical defect in mloop_nativeinfo; iConType=', iConType
                call util_error(ERROR%STOP_ALL, error_message)
            case (CONTYPE%PRO_DNA)
                coef_go(icon) = factor_go(icon)*inpro%cgo1210
            case (CONTYPE%PRO_DNA2)
                coef_go(icon) = factor_go(icon)*indna2%cgo1210
            case (CONTYPE%PRO_LIG)
                coef_go(icon) = factor_go(icon)*inpro%cgo1210
            case default
                write (error_message, *) 'Error: logical defect in mloop_nativeinfo; iConType=', iConType
                call util_error(ERROR%STOP_ALL, error_message)
            endselect
        endif

        ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
    end do

    ! ------------------------------------------------------------
    ! calc coef_rna_bp
    do ibp = 1, nrna_bp
        if (.not. inmisc%flg_coef_from_ninfo) then
            if (nhb_bp(ibp) == 2) then
                coef_rna_bp(ibp) = factor_rna_bp(ibp)*inrna%cbp1210_HB2
            elseif (nhb_bp(ibp) > 2) then
                coef_rna_bp(ibp) = factor_rna_bp(ibp)*inrna%cbp1210_HB3
            else
                write (error_message, *) &
                    'Error: invalid number of hydrogen bond in mloop_nativeinfo; ibp=', ibp
                call util_error(ERROR%STOP_ALL, error_message)
            endif
        endif

        iunit = irna_bp2unit(1, ibp)
        junit = irna_bp2unit(2, ibp)
        ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
    enddo

    ! ------------------------------------------------------------
    ! calc coef_rna_st
    do ist = 1, nrna_st
        if (.not. inmisc%flg_coef_from_ninfo) then
            coef_rna_st(ist) = factor_rna_st(ist)*inrna%cst1210
        endif

        !!! RNA stack intaractions are not included into Q-score calculation
        !iunit = irna_st2unit(1, ist)
        !junit = irna_st2unit(2, ist)
        !ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
    enddo

    ! ---------------------------------------------------------------
    ! calc coef_bd, coef_ba, coef_dih
    if (.not. inmisc%flg_coef_from_ninfo) then
        do ibd = 1, nbd
            iunit = imp2unit(ibd2mp(1, ibd))

            !     factor_bd(ibd) = inmisc%factor_local_unit(iunit, iunit) * factor_bd(ibd)
            if (iclass_unit(iunit) == CLASS%PRO) then
                coef_bd(1, ibd) = factor_bd(ibd)*inpro%cbd
                coef_bd(2, ibd) = 0.0

            else if (iclass_unit(iunit) == CLASS%DNA) then
                coef_bd(1, ibd) = factor_bd(ibd)*indna%cbd_dna
                coef_bd(2, ibd) = factor_bd(ibd)*indna%cbd2_dna

            else if (iclass_unit(iunit) == CLASS%RNA) then
                imp1 = ibd2mp(1, ibd)
                imp2 = ibd2mp(2, ibd)
                itype1 = imp2type(imp1)
                itype2 = imp2type(imp2)
                if (itype1 == MPTYPE%RNA_PHOS .AND. itype2 == MPTYPE%RNA_SUGAR) then
                    coef_bd(1, ibd) = factor_bd(ibd)*inrna%cbd_PS
                else if (itype1 == MPTYPE%RNA_SUGAR .AND. itype2 == MPTYPE%RNA_BASE) then
                    if (rna_base_type(imp2) == 'R') then
                        coef_bd(1, ibd) = factor_bd(ibd)*inrna%cbd_SR
                    else if (rna_base_type(imp2) == 'Y') then
                        coef_bd(1, ibd) = factor_bd(ibd)*inrna%cbd_SY
                    else
                        error_message = 'Error: invalid rna_base_type, in mloop_nativeinfo'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                else if (itype1 == MPTYPE%RNA_SUGAR .AND. itype2 == MPTYPE%RNA_PHOS) then
                    coef_bd(1, ibd) = factor_bd(ibd)*inrna%cbd_SP
                else
                    error_message = 'Error: bond in mloop_nativeinfo'
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                coef_bd(2, ibd) = 0.0e0_PREC

            else if (iclass_unit(iunit) == CLASS%LIG) then
                coef_bd(1, ibd) = factor_bd(ibd)*inligand%cbd
                coef_bd(2, ibd) = 0.0e0_PREC

            end if
        end do !ibd

        do iba = 1, nba
            iunit = imp2unit(iba2mp(1, iba))

            !     factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit) * factor_ba(iba)
            if (iclass_unit(iunit) == CLASS%PRO) then
                coef_ba(1, iba) = factor_ba(iba)*inpro%cba
                coef_ba(2, iba) = 0.0

            else if (iclass_unit(iunit) == CLASS%DNA) then
                coef_ba(1, iba) = factor_ba(iba)*indna%cba_dna
                coef_ba(2, iba) = 0.0

            else if (iclass_unit(iunit) == CLASS%RNA) then
                imp1 = iba2mp(1, iba)
                imp2 = iba2mp(2, iba)
                imp3 = iba2mp(3, iba)
                itype1 = imp2type(imp1)
                itype2 = imp2type(imp2)
                itype3 = imp2type(imp3)
                if (itype1 == MPTYPE%RNA_PHOS .AND. &
                    itype2 == MPTYPE%RNA_SUGAR .AND. &
                    itype3 == MPTYPE%RNA_BASE) then
                    if (rna_base_type(imp3) == 'R') then
                        coef_ba(1, iba) = factor_ba(iba)*inrna%cba_PSR
                    else if (rna_base_type(imp3) == 'Y') then
                        coef_ba(1, iba) = factor_ba(iba)*inrna%cba_PSY
                    else
                        error_message = 'Error: invalid rna_base_type, in mloop_nativeinfo'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                else if (itype1 == MPTYPE%RNA_PHOS .AND. &
                         itype2 == MPTYPE%RNA_SUGAR .AND. &
                         itype3 == MPTYPE%RNA_PHOS) then
                    coef_ba(1, iba) = factor_ba(iba)*inrna%cba_PSP
                else if (itype1 == MPTYPE%RNA_BASE .AND. &
                         itype2 == MPTYPE%RNA_SUGAR .AND. &
                         itype3 == MPTYPE%RNA_PHOS) then
                    if (rna_base_type(imp1) == 'R') then
                        coef_ba(1, iba) = factor_ba(iba)*inrna%cba_RSP
                    else if (rna_base_type(imp1) == 'Y') then
                        coef_ba(1, iba) = factor_ba(iba)*inrna%cba_YSP
                    else
                        error_message = 'Error: invalid rna_base_type, in mloop_nativeinfo'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                else if (itype1 == MPTYPE%RNA_SUGAR .AND. &
                         itype2 == MPTYPE%RNA_PHOS .AND. &
                         itype3 == MPTYPE%RNA_SUGAR) then
                    coef_ba(1, iba) = factor_ba(iba)*inrna%cba_SPS
                    !        else if (itype1 == MPTYPE%RNA_BASE  .AND. &
                    !                 itype2 == MPTYPE%RNA_SUGAR .AND. &
                    !                 itype3 == MPTYPE%RNA_BASE  ) then
                    !           coef_ba(1, iba) = factor_ba(iba) * inrna%cba_BSB
                else
                    error_message = 'Error: bond-angle in mloop_nativeinfo'
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                coef_ba(2, iba) = 0.0e0_PREC

            else if (iclass_unit(iunit) == CLASS%LIG) then
                coef_ba(1, iba) = factor_ba(iba)*inligand%cba
                coef_ba(2, iba) = 0.0e0_PREC

            end if
        end do !iba

        do idih = 1, ndih
            iunit = imp2unit(idih2mp(1, idih))

            !     factor_dih(idih) = inmisc%factor_local_unit(iunit, iunit) * factor_dih(idih)
            if (iclass_unit(iunit) == CLASS%PRO) then
                coef_dih(1, idih) = factor_dih(idih)*inpro%cdih_1
                coef_dih(2, idih) = factor_dih(idih)*inpro%cdih_3

            else if (iclass_unit(iunit) == CLASS%DNA) then
                coef_dih(1, idih) = factor_dih(idih)*indna%cdih_1_dna
                coef_dih(2, idih) = factor_dih(idih)*indna%cdih_3_dna

            else if (iclass_unit(iunit) == CLASS%RNA) then
                imp1 = idih2mp(1, idih)
                imp2 = idih2mp(2, idih)
                imp3 = idih2mp(3, idih)
                imp4 = idih2mp(4, idih)
                itype1 = imp2type(imp1)
                itype2 = imp2type(imp2)
                itype3 = imp2type(imp3)
                itype4 = imp2type(imp4)
                if (itype1 == MPTYPE%RNA_PHOS .AND. &
                    itype2 == MPTYPE%RNA_SUGAR .AND. &
                    itype3 == MPTYPE%RNA_PHOS .AND. &
                    itype4 == MPTYPE%RNA_SUGAR) then
                    coef_dih(1, idih) = factor_dih(idih)*inrna%cdih_1_PSPS
                    coef_dih(2, idih) = factor_dih(idih)*inrna%cdih_3_PSPS
                else if (itype1 == MPTYPE%RNA_SUGAR .AND. &
                         itype2 == MPTYPE%RNA_PHOS .AND. &
                         itype3 == MPTYPE%RNA_SUGAR .AND. &
                         itype4 == MPTYPE%RNA_BASE) then
                    if (rna_base_type(imp4) == 'R') then
                        coef_dih(1, idih) = factor_dih(idih)*inrna%cdih_1_SPSR
                        coef_dih(2, idih) = factor_dih(idih)*inrna%cdih_3_SPSR
                    else if (rna_base_type(imp4) == 'Y') then
                        coef_dih(1, idih) = factor_dih(idih)*inrna%cdih_1_SPSY
                        coef_dih(2, idih) = factor_dih(idih)*inrna%cdih_3_SPSY
                    else
                        error_message = 'Error: invalid rna_base_type, in mloop_nativeinfo'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                else if (itype1 == MPTYPE%RNA_SUGAR .AND. &
                         itype2 == MPTYPE%RNA_PHOS .AND. &
                         itype3 == MPTYPE%RNA_SUGAR .AND. &
                         itype4 == MPTYPE%RNA_PHOS) then
                    coef_dih(1, idih) = factor_dih(idih)*inrna%cdih_1_SPSP
                    coef_dih(2, idih) = factor_dih(idih)*inrna%cdih_3_SPSP
                else if (itype1 == MPTYPE%RNA_BASE .AND. &
                         itype2 == MPTYPE%RNA_SUGAR .AND. &
                         itype3 == MPTYPE%RNA_PHOS .AND. &
                         itype4 == MPTYPE%RNA_SUGAR) then
                    if (rna_base_type(imp1) == 'R') then
                        coef_dih(1, idih) = factor_dih(idih)*inrna%cdih_1_RSPS
                        coef_dih(2, idih) = factor_dih(idih)*inrna%cdih_3_RSPS
                    else if (rna_base_type(imp1) == 'Y') then
                        coef_dih(1, idih) = factor_dih(idih)*inrna%cdih_1_YSPS
                        coef_dih(2, idih) = factor_dih(idih)*inrna%cdih_3_YSPS
                    else
                        error_message = 'Error: invalid rna_base_type, in mloop_nativeinfo'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                else
                    error_message = 'Error: dihedral-angle in mloop_nativeinfo'
                    call util_error(ERROR%STOP_ALL, error_message)
                endif

            else if (iclass_unit(iunit) == CLASS%LIG) then
                coef_dih(1, idih) = factor_dih(idih)*inligand%cdih
                coef_dih(2, idih) = 0.0e0_PREC

            end if
        end do !idih
    endif

    ! ---------------------------------------------------------------------
    !Re-making iallcon2unit
    ! ---------------------------------------------------------------------
    if (allocated(iallcon2unit)) then
        deallocate (iallcon2unit)
    endif
    allocate (iallcon2unit(2, ncon + nmorse + nrna_bp), stat=ier)
    if (ier /= 0) then
        write (error_message, *) 'failed in memory allocation at mloop_nativeinfo'
        call util_error(ERROR%STOP_ALL, error_message)
    endif
    if (ncon /= 0) then
        iallcon2unit(1:2, 1:ncon) = icon2unit(1:2, 1:ncon)
    endif
    if (nmorse /= 0) then
        iallcon2unit(1:2, ncon + 1:ncon + nmorse) = imorse2unit(1:2, 1:nmorse)
    endif
    if (nrna_bp /= 0) then
        iallcon2unit(1:2, ncon + nmorse + 1:ncon + nmorse + nrna_bp) = irna_bp2unit(1:2, 1:nrna_bp)
    endif

end subroutine mloop_nativeinfo
