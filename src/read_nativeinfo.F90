! read_nativeinfo
!> @brief Reading nativeinfo file

! **********************************************************************
! This subroutine is for reading intra coordinates from the
! file defined by lun
! **********************************************************************
subroutine read_nativeinfo(lun, i_ninfo_type, iunit, junit)

    use const_maxsize
    use const_physical
    use const_index
    use var_setp, only: inmisc, indna, inpdpwm, inpdns
    use var_struct, only: nmp_all, lunit2mp, imp2unit, iclass_unit, imp2type, &
        nbd, ibd2mp, bd_nat, factor_bd, coef_bd, &
        nba, iba2mp, ba_nat, factor_ba, coef_ba, &
        ndih, idih2mp, dih_nat, dih_sin_nat, dih_cos_nat, factor_dih, coef_dih, &
        ncon, icon2mp, icon2unit, go_nat, go_nat2, factor_go, icon_dummy_mgo, coef_go, &
        nrna_bp, irna_bp2mp, irna_bp2unit, rna_bp_nat, rna_bp_nat2, &
        factor_rna_bp, irna_bp_dummy_mgo, nhb_bp, coef_rna_bp, &
        nrna_st, irna_st2mp, irna_st2unit, rna_st_nat, rna_st_nat2, &
        factor_rna_st, irna_st_dummy_mgo, coef_rna_st, &
        ibd2type, iba2type, idih2type, icon2type, &
        coef_aicg13_gauss, wid_aicg13_gauss, aicg13_nat, factor_aicg13, & ! AICG2
        coef_aicg14_gauss, wid_aicg14_gauss, aicg14_nat, factor_aicg14, & ! AICG2
        coef_dih_gauss, wid_dih_gauss, pro_mp2pdpwm_r0_cutoff, & ! AICG2
        ipdpwm2pro_mp, pro_mp2pdpwm_para, pdpwm_para_r0, pdpwm_para_a530, & ! PDPWM
        pdpwm_para_aNC0, pdpwm_para_aSB0, pdpwm_para_pwm_A, pdpwm_para_pwm_C, & ! PDPWM
        pdpwm_para_pwm_G, pdpwm_para_pwm_T, npdpwm_pro, npdpwm_para, & ! PDPWM
        pdpwm_mp_N_ca, pdpwm_mp_C_ca, pdpwm_mp_5_base, pdpwm_mp_3_base, & ! PDPWM
        pdpwm_local_coef, pdpwm_local_shift, & ! PDPWM
        pro_mp2pdns_r0_cutoff, & ! PDNNNNNNNNNS
        ipdns2pro_mp, pro_mp2pdns_para, pdns_para_r0, & ! PDNNNNNNNNNS
        pdns_para_aNC0, pdns_para_aSB0, pdns_coef, & ! PDNNNNNNNNNS
        npdns_pro, npdns_para, & ! PDNNNNNNNNNS
        pdns_mp_N_ca, pdns_mp_C_ca, & ! PDNNNNNNNNNS
        idtrna_st2mp, dtrna_st_nat, coef_dtrna_st, ndtrna_st, idtrna_st2nn, &
        idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, ndtrna_hb

    implicit none

    ! ---------------------------------------------------------------------
    integer, intent(in) :: lun, i_ninfo_type, iunit, junit
    ! intent(out) :: nbd, ibd2mp, bd_nat, factor_bd, &
    !     nba, iba2mp, ba_nat, factor_ba, &
    !     ndih, idih2mp, dih_nat, factor_dih, &
    !     ncon, icon2mp, icon2unit, go_nat, factor_go

    ! ---------------------------------------------------------------------
    ! local variables
    integer :: input_status
    integer :: ii, jj, imp1, imp2, imp3, imp4, iunit1, iunit2, imp_tmp
    integer :: imp1un, imp2un, imp3un, imp4un, kunit1, kunit2, iunit_tmp
    integer :: ibd, iba, idih, icon
    integer :: iba_aicg2, idih_aicg2, idih_aicg2p
    integer :: ibd_read, iba_read, idih_read
    integer :: icon_read, idummy
    integer :: ibp, ibp_read, nHB
    integer :: ist, ist_read
    integer :: ibsdist
    integer :: ihbdist, ihb_read
    integer :: itype
    integer :: iread2ba(MXBA), iread2dih(MXDIH)
    integer :: ipdpwm_pro, ipdpwm_para
    integer :: ipdpwm_read, ipdpwm_pro_tmp_id
    integer :: ipdpwm_flag_1
    integer :: ipdpwm_i1, ipdpwm_i2
    integer :: ipdns_pro, ipdns_para
    integer :: ipdns_read, ipdns_pro_tmp_id
    integer :: ipdns_flag_1
    integer :: ipdns_i1, ipdns_i2
    real(PREC) :: bl, ba, dih, go, factor, correct, coef, coef3
    real(PREC) :: aicg13, aicg14, wid ! AICG2
    real(PREC) :: dist, energy0
    real(PREC) :: pdpwm_r0, pdpwm_a53, pdpwm_aNC, pdpwm_aSB, pdpwm_rtmp
    real(PREC) :: pdpwm_pwm_A, pdpwm_pwm_C, pdpwm_pwm_G, pdpwm_pwm_T
    real(PREC) :: pdns_r0, pdns_aNC, pdns_aSB, pdns_rtmp
    real(PREC) :: pdns_coef_local
    character(256) :: cline, cline_head
    character(CARRAY_MSG_ERROR) :: error_message
    character(2) :: ctype2
    character(3) :: ctype3
    character(4) :: ctype4

    integer ::ifunc_nn2id
    ! ---------------------------------------------------------------------

    ibd = nbd
    iba = nba
    idih = ndih
    iba_aicg2 = nba
    idih_aicg2 = ndih
    idih_aicg2p = ndih
    icon = ncon
    ibp = nrna_bp
    ist = nrna_st
    ibsdist = ndtrna_st
    ihbdist = ndtrna_hb
    ii = lunit2mp(1, iunit) - 1
    jj = lunit2mp(1, junit) - 1
    iread2ba(:) = -1
    iread2dih(:) = -1
    ipdpwm_pro = npdpwm_pro
    ipdpwm_para = npdpwm_para
    ipdns_pro = npdns_pro
    ipdns_para = npdns_para

    ! ---------------------------------------------------------------------
    ! reading intra coodinates
    ! ---------------------------------------------------------------------
    do
        read (lun, '(a256)', iostat=input_status) cline
        if (input_status < 0) then
            exit
        else if (input_status > 0) then
            error_message = 'Error: cannot read intra coordinate in read_nativeinfo'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! ------------------------------------------------------------------
        ! read the bond length
        if (cline(1:4) == 'bond') then
            ctype2 = '  '
            read (cline, *, iostat=input_status) &
                cline_head, ibd_read, iunit1, iunit2, &
                imp1, imp2, imp1un, imp2un, &
                bl, factor, correct, coef, ctype2
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            ibd = ibd + 1
            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                ibd2mp(1, ibd) = imp1 + ii
                ibd2mp(2, ibd) = imp2 + ii
                kunit1 = iunit1
            else ! NATIVEINFO%ONE_BY_ONE
                ibd2mp(1, ibd) = imp1un + ii
                ibd2mp(2, ibd) = imp2un + ii
                kunit1 = iunit
            end if
            bd_nat(ibd) = bl
            factor_bd(ibd) = factor
            if (inmisc%flg_coef_from_ninfo) then
                coef_bd(1, ibd) = coef
                if (iclass_unit(kunit1) == CLASS%DNA) then
                    coef_bd(2, ibd) = factor_bd(ibd)*indna%cbd2_dna
                endif
                if (iclass_unit(kunit1) == CLASS%DNA2) then
                    coef_bd(2, ibd) = 100.0e0_PREC*coef_bd(1, ibd)
                endif
            endif

            if (ctype2 /= '  ') then
                ibd2type(ibd) = str2bondtype(ctype2)
            endif
        end if

        ! ------------------------------------------------------------------
        ! read the bond angle
        if (cline(1:4) == 'angl') then
            ctype3 = '   '
            read (cline, *, iostat=input_status) &
                cline_head, iba_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp1un, imp2un, imp3un, &
                ba, factor, correct, coef, ctype3
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            iba = iba + 1
            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                if (iba_read <= 0 .or. iba_read > MXBA) then
                    write (error_message, *) 'Invalid iba in angl: ', iba, iba_read
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
                iread2ba(iba_read) = iba
                iba2mp(1, iba) = imp1 + ii
                iba2mp(2, iba) = imp2 + ii
                iba2mp(3, iba) = imp3 + ii
            else ! NATIVEINFO%ONE_BY_ONE
                iba2mp(1, iba) = imp1un + ii
                iba2mp(2, iba) = imp2un + ii
                iba2mp(3, iba) = imp3un + ii
            end if
            ba_nat(iba) = (ba/180.0e0_PREC)*F_PI
            factor_ba(iba) = factor
            if (inmisc%flg_coef_from_ninfo) then
                coef_ba(1, iba) = coef
            endif

            if (ctype3 /= '   ') then
                iba2type(iba) = str2angletype(ctype3)
            endif
        end if

        ! ------------------------------------------------------------------
        !  read the dihedral angle
        if (cline(1:4) == 'dihd') then
            ctype4 = '    '
            read (cline, *, iostat=input_status) &
                cline_head, idih_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
                dih, factor, correct, coef, coef3, ctype4
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            idih = idih + 1
            if (ctype4 /= '    ') then
                itype = str2dihtype(ctype4)
                if (itype == DIHTYPE%VOID) then
                    write (error_message, *) 'Warning: Dihedral angle ', ctype4, &
                        ' is not currently supported.'
                    call util_error(ERROR%WARN_ALL, error_message)
                    idih = idih - 1
                    cycle
                endif
                idih2type(idih) = itype
            endif

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                if (idih_read <= 0 .or. idih_read > MXDIH) then
                    write (error_message, *) 'Invalid idih in dihd: ', idih, idih_read
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
                iread2dih(idih_read) = idih
                idih2mp(1, idih) = imp1 + ii
                idih2mp(2, idih) = imp2 + ii
                idih2mp(3, idih) = imp3 + ii
                idih2mp(4, idih) = imp4 + ii
                kunit1 = iunit1
            else ! NATIVEINFO%ONE_BY_ONE
                idih2mp(1, idih) = imp1un + ii
                idih2mp(2, idih) = imp2un + ii
                idih2mp(3, idih) = imp3un + ii
                idih2mp(4, idih) = imp4un + ii
                kunit1 = iunit
            end if
            dih_nat(idih) = (dih/180.0e0_PREC)*F_PI
            dih_sin_nat(idih) = sin(dih_nat(idih))
            dih_cos_nat(idih) = cos(dih_nat(idih))
            factor_dih(idih) = factor
            if (inmisc%flg_coef_from_ninfo) then

                if (iclass_unit(kunit1) == CLASS%DNA2) then
                    coef_dih_gauss(idih) = coef
                    wid_dih_gauss(idih) = coef3
                else if (iclass_unit(kunit1) == CLASS%DNA) then
                    coef_dih(1, idih) = coef
                    coef_dih(2, idih) = coef3
                else
                    coef_dih(1, idih) = coef
                    if (inmisc%i_triple_angle_term == 0) then
                        coef_dih(2, idih) = 0.0e0_PREC
                    elseif (inmisc%i_triple_angle_term == 1) then
                        coef_dih(2, idih) = 0.5e0_PREC*coef
                    elseif (inmisc%i_triple_angle_term == 2) then
                        coef_dih(2, idih) = coef3
                    else
                        error_message = 'Error: invalid value for i_triple_angle_term in read_nativeinfo.F90'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                end if
            endif
        end if

        ! ------------------------------------------------------------------
        ! read contact
        if (cline(1:7) == 'contact') then
            ctype3 = '   '
            read (cline, *, iostat=input_status) &
                cline_head, icon_read, iunit1, iunit2, &
                imp1, imp2, imp1un, imp2un, &
                go, factor, idummy, coef, ctype3
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            icon = icon + 1
            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + jj
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + jj
            end if
            icon2unit(1, icon) = imp2unit(imp1)
            icon2unit(2, icon) = imp2unit(imp2)
            icon2mp(1, icon) = imp1
            icon2mp(2, icon) = imp2
            go_nat(icon) = go
            go_nat2(icon) = go**2
            factor_go(icon) = factor
            icon_dummy_mgo(icon) = idummy
            if (inmisc%flg_coef_from_ninfo) then
                coef_go(icon) = coef
            endif

            if (ctype3 /= '   ') then
                icon2type(icon) = str2gotype(ctype3)
            endif
        end if

        ! ------------------------------------------------------------------
        ! read basepair
        if (cline(1:8) == 'basepair') then
            ctype3 = '   '
            nHB = 0
            read (cline, *, iostat=input_status) &
                cline_head, ibp_read, iunit1, iunit2, &
                imp1, imp2, imp1un, imp2un, &
                dist, factor, idummy, coef, ctype3, nHB
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            ibp = ibp + 1
            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + jj
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + jj
            end if
            irna_bp2unit(1, ibp) = imp2unit(imp1)
            irna_bp2unit(2, ibp) = imp2unit(imp2)
            irna_bp2mp(1, ibp) = imp1
            irna_bp2mp(2, ibp) = imp2
            rna_bp_nat(ibp) = dist
            rna_bp_nat2(ibp) = dist**2
            factor_rna_bp(ibp) = factor
            irna_bp_dummy_mgo(ibp) = idummy
            if (inmisc%flg_coef_from_ninfo) then
                coef_rna_bp(ibp) = coef
            endif

            if (ctype3 /= '   ') then
                if (nHB < 2) then
                    write (error_message, *) 'Error: in read_nativeinfo, invalid nHB', nHB
                    call util_error(ERROR%STOP_ALL, error_message)
                else
                    nhb_bp(ibp) = nHB
                endif
            endif
        end if

        ! ------------------------------------------------------------------
        ! read basestack
        if (cline(1:9) == 'basestack') then
            ctype3 = '   '
            read (cline, *, iostat=input_status) &
                cline_head, ist_read, iunit1, iunit2, &
                imp1, imp2, imp1un, imp2un, &
                dist, factor, idummy, coef, ctype3
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + jj
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + jj
            end if

            ist = ist + 1
            irna_st2unit(1, ist) = imp2unit(imp1)
            irna_st2unit(2, ist) = imp2unit(imp2)
            irna_st2mp(1, ist) = imp1
            irna_st2mp(2, ist) = imp2
            rna_st_nat(ist) = dist
            rna_st_nat2(ist) = dist**2
            factor_rna_st(ist) = factor
            irna_st_dummy_mgo(ist) = idummy
            if (inmisc%flg_coef_from_ninfo) then
                coef_rna_st(ist) = coef
            endif
        end if

        ! ------------------------------------------------------------------
        ! read basestack of DT model
        if (cline(1:7) == 'bs-dist') then
            ctype3 = '   '
            read (cline, *, iostat=input_status) &
                cline_head, ist_read, iunit1, iunit2, &
                imp1, imp2, imp1un, imp2un, &
                energy0, dist, coef, ctype3
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + jj
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + jj
            end if

            if (imp1 > imp2) then
                imp_tmp = imp1
                imp1 = imp2
                imp2 = imp_tmp
            endif

            if (imp2 /= imp1 + 3) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            ibsdist = ibsdist + 1
            !!!! Caution: read this code carefully!
            !!!! The usage of ibsdist is unusual.
            !!!! Below, index is ist_read, but not ibsdist.

            idtrna_st2mp(1, ist_read) = imp1      ! B1
            idtrna_st2mp(2, ist_read) = imp2      ! B1
            idtrna_st2mp(3, ist_read) = imp1 - 2  ! P1
            idtrna_st2mp(4, ist_read) = imp1 - 1  ! S1
            idtrna_st2mp(5, ist_read) = imp1 + 1  ! P2
            idtrna_st2mp(6, ist_read) = imp1 + 2  ! S2
            idtrna_st2mp(7, ist_read) = imp1 + 4  ! P3

            dtrna_st_nat(1, ist_read) = dist
            if (inmisc%flg_coef_from_ninfo) then
                coef_dtrna_st(1, ist_read, :) = coef
                coef_dtrna_st(0, ist_read, :) = energy0
            endif

            ctype2(1:1) = ctype3(1:1)
            ctype2(2:2) = ctype3(3:3)
            idtrna_st2nn(ist_read) = ifunc_nn2id(ctype2)
        end if

        ! ------------------------------------------------------------------
        ! read hydrogen bond of DT model
        if (cline(1:7) == 'hb-dist') then
            read (cline, *, iostat=input_status) &
                cline_head, ihb_read, iunit1, iunit2, &
                imp1, imp2, imp1un, imp2un, &
                energy0, dist, coef
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + jj
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + jj
            end if

            ihbdist = ihbdist + 1
            if (imp1 < imp2) then
                idtrna_hb2mp(1, ihb_read) = imp1
                idtrna_hb2mp(2, ihb_read) = imp2
            else
                idtrna_hb2mp(1, ihb_read) = imp2
                idtrna_hb2mp(2, ihb_read) = imp1
            endif
            dtrna_hb_nat(1, ihb_read) = dist
            if (inmisc%flg_coef_from_ninfo) then
                coef_dtrna_hb(1, ihb_read) = coef
                coef_dtrna_hb(0, ihb_read) = energy0
            endif
        end if

        ! ------------------------------------------------------------------
        ! read the protein-DNA sequence-specific interaction information
        if (cline(1:5) == 'pdpwm') then
            if ((inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PRO_DNA_PWM) .eqv. .false.) .and. &
                (inmisc%flag_nlocal_unit(junit, iunit, INTERACT%PRO_DNA_PWM) .eqv. .false.)) then
                cycle
            end if
            read (cline, *, iostat=input_status) &
                cline_head, ipdpwm_read, iunit1, imp1, imp1un, &
                pdpwm_r0, pdpwm_aSB, pdpwm_a53, pdpwm_aNC, &
                pdpwm_pwm_A, pdpwm_pwm_C, pdpwm_pwm_G, pdpwm_pwm_T, &
                factor, correct
            if (input_status > 0) then
                error_message = 'read error (in pdpwm block) =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            if (ipdpwm_read <= 0 .or. ipdpwm_read > MXPDPWM) then
                write (error_message, *) 'Invalid ipdpwm in ninfo file (PDPWM block): ', ipdpwm_read
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (ipdpwm_read <= ipdpwm_para) then
                cycle
            end if
            ipdpwm_para = ipdpwm_para + 1
            if (ipdpwm_para > MXPDPWM) then
                write (error_message, *) 'Too many PDPWM interactions: ', ipdpwm_para
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                if (iclass_unit(iunit) == CLASS%PRO) then
                    if (iunit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDPWM: ', iunit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdpwm_pro_tmp_id = imp1 + ii
                else if (iclass_unit(junit) == CLASS%PRO) then
                    if (junit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDPWM: ', junit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdpwm_pro_tmp_id = imp1 + jj
                else
                    write (error_message, *) 'At least one chain should be protein in PDPWM: ', cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            else ! NATIVEINFO%ONE_BY_ONE
                if (iclass_unit(iunit) == CLASS%PRO) then
                    if (iunit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDPWM: ', iunit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdpwm_pro_tmp_id = imp1un + ii
                else if (iclass_unit(junit) == CLASS%PRO) then
                    if (junit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDPWM: ', junit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdpwm_pro_tmp_id = imp1un + jj
                else
                    write (error_message, *) 'At least one chain should be protein in PDPWM: ', cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            end if

            ! -------------------- add ipdpwm_pro_tmp_id to PDPWM_pro list --------------------
            ipdpwm_flag_1 = 0
            ipdpwm_i1 = 0
            do ipdpwm_i1 = 1, ipdpwm_pro
                if (ipdpwm2pro_mp(ipdpwm_i1) == ipdpwm_pro_tmp_id) then
                    ipdpwm_flag_1 = 1  ! if ipdpwm_pro_tmp_id is already in the list ipdpwm2pro_mp
                    exit
                end if
                if (ipdpwm2pro_mp(ipdpwm_i1) == -1) then
                    error_message = ' OH NO! Something error (in pdpwm pro local list) =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            end do
            if (ipdpwm_flag_1 == 0) then
                ipdpwm_pro = ipdpwm_pro + 1
                ipdpwm2pro_mp(ipdpwm_pro) = ipdpwm_pro_tmp_id
            end if
            ! -------------------- add pdpwm_para_id to current aa residue --------------------
            if (pro_mp2pdpwm_para(MXMPPDPWM, ipdpwm_pro_tmp_id) /= 0) then
                error_message = ' OH NO! MXMPPDPWM is not large enough! =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            do ipdpwm_i2 = 1, MXMPPDPWM
                if (pro_mp2pdpwm_para(ipdpwm_i2, ipdpwm_pro_tmp_id) == ipdpwm_para) then
                    error_message = ' OH NO! Conflicts in pro_mp2pdpwm_para! =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
                if (pro_mp2pdpwm_para(ipdpwm_i2, ipdpwm_pro_tmp_id) == 0) then
                    pro_mp2pdpwm_para(ipdpwm_i2, ipdpwm_pro_tmp_id) = ipdpwm_para
                    exit
                end if
            end do
            pdpwm_para_r0(ipdpwm_para) = pdpwm_r0
            pdpwm_para_aNC0(ipdpwm_para) = pdpwm_aNC*F_PI/1.800e2_PREC
            pdpwm_para_a530(ipdpwm_para) = pdpwm_a53*F_PI/1.800e2_PREC
            pdpwm_para_aSB0(ipdpwm_para) = pdpwm_aSB*F_PI/1.800e2_PREC
            pdpwm_para_pwm_A(ipdpwm_para) = pdpwm_pwm_A*inpdpwm%energy_unit_pro_dna_pwm
            pdpwm_para_pwm_C(ipdpwm_para) = pdpwm_pwm_C*inpdpwm%energy_unit_pro_dna_pwm
            pdpwm_para_pwm_G(ipdpwm_para) = pdpwm_pwm_G*inpdpwm%energy_unit_pro_dna_pwm
            pdpwm_para_pwm_T(ipdpwm_para) = pdpwm_pwm_T*inpdpwm%energy_unit_pro_dna_pwm
            pdpwm_local_coef(ipdpwm_para) = factor
            pdpwm_local_shift(ipdpwm_para) = correct
        end if

        ! ------------------------------------------------------------------
        ! read the protein-DNA sequence-non-specific interaction information
        if (cline(1:4) == 'pdns') then
            if ((inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PRO_DNA_NONSPEC) .eqv. .false.) .and. &
                (inmisc%flag_nlocal_unit(junit, iunit, INTERACT%PRO_DNA_NONSPEC) .eqv. .false.)) then
                cycle
            end if
            read (cline, *, iostat=input_status) &
                cline_head, ipdns_read, iunit1, imp1, imp1un, &
                pdns_r0, pdns_aNC, pdns_aSB, &
                pdns_coef_local
            if (input_status > 0) then
                error_message = 'read error (in pdns block) =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            if (ipdns_read <= 0 .or. ipdns_read > MXPDNS) then
                write (error_message, *) 'Invalid ipdns in ninfo file (PDNS block): ', ipdns_read
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (ipdns_read <= ipdns_para) then
                cycle
            end if
            ipdns_para = ipdns_para + 1
            if (ipdns_para > MXPDNS) then
                write (error_message, *) 'Too many PDNS interactions: ', ipdns_para
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                if (iclass_unit(iunit) == CLASS%PRO) then
                    if (iunit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDNS: ', iunit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdns_pro_tmp_id = imp1 + ii
                else if (iclass_unit(junit) == CLASS%PRO) then
                    if (junit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDNS: ', junit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdns_pro_tmp_id = imp1 + jj
                else
                    write (error_message, *) 'At least one chain should be protein in PDNS: ', cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            else ! NATIVEINFO%ONE_BY_ONE
                if (iclass_unit(iunit) == CLASS%PRO) then
                    if (iunit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDNS: ', iunit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdns_pro_tmp_id = imp1un + ii
                else if (iclass_unit(junit) == CLASS%PRO) then
                    if (junit /= iunit1) then
                        write (error_message, *) 'Wrong protein chain in PDNS: ', junit, iunit1
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    ipdns_pro_tmp_id = imp1un + jj
                else
                    write (error_message, *) 'At least one chain should be protein in PDNS: ', cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            end if

            ! -------------------- add ipdns_pro_tmp_id to PDNS_pro list --------------------
            ipdns_flag_1 = 0
            ipdns_i1 = 0
            do ipdns_i1 = 1, ipdns_pro
                if (ipdns2pro_mp(ipdns_i1) == ipdns_pro_tmp_id) then
                    ipdns_flag_1 = 1  ! if ipdns_pro_tmp_id is already in the list ipdns2pro_mp
                    exit
                end if
                if (ipdns2pro_mp(ipdns_i1) == -1) then
                    error_message = ' OH NO! Something error (in pdns pro local list) =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            end do
            if (ipdns_flag_1 == 0) then
                ipdns_pro = ipdns_pro + 1
                ipdns2pro_mp(ipdns_pro) = ipdns_pro_tmp_id
            end if
            ! -------------------- add pdns_para_id to current aa residue --------------------
            if (pro_mp2pdns_para(MXMPPDNS, ipdns_pro_tmp_id) /= 0) then
                error_message = ' OH NO! MXMPPDNS is not large enough! =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            do ipdns_i2 = 1, MXMPPDNS
                if (pro_mp2pdns_para(ipdns_i2, ipdns_pro_tmp_id) == ipdns_para) then
                    error_message = ' OH NO! Conflicts in pro_mp2pdns_para! =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
                if (pro_mp2pdns_para(ipdns_i2, ipdns_pro_tmp_id) == 0) then
                    pro_mp2pdns_para(ipdns_i2, ipdns_pro_tmp_id) = ipdns_para
                    exit
                end if
            end do
            pdns_para_r0(ipdns_para) = pdns_r0
            pdns_para_aNC0(ipdns_para) = pdns_aNC*F_PI/1.800e2_PREC
            pdns_para_aSB0(ipdns_para) = pdns_aSB*F_PI/1.800e2_PREC
            pdns_coef(ipdns_para) = pdns_coef_local*inpdns%energy_unit_pro_dna_nonspec
        end if

    end do

    ! ---------------------------------------------------------------------
    ! rereading nativeinfo file
    ! ---------------------------------------------------------------------
    rewind (lun)
    do
        read (lun, '(a256)', iostat=input_status) cline
        if (input_status < 0) then
            exit
        else if (input_status > 0) then
            error_message = 'Error: cannot read intra coordinate in read_nativeinfo'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! ------------------------------------------------------------------
        ! read the aicg13 with L_AICG2 or L_AICG2_PLUS
        if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
            inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then

            if (cline(1:6) == 'aicg13') then
                ctype3 = '   '
                read (cline, *, iostat=input_status) &
                    cline_head, iba_read, iunit1, iunit2, &
                    imp1, imp2, imp3, imp1un, imp2un, imp3un, &
                    aicg13, factor, correct, coef, wid, ctype3
                if (input_status > 0) then
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                end if

                if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                    if (iba_read <= 0 .or. iba_read > MXBA) then
                        write (error_message, *) 'Invalid iba in aicg13', iba_read
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    iba_aicg2 = iread2ba(iba_read)
                    if (iba_aicg2 <= 0 .or. iba_aicg2 > MXBA) then
                        write (error_message, *) 'Invalid iba_aicg2 in aicg13', iba_read, iba_aicg2
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if
                else
                    iba_aicg2 = iba_aicg2 + 1
                end if

                if (iba2mp(1, iba_aicg2) /= imp1 .or. iba2mp(2, iba_aicg2) /= imp2 .or. iba2mp(3, iba_aicg2) /= imp3) then
           write (error_message, *) "Invalid aicg13 parameters in nativeinfo file ", iba_read, iba2mp(1, iba_aicg2), imp1, iba2mp(2, iba_aicg2), imp2, iba2mp(3, iba_aicg2), imp3
                    call util_error(ERROR%STOP_ALL, error_message)
                end if

                aicg13_nat(iba_aicg2) = aicg13
                factor_aicg13(iba_aicg2) = factor
                wid_aicg13_gauss(iba_aicg2) = wid
                if (inmisc%flg_coef_from_ninfo) then
                    coef_aicg13_gauss(iba_aicg2) = coef
                endif
            end if
        end if

        ! ------------------------------------------------------------------
        !  read the aicg14 with L_AICG2
        if (inmisc%force_flag_local(LINTERACT%L_AICG2)) then
        if (cline(1:6) == 'aicg14') then
            ctype4 = '    '
            read (cline, *, iostat=input_status) &
                cline_head, idih_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
                aicg14, factor, correct, coef, wid, ctype4
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                if (idih_read <= 0 .or. idih_read > MXBA) then
                    write (error_message, *) 'Invalid idih in aicg14', idih_read
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
                idih_aicg2 = iread2dih(idih_read)
                if (idih_aicg2 <= 0 .or. idih_aicg2 > MXDIH) then
                    write (error_message, *) 'Invalid idih_aicg2 in aicg14', idih_read, idih_aicg2
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            else
                idih_aicg2 = idih_aicg2 + 1
            end if

        if(idih2mp(1, idih_aicg2) /= imp1 .or. idih2mp(2, idih_aicg2) /= imp2 .or. idih2mp(3, idih_aicg2) /= imp3 .or. idih2mp(4, idih_aicg2) /= imp4) then
           write (error_message, *) "Invalid aicg14 parameters in nativeinfo file ", idih_read, idih2mp(1, idih_aicg2), imp1, idih2mp(2, idih_aicg2), imp2, idih2mp(3, idih_aicg2), imp3, idih2mp(4, idih_aicg2), imp4
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            aicg14_nat(idih_aicg2) = aicg14
            factor_aicg14(idih_aicg2) = factor
            wid_aicg14_gauss(idih_aicg2) = wid
            if (inmisc%flg_coef_from_ninfo) then
                coef_aicg14_gauss(idih_aicg2) = coef
            endif
        end if
        end if

        ! ------------------------------------------------------------------
        !  read the aicg14 with L_AICG2_PLUS
        if (inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
        if (cline(1:7) == 'aicgdih') then
            ctype4 = '    '
            read (cline, *, iostat=input_status) &
                cline_head, idih_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
                dih, factor, correct, coef, wid, ctype4
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                if (idih_read <= 0 .or. idih_read > MXBA) then
                    write (error_message, *) 'Invalid idih in aicgdih', idih_read
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
                idih_aicg2p = iread2dih(idih_read)
                if (idih_aicg2p <= 0 .or. idih_aicg2p > MXDIH) then
                    write (error_message, *) 'Invalid idih_aicg2p in aicgdih', idih_read, idih_aicg2p
                    call util_error(ERROR%STOP_ALL, error_message)
                end if
            else
                idih_aicg2p = idih_aicg2p + 1
            end if

        if(idih2mp(1, idih_aicg2p) /= imp1 .or. idih2mp(2, idih_aicg2p) /= imp2 .or. idih2mp(3, idih_aicg2p) /= imp3 .or. idih2mp(4, idih_aicg2p) /= imp4) then
           write (error_message, *) "Invalid aicgdih parameters in nativeinfo file ", idih_aicg2p, idih2mp(1, idih_aicg2p), imp1, idih2mp(2, idih_aicg2p), imp2, idih2mp(3, idih_aicg2p), imp3, idih2mp(4, idih_aicg2p), imp4
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                kunit1 = iunit1
                kunit2 = iunit2
            else ! NATIVEINFO%ONE_BY_ONE
                kunit1 = iunit
                kunit2 = junit
            end if
            if (inmisc%flag_local_unit(kunit1, kunit2, LINTERACT%L_AICG2_PLUS)) then
                dih_nat(idih_aicg2p) = (dih/180.0e0_PREC)*F_PI
                factor_aicg14(idih_aicg2p) = factor
                wid_dih_gauss(idih_aicg2p) = wid
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dih_gauss(idih_aicg2p) = coef
                endif
            end if
        end if
        end if
    end do

    ! ---------------------------------------------------------------------
    nbd = ibd
    nba = iba
    ndih = idih
    ncon = icon
    nrna_bp = ibp
    nrna_st = ist
    ndtrna_st = ibsdist
    ndtrna_hb = ihbdist
    npdpwm_pro = ipdpwm_pro
    npdpwm_para = ipdpwm_para
    npdns_pro = ipdns_pro
    npdns_para = ipdns_para

    ! ---------------------------------------------------------------------
    ! check the input
    do icon = 1, ncon
        if (icon2mp(1, icon) == 0 .or. icon2mp(2, icon) == 0) then
            error_message = 'Error: at contact in read_nativeinfo'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
    end do

    ! ------------------------------------------------------------
    ! Assign distance cutoff to Protein-DNA sequence-specific intaractions
    ! ------------------------------------------------------------
    if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PRO_DNA_PWM) .or. &
        inmisc%flag_nlocal_unit(junit, iunit, INTERACT%PRO_DNA_PWM)) then
        if ((iclass_unit(iunit) == CLASS%DNA2 .and. iclass_unit(junit) == CLASS%PRO) &
            .or. (iclass_unit(iunit) == CLASS%PRO .and. iclass_unit(junit) == CLASS%DNA2)) then
            do ipdpwm_i1 = 1, npdpwm_pro
                ipdpwm_pro_tmp_id = ipdpwm2pro_mp(ipdpwm_i1)
                pdpwm_rtmp = 0
                do ipdpwm_i2 = 1, MXMPPDPWM
                    imp_tmp = pro_mp2pdpwm_para(ipdpwm_i2, ipdpwm_pro_tmp_id)
                    if (imp_tmp == 0) then
                        if (ipdpwm_i2 == 1) then
                            write (error_message, *) ' BUG! WTF! Error in PDPWM params: aa ', ipdpwm_pro_tmp_id
                            call util_error(ERROR%STOP_ALL, error_message)
                        end if
                        exit
                    end if
                    if (pdpwm_para_r0(imp_tmp) > pdpwm_rtmp) then
                        pdpwm_rtmp = pdpwm_para_r0(imp_tmp)
                    end if
                end do
                pro_mp2pdpwm_r0_cutoff(ipdpwm_pro_tmp_id) = (pdpwm_rtmp + 5.0e0_PREC)**2.0e0_PREC

                ! ============================== assign N-, C- C_alphas ==============================
                ipdpwm_i2 = ipdpwm_pro_tmp_id - 1
                if (ipdpwm_i2 < 1) then
                    pdpwm_mp_N_ca(ipdpwm_pro_tmp_id) = ipdpwm_pro_tmp_id
                else if (imp2unit(ipdpwm_i2) /= imp2unit(ipdpwm_pro_tmp_id)) then
                    pdpwm_mp_N_ca(ipdpwm_pro_tmp_id) = ipdpwm_pro_tmp_id
                else
                    pdpwm_mp_N_ca(ipdpwm_pro_tmp_id) = ipdpwm_i2
                end if
                ipdpwm_i2 = ipdpwm_pro_tmp_id + 1
                if (ipdpwm_i2 > nmp_all) then
                    pdpwm_mp_C_ca(ipdpwm_pro_tmp_id) = ipdpwm_pro_tmp_id
                else if (imp2unit(ipdpwm_i2) /= imp2unit(ipdpwm_pro_tmp_id)) then
                    pdpwm_mp_C_ca(ipdpwm_pro_tmp_id) = ipdpwm_pro_tmp_id
                else
                    pdpwm_mp_C_ca(ipdpwm_pro_tmp_id) = ipdpwm_i2
                end if
            end do

            ! ============================== assign 5'-, 3'- base ==============================
            if (iclass_unit(iunit) == CLASS%DNA2) then
                iunit_tmp = iunit
            else
                iunit_tmp = junit
            end if
            do ipdpwm_i1 = lunit2mp(1, iunit_tmp), lunit2mp(2, iunit_tmp)
                if (imp2type(ipdpwm_i1) /= MPTYPE%DNA2_BASE) cycle
                ipdpwm_i2 = ipdpwm_i1 - 3
                if (ipdpwm_i2 < 1) then
                    pdpwm_mp_5_base(ipdpwm_i1) = ipdpwm_i1
                else if (imp2unit(ipdpwm_i2) /= iunit_tmp) then
                    pdpwm_mp_5_base(ipdpwm_i1) = ipdpwm_i1
                else
                    pdpwm_mp_5_base(ipdpwm_i1) = ipdpwm_i2
                end if
                ipdpwm_i2 = ipdpwm_i1 + 3
                if (ipdpwm_i2 > nmp_all) then
                    pdpwm_mp_3_base(ipdpwm_i1) = ipdpwm_i1
                else if (imp2unit(ipdpwm_i2) /= imp2unit(ipdpwm_i1)) then
                    pdpwm_mp_3_base(ipdpwm_i1) = ipdpwm_i1
                else
                    pdpwm_mp_3_base(ipdpwm_i1) = ipdpwm_i2
                end if
            end do
        else
            error_message = 'Error: What the fucking are you doing?! Wrong pro-DNA interaction in inp!'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
    end if

    ! ------------------------------------------------------------
    ! Assign distance cutoff to Protein-DNA sequence-non-specific intaractions
    ! ------------------------------------------------------------
    if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PRO_DNA_NONSPEC) .or. &
        inmisc%flag_nlocal_unit(junit, iunit, INTERACT%PRO_DNA_NONSPEC)) then
        if ((iclass_unit(iunit) == CLASS%DNA2 .and. iclass_unit(junit) == CLASS%PRO) &
            .or. (iclass_unit(iunit) == CLASS%PRO .and. iclass_unit(junit) == CLASS%DNA2)) then
            do ipdns_i1 = 1, npdns_pro
                ipdns_pro_tmp_id = ipdns2pro_mp(ipdns_i1)
                pdns_rtmp = 0
                do ipdns_i2 = 1, MXMPPDNS
                    imp_tmp = pro_mp2pdns_para(ipdns_i2, ipdns_pro_tmp_id)
                    if (imp_tmp == 0) then
                        if (ipdns_i2 == 1) then
                            write (error_message, *) ' BUG! WTF! Error in PDNS params: aa ', ipdns_pro_tmp_id
                            call util_error(ERROR%STOP_ALL, error_message)
                        end if
                        exit
                    end if
                    if (pdns_para_r0(imp_tmp) > pdns_rtmp) then
                        pdns_rtmp = pdns_para_r0(imp_tmp)
                    end if
                end do
                pro_mp2pdns_r0_cutoff(ipdns_pro_tmp_id) = (pdns_rtmp + 5.0e0_PREC)**2.0e0_PREC

                ! ============================== assign N-, C- C_alphas ==============================
                ipdns_i2 = ipdns_pro_tmp_id - 1
                if (ipdns_i2 < 1) then
                    pdns_mp_N_ca(ipdns_pro_tmp_id) = ipdns_pro_tmp_id
                else if (imp2unit(ipdns_i2) /= imp2unit(ipdns_pro_tmp_id)) then
                    pdns_mp_N_ca(ipdns_pro_tmp_id) = ipdns_pro_tmp_id
                else
                    pdns_mp_N_ca(ipdns_pro_tmp_id) = ipdns_i2
                end if
                ipdns_i2 = ipdns_pro_tmp_id + 1
                if (ipdns_i2 > nmp_all) then
                    pdns_mp_C_ca(ipdns_pro_tmp_id) = ipdns_pro_tmp_id
                else if (imp2unit(ipdns_i2) /= imp2unit(ipdns_pro_tmp_id)) then
                    pdns_mp_C_ca(ipdns_pro_tmp_id) = ipdns_pro_tmp_id
                else
                    pdns_mp_C_ca(ipdns_pro_tmp_id) = ipdns_i2
                end if
            end do
        else
            error_message = 'Error: What the fucking are you doing?! Wrong pro-DNA interaction in inp!'
            call util_error(ERROR%STOP_ALL, error_message)
        end if
    end if

    ! ------------------------------------------------------------------
    ! read again for bs-dihd, hb-angl, and hb-dihd
    ! ------------------------------------------------------------------
    rewind (lun)
    do
        read (lun, '(a256)', iostat=input_status) cline
        if (input_status < 0) then
            exit
        else if (input_status > 0) then
            error_message = 'Error: cannot read intra coordinate in read_nativeinfo'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! ------------------------------------------------------------------
        ! read basestack dihedral of DT model
        if (cline(1:7) == 'bs-dihd') then
            ctype4 = '    '
            read (cline, *, iostat=input_status) &
                cline_head, ist_read, idih_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
                dih, coef, ctype4
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            dih = dih*F_PI/180.0

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + ii
                imp3 = imp3 + ii
                imp4 = imp4 + ii
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + ii
                imp3 = imp3un + ii
                imp4 = imp4un + ii
            end if

            if (ctype4 == 'PSPS') then
                if ((imp1 == idtrna_st2mp(3, ist_read) .and. &
                     imp2 == idtrna_st2mp(4, ist_read) .and. &
                     imp3 == idtrna_st2mp(5, ist_read) .and. &
                     imp4 == idtrna_st2mp(6, ist_read)) .or. &
                    (imp4 == idtrna_st2mp(3, ist_read) .and. &
                     imp3 == idtrna_st2mp(4, ist_read) .and. &
                     imp2 == idtrna_st2mp(5, ist_read) .and. &
                     imp1 == idtrna_st2mp(6, ist_read))) then
                    continue
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif

                dtrna_st_nat(2, ist_read) = dih
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_st(2, ist_read, :) = coef
                endif

            else if (ctype4 == 'SPSP') then
                if ((imp1 == idtrna_st2mp(4, ist_read) .and. &
                     imp2 == idtrna_st2mp(5, ist_read) .and. &
                     imp3 == idtrna_st2mp(6, ist_read) .and. &
                     imp4 == idtrna_st2mp(7, ist_read)) .or. &
                    (imp4 == idtrna_st2mp(4, ist_read) .and. &
                     imp3 == idtrna_st2mp(5, ist_read) .and. &
                     imp2 == idtrna_st2mp(6, ist_read) .and. &
                     imp1 == idtrna_st2mp(7, ist_read))) then
                    continue
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif

                dtrna_st_nat(3, ist_read) = dih
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_st(3, ist_read, :) = coef
                endif

            else
                error_message = 'Error: unknown type of bs-dihd in read_nativeinfo'
                call util_error(ERROR%STOP_ALL, error_message)
            endif
        end if

        ! ------------------------------------------------------------------
        ! read hydrogen-bond angle of DT model
        if (cline(1:7) == 'hb-angl') then
            read (cline, *, iostat=input_status) &
                cline_head, ihb_read, iba_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp1un, imp2un, imp3un, &
                ba, coef
            if (input_status > 0) then
                write (*, *) input_status
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            ba = ba*F_PI/180.0

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + ii
                imp3 = imp3 + ii
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + ii
                imp3 = imp3un + ii
            end if

            if (imp1 > imp3) then
                imp_tmp = imp1
                imp1 = imp3
                imp3 = imp_tmp
            endif

            if (imp2 == idtrna_hb2mp(1, ihb_read) .and. imp3 == idtrna_hb2mp(2, ihb_read)) then
                if (imp1 == idtrna_hb2mp(3, ihb_read)) then
                    continue
                else if (idtrna_hb2mp(3, ihb_read) == 0) then
                    idtrna_hb2mp(3, ihb_read) = imp1
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                dtrna_hb_nat(2, ihb_read) = ba
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_hb(2, ihb_read) = coef
                endif
            else if (imp1 == idtrna_hb2mp(1, ihb_read) .and. imp2 == idtrna_hb2mp(2, ihb_read)) then
                if (imp3 == idtrna_hb2mp(4, ihb_read)) then
                    continue
                else if (idtrna_hb2mp(4, ihb_read) == 0) then
                    idtrna_hb2mp(4, ihb_read) = imp3
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                dtrna_hb_nat(3, ihb_read) = ba
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_hb(3, ihb_read) = coef
                endif
            else
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            endif
        end if

        ! ------------------------------------------------------------------
        ! read hydrogen-bond dihedral of DT model
        if (cline(1:7) == 'hb-dihd') then
            read (cline, *, iostat=input_status) &
                cline_head, ihb_read, idih_read, iunit1, iunit2, &
                imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
                dih, coef
            if (input_status > 0) then
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            dih = dih*F_PI/180.0

            if (i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
                imp1 = imp1 + ii
                imp2 = imp2 + ii
                imp3 = imp3 + ii
                imp4 = imp4 + ii
            else ! NATIVEINFO%ONE_BY_ONE
                imp1 = imp1un + ii
                imp2 = imp2un + ii
                imp3 = imp3un + ii
                imp4 = imp4un + ii
            end if

            if (imp1 > imp4) then
                imp_tmp = imp1
                imp1 = imp4
                imp4 = imp_tmp
                imp_tmp = imp2
                imp2 = imp3
                imp3 = imp_tmp
            endif

            if (imp2 == idtrna_hb2mp(1, ihb_read) .and. imp3 == idtrna_hb2mp(2, ihb_read)) then
                if (imp1 == idtrna_hb2mp(3, ihb_read)) then
                    continue
                else if (idtrna_hb2mp(3, ihb_read) == 0) then
                    idtrna_hb2mp(3, ihb_read) = imp1
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                if (imp4 == idtrna_hb2mp(4, ihb_read)) then
                    continue
                else if (idtrna_hb2mp(4, ihb_read) == 0) then
                    idtrna_hb2mp(4, ihb_read) = imp4
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                dtrna_hb_nat(4, ihb_read) = dih
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_hb(4, ihb_read) = coef
                endif

            else if (imp3 == idtrna_hb2mp(1, ihb_read) .and. imp4 == idtrna_hb2mp(2, ihb_read)) then
                if (imp2 == idtrna_hb2mp(3, ihb_read)) then
                    continue
                else if (idtrna_hb2mp(3, ihb_read) == 0) then
                    idtrna_hb2mp(3, ihb_read) = imp2
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                idtrna_hb2mp(5, ihb_read) = imp1
                dtrna_hb_nat(5, ihb_read) = dih
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_hb(5, ihb_read) = coef
                endif

            else if (imp1 == idtrna_hb2mp(1, ihb_read) .and. imp2 == idtrna_hb2mp(2, ihb_read)) then
                if (imp3 == idtrna_hb2mp(4, ihb_read)) then
                    continue
                else if (idtrna_hb2mp(4, ihb_read) == 0) then
                    idtrna_hb2mp(4, ihb_read) = imp3
                else
                    error_message = 'read error =>'//cline
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                idtrna_hb2mp(6, ihb_read) = imp4
                dtrna_hb_nat(6, ihb_read) = dih
                if (inmisc%flg_coef_from_ninfo) then
                    coef_dtrna_hb(6, ihb_read) = coef
                endif
            else
                error_message = 'read error =>'//cline
                call util_error(ERROR%STOP_ALL, error_message)
            endif
        endif
    enddo

    ! ---------------------------------------------------------------------

!############################################################################
contains
    integer function str2bondtype(c2)
        use const_index
        implicit none
        character(2), intent(in) :: c2
        if (c2 == 'pp') then
            str2bondtype = BDTYPE%PRO
        else if (c2 == 'PS') then
            str2bondtype = BDTYPE%RNA_PS
        else if (c2 == 'SP') then
            str2bondtype = BDTYPE%RNA_SP
        else if (c2 == 'SA' .OR. c2 == 'SG' .OR. c2 == 'SR') then
            str2bondtype = BDTYPE%RNA_SR
        else if (c2 == 'SU' .OR. c2 == 'SC' .OR. c2 == 'SY') then
            str2bondtype = BDTYPE%RNA_SY
        else
            error_message = 'Error: in read_nativeinfo, unknown bondtype'//c2
            call util_error(ERROR%STOP_ALL, error_message)
        endif
    endfunction str2bondtype

    integer function str2angletype(c3)
        use const_index
        implicit none
        character(3), intent(in) :: c3
        if (c3 == 'ppp') then
            str2angletype = BATYPE%PRO
        else if (c3 == 'PSP') then
            str2angletype = BATYPE%RNA_PSP
        else if (c3 == 'SPS') then
            str2angletype = BATYPE%RNA_SPS
        else if (c3 == 'ASP' .OR. c3 == 'GSP' .OR. c3 == 'RSP') then
            str2angletype = BATYPE%RNA_RSP
        else if (c3 == 'USP' .OR. c3 == 'CSP' .OR. c3 == 'YSP') then
            str2angletype = BATYPE%RNA_YSP
        else if (c3 == 'PSA' .OR. c3 == 'PSG' .OR. c3 == 'PSR') then
            str2angletype = BATYPE%RNA_PSR
        else if (c3 == 'PSU' .OR. c3 == 'PSC' .OR. c3 == 'PSY') then
            str2angletype = BATYPE%RNA_PSY
        else
            error_message = 'Error: in read_nativeinfo, unknown angletype'//c3
            call util_error(ERROR%STOP_ALL, error_message)
        endif
    endfunction str2angletype

    integer function str2dihtype(c4)
        use const_index
        implicit none
        character(4), intent(in) :: c4
        if (c4 == 'pppp') then
            str2dihtype = DIHTYPE%PRO
        else if (c4 == 'PSPS') then
            str2dihtype = DIHTYPE%RNA_PSPS
        else if (c4 == 'SPSP') then
            str2dihtype = DIHTYPE%RNA_SPSP
        else if (c4 == 'ASPS' .OR. c4 == 'GSPS' .OR. c4 == 'RSPS') then
            str2dihtype = DIHTYPE%RNA_RSPS
        else if (c4 == 'USPS' .OR. c4 == 'CSPS' .OR. c4 == 'YSPS') then
            str2dihtype = DIHTYPE%RNA_YSPS
        else if (c4 == 'SPSA' .OR. c4 == 'SPSG' .OR. c4 == 'SPSR') then
            str2dihtype = DIHTYPE%RNA_SPSR
        else if (c4 == 'SPSU' .OR. c4 == 'SPSC' .OR. c4 == 'SPSY') then
            str2dihtype = DIHTYPE%RNA_SPSY
        else if (c4 == 'N2P1') then
            str2dihtype = DIHTYPE%DNA_PER1
        else if (c4 == 'N2P2') then
            str2dihtype = DIHTYPE%DNA_PER2
        else if (c4 == 'N2P3') then
            str2dihtype = DIHTYPE%DNA_PER3
        else if (c4 == 'ASSA' .OR. c4 == 'ASSG' .OR. c4 == 'ASSR' .OR. &
                 c4 == 'GSSA' .OR. c4 == 'GSSG' .OR. c4 == 'GSSR' .OR. &
                 c4 == 'RSSA' .OR. c4 == 'RSSG' .OR. c4 == 'RSSR' .OR. &
                 c4 == 'ASSU' .OR. c4 == 'ASSC' .OR. c4 == 'ASSY' .OR. &
                 c4 == 'GSSU' .OR. c4 == 'GSSC' .OR. c4 == 'GSSY' .OR. &
                 c4 == 'RSSU' .OR. c4 == 'RSSC' .OR. c4 == 'RSSY' .OR. &
                 c4 == 'USSA' .OR. c4 == 'USSG' .OR. c4 == 'USSR' .OR. &
                 c4 == 'CSSA' .OR. c4 == 'CSSG' .OR. c4 == 'CSSR' .OR. &
                 c4 == 'YSSA' .OR. c4 == 'YSSG' .OR. c4 == 'YSSR' .OR. &
                 c4 == 'USSU' .OR. c4 == 'USSC' .OR. c4 == 'USSY' .OR. &
                 c4 == 'CSSU' .OR. c4 == 'CSSC' .OR. c4 == 'CSSY' .OR. &
                 c4 == 'YSSU' .OR. c4 == 'YSSC' .OR. c4 == 'YSSY') then
            ! These are "stack dihedral" supported in previous version.
            str2dihtype = DIHTYPE%VOID
        else
            error_message = 'Error: in read_nativeinfo, unknown dihtype'//c4
            call util_error(ERROR%STOP_ALL, error_message)
        endif
    endfunction str2dihtype

    integer function str2gotype(c3)
        use const_index
        implicit none
        character(3), intent(in) :: c3

        if (c3 == 'p-p') then
            str2gotype = CONTYPE%PRO_PRO
        else if (c3 == 'p-P') then
            str2gotype = CONTYPE%PRO_RP
        else if (c3 == 'p-S') then
            str2gotype = CONTYPE%PRO_RS
        else if (c3 == 'p-A' .OR. c3 == 'p-G' .OR. c3 == 'p-R' .OR. &
                 c3 == 'p-U' .OR. c3 == 'p-C' .OR. c3 == 'p-Y') then
            str2gotype = CONTYPE%PRO_RB
        else if (c3 == 'P-P') then
            str2gotype = CONTYPE%RP_RP
        else if (c3 == 'P-S') then
            str2gotype = CONTYPE%RP_RS
        else if (c3 == 'P-A' .OR. c3 == 'P-G' .OR. c3 == 'P-R' .OR. &
                 c3 == 'P-U' .OR. c3 == 'P-C' .OR. c3 == 'P-Y') then
            str2gotype = CONTYPE%RP_RB
        else if (c3 == 'S-S') then
            str2gotype = CONTYPE%RS_RS
        else if (c3 == 'S-A' .OR. c3 == 'S-G' .OR. c3 == 'S-R' .OR. &
                 c3 == 'S-U' .OR. c3 == 'S-C' .OR. c3 == 'S-Y') then
            str2gotype = CONTYPE%RS_RB
        else if (c3 == 'A-A' .OR. c3 == 'A-G' .OR. c3 == 'A-R' .OR. &
                 c3 == 'A-U' .OR. c3 == 'A-C' .OR. c3 == 'A-Y' .OR. &
                 c3 == 'U-G' .OR. c3 == 'U-R' .OR. &
                 c3 == 'U-U' .OR. c3 == 'U-C' .OR. c3 == 'U-Y' .OR. &
                 c3 == 'G-G' .OR. c3 == 'G-R' .OR. &
                 c3 == 'G-C' .OR. c3 == 'G-Y' .OR. &
                 c3 == 'C-R' .OR. c3 == 'C-C' .OR. c3 == 'C-Y' .OR. &
                 c3 == 'R-R' .OR. c3 == 'R-Y' .OR. c3 == 'Y-Y') then
            str2gotype = CONTYPE%RB_RB
        else
            error_message = 'Error: in read_nativeinfo, unknown contact type'//c3
            call util_error(ERROR%STOP_ALL, error_message)
        endif
    endfunction str2gotype

end subroutine read_nativeinfo
