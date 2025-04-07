! simu_neighbor_assign
!> @brief Assignment and sorting of neighborling list for each energy type

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_assign(irep, num_pairs, pairs)

    use if_neighbor
    use const_maxsize
    use const_index
    use const_physical
    use var_inp, only: inperi
    use var_neighbor_list, only: push_pair, clear_list, &
        exv_list, &
        bp_at_list, bp_gc_list, bp_mismatch_list, &
        exv_dna_list, exv_dna2_list, exv_wca_list, &
        pro_dna_nonspec_list, pro_dna_pwm_list, &
        exv_ion_list, hyd_ion_list, sasa_list
    use var_setp, only: inpro, inmisc, inrna, indtrna
    use var_struct, only: nunit_real, iontype_mp, &
        pxyz_mp_rep, cmp2seq, imp2unit, lmp2con, icon2mp, coef_go, lmp2stack, &
        istack2mp, imp2type, iclass_unit, nmp_real, &
        ires_mp, lmp2morse, lmp2rna_bp, lmp2rna_st, imorse2mp, &
        irna_bp2mp, irna_st2mp, coef_morse_a, coef_morse_fD, coef_rna_bp, &
        coef_rna_bp_a, coef_rna_bp_fD, coef_rna_st, coef_rna_st_a, &
        coef_rna_st_fD, cmp2atom, pro_mp2pdpwm_para, pro_mp2pdns_para
    use time
    use mpiconst

    implicit none

    ! -------------------------------------------------------------------
    integer, intent(in) :: irep
    integer :: num_pairs(0:nthreads - 1)
    integer :: pairs(2, MXMPNEIGHBOR*nmp_real/nthreads, 0:nthreads - 1)

    ! -------------------------------------------------------------------
    ! local variables
    integer :: n
    integer :: inum, imp, jmp, kmp, jimp
    integer :: imp1, imp2, imirror
    integer :: iunit, junit
    integer :: icon, istack
    integer :: istart
    integer :: i_exvol, i_dna
    logical :: flag_exv_wca
    real(PREC) :: vx(3)
    integer :: bp_dna2
    logical :: is_non_local, is_wc_bp

    type calc_type
        integer :: GO
        integer :: MORSE
        integer :: EXV
        integer :: DNA
        integer :: DNA2
        integer :: PAIR_RNA
        integer :: STACK_RNA
        integer :: ION_HYD
        integer :: ION_EXV
        integer :: AICG1
        integer :: AICG2
        integer :: SASA
        integer :: EXV_WCA
        integer :: AICG2_PLUS
        integer :: PRO_DNA_PWM
        integer :: PRO_DNA_NONSPEC
        integer :: MAX
    endtype calc_type
    type(calc_type), parameter :: CALC = calc_type(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 16)
    integer :: icalc(CALC%MAX, nunit_real, nunit_real)

    character(CARRAY_MSG_ERROR) :: error_message

    integer :: imp_l

    ! -------------------------------------------------------------------

    call clear_list(exv_list(irep))
    call clear_list(bp_at_list(irep))
    call clear_list(bp_gc_list(irep))
    call clear_list(bp_mismatch_list(irep))
    call clear_list(exv_dna_list(irep))
    call clear_list(exv_dna2_list(irep))

    if (inmisc%force_flag(INTERACT%EXV_WCA)) then
        call clear_list(exv_wca_list(irep))
    end if

    if (inmisc%force_flag(INTERACT%PRO_DNA_NONSPEC)) then
        call clear_list(pro_dna_nonspec_list(irep))
    end if

    if (inmisc%force_flag(INTERACT%PRO_DNA_PWM)) then
        call clear_list(pro_dna_pwm_list(irep))
    end if

    if (inmisc%class_flag(CLASS%ION)) then
        call clear_list(exv_ion_list(irep))
        call clear_list(hyd_ion_list(irep))
    end if

    if (inmisc%force_flag(INTERACT%SASA)) then
        call clear_list(sasa_list(irep))
    end if

    ! --------------------------------------------------------------------
    icalc(1:CALC%MAX, 1:nunit_real, 1:nunit_real) = 0

    do iunit = 1, nunit_real
        do junit = iunit, nunit_real
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%GO)) then
                icalc(CALC%GO, iunit, junit) = 1
            end if
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%MORSE)) then
                icalc(CALC%MORSE, iunit, junit) = 1
            end if
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV)) then
                icalc(CALC%EXV, iunit, junit) = 1
            end if
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA)) then
                icalc(CALC%DNA, iunit, junit) = 1
            end if
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA2) .OR. &
                inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA2C)) then
                icalc(CALC%DNA2, iunit, junit) = 1
            end if
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PAIR_RNA)) then
                icalc(CALC%PAIR_RNA, iunit, junit) = 1
            endif
            if ((iunit == junit) .AND. (iclass_unit(iunit) == CLASS%RNA)) then
                icalc(CALC%STACK_RNA, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ION_HYD)) then
                icalc(CALC%ION_HYD, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ION_EXV)) then
                icalc(CALC%ION_EXV, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG1)) then
                icalc(CALC%AICG1, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2)) then
                icalc(CALC%AICG2, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2_PLUS)) then
                icalc(CALC%AICG2_PLUS, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%SASA)) then
                icalc(CALC%SASA, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PRO_DNA_PWM)) then
                icalc(CALC%PRO_DNA_PWM, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PRO_DNA_NONSPEC)) then
                icalc(CALC%PRO_DNA_NONSPEC, iunit, junit) = 1
            endif
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV_WCA)) then
                icalc(CALC%EXV_WCA, iunit, junit) = 1
            endif
        end do
    end do

    ! --------------------------------------------------------------------
    do n = 0, nthreads - 1
        do imp_l = 1, num_pairs(n)
            imp = pairs(1, imp_l, n)
            jmp = pairs(2, imp_l, n)

            iunit = imp2unit(imp)
            junit = imp2unit(jmp)

            i_exvol = 1

            ! -----------------------------------------------------------------
            ! go
            if (icalc(CALC%GO, iunit, junit) == 1 .OR. &
                icalc(CALC%AICG1, iunit, junit) == 1 .OR. &
                icalc(CALC%AICG2, iunit, junit) == 1 .OR. &
                icalc(CALC%AICG2_PLUS, iunit, junit) == 1) then
                if (imp == 1) then
                    istart = 1
                else
                    istart = lmp2con(imp - 1) + 1
                end if
                do icon = istart, lmp2con(imp)
                    kmp = icon2mp(2, icon)

                    if (jmp > kmp) cycle
                    if (jmp < kmp) exit

                    if (coef_go(icon) > ZERO_JUDGE) then
                        i_exvol = 0
                    end if
                end do
            end if

            ! -----------------------------------------------------------------
            ! morse
            if (icalc(CALC%MORSE, iunit, junit) == 1) then
                if (imp == 1) then
                    istart = 1
                else
                    istart = lmp2morse(imp - 1) + 1
                end if
                do icon = istart, lmp2morse(imp)
                    kmp = imorse2mp(2, icon)

                    if (jmp > kmp) cycle
                    if (jmp < kmp) exit

                    if (coef_morse_a(icon) > ZERO_JUDGE .and. coef_morse_fD(icon) > ZERO_JUDGE) then
                        i_exvol = 0
                    end if
                end do
            end if

            ! -----------------------------------------------------------------
            ! RNA base pairing
            if (icalc(CALC%PAIR_RNA, iunit, junit) == 1) then
                if (imp == 1) then
                    istart = 1
                else
                    istart = lmp2rna_bp(imp - 1) + 1
                end if
                do icon = istart, lmp2rna_bp(imp)
                    kmp = irna_bp2mp(2, icon)

                    if (jmp > kmp) cycle
                    if (jmp < kmp) exit

                    ! if (coef_rna_bp(icon) > ZERO_JUDGE .or. &
                    !     (coef_rna_bp_a(icon) > ZERO_JUDGE .and. coef_rna_bp_fD(icon) > ZERO_JUDGE)) then
                    !     i_exvol = 0
                    ! end if
                    i_exvol = 0
                end do
            end if

            ! -----------------------------------------------------------------
            ! RNA base stacking
            if (icalc(CALC%STACK_RNA, iunit, junit) == 1) then
                if (imp == 1) then
                    istart = 1
                else
                    istart = lmp2rna_st(imp - 1) + 1
                end if
                do icon = istart, lmp2rna_st(imp)
                    kmp = irna_st2mp(2, icon)

                    if (jmp > kmp) cycle
                    if (jmp < kmp) exit

                    if (coef_rna_st(icon) > ZERO_JUDGE .or. &
                        (coef_rna_st_a(icon) > ZERO_JUDGE .and. coef_rna_st_fD(icon) > ZERO_JUDGE)) then
                        i_exvol = 0
                    end if
                end do
            end if

            ! -----------------------------------------------------------------
            ! exvol
            if (icalc(CALC%EXV, iunit, junit) == 1) then
                if (iunit == junit) then
                    if (iclass_unit(iunit) == CLASS%RNA) then
                        select case (imp2type(imp))
                        case (MPTYPE%RNA_PHOS)
                            if (jmp < imp + inrna%n_sep_nlocal_P) then
                                i_exvol = 0
                            end if
                        case (MPTYPE%RNA_SUGAR)
                            if (jmp < imp + inrna%n_sep_nlocal_S) then
                                i_exvol = 0
                            end if
                        case (MPTYPE%RNA_BASE)
                            if (jmp < imp + inrna%n_sep_nlocal_B) then
                                i_exvol = 0
                            end if
                        case default
                            error_message = 'Error: logical defect in simu_neighbor_assign'
                            call util_error(ERROR%STOP_ALL, error_message)
                        endselect

                    else if (iclass_unit(iunit) == CLASS%LIG) then
                        if (ires_mp(imp) == ires_mp(jmp)) then
                            i_exvol = 0
                        end if

                    else
                        if (jmp < imp + inpro%n_sep_nlocal) then
                            i_exvol = 0
                        end if

                    endif
                end if

                if (i_exvol == 1) then
                    call push_pair(exv_list(irep), imp, jmp)
                end if
            end if

            ! -----------------------------------------------------------------
            ! DNA2-DNA2 (for 3SPN.2 model)
            if (icalc(CALC%DNA2, iunit, junit) == 1) then
                is_non_local = .true.

                if (iunit == junit) then
                    is_non_local = abs(jmp - imp) > 4
                end if

                is_wc_bp = .false.

                if (is_non_local) then
                    ! check if the particle pair is W.C. base pair
                    if (imp2type(imp) == MPTYPE%DNA2_BASE .and. imp2type(jmp) == MPTYPE%DNA2_BASE) then
                        call seq2bp(cmp2seq(imp), cmp2seq(jmp), bp_dna2, is_wc_bp)
                    end if

                    if (.not. is_wc_bp) then
                        call push_pair(exv_dna2_list(irep), imp, jmp)
                    end if
                end if
            end if

            ! -----------------------------------------------------------------
            ! excluded volume (DTRNA)
            if (icalc(CALC%EXV_WCA, iunit, junit) == 1) then
                flag_exv_wca = .true.
                if (iunit == junit) then
                    if (iclass_unit(iunit) == CLASS%RNA) then
                        select case (imp2type(imp))
                        case (MPTYPE%RNA_PHOS)
                            flag_exv_wca = (jmp >= imp + indtrna%n_sep_nlocal_P)
                        case (MPTYPE%RNA_SUGAR)
                            flag_exv_wca = (jmp >= imp + indtrna%n_sep_nlocal_S)
                        case (MPTYPE%RNA_BASE)
                            flag_exv_wca = (jmp >= imp + indtrna%n_sep_nlocal_B)
                        case default
                            error_message = 'Error: logical defect in simu_neighbor_assign'
                            call util_error(ERROR%STOP_ALL, error_message)
                        endselect
                    else
                        error_message = 'Error: EXV_WCA is only for RNA currently, in simu_neighbor_assign'
                        call util_error(ERROR%STOP_ALL, error_message)
                    endif
                end if

                if (flag_exv_wca) then
                    call push_pair(exv_wca_list(irep), imp, jmp)
                end if
            end if

            ! -----------------------------------------------------------------
            ! DNA-DNA
            i_dna = 1
            if (icalc(CALC%DNA, iunit, junit) == 1) then
                if (iunit == junit) then
                    if (imp == 1) then
                        istart = 1
                    else
                        istart = lmp2stack(imp - 1) + 1
                    end if
                    do istack = istart, lmp2stack(imp)
                        kmp = istack2mp(2, istack)

                        if (jmp > kmp) cycle
                        if (jmp == kmp) then
                            i_dna = 0
                        end if
                        exit
                    end do
                end if

                if (i_dna == 1) then
                    if (imp2type(imp) == MPTYPE%DNA_BASE .and. imp2type(jmp) == MPTYPE%DNA_BASE) then
                        if ((cmp2seq(imp) == ' DA' .and. cmp2seq(jmp) == ' DT') .or. &
                            (cmp2seq(imp) == ' DT' .and. cmp2seq(jmp) == ' DA')) then
                            call push_pair(bp_at_list(irep), imp, jmp)
                            i_dna = 0

                        else if ((cmp2seq(imp) == ' DG' .and. cmp2seq(jmp) == ' DC') .or. &
                                 (cmp2seq(imp) == ' DC' .and. cmp2seq(jmp) == ' DG')) then
                            call push_pair(bp_gc_list(irep), imp, jmp)
                            i_dna = 0

                        else
                            call push_pair(bp_mismatch_list(irep), imp, jmp)
                            i_dna = 0
                        end if
                    end if

                    if (iunit == junit) then
                        jimp = jmp - imp
                        if (jimp <= 2) then
                            if (imp2type(imp) == MPTYPE%DNA_SUGAR) then
                                i_dna = 0
                            else if (imp2type(imp) == MPTYPE%DNA_BASE) then
                                if (jimp <= 0) then
                                    i_dna = 0
                                end if
                            else
                                if (jimp <= 1) then
                                    i_dna = 0
                                end if
                            end if
                        end if
                    end if

                    ! skip phosphate-phosphate interaction when using explicit ion model
                    if (icalc(CALC%ION_HYD, iunit, junit) == 1 .and. iontype_mp(imp) /= 0 .and. iontype_mp(imp) /= 0) then
                        i_dna = 0
                    end if

                    if (i_dna == 1) then
                        call push_pair(exv_dna_list(irep), imp, jmp)
                    end if
                end if
            end if

            ! -----------------------------------------------------------------
            ! Ion
            if (icalc(CALC%ION_HYD, iunit, junit) == 1 .or. icalc(CALC%ION_EXV, iunit, junit) == 1) then
                if (iontype_mp(imp) /= 0 .and. iontype_mp(jmp) /= 0 .and. icalc(CALC%ION_HYD, iunit, junit) == 1) then
                    call push_pair(hyd_ion_list(irep), imp, jmp)
                else if (iontype_mp(imp) /= 0 .or. iontype_mp(jmp) /= 0) then
                    if (iontype_mp(imp) /= 0 .or. iontype_mp(jmp) /= IONTYPE%P) then
                        if (iontype_mp(imp) /= IONTYPE%P .or. iontype_mp(jmp) /= 0) then
                            call push_pair(exv_ion_list(irep), imp, jmp)
                        end if
                    end if
                end if
            end if

            ! -----------------------------------------------------------------
            ! SASA
            if (icalc(CALC%SASA, iunit, junit) == 1) then
                if (iclass_unit(iunit) == CLASS%PRO .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ') then
                    if (iclass_unit(junit) == CLASS%PRO .or. cmp2atom(jmp) == ' P  ' .or. cmp2atom(jmp) == ' O  ') then
                        call push_pair(sasa_list(irep), imp, jmp); 
                    end if
                end if
            end if
            ! ------------------------------------------------------------------

            ! -----------------------------------------------------------------
            ! PDPWM
            if (icalc(CALC%PRO_DNA_PWM, iunit, junit) == 1) then
                ! ========== sequence: (aa_index, dna_index, E_TYPE%PRO_DNA_PWM) ==========
                if (imp2type(imp) == MPTYPE%DNA2_BASE .and. pro_mp2pdpwm_para(1, jmp) > 0) then
                    call push_pair(pro_dna_pwm_list(irep), jmp, imp)
                else if (imp2type(jmp) == MPTYPE%DNA2_BASE .and. pro_mp2pdpwm_para(1, imp) > 0) then
                    call push_pair(pro_dna_pwm_list(irep), imp, jmp)
                end if
            end if
            ! -----------------------------------------------------------------

            ! -----------------------------------------------------------------
            ! PDNS
            if (icalc(CALC%PRO_DNA_NONSPEC, iunit, junit) == 1) then
                ! ========== sequence: (aa_index, dna_index, E_TYPE%PRO_DNA_NONSPEC) ==========
                if (imp2type(imp) == MPTYPE%DNA2_PHOS .and. pro_mp2pdns_para(1, jmp) > 0) then
                    call push_pair(pro_dna_nonspec_list(irep), jmp, imp)
                else if (imp2type(jmp) == MPTYPE%DNA2_PHOS .and. pro_mp2pdns_para(1, imp) > 0) then
                    call push_pair(pro_dna_nonspec_list(irep), imp, jmp)
                end if
            end if
            ! -----------------------------------------------------------------
        end do
    end do

    ! --------------------------------------------------------------------

    if (inperi%i_periodic == 1) then
        call apply_pbc(exv_list, irep)
        call apply_pbc(bp_at_list, irep)
        call apply_pbc(bp_gc_list, irep)
        call apply_pbc(bp_mismatch_list, irep)
        call apply_pbc(exv_dna_list, irep)
        call apply_pbc(exv_dna2_list, irep)

        if (inmisc%force_flag(INTERACT%EXV_WCA)) then
            call apply_pbc(exv_wca_list, irep)
        end if

        if (inmisc%force_flag(INTERACT%PRO_DNA_NONSPEC)) then
            call apply_pbc(pro_dna_nonspec_list, irep)
        end if

        if (inmisc%force_flag(INTERACT%PRO_DNA_PWM)) then
            call apply_pbc(pro_dna_pwm_list, irep)
        end if

        if (inmisc%class_flag(CLASS%ION)) then
            call apply_pbc(exv_ion_list, irep)
            call apply_pbc(hyd_ion_list, irep)
        end if

        if (inmisc%force_flag(INTERACT%SASA)) then
            call apply_pbc(sasa_list, irep)
        end if
    end if

contains

    subroutine apply_pbc(list, irep)
        use var_neighbor_list, only: neighbor_list_type
        implicit none
        type(neighbor_list_type), allocatable, intent(inout) :: list(:)
        integer, intent(in) :: irep
        integer :: i, imp, jmp, imirror
        real(PREC) :: vec(3)

        do i = 1, list(irep)%num_pairs
            imp = list(irep)%pairs(1, i)
            jmp = list(irep)%pairs(2, i)

            vec(:) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
            call util_pbneighbor(vec, imirror)

            list(irep)%pairs(3, i) = imirror
        end do
    end subroutine

end subroutine simu_neighbor_assign
