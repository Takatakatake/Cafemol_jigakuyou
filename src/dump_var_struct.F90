subroutine dump_var_struct(lunout)

    use const_maxsize
    use var_neighbor_list, only: ele_list, &
        solv_list, exv_list, bp_at_list, bp_gc_list, bp_mismatch_list, &
        exv_dna_list, exv_dna2_list, exv_wca_list, &
        pro_dna_nonspec_list, pro_dna_pwm_list, &
        exv_ion_list, hyd_ion_list, sasa_list, hp_list
    use var_struct

    integer :: i, j, k
    integer :: ni, nj, nk
    integer, intent(in) :: lunout

    write (lunout, '(a)') ''
    write (lunout, '(a)') '### var_struct'
    write (lunout, *) 'nunit_real,', nunit_real
    write (lunout, *) 'nunit_all,', nunit_all
    write (lunout, *) 'nmp_real,', nmp_real
    write (lunout, *) 'nmp_all,', nmp_all
    write (lunout, *) 'nres,', nres
    do i = 1, MXUNIT
        write (lunout, *) 'lunit2mp(:,', i, '),', lunit2mp(1, i), lunit2mp(2, i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'ires_mp(', i, '),', ires_mp(i)
    enddo
    do i = 1, MXUNIT
        write (lunout, *) 'iclass_unit(', i, '),', iclass_unit(i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'iclass_mp(', i, '),', iclass_mp(i)
    enddo

! xyz_mp_rep
    write (lunout, *) '# xyz_mp_rep'
    if (allocated(xyz_mp_rep)) then
        nj = ubound(xyz_mp_rep, 2)
        nk = ubound(xyz_mp_rep, 3)
        do k = 1, nk
            do j = 1, nj
                write (lunout, *) 'xyz_mp_rep(:,', j, ',', k, '),'
                write (lunout, *) xyz_mp_rep(:, j, k)
            enddo
        enddo
    else
        write (lunout, *) '#xyz_mp_rep(:,:,:) is not allocated.'
    endif

    do i = 1, MXMP
        write (lunout, *) 'cmass_mp(', i, '),', cmass_mp(i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'fric_mp(', i, '),', fric_mp(i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'cmp2seq(', i, '),', cmp2seq(i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'cmp2atom(', i, '),', cmp2atom(i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'imp2type(', i, '),', imp2type(i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'imp2unit(', i, '),', imp2unit(i)
    enddo

! secondary structure
    do i = 1, MXMP
        write (lunout, *) 'istype_mp(', i, '),', istype_mp(i)
    enddo

    ! ----------------------------------------------------------------
    ! protein structure
    ! bond length
    write (lunout, *) ''
    write (lunout, *) '# bond length'
    write (lunout, *) 'nbd,', nbd
    do i = 1, MXMPBD*nmp_all
        write (lunout, *) 'ibd2mp(:,', i, '),', ibd2mp(1, i), ibd2mp(2, i)
    enddo
    do i = 1, MXMPBD*nmp_all
        write (lunout, *) 'bd_nat(', i, '),', bd_nat(i)
    enddo
    do i = 1, MXMPBD*nmp_all
        write (lunout, *) 'factor_bd(', i, '),', factor_bd(i)
    enddo
    do i = 1, MXMPBD*nmp_all
        write (lunout, *) 'coef_bd(:,', i, '),', coef_bd(1, i), coef_bd(2, i)
    enddo
    do i = 1, MXMPBD*nmp_all
        write (lunout, *) 'correct_bd_mgo(', i, '),', correct_bd_mgo(i)
    enddo

    ! bond angle
    write (lunout, *) ''
    write (lunout, *) '# bond angle'
    write (lunout, *) 'nba,', nba
    do i = 1, MXMPBA*nmp_all
        write (lunout, *) 'iba2mp(:,', i, '),', iba2mp(1, i), iba2mp(2, i), iba2mp(3, i)
    enddo
    do i = 1, MXMPBA*nmp_all
        write (lunout, *) 'ba_nat(', i, '),', ba_nat(i)
    enddo
    do i = 1, MXMPBA*nmp_all
        write (lunout, *) 'factor_ba(', i, '),', factor_ba(i)
    enddo
    do i = 1, MXMPBA*nmp_all
        write (lunout, *) 'coef_ba(:,', i, '),', coef_ba(1, i), coef_ba(2, i)
    enddo
    do i = 1, MXMPBA*nmp_all
        write (lunout, *) 'correct_ba_mgo(', i, '),', correct_ba_mgo(i)
    enddo

    ! dihedral angle
    write (lunout, *) ''
    write (lunout, *) '# dihedral angle'
    write (lunout, *) 'ndih,', ndih
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'idih2mp(:,', i, '),', idih2mp(1, i), idih2mp(2, i), idih2mp(3, i), idih2mp(4, i)
    enddo
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'dih_nat(', i, '),', dih_nat(i)
    enddo
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'factor_dih(', i, '),', factor_dih(i)
    enddo
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'coef_dih(:,', i, '),', coef_dih(1, i), coef_dih(2, i)
    enddo
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'dih_sin_nat(', i, '),', dih_sin_nat(i)
    enddo
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'dih_cos_nat(', i, '),', dih_cos_nat(i)
    enddo
    do i = 1, MXMPDIH*nmp_all
        write (lunout, *) 'correct_dih_mgo(', i, '),', correct_dih_mgo(i)
    enddo

    ! go
    write (lunout, *) ''
    write (lunout, *) '# go'
    write (lunout, *) 'ncon,', ncon
    do i = 1, MXCON
        write (lunout, *) 'icon2mp(:,', i, '),', icon2mp(1, i), icon2mp(2, i)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'lmp2con(', i, '),', lmp2con(i)
    enddo
    do i = 1, MXCON
        write (lunout, *) 'icon2unit(:,', i, '),', icon2unit(1, i), icon2unit(2, i)
    enddo
    do i = 1, MXCON
        write (lunout, *) 'icon_dummy_mgo(', i, '),', icon_dummy_mgo(i)
    enddo
    do i = 1, MXCON
        write (lunout, *) 'go_nat(', i, '),', go_nat(i)
    enddo
    do i = 1, MXCON
        write (lunout, *) 'go_nat2(', i, '),', go_nat2(i)
    enddo
    do i = 1, MXCON
        write (lunout, *) 'factor_go(', i, '),', factor_go(i)
    enddo
    do i = 1, MXCON
        write (lunout, *) 'coef_go(', i, '),', coef_go(i)
    enddo

    ! qscore
    write (lunout, *) '# qsocre'
    do i = 1, MXUNIT
        do j = 1, MXUNIT
            write (lunout, *) 'ncon_unit(', i, ',', j, '),', ncon_unit(i, j)
        enddo
    enddo

    ! base stacking DNA
    write (lunout, *) '# base stacking DNA'
    write (lunout, *) 'nstack,', nstack
    do i = 1, MXSTACK
        write (lunout, *) 'istack2mp(:,', i, '),', istack2mp(1, MXSTACK), istack2mp(2, MXSTACK)
    enddo
    do i = 1, MXMP
        write (lunout, *) 'lmp2stack(', i, '),', lmp2stack(i)
    enddo
    do i = 1, MXSTACK
        write (lunout, *) 'stack_nat(', i, '),', stack_nat(i)
    enddo

    ! electrostatic DNA
    write (lunout, *) '# electrostatic'
    write (lunout, *) 'ncharge,', ncharge
    do i = 1, MXCHARGE
        write (lunout, *) 'icharge2mp(', i, '),', icharge2mp(i)
    enddo
    do i = 1, MXCHARGE
        nj = ubound(coef_charge, 2)
        do j = 1, nj
            write (lunout, *) 'coef_charge(', i, ',', j, '),', coef_charge(i, j)
        enddo
    enddo

    call dump_neighbor_list(lunout, ele_list)

    ! solvation DNA
    write (lunout, *) '# solvation DNA'
    call dump_neighbor_list(lunout, solv_list)

    ! ----------------------------------------------------------------
    ! neighbor_list
    write (lunout, *) '# neighbor list for EXV'
    call dump_neighbor_list(lunout, exv_list)

    write (lunout, *) '# neighbor list for BP_AT'
    call dump_neighbor_list(lunout, bp_at_list)

    write (lunout, *) '# neighbor list for BP_GC'
    call dump_neighbor_list(lunout, bp_gc_list)

    write (lunout, *) '# neighbor list for MBP'
    call dump_neighbor_list(lunout, bp_mismatch_list)

    write (lunout, *) '# neighbor list for EXV_DNA'
    call dump_neighbor_list(lunout, exv_dna_list)

    write (lunout, *) '# neighbor list for EXV_DNA2'
    call dump_neighbor_list(lunout, exv_dna2_list)

    write (lunout, *) '# neighbor list for EXV_WCA'
    call dump_neighbor_list(lunout, exv_wca_list)

    write (lunout, *) '# neighbor list for PRO_DNA_NONSPEC'
    call dump_neighbor_list(lunout, pro_dna_nonspec_list)

    write (lunout, *) '# neighbor list for PRO_DNA_PWM'
    call dump_neighbor_list(lunout, pro_dna_pwm_list)

    write (lunout, *) '# neighbor list for EXV_ION'
    call dump_neighbor_list(lunout, exv_ion_list)

    write (lunout, *) '# neighbor list for HYD_ION'
    call dump_neighbor_list(lunout, hyd_ion_list)

    write (lunout, *) '# neighbor list for SASA'
    call dump_neighbor_list(lunout, sasa_list)

    ! hydrophobic and binding site
    write (lunout, *) '# hydrophobic'
    call dump_neighbor_list(lunout, hp_list)

    write (lunout, *) 'nhp,', nhp
    do i = 1, MXHP
        write (lunout, *) 'ihp2mp(', i, '),', ihp2mp(i)
    enddo
    do i = 1, MXUNIT
        write (lunout, *) 'lunit2hp(1:2,', i, '),', lunit2hp(1:2, i)
    enddo
    do i = 1, MXHP
        write (lunout, *) 'ncoor_hp(', i, '),', ncoor_hp(i)
    enddo
    do i = 1, MXHP
        write (lunout, *) 'ncoor_max_hp(', i, '),', ncoor_max_hp(i)
    enddo
    do i = 1, MXHP
        write (lunout, *) 'coef_aa_hp(', i, '),', coef_aa_hp(i)
    enddo
contains

    subroutine dump_neighbor_list(lunout, list)
        use var_neighbor_list, only: neighbor_list_type
        implicit none

        integer, intent(in) :: lunout
        type(neighbor_list_type), allocatable, intent(in) :: list(:)
        integer :: i, j

        if (allocated(list)) then
            do i = 1, ubound(list, 1)
                write (lunout, '("replica(", i0, ")")') i
                write (lunout, '(2x, "num_pairs = ", i0)') list(i)%num_pairs
                do j = 1, list(i)%num_pairs
                    write (lunout, '(2x, "pairs(", i0, ") = ", 3i0)') j, list(i)%pairs(:, j)
                    write (lunout, '(2x, "coefs(", i0, ") = ", 10f10.3)') j, list(i)%coefs(:, j)
                end do
            end do
        else
            write (lunout, *) 'Not allocated'
        end if
    end subroutine

endsubroutine dump_var_struct
