subroutine dump_var_struct_neighbor(lunout)

    use const_maxsize
    use var_neighbor_list, only: ele_list, &
        solv_list, exv_list, bp_at_list, bp_gc_list, bp_mismatch_list, &
        exv_dna_list, exv_dna2_list, exv_wca_list, &
        pro_dna_nonspec_list, pro_dna_pwm_list, &
        exv_ion_list, hyd_ion_list, sasa_list, hp_list
    use var_struct

    integer :: i, j, k
    integer :: ni, nj, nk

    write (lunout, *) ''
    write (lunout, *) '### var_struct neighbor list'

    ! electrostatic DNA
    write (lunout, *) '# neighbor list of electrostatic'

    call dump_neighbor_list(lunout, ele_list)

    ! solvation DNA
    write (lunout, *) '# neighbor list of solvation DNA'
    call dump_neighbor_list(lunout, solv_list)

    ! ----------------------------------------------------------------
    ! neighbor_list
    call dump_neighbor_list(lunout, hp_list)

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

    ! -------------------------------------------------------------------
    ! hydrophobic
    write (lunout, *) '# hydrophobic'
    call dump_neighbor_list(lunout, hp_list)

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

endsubroutine dump_var_struct_neighbor
