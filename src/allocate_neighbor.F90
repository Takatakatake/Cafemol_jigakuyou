! allocate_neighbor
!> @brief Allocate/Deallocate arrays of neighborling list

subroutine allocate_neighbor()

    use const_maxsize
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: allocate_neighbor_list, &
        ele_gb_list, pnl_gb_list, &
        exv_list, &
        bp_at_list, bp_gc_list, bp_mismatch_list, &
        exv_dna_list, exv_dna2_list, exv_wca_list, &
        pro_dna_nonspec_list, pro_dna_pwm_list, &
        exv_ion_list, hyd_ion_list, &
        sasa_list, ele_list, &
        solv_list, hp_list, &
        lele_k, iele2charge_k
    use var_setp, only: inmisc, inele
    use var_replica, only: n_replica_mpi
    use var_struct, only: ncharge, nmp_all, nhp

#ifdef MPI_PAR2
    use mpiconst
#endif

    implicit none

    ! -------------------------------------------------------------------
    ! local variables
    integer :: ier
    integer :: ncharge_mpi
    integer :: n_index
    integer :: neighbor_list_capacity
    logical :: flg_error
    character(CARRAY_MSG_ERROR) :: error_message

    ! -------------------------------------------------------------------
    ncharge_mpi = ncharge
#ifdef MPI_PAR2
    ncharge_mpi = ncharge_l
#endif

    n_index = 2 + inperi%n_mirror_index

    ! check
    flg_error = .false.
    if (allocated(exv_list)) flg_error = .true.
    if (allocated(bp_at_list)) flg_error = .true.
    if (allocated(bp_gc_list)) flg_error = .true.
    if (allocated(bp_mismatch_list)) flg_error = .true.
    if (allocated(exv_dna_list)) flg_error = .true.
    if (allocated(exv_dna2_list)) flg_error = .true.
    if (allocated(exv_wca_list)) flg_error = .true.
    if (allocated(pro_dna_nonspec_list)) flg_error = .true.
    if (allocated(pro_dna_pwm_list)) flg_error = .true.
    if (allocated(exv_ion_list)) flg_error = .true.
    if (allocated(hyd_ion_list)) flg_error = .true.
    if (allocated(sasa_list)) flg_error = .true.
    if (allocated(ele_gb_list)) flg_error = .true.
    if (allocated(pnl_gb_list)) flg_error = .true.
    if (allocated(ele_list)) flg_error = .true.
    if (allocated(lele_k)) flg_error = .true.
    if (allocated(iele2charge_k)) flg_error = .true.
    if (allocated(solv_list)) flg_error = .true.
    if (allocated(hp_list)) flg_error = .true.

    if (flg_error) then
        error_message = 'defect at allocate_neighbor, PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
    endif

    !-----------
    ! allocate
    !-----------
    error_message = 'failed in memory allocation at allocate_neighbor, PROGRAM STOP'

#ifdef MPI_PAR2
#ifdef SHARE_NEIGH_PNL
    neighbor_list_capacity = MXMPNEIGHBOR*nmp_all
#else
    neighbor_list_capacity = MXMPNEIGHBOR*nmp_all/npar_mpi + 1
#endif
#else
    neighbor_list_capacity = MXMPNEIGHBOR*nmp_all
#endif

    call allocate_neighbor_list(exv_list, neighbor_list_capacity, stat=ier)
    if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

    call allocate_neighbor_list(bp_at_list, neighbor_list_capacity, stat=ier)
    if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

    call allocate_neighbor_list(bp_gc_list, neighbor_list_capacity, stat=ier)
    if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

    call allocate_neighbor_list(bp_mismatch_list, neighbor_list_capacity, stat=ier)
    if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

    call allocate_neighbor_list(exv_dna_list, neighbor_list_capacity, stat=ier)
    if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

    call allocate_neighbor_list(exv_dna2_list, neighbor_list_capacity, stat=ier)
    if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

    if (inmisc%force_flag(INTERACT%EXV_WCA)) then
        call allocate_neighbor_list(exv_wca_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    end if

    if (inmisc%force_flag(INTERACT%PRO_DNA_NONSPEC)) then
        call allocate_neighbor_list(pro_dna_nonspec_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    end if

    if (inmisc%force_flag(INTERACT%PRO_DNA_PWM)) then
        call allocate_neighbor_list(pro_dna_pwm_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    end if

    if (inmisc%class_flag(CLASS%ION)) then
        call allocate_neighbor_list(exv_ion_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

        call allocate_neighbor_list(hyd_ion_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    end if

    if (inmisc%force_flag(INTERACT%SASA)) then
        call allocate_neighbor_list(sasa_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    end if

    ! for generalized born
    if (inmisc%force_flag(INTERACT%GB)) then
#if defined(MPI_PAR2) && !defined(SHARE_NEIGH_PNL)
        neighbor_list_capacity = MXMPNEIGHBOR*nmp_all/npar_mpi + 1
#else
        neighbor_list_capacity = MXMPNEIGHBOR*nmp_all
#endif
        call allocate_neighbor_list(ele_gb_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

        call allocate_neighbor_list(pnl_gb_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    endif

    ! for electrostatic interaction
    if (inmisc%force_flag(INTERACT%ELE) .or. inmisc%force_flag(INTERACT%GB)) then
#if defined(MPI_PAR2) && !defined(SHARE_NEIGH)
        neighbor_list_capacity = MXMPELE*ncharge/npar_mpi + 1
#else
        neighbor_list_capacity = MXMPELE*ncharge
#endif
        call allocate_neighbor_list(ele_list, neighbor_list_capacity, 1, ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

        if (inele%i_calc_method == 0) then
        else if (inele%i_calc_method == 1) then
            allocate (lele_k(ncharge, n_replica_mpi))
            if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
            lele_k(:, :) = 0.0e0_PREC

            allocate (iele2charge_k(ncharge, ncharge_mpi, n_replica_mpi))
            if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
            iele2charge_k(:, :, :) = 0
        end if
    endif

    ! for DNA
    if (inmisc%force_flag(INTERACT%DNA)) then
#if defined(MPI_PAR2) && !defined(SHARE_NEIGH_SOLV)
        neighbor_list_capacity = MXMPSOLV*nmp_all/npar_mpi + 1
#else
        neighbor_list_capacity = MXMPSOLV*nmp_all
#endif
        call allocate_neighbor_list(solv_list, neighbor_list_capacity, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    endif

    ! for hydrophobic interaction
    if (inmisc%force_flag(INTERACT%HP)) then
#if defined(MPI_PAR2) && !defined(SHARE_NEIGH_HP)
        neighbor_list_capacity = MXMPHP*nhp/npar_mpi + 1
#else
        neighbor_list_capacity = MXMPHP*nhp
#endif
        call allocate_neighbor_list(hp_list, neighbor_list_capacity, 2, stat=ier)
        if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
    end if

end subroutine allocate_neighbor

!######################################################################################

subroutine deallocate_neighbor

    use var_neighbor_list, only: deallocate_neighbor_list, &
        ele_gb_list, pnl_gb_list, &
        exv_list, &
        bp_at_list, bp_gc_list, bp_mismatch_list, &
        exv_dna_list, exv_dna2_list, exv_wca_list, &
        pro_dna_nonspec_list, pro_dna_pwm_list, &
        exv_ion_list, hyd_ion_list, sasa_list, ele_list, &
        solv_list, hp_list, &
        lele_k, iele2charge_k

    if (allocated(exv_list)) call deallocate_neighbor_list(exv_list)
    if (allocated(bp_at_list)) call deallocate_neighbor_list(bp_at_list)
    if (allocated(bp_gc_list)) call deallocate_neighbor_list(bp_gc_list)
    if (allocated(bp_mismatch_list)) call deallocate_neighbor_list(bp_mismatch_list)
    if (allocated(exv_dna_list)) call deallocate_neighbor_list(exv_dna_list)
    if (allocated(exv_dna2_list)) call deallocate_neighbor_list(exv_dna2_list)
    if (allocated(exv_wca_list)) call deallocate_neighbor_list(exv_wca_list)
    if (allocated(pro_dna_nonspec_list)) call deallocate_neighbor_list(pro_dna_nonspec_list)
    if (allocated(pro_dna_pwm_list)) call deallocate_neighbor_list(pro_dna_pwm_list)
    if (allocated(exv_ion_list)) call deallocate_neighbor_list(exv_ion_list)
    if (allocated(hyd_ion_list)) call deallocate_neighbor_list(hyd_ion_list)
    if (allocated(sasa_list)) call deallocate_neighbor_list(sasa_list)
    if (allocated(ele_gb_list)) call deallocate_neighbor_list(ele_gb_list)
    if (allocated(pnl_gb_list)) call deallocate_neighbor_list(pnl_gb_list)
    if (allocated(ele_list)) call deallocate_neighbor_list(ele_list)
    if (allocated(lele_k)) deallocate (lele_k)
    if (allocated(iele2charge_k)) deallocate (iele2charge_k)
    if (allocated(solv_list)) call deallocate_neighbor_list(solv_list)
    if (allocated(hp_list)) call deallocate_neighbor_list(hp_list)

endsubroutine deallocate_neighbor
