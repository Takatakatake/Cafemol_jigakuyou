! simu_energy_allrep
!> @brief This subroutine is to calculate the replica energy.

!#define MEM_ALLOC

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! ************************************************************************
! subroutine for the energy
! ************************************************************************
subroutine simu_energy_allrep(pnle_unit, &
                              pnlet, &
                              velo_mp, &
                              replica_energy, &
                              flg_replica, &  ! flag for calculating replica_energy
                              tempk)
    use const_physical
    use const_maxsize
    use const_index
    use var_inp, only: outfile
    use var_neighbor_list, only: ele_list
    use var_setp, only: inmisc, inele, inwind
    use var_struct, only: nunit_all, lmp2charge, &
        ncharge, coef_charge
    use var_simu, only: qscore_unit, qscore
    use var_replica, only: n_replica_all, n_replica_mpi, irep2grep, &
        lab2val, rep2lab, flg_rep, get_pair, inrep
    use time
#ifdef MPI_PAR3
    use mpiconst
#endif
    implicit none
! --------------------------------------------------------------------
    real(PREC), intent(out) :: pnle_unit(:, :, :, :)  ! (nunit_all, nunit_all, E_TYPE%MAX, replica)
    real(PREC), intent(out) :: pnlet(:, :)          ! (E_TYPE%MAX,replica)
    real(PREC), intent(in)  :: velo_mp(:, :, :)      ! (3, nmp_real, replica)
    real(PREC), intent(out) :: replica_energy(:, :)
    logical, intent(in)  :: flg_replica         ! flag for calculating replica_energy
    real(PREC), intent(in)  :: tempk

    ! --------------------------------------------------------------------
    ! local variables
    integer :: irep, grep
    integer :: iele, icharge, jcharge
    integer :: label_own, label_opp
    real(PREC) :: ionic_strength, tempk_rep, pullforce
    integer :: ksta, kend
#ifdef MPI_PAR2
#ifdef SHARE_NEIGH
    integer :: klen
    real(PREC) :: coef_ele_save(MXMPELE*ncharge)
#else
    real(PREC) :: coef_ele_save(MXMPELE*ncharge/npar_mpi + 1)
#endif
#else
    real(PREC) :: coef_ele_save(MXMPELE*ncharge)
#endif
    character(CARRAY_MSG_ERROR), parameter :: msg_er_allocate = &
                                              'failed in memory allocation at mloop_simulator, PROGRAM STOP'
    character(CARRAY_MSG_ERROR), parameter :: msg_er_deallocate = &
                                              'failed in memory deallocation at mloop_simulator, PROGRAM STOP'

    ! for energy_replica
    real(PREC) :: new_pnlet(E_TYPE%MAX, n_replica_mpi)
    real(PREC) :: new_pnle_unit(nunit_all, nunit_all, E_TYPE%MAX, n_replica_mpi)

#ifdef _DUMP_COMMON
    integer :: i_dump, j_dump, k_dump, l_dump
    integer :: ni_dump, nj_dump, nk_dump, nl_dump
    integer :: istep = -999
    integer :: lundump
#endif

    interface
        subroutine simu_energy(irep, velo_mp, pnlet, pnle_unit)
            use const_maxsize
            implicit none
            integer, intent(in)  :: irep
            real(PREC), intent(in)  :: velo_mp(:, :)      ! (3, nmp_real)
            real(PREC), intent(out) :: pnlet(:)          ! (E_TYPE%MAX)
            real(PREC), intent(out) :: pnle_unit(:, :, :)  ! (nunit_all, nunit_all, E_TYPE%MAX)
        endsubroutine simu_energy
    endinterface

! ------------------------------------------------------------------------
! zero clear
    pnlet(:, :) = 0.0e0_PREC
    pnle_unit(:, :, :, :) = 0.0e0_PREC

#ifdef _DUMP_COMMON
    lundump = outfile%dump
#endif
#ifdef _DUMP_ENERGY
    ! about energy
    write (lundump, *) '# energy  step=', istep
!  real(PREC),allocatable :: pnlet(:,:)          ! (E_TYPE%MAX, replica)
    ni_dump = ubound(pnlet, 1)
    nj_dump = ubound(pnlet, 2)
    do j_dump = 1, nj_dump
        do i_dump = 1, ni_dump
            write (lundump, *) 'pnlet(,', i_dump, ',', j_dump, '),', pnlet(i_dump, j_dump)
        enddo
    enddo

!  real(PREC),allocatable :: pnle_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
    ni_dump = ubound(pnle_unit, 1)
    nj_dump = ubound(pnle_unit, 2)
    nk_dump = ubound(pnle_unit, 3)
    nl_dump = ubound(pnle_unit, 4)
    do l_dump = 1, nl_dump
        do k_dump = 1, nk_dump
            do j_dump = 1, nj_dump
                do i_dump = 1, ni_dump
                    write (lundump, *) 'pnle_unit(', i_dump, ',', j_dump, ',', k_dump, ',', l_dump, '),', &
                        pnle_unit(i_dump, j_dump, k_dump, l_dump)
                enddo
            enddo
        enddo
    enddo
#endif

#ifdef _DEBUG
    if ((size(velo_mp, 3) /= n_replica_mpi) &
        .OR. (size(pnlet, 2) /= n_replica_mpi) .OR. (size(pnle_unit, 4) /= n_replica_mpi) &
        .OR. (size(qscore, 1) /= n_replica_mpi) .OR. (size(qscore_unit, 3) /= n_replica_mpi)) &
        call util_error(ERROR%STOP_ALL, 'defect at simu_energy')
#endif

    do irep = 1, n_replica_mpi

        call simu_energy(irep, velo_mp(:, :, irep), &
                         pnlet(:, irep), pnle_unit(:, :, :, irep))

    enddo

    if (.not. flg_replica) then
        return                 ! <<====== If not replica simulation, return here!
    endif

    !#############################################################################

    TIME_S(tm_energy_replica)
    do irep = 1, n_replica_mpi

        grep = irep2grep(irep)
        replica_energy(1, grep) = pnlet(E_TYPE%TOTAL, irep)

        label_own = rep2lab(grep)
        label_opp = get_pair(label_own)

        if (label_opp == 0 .OR. label_opp > n_replica_all) then
            ! In this case, exchange does not occur.
            replica_energy(2, grep) = 0.0e0_PREC

        else

            !#########################################
            ! Set replica variables of counterpart
            !#########################################
            if (flg_rep(REPTYPE%TEMP)) then
                tempk_rep = lab2val(label_opp, REPTYPE%TEMP)
            else
                tempk_rep = tempk
            endif

            if (flg_rep(REPTYPE%ION)) then
                ionic_strength = lab2val(label_opp, REPTYPE%ION)
            else
                ionic_strength = inele%ionic_strength
            endif

            if (flg_rep(REPTYPE%WIND)) then
                inwind%iwind(grep) = int(lab2val(label_opp, REPTYPE%WIND))
            endif

            if (flg_rep(REPTYPE%WINZ)) then
                inwind%iwinz(grep) = int(lab2val(label_opp, REPTYPE%WINZ))
            endif

            if (flg_rep(REPTYPE%PULL)) then
                pullforce = lab2val(label_opp, REPTYPE%PULL)
                inmisc%pull_unravel_xyz(:, 1:inrep%n_pull, grep) = pullforce*inrep%pull_direction(:, 1:inrep%n_pull)
            endif

            !#############################################
            ! Re-calculate replica-dependent coefficients
            !#############################################
            if (inmisc%force_flag(INTERACT%ELE) .or. inmisc%force_flag(INTERACT%GB)) then
                call simu_ele_set(grep, tempk_rep, ionic_strength)

                if (inele%i_diele /= 0 .and. inele%i_calc_method == 0) then
                    coef_ele_save(:) = ele_list(irep)%coefs(1, :)
#if defined(MPI_PAR3) && defined(SHARE_NEIGH)
                    klen = (ele_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
                    ksta = 1 + klen*local_rank_mpi
                    kend = min(ksta + klen - 1, ele_list(irep)%num_pairs)
#else
                    ksta = 1
                    kend = ele_list(irep)%num_pairs
#endif
                    do iele = ksta, kend
                        icharge = lmp2charge(ele_list(irep)%pairs(1, iele))
                        jcharge = lmp2charge(ele_list(irep)%pairs(2, iele))
                        ele_list(irep)%coefs(1, iele) = coef_charge(icharge, grep) &
                                                        *coef_charge(jcharge, grep)*inele%coef(grep)
                    enddo
                endif
            endif
            if (inmisc%force_flag(INTERACT%DNA)) then
                call simu_solv_set(grep, ionic_strength)
            endif
            if (inmisc%force_flag(INTERACT%DTRNA)) then
                call simu_set_dtrna(grep, tempk_rep)
            endif

            !#############################################
            ! Calculate energy
            !#############################################
            call simu_energy(irep, velo_mp(:, :, irep), &
                             new_pnlet(:, irep), new_pnle_unit(:, :, :, irep))
            replica_energy(2, grep) = new_pnlet(E_TYPE%TOTAL, irep)

            !#########################################
            ! Set back replica variables
            !#########################################
            if (flg_rep(REPTYPE%TEMP)) then
                tempk_rep = lab2val(label_own, REPTYPE%TEMP)
            else
                tempk_rep = tempk
            endif

            if (flg_rep(REPTYPE%ION)) then
                ionic_strength = lab2val(label_own, REPTYPE%ION)
            else
                ionic_strength = inele%ionic_strength
            endif

            if (flg_rep(REPTYPE%WIND)) then
                inwind%iwind(grep) = int(lab2val(label_own, REPTYPE%WIND))
            endif

            if (flg_rep(REPTYPE%WINZ)) then
                inwind%iwinz(grep) = int(lab2val(label_own, REPTYPE%WINZ))
            endif

            if (flg_rep(REPTYPE%PULL)) then
                pullforce = lab2val(label_own, REPTYPE%PULL)
                inmisc%pull_unravel_xyz(:, 1:inrep%n_pull, grep) = pullforce*inrep%pull_direction(:, 1:inrep%n_pull)
            endif

            !#############################################
            ! Set back original coefficients
            !#############################################
            if (inmisc%force_flag(INTERACT%ELE) .or. inmisc%force_flag(INTERACT%GB)) then
                call simu_ele_set(grep, tempk_rep, ionic_strength)

                if (inele%i_diele /= 0 .and. inele%i_calc_method == 0) then
                    ele_list(irep)%coefs(1, :) = coef_ele_save(:)
                endif
            endif
            if (inmisc%force_flag(INTERACT%DNA)) then
                call simu_solv_set(grep, ionic_strength)
            endif
            if (inmisc%force_flag(INTERACT%DTRNA)) then
                call simu_set_dtrna(grep, tempk_rep)
            endif

        endif

    enddo  ! irep
    TIME_E(tm_energy_replica)

end subroutine simu_energy_allrep
