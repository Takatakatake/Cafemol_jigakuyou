! simu_neighbor_list_hp
!> @brief This subroutine is to make neighborlist especially for the hydrophobic interaction.

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
! neighbor list for hydrohpobic interaction
subroutine simu_neighbor_list_hp(irep)

    use const_maxsize
    use const_index
    use var_neighbor_list, only: hp_list, clear_list, push_pair
    use var_setp, only: inmisc, inhp
    use var_struct, only: nunit_real, xyz_mp_rep, &
        imp2unit, nhp, ihp2mp, lunit2hp, cmp2seq, iclass_mp
    use time
#ifdef MPI_PAR2
    use mpiconst
    use var_neighbor_list, only: neighbor_list_type, &
        reserve_neighbor_list, dealloc_neighbor_list
#endif

    implicit none

    integer, intent(in) :: irep
    ! -------------------------------------------------------------------
    ! local variables
    integer :: imp, jmp, iunit, junit
    integer :: ihp, jhp, itype1, itype2
    integer :: icalc(MXUNIT, MXUNIT)
    real(PREC) :: dist2, rneighbor2_hp, v21(3)
    integer :: ksta, kend

#ifdef MPI_PAR2
    integer :: klen
#ifdef SHARE_NEIGH_HP
    type(neighbor_list_type) :: local_hp_list

    integer :: n_index
    integer :: n
    integer :: num_all_pairs(0:npar_mpi - 1)
    integer :: disp(0:npar_mpi - 1)
#endif
#endif

    integer :: ifunc_seq2id

    ! -------------------------------------------------------------------
    if (.not. inmisc%force_flag(INTERACT%HP)) then
        return
    end if

    ! -------------------------------------------------------------------
    icalc(1:nunit_real, 1:nunit_real) = 0
    do iunit = 1, nunit_real
        do junit = iunit, nunit_real
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%HP)) then
                icalc(iunit, junit) = 1
                icalc(junit, iunit) = 1
            end if
        end do
    end do

    ! --------------------------------------------------------------------

    call clear_list(hp_list(irep))

#ifdef MPI_PAR2
#ifdef SHARE_NEIGH_HP
    call reserve_neighbor_list(local_hp_list, hp_list(irep), ierr)
#endif
    klen = (nhp - 1 + npar_mpi)/npar_mpi
    ksta = 1 + klen*local_rank_mpi
    kend = min(ksta + klen - 1, nhp)
#else
    ksta = 1
    kend = nhp
#endif

    do ihp = ksta, kend
        imp = ihp2mp(ihp)
        iunit = imp2unit(imp)

        jhp = 1

        do while (jhp <= nhp)
            if (jhp == ihp) then
                jhp = jhp + 1
                cycle
            endif

            jmp = ihp2mp(jhp)
            junit = imp2unit(jmp)

            if (icalc(iunit, junit) == 1) then
                if (iclass_mp(imp) /= CLASS%PRO) then
                    itype1 = 21
                else
                    itype1 = ifunc_seq2id(cmp2seq(imp))
                end if
                if (iclass_mp(jmp) /= CLASS%PRO) then
                    itype2 = 21
                else
                    itype2 = ifunc_seq2id(cmp2seq(jmp))
                end if
                rneighbor2_hp = inhp%cutoffdmax_para_hp(itype1, itype2)* &
                                inhp%cutoffdmax_para_hp(itype1, itype2)

                v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
                dist2 = dot_product(v21, v21)

                if (dist2 >= rneighbor2_hp) then
                    jhp = jhp + 1
                    cycle
                end if
#if defined(MPI_PAR2) && defined(SHARE_NEIGH_HP)
                call push_pair(local_hp_list, ihp, jhp)
                local_hp_list%coefs(1:2, local_hp_list%num_pairs) = &
                    [inhp%cutoffdmin_para_hp(itype1, itype2), &
                     inhp%cutoffdmax_para_hp(itype1, itype2)]
#else
                call push_pair(hp_list(irep), ihp, jhp)
                hp_list(irep)%coefs(1:2, hp_list(irep)%num_pairs) = &
                    [inhp%cutoffdmin_para_hp(itype1, itype2), &
                     inhp%cutoffdmax_para_hp(itype1, itype2)]
#endif
            else
                jhp = lunit2hp(2, junit)
            end if

            jhp = jhp + 1
        end do
    end do

#if defined(MPI_PAR2) && defined(SHARE_NEIGH_HP)
    TIME_S(tmc_neighbor)

    call mpi_allgather(local_hp_list%num_pairs, 1, MPI_INTEGER, &
                       num_all_pairs, 1, MPI_INTEGER, &
                       mpi_comm_local, ierr)

    hp_list(irep)%num_pairs = sum(num_all_pairs(0:npar_mpi - 1))

    disp(0) = 0
    do n = 1, npar_mpi - 1
        disp(n) = disp(n - 1) + num_all_pairs(n - 1)
    end do

    n_index = 2 + inperi%n_mirror_index
    call mpi_allgatherv(local_hp_list%pairs, n_index*local_hp_list%num_pairs, MPI_INTEGER, &
                        hp_list(irep)%pairs, n_index*num_all_pairs, n_index*disp, MPI_INTEGER, &
                        mpi_comm_local, ierr)

    call mpi_allgatherv(local_hp_list%coefs, 2*local_hp_list%num_pairs, PREC_MPI, &
                        hp_list(irep)%coefs, 2*num_all_pairs, 2*disp, PREC_MPI, &
                        mpi_comm_local, ierr)

    TIME_E(tmc_neighbor)
#endif

end subroutine simu_neighbor_list_hp
