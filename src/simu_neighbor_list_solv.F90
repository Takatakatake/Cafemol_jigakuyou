! simu_neighbor_list_solv
!> @brief Construct neighborling list for DNA solvation interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_list_solv(irep)

    use const_maxsize
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: solv_list, clear_list, push_pair
    use var_setp, only: inmisc, indna
    use var_struct, only: nmp_real, nunit_real, lunit2mp, xyz_mp_rep, &
        pxyz_mp_rep, imp2unit, imp2type, nmp_all
    use time
#ifdef MPI_PAR
    use mpiconst
#ifdef SHARE_NEIGH
    use var_neighbor_list, only: neighbor_list_type, &
        reserve_neighbor_list, dealloc_neighbor_list
#endif
#endif

    implicit none

    ! ---------------------------------------------------------------------
    integer, intent(in) :: irep

    ! ---------------------------------------------------------------------
    ! local variables
    integer :: from, to, i
    integer :: imp, jmp, iunit, junit
    integer :: imirror
    integer :: icalc(MXUNIT, MXUNIT)
    real(PREC) :: v21(3)
    real(PREC) :: rneighbor2_solv
    character(CARRAY_MSG_ERROR) :: error_message

#if defined(MPI_PAR2) && defined(SHARE_NEIGH_SOLV)
    type(neighbor_list_type) :: local_solv_list
    integer :: n, n_index
    integer :: num_all_pairs(0:npar_mpi - 1)
    integer :: disp(0:npar_mpi - 1)
#endif

    ! -------------------------------------------------------------------
    if (.not. inmisc%force_flag(INTERACT%DNA)) then
        return
    end if

    call clear_list(solv_list(irep))

    ! -------------------------------------------------------------------
    ! calc isolv2mp
    rneighbor2_solv = (1.2*(indna%cutoff_solv_dna*indna%cralpha_solv_dna &
                            + indna%cutoff_solv_dna))**2

    icalc(1:nunit_real, 1:nunit_real) = 0
    do iunit = 1, nunit_real
        do junit = iunit, nunit_real
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA)) then
                icalc(iunit, junit) = 1
            end if
        end do
    end do

    ! --------------------------------------------------------------------

#ifdef MPI_PAR2
    from = 1
    to = nmp_l
#else
    from = 1
    to = nmp_real - 1
#endif

    do i = from, to
#ifdef MPI_PAR2
        imp = imp_l2g(i)
#else
        imp = i
#endif
        if (imp2type(imp) /= MPTYPE%DNA_SUGAR) then
            cycle
        end if

        iunit = imp2unit(imp)
        if (iunit >= nunit_real) then
            exit
        end if

        jmp = lunit2mp(1, iunit + 1)
        do while (jmp <= nmp_real)
            if (imp2type(jmp) /= MPTYPE%DNA_SUGAR) then
                jmp = jmp + 1
                cycle
            end if

            junit = imp2unit(jmp)
            if (icalc(iunit, junit) /= 1) then
                jmp = lunit2mp(2, junit) + 1
                cycle
            end if

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
            else
                v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
                call util_pbneighbor(v21, imirror)
            end if

            if (dot_product(v21, v21) < rneighbor2_solv) then
#if defined(MPI_PAR2) && defined(SHARE_NEIGH_SOLV)
                call push_pair(local_solv_list, imp, jmp)
                if (inperi%i_periodic == 1) then
                    local_solv_list%pairs(3, local_solv_list%num_pairs) = imirror
                end if
#else
                call push_pair(solv_list(irep), imp, jmp)
                if (inperi%i_periodic == 1) then
                    solv_list(irep)%pairs(3, solv_list(irep)%num_pairs) = imirror
                end if
#endif
            end if

            jmp = jmp + 1
        end do
    end do

#if defined(MPI_PAR2) && defined(SHARE_NEIGH_SOLV)
    n_index = 2 + inperi%n_mirror_index

    TIME_S(tmc_neighbor)

    call mpi_allgather(local_solv_list%num_pairs, 1, MPI_INTEGER, &
                       num_all_pairs, 1, MPI_INTEGER, &
                       mpi_comm_local, ierr)

    solv_list(irep)%num_pairs = sum(num_all_pairs(0:npar_mpi - 1))

    disp(0) = 0
    do n = 1, npar_mpi - 1
        disp(n) = disp(n - 1) + num_all_pairs(n - 1)
    end do

    call mpi_allgatherv(local_solv_list%pairs, n_index*local_solv_list%num_pairs, MPI_INTEGER, &
                        solv_list(irep)%pairs(1, 1), n_index*num_all_pairs, n_index*disp, MPI_INTEGER, &
                        mpi_comm_local, ierr)

    TIME_E(tmc_neighbor)
#endif

    if (solv_list(irep)%num_pairs > (MXMPSOLV*nmp_all)) then
        error_message = 'Error: too many pairs in simu_neighbor_list_solv'
        call util_error(ERROR%STOP_ALL, error_message)
    end if

end subroutine simu_neighbor_list_solv
