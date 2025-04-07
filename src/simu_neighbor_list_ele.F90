! simu_neighbor_list_ele
!> @brief Construct a neighbor list for electrostatic interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_list_ele(jrep)

    use const_maxsize
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: ele_list, clear_list, push_pair
    use var_setp, only: inmisc, inele
    use var_struct, only: nunit_real, xyz_mp_rep, pxyz_mp_rep, &
        imp2unit, iclass_unit, &
        ncharge, icharge2mp, coef_charge
    use var_replica, only: irep2grep
    use time
#ifdef MPI_PAR2
    use mpiconst
#ifdef SHARE_NEIGH
    use var_neighbor_list, only: neighbor_list_type, &
        reserve_neighbor_list, dealloc_neighbor_list
#endif
#endif

    implicit none

    ! -------------------------------------------------------------------
    integer, intent(in) :: jrep

    ! -------------------------------------------------------------------
    ! local variables
    integer :: imp, jmp, iunit, junit, irep, grep
    integer :: icharge, jcharge, imirror, n_index
    integer :: icalc(MXUNIT, MXUNIT)
    real(PREC) :: dist2, rneighbor2_ele, v21(3)
    character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR2
    integer :: icharge_l

#ifdef SHARE_NEIGH
    type(neighbor_list_type) :: local_ele_list
    integer :: ierr

    integer :: n
    integer :: num_all_pairs(0:npar_mpi - 1)
    integer :: disp(0:npar_mpi - 1)
#endif

#endif

    ! -------------------------------------------------------------------
    n_index = 2 + inperi%n_mirror_index

    if (inmisc%force_flag(INTERACT%ELE) .or. inmisc%force_flag(INTERACT%GB)) then
        continue
    else
        return
    end if

    call clear_list(ele_list(jrep))

    ! -------------------------------------------------------------------
    icalc(1:nunit_real, 1:nunit_real) = 0
    do iunit = 1, nunit_real
        do junit = iunit, nunit_real
!        if(mod(inmisc%itype_nlocal(iunit, junit), INTERACT%ELE) == 0) then
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ELE) .or. &
                inmisc%flag_nlocal_unit(iunit, junit, INTERACT%GB)) then
                icalc(iunit, junit) = 1
            end if
        end do
    end do

    ! --------------------------------------------------------------------

    irep = jrep   ! to avoid intel compiler internal error.(@dan)

    grep = irep2grep(irep)
    rneighbor2_ele = (1.2*inele%cutoff_ele*inele%cdist(grep))**2

#if defined(MPI_PAR2) && defined(SHARE_NEIGH)
    call reserve_neighbor_list(local_ele_list, ele_list(irep), ierr)
#endif

#ifdef MPI_PAR2
    do icharge_l = 1, ncharge_l
        icharge = icharge_l2g(icharge_l)
#else
        do icharge = 1, ncharge - 1
#endif

            imp = icharge2mp(icharge)
            iunit = imp2unit(imp)

            jcharge = icharge + 1

            if (jcharge > ncharge) cycle

            if (iclass_unit(iunit) == CLASS%DNA .or. iclass_unit(iunit) == CLASS%DNA2) then
                jmp = icharge2mp(jcharge)
                junit = imp2unit(jmp)
                if (iunit == junit) then
                    jcharge = jcharge + 1
                end if
            end if

            do while (jcharge <= ncharge)
                jmp = icharge2mp(jcharge)
                junit = imp2unit(jmp)

                if (icalc(iunit, junit) == 1) then

                    if (inperi%i_periodic == 0) then
                        v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
                    else
                        v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
                        call util_pbneighbor(v21, imirror)
                    end if

                    dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

                    if (dist2 < rneighbor2_ele) then
#if defined(MPI_PAR2) && defined(SHARE_NEIGH)
                        call push_pair(local_ele_list, imp, jmp)
                        if (inperi%i_periodic == 1) then
                            local_ele_list%pairs(3, local_ele_list%num_pairs) = imirror
                        end if
                        local_ele_list%coefs(1, local_ele_list%num_pairs) = &
                            coef_charge(icharge, grep)*coef_charge(jcharge, grep)*inele%coef(grep)
#else
                        call push_pair(ele_list(irep), imp, jmp)
                        if (inperi%i_periodic == 1) then
                            ele_list(irep)%pairs(3, ele_list(irep)%num_pairs) = imirror
                        end if
                        ele_list(irep)%coefs(1, ele_list(irep)%num_pairs) = &
                            coef_charge(icharge, grep)*coef_charge(jcharge, grep)*inele%coef(grep)
#endif
                    end if
                end if

                jcharge = jcharge + 1
            end do
        end do

#if defined(MPI_PAR2) && defined(SHARE_NEIGH)
        TIME_S(tmc_neighbor)

        call mpi_allgather(local_ele_list%num_pairs, 1, MPI_INTEGER, &
                           num_all_pairs, 1, MPI_INTEGER, &
                           mpi_comm_local, ierr)

        ele_list(irep)%num_pairs = sum(num_all_pairs(0:npar_mpi - 1))

        disp(0) = 0
        do n = 1, npar_mpi - 1
            disp(n) = disp(n - 1) + num_all_pairs(n - 1)
        end do

        call mpi_allgatherv(local_ele_list%pairs, n_index*local_ele_list%num_pairs, MPI_INTEGER, &
                            ele_list(irep)%pairs, n_index*num_all_pairs, n_index*disp, MPI_INTEGER, &
                            mpi_comm_local, ierr)

        call mpi_allgatherv(local_ele_list%coefs, local_ele_list%num_pairs, PREC_MPI, &
                            ele_list(irep)%coefs, num_all_pairs, disp, PREC_MPI, &
                            mpi_comm_local, ierr)

        TIME_E(tmc_neighbor)

        call dealloc_neighbor_list(local_ele_list, MXMPELE*ncharge, 1, ierr)
#endif

        if (ele_list(irep)%num_pairs > (MXMPELE*ncharge)) then
            error_message = 'Error: too big nele in simu_neighbor_list_ele'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

! enddo ! irep

! DBG
#ifdef _DEBUG
        write (6, *) 'simu_neighbor_list_ele, irep, lele(irep) ', irep, ele_list(irep)%num_pairs
        call flush(6)
#endif

!  write (*, *) ncharge, lele(irep), sqrt(rneighbor2_ele)

    end subroutine simu_neighbor_list_ele
