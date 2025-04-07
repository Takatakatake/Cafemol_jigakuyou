! simu_neighbor_list_gb
!> @brief Construct a neighbor list for gb interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_list_gb(jrep)

    use const_maxsize
    use const_index
    use var_neighbor_list, only: ele_gb_list, pnl_gb_list, clear_list, push_pair
    use var_setp, only: inmisc, ingb
    use var_struct, only: nunit_real, pxyz_mp_rep, &
        imp2unit, nmp_all
    use var_replica, only: irep2grep
    use time
#ifdef MPI_PAR2
    use mpiconst
#endif

    implicit none

    ! -------------------------------------------------------------------
    integer, intent(in) :: jrep

    ! -------------------------------------------------------------------
    ! local variables
    integer :: imp, jmp, iunit, junit, irep, grep
    integer :: icalc(MXUNIT, MXUNIT)
    real(PREC) :: dist2, rneighbor2_gb, v21(3)
    character(CARRAY_MSG_ERROR) :: error_message

    ! -------------------------------------------------------------------

    call clear_list(ele_gb_list(jrep))
    call clear_list(pnl_gb_list(jrep))

    if (.not. inmisc%force_flag(INTERACT%GB)) then
        return
    end if

    ! -------------------------------------------------------------------
    icalc(1:nunit_real, 1:nunit_real) = 0
    do iunit = 1, nunit_real
        do junit = iunit, nunit_real
            if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%GB)) then
                icalc(iunit, junit) = 1
            end if
        end do
    end do

    ! --------------------------------------------------------------------
    irep = jrep
    grep = irep2grep(irep)
    rneighbor2_gb = ingb%cutoff_ele*ingb%cutoff_ele

    do imp = 1, nmp_all - 1

        iunit = imp2unit(imp)

        jmp = imp + 1

        if (jmp > nmp_all) cycle

        do while (jmp <= nmp_all)
            junit = imp2unit(jmp)
            v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
            dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
            if ((junit == iunit) .and. (jmp > imp + 2) .and. (dist2 < rneighbor2_gb)) then
                call push_pair(ele_gb_list(irep), imp, jmp)
            else if ((junit /= iunit) .and. (dist2 < rneighbor2_gb)) then
                call push_pair(ele_gb_list(irep), imp, jmp)
            end if
            jmp = jmp + 1
        end do
    end do

    if (ele_gb_list(irep)%num_pairs > MXMPNEIGHBOR*nmp_all) then
        error_message = 'Error: too big npnl in simu_neighbor_list_gb'
        call util_error(ERROR%STOP_ALL, error_message)
    end if

    ! --------------------------------------------------------------------
    irep = jrep
    grep = irep2grep(irep)
    rneighbor2_gb = ingb%cutoff_gb*ingb%cutoff_gb

    do imp = 1, nmp_all - 1

        iunit = imp2unit(imp)

        jmp = imp + 1

        if (jmp > nmp_all) cycle

        do while (jmp <= nmp_all)
            junit = imp2unit(jmp)
            v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
            dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
            if ((junit == iunit) .and. (jmp > imp + 2) .and. (dist2 < rneighbor2_gb)) then
                call push_pair(pnl_gb_list(irep), imp, jmp)
            else if ((junit /= iunit) .and. (dist2 < rneighbor2_gb)) then
                call push_pair(pnl_gb_list(irep), imp, jmp)
            end if
            jmp = jmp + 1
        end do
    end do

    if (pnl_gb_list(irep)%num_pairs > MXMPNEIGHBOR*nmp_all) then
        error_message = 'Error: too big npnl in simu_neighbor_list_gb'
        call util_error(ERROR%STOP_ALL, error_message)
    end if

end subroutine simu_neighbor_list_gb
