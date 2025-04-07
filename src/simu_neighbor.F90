! simu_neighbor
!> @brief This subroutine is to make neighborlist for the whole of non-local interaction.

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor(irep)

    use const_maxsize
    use const_index
    use var_inp, only: inperi
    use var_setp, only: inmisc, inele
    use var_struct, only: nmp_real

    use time
    use mpiconst
    implicit none

    integer, intent(in) :: irep
    ! -------------------------------------------------------------------
    ! local variables
    integer :: num_pairs(0:nthreads - 1)
    integer :: pairs(2, MXMPNEIGHBOR*nmp_real/nthreads, 0:nthreads - 1)

    ! -------------------------------------------------------------------

    if (inperi%i_periodic == 1) then
        call util_periodic(irep)
    end if

    TIME_S(tm_neighbor_pnl)
    call simu_neighbor_list(irep, num_pairs, pairs)
    call simu_neighbor_assign(irep, num_pairs, pairs)
    TIME_E(tm_neighbor_pnl)

    TIME_S(tm_neighbor_ele)
    if (inmisc%force_flag(INTERACT%ELE) .or. &
        inmisc%force_flag(INTERACT%GB)) then
        if (inele%i_calc_method == 0) then
            call simu_neighbor_list_ele(irep)
        else if (inele%i_calc_method == 1) then
            call simu_neighbor_list_ele2(irep)
        end if
    endif
    TIME_E(tm_neighbor_ele)

    TIME_S(tm_neighbor_gb)
    if (inmisc%force_flag(INTERACT%GB)) then
        call simu_neighbor_list_gb(irep)
    endif
    TIME_E(tm_neighbor_gb)

    TIME_S(tm_neighbor_solv)
    if (inmisc%force_flag(INTERACT%DNA)) then
        call simu_neighbor_list_solv(irep)
    endif
    TIME_E(tm_neighbor_solv)

    TIME_S(tm_neighbor_hp)
    if (inmisc%force_flag(INTERACT%HP)) then
        call simu_neighbor_list_hp(irep)
    endif
    TIME_E(tm_neighbor_hp)

end subroutine simu_neighbor
