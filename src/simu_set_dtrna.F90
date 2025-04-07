! simu_set_dtrna
!> @brief Calculate stack parameter U0 in DT-RNA model

! ***********************************************************************
subroutine simu_set_dtrna(grep, tempk)

    use const_maxsize
    use const_physical
    use const_index
    use var_setp, only: indtrna, inmisc
    use var_struct, only: idtrna_st2nn, ndtrna_st, coef_dtrna_st

    implicit none
    ! ----------------------------------------------------------------------
    integer :: ist
    integer :: inn
    real(PREC) :: h, s, Tm
    integer, intent(in) :: grep
    real(PREC), intent(in) :: tempk

    ! -----------------------------------------------------------------------

    do ist = 1, ndtrna_st
        inn = idtrna_st2nn(ist)
        h = indtrna%cst_h(inn)
        s = indtrna%cst_s(inn)
        Tm = indtrna%cst_Tm(inn)
        if (inmisc%i_temp_independent == 0) then
            coef_dtrna_st(0, ist, grep) = -h + BOLTZC*(tempk - Tm)*s
        else if (inmisc%i_temp_independent > 0) then
            coef_dtrna_st(0, ist, grep) = -h - BOLTZC*Tm*s
        endif
    enddo

end subroutine simu_set_dtrna
