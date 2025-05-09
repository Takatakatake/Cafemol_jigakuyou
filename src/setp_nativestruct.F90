! setp_nativestruct
!> @brief This subroutine is to set up the information for native structure that is related to &
!>        not only local-structure but also nonlocal structure. &
!>        (This subroutine makes the neighborlist.)

! ***********************************************************************
subroutine setp_nativestruct(xyz_mp_init &  ! [i ]
                             , iatomnum &  ! [i ]
                             , xyz &  ! [i ]
                             , cname_ha &  ! [i ]
                             )

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: i_go_native_read_style
    use var_struct, only: nmp_all
    use var_setp, only: inmisc

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! --------------------------------------------------------------------
    real(PREC), intent(in) :: xyz_mp_init(SPACE_DIM, MXMP)
    integer, intent(in) :: iatomnum(MXMP)
    real(PREC), intent(in) :: xyz(SPACE_DIM, MXATOM_MP, MXMP)
    character(4), intent(in):: cname_ha(MXATOM_MP, MXMP)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: lmp2neigh(MXMP)
    integer :: ineigh2mp(MXMPNEIGHBOR*nmp_all)
    real(PREC) :: xyz_dna(SPACE_DIM, MXMP)

    ! --------------------------------------------------------------------
#ifdef _DEBUG
    write (*, *) '##### start setp_nativestruct'
#endif

    if (inmisc%class_flag(CLASS%DNA)) then
        call setp_make_dna(xyz_dna)
        call setp_native_stack_dna(xyz_dna)
    endif

    if (i_go_native_read_style /= NATIVEREAD%INFO) then
        ! secondary structure
        call setp_native_secstruct()

        ! local
        call setp_native_bond(xyz_mp_init, xyz_dna)
        call setp_native_bangle(xyz_mp_init, xyz_dna)
        call setp_native_dih(xyz_mp_init, xyz_dna)

        ! make neighborlist
        call setp_native_neighbor_list(xyz_mp_init, ineigh2mp, lmp2neigh)

        ! nonlocal
        call setp_native_go(xyz_mp_init, iatomnum, xyz, ineigh2mp, lmp2neigh, cname_ha)

    endif

    !if (inmisc%force_flag(INTERACT%DTRNA)) then
    !   call setp_native_dtrna(xyz_mp_init)
    !endif

    ! DNA (3SPN.2 model)
    call setp_native_dna2()

    if (inmisc%force_flag(INTERACT%ELE) .or. inmisc%force_flag(INTERACT%GB)) then
        call setp_native_charge()
    endif

    !set up variables for hydrohpobic interaction
    call setp_native_hp()

    !set up variables for sasa
    if (inmisc%force_flag(INTERACT%SASA)) then
        call setp_native_sasa()  !sasa
    end if

#ifdef _DEBUG
    write (*, *) '##### end setp_nativestruct'
#endif

end subroutine setp_nativestruct
