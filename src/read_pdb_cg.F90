!read_pdb_cg
!> @brief Reads structure information about proteins from PDB files.

subroutine read_pdb_cg(lun, & ![i ] target I/O unit
                       nunit, & ![io] the number of unit
                       nmp, & ![io] the number of mp(mass point)
                       nres, & ![io] the number of residue
                       lunit2mp, & ![ o] correspondence list (unit -> mp  )
                       ires_mp, & ![ o] residue index       (mp   -> res )
                       xyz_mp, & ![ o] coordinate of mp    (mp   -> xyz )
                       cmp2seq, & ![ o] correspondence list (mp   -> seq )
                       cmp2atom, & ![ o] correspondence list (mp   -> atom)
                       imp2type, & ![ o] correspondence list (mp   -> type)
                       imp2base)   ![ o] correspondence list (mp   -> base)

    use const_maxsize
    use const_index
    implicit none

    ! ---------------------------------------------------------------------
    integer, intent(in)    :: lun
    integer, intent(inout) :: nmp, nunit, nres
    integer, intent(out)   :: lunit2mp(2, MXUNIT)
    integer, intent(out)   :: ires_mp(MXMP)
    real(PREC), intent(out)   :: xyz_mp(3, MXMP)
    character(3), intent(out)   :: cmp2seq(MXMP)
    character(4), intent(out)   :: cmp2atom(MXMP)
    integer, intent(out)   :: imp2type(MXMP)
    integer, intent(out)   :: imp2base(MXMP)

    ! ---------------------------------------------------------------------
    ! local variables
    integer :: input_status
    integer :: imp, iatom, iunit, ires
    integer :: iresnum, iresnum_save
    integer :: iclass, iclass_save
    integer :: itype, ibase
    logical :: flg_reading
    real(PREC) :: x, y, z
    real(PREC) :: tempfactor, occupancy  ! These variables are just read, not used.
    character(72) :: char72
    character(1) :: multistruct
    character(4) :: nameofatom
    character(6) :: nameid
    character(3) :: nameofmp
    character(2) :: chainid, chainid_save
    character(CARRAY_MSG_ERROR) :: error_message
#ifdef _DEBUG
    integer :: line
    write (*, *) '#### start read_pdb_cg'
#endif

    ! ---------------------------------------------------------------------
    flg_reading = .false.
    imp = nmp
    iunit = nunit - 1
    ires = nres
#ifdef _DEBUG
    line = 0
    write (*, *) 'imp = ', imp
    write (*, *) 'iunit = ', iunit
    write (*, *) 'ires = ', ires
#endif
    ! ---------------------------------------------------------------------
    rewind (lun)

    do
        read (lun, '(a72)', iostat=input_status) char72
        if (input_status < 0) then
            exit
        else if (input_status > 0) then
            error_message = 'Error: input error in read_pdb_cg'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef _DEBUG
        line = line + 1
#endif
        if (flg_reading .and. (char72(1:3) == 'TER' .OR. char72(1:2) == '>>')) then
            lunit2mp(2, iunit) = imp
            flg_reading = .false.
        end if

        if (char72(1:4) == 'ATOM') then
            read (char72, '(a6, i5, 1x, a4, a1, a3, a2, i4, 4x, 3f8.3, 2f6.2)', &
                  iostat=input_status) &
                nameid, iatom, nameofatom, multistruct, nameofmp, &
                chainid, iresnum, x, y, z, occupancy, tempfactor
            if (input_status > 0) then
                error_message = 'Error: cannot read pdb file in read_pdb_cg'
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            call sub_identify_class_and_type(nameofatom, nameofmp, iclass, itype, ibase)

            if (.not. flg_reading) then
                iunit = iunit + 1
                ires = ires + 1
                chainid_save = chainid
                iresnum_save = iresnum
                iclass_save = iclass
                flg_reading = .true.

            else
                if (chainid /= chainid_save) then
                    error_message = 'Error: chainid /= chainid_save in read_pdb_cg'
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                if (iclass /= iclass_save) then
                    error_message = 'Error: iclass /= iclass_save in read_pdb_cg'
                    call util_error(ERROR%STOP_ALL, error_message)
                endif
                if (iresnum /= iresnum_save) then
                    ires = ires + 1
                    iresnum_save = iresnum
                endif
            end if

            if (multistruct /= " " .and. multistruct /= "A") then
                cycle
            end if

            imp = imp + 1
            ires_mp(imp) = ires
            imp2type(imp) = itype
            imp2base(imp) = ibase
            cmp2seq(imp) = nameofmp
            cmp2atom(imp) = nameofatom
            xyz_mp(1, imp) = x
            xyz_mp(2, imp) = y
            xyz_mp(3, imp) = z

        end if
    end do

    if (flg_reading) then
        lunit2mp(2, iunit) = imp
    endif

    nunit = iunit
    nmp = imp
    nres = ires

#ifdef _DEBUG
    write (*, *) 'nmp = ', nmp
    write (*, *) 'nunit = ', nunit
    write (*, *) 'nres = ', nres
    write (*, *) '#### end read_pdb_cg'
#endif

contains

    subroutine sub_identify_class_and_type(c_atom, c_res, i_class, i_type, i_base)
        use const_index
        implicit none
        character(4), intent(in) :: c_atom
        character(3), intent(in) :: c_res
        integer, intent(out) :: i_class
        integer, intent(out) :: i_type, i_base

        i_class = CLASS%VOID
        i_type = MPTYPE%VOID
        i_base = BASETYPE%VOID

        ! Protein
        if (c_atom == ' CA ') then
            i_class = CLASS%PRO
            i_type = MPTYPE%PRO

            ! DNA
        elseif (c_res == ' DA' .OR. c_res == ' DT' .OR. &
                c_res == ' DG' .OR. c_res == ' DC') then
            i_class = CLASS%DNA

            if (c_atom == ' S  ') then
                i_type = MPTYPE%DNA_SUGAR
            elseif (c_atom == ' N  ') then
                i_type = MPTYPE%DNA_BASE
            elseif (c_atom == ' O  ') then
                i_type = MPTYPE%DNA_PHOS
            endif

            ! DNA2
        elseif (c_res == 'DA ' .OR. c_res == 'DT ' .OR. &
                c_res == 'DG ' .OR. c_res == 'DC ') then
            i_class = CLASS%DNA2

            if (c_atom == 'DS  ') then
                i_type = MPTYPE%DNA2_SUGAR
                i_base = BASETYPE%S
            elseif (c_atom == 'DB  ') then
                i_type = MPTYPE%DNA2_BASE
                if (c_res == 'DA ') then
                    i_base = BASETYPE%A
                else if (c_res == 'DT ') then
                    i_base = BASETYPE%T
                else if (c_res == 'DG ') then
                    i_base = BASETYPE%G
                else if (c_res == 'DC ') then
                    i_base = BASETYPE%C
                end if
            elseif (c_atom == 'DP  ') then
                i_type = MPTYPE%DNA2_PHOS
                i_base = BASETYPE%P
            endif

            ! RNA
        elseif (c_res == 'RA ' .OR. c_res == 'RU ' .OR. c_res == 'RG ' .OR. &
                c_res == 'RC ' .OR. c_res == 'OMG' .OR. c_res == '2MG' .OR. &
                c_res == 'M2G' .OR. c_res == 'YYG' .OR. c_res == '1MG' .OR. &
                c_res == '7MG' .OR. c_res == 'OMU' .OR. c_res == 'PSU' .OR. &
                c_res == 'H2U' .OR. c_res == 'OMC' .OR. c_res == '5MC' .OR. &
                c_res == '5MU' .OR. c_res == '1MA' .OR. c_res == 'UR3' .OR. &
                c_res == 'MA6' .OR. c_res == '4OC' .OR. c_res == '  I' .OR. &
                c_res == 'PYY' .OR. c_res == 'CB2' .OR. c_res == 'AET' .OR. &
                c_res == 'QUO' .OR. c_res == 'FHU' .OR. c_res == 'T6A' .OR. &
                c_res == 'RIA' .OR. c_res == 'MIA' .OR. c_res == 'PRF') then

            i_class = CLASS%RNA

            if (c_atom == ' S  ') then
                i_type = MPTYPE%RNA_SUGAR
            elseif (c_atom == ' P  ') then
                i_type = MPTYPE%RNA_PHOS
            elseif (c_atom == ' Ab ' .OR. c_atom == ' Ub ' .OR. c_atom == ' Rb ' .OR. &
                    c_atom == ' Gb ' .OR. c_atom == ' Cb ' .OR. c_atom == ' Yb ') then
                i_type = MPTYPE%RNA_BASE
            endif

        endif

        ! Error
        if (i_class == CLASS%VOID) then
            error_message = 'Error: unknown class in read_pdb_cg; '//c_res//','//c_atom
            call util_error(ERROR%STOP_ALL, error_message)
        endif
        if (i_type == MPTYPE%VOID) then
            error_message = 'Error: unknown mp-type in read_pdb_cg; '//c_res//','//c_atom
            call util_error(ERROR%STOP_ALL, error_message)
        endif

    endsubroutine sub_identify_class_and_type

end subroutine read_pdb_cg
