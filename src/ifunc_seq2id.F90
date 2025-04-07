! ifunc_seq2id
!> @brief This function converts 3-letter for the amino-acid name into the identification number (id).

!********************************************************************
integer function ifunc_seq2id(name)

    use const_maxsize
    use const_index

    implicit none
    ! ----------------------------------------------------------------
    character(3), intent(in) :: name

    ! ----------------------------------------------------------------
    ! local variables
    character(CARRAY_MSG_ERROR) :: error_message

    ifunc_seq2id = 1
    if (name == 'ALA') then
        ifunc_seq2id = 1
    else if (name == 'ARG') then
        ifunc_seq2id = 2
    else if (name == 'ASN') then
        ifunc_seq2id = 3
    else if (name == 'ASP') then
        ifunc_seq2id = 4
    else if (name == 'CYS') then
        ifunc_seq2id = 5
    else if (name == 'GLN') then
        ifunc_seq2id = 6
    else if (name == 'GLU') then
        ifunc_seq2id = 7
    else if (name == 'GLY') then
        ifunc_seq2id = 8
    else if (name == 'HIS' .or. name == 'HIE' .or. name == 'HID' .or. name == 'HIP' .or. &
             name == 'HSE' .or. name == 'HSD' .or. name == 'HSP') then
        ifunc_seq2id = 9
    else if (name == 'ILE') then
        ifunc_seq2id = 10
    else if (name == 'LEU') then
        ifunc_seq2id = 11
    else if (name == 'LYS') then
        ifunc_seq2id = 12
    else if (name == 'MET') then
        ifunc_seq2id = 13
    else if (name == 'PHE') then
        ifunc_seq2id = 14
    else if (name == 'PRO') then
        ifunc_seq2id = 15
    else if (name == 'SER') then
        ifunc_seq2id = 16
    else if (name == 'THR') then
        ifunc_seq2id = 17
    else if (name == 'TRP') then
        ifunc_seq2id = 18
    else if (name == 'TYR') then
        ifunc_seq2id = 19
    else if (name == 'VAL') then
        ifunc_seq2id = 20
    else if (name == 'OTH') then
        ifunc_seq2id = 21
    else if (name == 'P  ') then
        ifunc_seq2id = 22

! added for sasa
!-----------------------------------------
    else if (name == ' RA' .OR. &
             name == 'RA ' .OR. &
             name == 'A  ' .OR. &
             name == 'RA5' .OR. &
             name == 'RA3' .OR. &
             name == '1MA' .OR. &
             name == 'MA6' .OR. &
             name == 'T6A' .OR. &
             name == 'RIA' .OR. &
             name == 'MIA') then
        ifunc_seq2id = 21

    else if (name == ' RC' .OR. &
             name == 'RC ' .OR. &
             name == '  C' .OR. &
             name == 'RC5' .OR. &  ! AMBER
             name == 'RC3' .OR. &  ! AMBER
             name == 'OMC' .OR. &
             name == '5MC' .OR. &
             name == '4OC') then
        ifunc_seq2id = 22

    else if (name == ' RG' .OR. &
             name == 'RG ' .OR. &
             name == '  G' .OR. &
             name == 'RG ' .OR. &
             name == 'RG5' .OR. &
             name == 'RG3' .OR. &
             name == 'OMG' .OR. &
             name == '2MG' .OR. &
             name == 'M2G' .OR. &
             name == 'YYG' .OR. &
             name == '1MG' .OR. &
             name == '7MG') then
        ifunc_seq2id = 23

    else if (name == ' RT' .OR. &
             name == 'RT ' .OR. &
             name == '  T') then
        ifunc_seq2id = 24

    else if (name == ' RU' .OR. &
             name == 'RU ' .OR. &
             name == 'U  ' .OR. &
             name == 'RU5' .OR. &  ! AMBER
             name == 'RU3' .OR. &  ! AMBER
             name == 'OMU' .OR. &
             name == 'PSU' .OR. &
             name == 'H2U' .OR. &
             name == '5MU') then
        ifunc_seq2id = 25

    else if (name == '  I' .OR. &
             name == 'RI ' .OR. &
             name == ' RI') then
        ifunc_seq2id = 26
    else if (name == ' DA' .OR. &
             name == 'DA ') then
        ifunc_seq2id = 27
    else if (name == ' DC' .OR. &
             name == 'DC ') then
        ifunc_seq2id = 28
    else if (name == ' DG' .OR. &
             name == 'DG ') then
        ifunc_seq2id = 29
    else if (name == ' DT' .OR. &
             name == 'DT ') then
        ifunc_seq2id = 30
    else if (name == ' DU' .OR. &
             name == 'DU ') then
        ifunc_seq2id = 31
    else if (name == ' DI' .OR. &
             name == 'DI ') then
        ifunc_seq2id = 32
!------------------------------------------

    else
        error_message = 'Error: in ifunc_seq2id there is no name such '//name
        call util_error(ERROR%STOP_ALL, error_message)
    end if

end function ifunc_seq2id
