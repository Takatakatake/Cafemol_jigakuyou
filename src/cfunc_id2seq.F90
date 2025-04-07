! cfunc_id2seq
!> @brief This function converts amino-acid identification number(id) into the 3-letter codes.

!********************************************************************
character(3) function cfunc_id2seq(id)

    use const_maxsize
    use const_index

    implicit none
    ! ----------------------------------------------------------------
    integer, intent(in) :: id

    ! ----------------------------------------------------------------
    ! local variables
    character(CARRAY_MSG_ERROR) :: error_message

    select case (id)
    case (CHEMICALTYPE%ALA)
        cfunc_id2seq = 'ALA'
    case (CHEMICALTYPE%ARG)
        cfunc_id2seq = 'ARG'
    case (CHEMICALTYPE%ASN)
        cfunc_id2seq = 'ASN'
    case (CHEMICALTYPE%ASP)
        cfunc_id2seq = 'ASP'
    case (CHEMICALTYPE%CYS)
        cfunc_id2seq = 'CYS'
    case (CHEMICALTYPE%GLN)
        cfunc_id2seq = 'GLN'
    case (CHEMICALTYPE%GLU)
        cfunc_id2seq = 'GLU'
    case (CHEMICALTYPE%GLY)
        cfunc_id2seq = 'GLY'
    case (CHEMICALTYPE%HIS)
        cfunc_id2seq = 'HIS'
    case (CHEMICALTYPE%ILE)
        cfunc_id2seq = 'ILE'
    case (CHEMICALTYPE%LEU)
        cfunc_id2seq = 'LEU'
    case (CHEMICALTYPE%LYS)
        cfunc_id2seq = 'LYS'
    case (CHEMICALTYPE%MET)
        cfunc_id2seq = 'MET'
    case (CHEMICALTYPE%PHE)
        cfunc_id2seq = 'PHE'
    case (CHEMICALTYPE%PRO)
        cfunc_id2seq = 'PRO'
    case (CHEMICALTYPE%SER)
        cfunc_id2seq = 'SER'
    case (CHEMICALTYPE%THR)
        cfunc_id2seq = 'THR'
    case (CHEMICALTYPE%TRP)
        cfunc_id2seq = 'TRP'
    case (CHEMICALTYPE%TYR)
        cfunc_id2seq = 'TYR'
    case (CHEMICALTYPE%VAL)
        cfunc_id2seq = 'VAL'

    case default
        write (error_message, '(a,i10)') 'Error: in cfunc_id2seq there is no id such ', id
        call util_error(ERROR%STOP_ALL, error_message)
    endselect

end function cfunc_id2seq
