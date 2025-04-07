! setp_exv
!> @brief store residue type specific excluded radius
!>        for protein and DNA2. For DNA1 and RNA, this routine gives error and stop the execlution.

!#define _DEBUG
! ***********************************************************************
subroutine setp_exv()
    use const_maxsize
    use const_index
    use const_physical
    use var_setp, only: inexv, inmisc
    use var_struct, only: iclass_unit, imp2type, imp2unit, &
        exv_radius_mp, cmp2seq, imp2type, nmp_all

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! -----------------------------------------------------------------
    ! local variables
    integer :: imp, iunit
    character(CARRAY_MSG_ERROR) :: error_message
    !real(PREC) :: dist2

#ifdef _DEBUG
    write (*, *) '#### start set up residue-type-specific'
#endif

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        ! -----------------------------------------------------------------
        do imp = 1, nmp_all
            iunit = imp2unit(imp)

            ! if the unit is protein or ligand
            if (iclass_unit(iunit) == CLASS%PRO .or. iclass_unit(iunit) == CLASS%LIG) then

                ! if the unit is CA model or ligand
                if (inmisc%i_use_atom_protein == USE_PRO%CA .or. iclass_unit(iunit) == CLASS%LIG) then

                    if (cmp2seq(imp) == 'ALA') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ALA)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'ARG') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ARG)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'ASN') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ASN)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'ASP') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ASP)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'CYS') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%CYS)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'GLN') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%GLN)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'GLU') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%GLU)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'GLY') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%GLY)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'HIS') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%HIS)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'ILE') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ILE)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'LEU') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%LEU)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'LYS') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%LYS)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'MET') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%MET)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'PHE') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%PHE)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'PRO') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%PRO)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'SER') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%SER)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'THR') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%THR)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'TRP') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%TRP)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'TYR') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%TYR)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'VAL') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%VAL)*0.5e0_PREC
                    else if (cmp2seq(imp) == 'OTH') then
                        exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%OTH)*0.5e0_PREC
                    else
                        write (error_message, *) 'Error: Unknown residue type for protein (or ligand) in setp_exv'
                        call util_error(ERROR%STOP_ALL, error_message)
                    end if

                else
           write(error_message,*) 'Error: Non C-alpha representation for Protein model does not support residue-type-specific excluded volume radii'
                    call util_error(ERROR%STOP_ALL, error_message)
                end if

            else if (iclass_unit(iunit) == CLASS%DNA2) then

                LABEL_ATOM:select case(imp2type(imp))
                case (MPTYPE%DNA2_PHOS) LABEL_ATOM
                exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DP)*0.5e0_PREC

                case (MPTYPE%DNA2_SUGAR) LABEL_ATOM
                exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DS)*0.5e0_PREC

                case (MPTYPE%DNA2_BASE) LABEL_ATOM
                if (cmp2seq(imp) == ' DA' .or. cmp2seq(imp) == 'DA ') then
                    exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DA)*0.5e0_PREC
                else if (cmp2seq(imp) == ' DT' .or. cmp2seq(imp) == 'DT ') then
                    exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DT)*0.5e0_PREC
                else if (cmp2seq(imp) == ' DC' .or. cmp2seq(imp) == 'DC ') then
                    exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DC)*0.5e0_PREC
                else if (cmp2seq(imp) == ' DG' .or. cmp2seq(imp) == 'DG ') then
                    exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DG)*0.5e0_PREC
                endif
                case default LABEL_ATOM
                write (error_message, *) 'Error: failed in setp_exv', &
                    'iunit=', iunit, 'imp=', imp
                call util_error(ERROR%STOP_ALL, error_message)
                end select LABEL_ATOM

            else
                write (error_message, *) 'Error: This model does not support residue-type-specific excluded volume radii'
                call util_error(ERROR%STOP_ALL, error_message)
            end if

        end do

#ifdef MPI_PAR
    end if
!call MPI_Bcast(inexv, inexv%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
    call MPI_Bcast(exv_radius_mp, MXMP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_exv
