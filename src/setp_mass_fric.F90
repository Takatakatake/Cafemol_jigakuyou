! setp_mass_fric
!> @brief This subroutine is to read the information from input-file and &
!>        set the mass and/or friction of either all, by chain, or by particles.

! *************************************************************************
subroutine setp_mass_fric()

    use const_maxsize
    use const_index
    use const_physical
    use var_inp, only: infile, outfile
    use var_setp, only: inmisc, inpara
    use var_struct, only: lunit2mp, cmass_mp, fric_mp, nmp_all

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! -----------------------------------------------------------------------
    ! local variables
    integer :: imp, icol, ichemical
    integer :: inunit(2), instate
    integer :: luninp, lunout
    integer :: iline, nlines, iequa, nequat
    real(PREC) :: rmass, fric, visc, radius, rmassunit2fric

    character(4) :: kfind
    character(12) :: char12
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MXCOLM) :: char02
    character(CARRAY_MSG_ERROR) :: error_message

    ! -----------------------------------------------------------------------
    luninp = infile%inp
    lunout = outfile%data

    ! -----------------------------------------------------------------------
    ! Mass unit
    if (inmisc%i_mass_unit == 0) then
        rmassunit2fric = sqrt(13.7e0_PREC)
    else if (inmisc%i_mass_unit == 1) then
        rmassunit2fric = 1.0e0_PREC
    endif

    ! -----------------------------------------------------------------------
    ! Default value
    cmass_mp(1:MXMP) = inpara%cmass(CHEMICALTYPE%UNKNOWN)
    fric_mp(1:MXMP) = inpara%fric_const

    ! -----------------------------------------------------------------------
    ! Particle-type dependent mass
    if (inmisc%i_mass == 1) then
        do imp = 1, nmp_all
            ichemical = imp2chemicaltype(imp)
            if (inpara%cmass(ichemical) <= 0) then
                write (error_message, *) 'Invalid mass: imp=', imp, ", ichemical=", ichemical, ", cmass=", inpara%cmass(ichemical)
                call util_error(ERROR%STOP_ALL, error_message)
            end if
            cmass_mp(imp) = inpara%cmass(ichemical)
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! Friction coefficients
    if (inmisc%i_fric == 1) then
        fric = rmassunit2fric*inpara%fric_scale*inpara%fric_typical_coef
        do imp = 1, nmp_all
            fric_mp(imp) = fric/cmass_mp(imp)
        end do

#ifdef MPI_PAR
        if (myrank == 0) then
#endif
            write (lunout, *) "Friction coefficients is dependent of mass (i_fric = 1)"
            do imp = 1, nmp_all
                write (lunout, *) '   friction coefficient of (', imp, '): ', fric_mp(imp)*cmass_mp(imp)
            enddo
#ifdef MPI_PAR
        end if
#endif

    end if

    ! -----------------------------------------------------------------------
    ! Friction coefficients by Stokes' law
    if (inmisc%i_fric == 11) then

        ! (Pa)(s) ---> (Da)(A^-1)(tau^-1)   (tau: cafemol time unit)
        visc = rmassunit2fric*inpara%viscosity*sqrt(1.0e3_PREC/JOUL2KCAL)*N_AVO*1.0e-20_PREC

        do imp = 1, nmp_all
            ichemical = imp2chemicaltype(imp)
            radius = inpara%stokes_radius(ichemical)
            fric_mp(imp) = 6.0e0_PREC*F_PI*visc*radius/cmass_mp(imp)
            if (ichemical <= 0 .or. radius <= 0 .or. fric_mp(imp) <= 0) then
      write (error_message, *) 'Invalid fric_mp: imp=', imp, ", ichemical=", ichemical, ", radius=", radius, ", fric=", fric_mp(imp)
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            !! Notice:
            !! This friction coefficient is divided by mass for a technical reason.
            !! In CafeMol, the Langevin equation contains a friction term
            !! represented as (mass x gamma x dr/dt). [See eq. 3.77 in the manual]
            !! Thus this mass-divided friction (gamma) will be multiplied back by mass
            !! in the velocity-Verlet procedure.
        enddo

#ifdef MPI_PAR
        if (myrank == 0) then
#endif
            write (lunout, *) "Friction coefficients by Stokes' law (i_fric = 11)"
            write (lunout, *) '   viscosity: ', inpara%viscosity, ' [(Pa)(s)]'
            write (lunout, *) '            = ', visc, ' [(Da)(A^-1)(tau^-1)]'
            do imp = 1, nmp_all
                write (lunout, *) '   friction coefficient of (', imp, '): ', fric_mp(imp)*cmass_mp(imp)
            enddo
#ifdef MPI_PAR
        end if
#endif
    endif

    ! -----------------------------------------------------------------------
    ! Re-define mass ans friction
    if (inmisc%i_redef_mass_fric == 1) then
#ifdef MPI_PAR
        if (myrank == 0) then
#endif
            ! -----------------------------------------------------------------------
            ! read from input file
            rewind (luninp)
            call ukoto_uiread2(luninp, lunout, 'mass_friction   ', kfind, &
                               CARRAY_MXLINE, nlines, cwkinp)
            if (kfind /= 'FIND') then
                error_message = 'Error: cannot find "mass_friction" field in the input file'
                call util_error(ERROR%STOP_ALL, error_message)
            end if

            do iline = 1, nlines
                call ukoto_uiequa2(cwkinp(iline), nequat, csides)
                do iequa = 1, nequat
                    char02 = csides(1, iequa)

                    cvalue = 'rmass_all'
                    if (char02(1:9) == cvalue) then
                        write (6, *) ' setp_mass_fric:rmass_all '

                        call ukoto_rvalue2(lunout, csides(1, iequa), &
                                           rmass, cvalue)
                        cmass_mp(1:MXMP) = rmass
                    end if

                    cvalue = 'fric_all'
                    if (char02(1:8) == cvalue) then
                        call ukoto_rvalue2(lunout, csides(1, iequa), &
                                           fric, cvalue)
                        fric_mp(1:MXMP) = fric
                    end if
                end do
            end do

            do iline = 1, nlines
                call ukoto_uiequa2(cwkinp(iline), nequat, csides)
                do iequa = 1, nequat
                    char02 = csides(1, iequa)

                    if (char02(1:10) == 'rmass_unit') then
!
                        write (6, *) ' setp_mass_fric:rmass_unit '
!
                        do icol = 1, CARRAY_MXCOLM
                            if (char02(icol:icol) == ')') exit
                        end do
                        read (char02(12:icol - 1), *) char12
                        call util_unitstate(char12, inunit, instate)

                        cvalue = char02(1:icol)
                        call ukoto_rvalue2(lunout, csides(1, iequa), &
                                           rmass, cvalue)
                        do imp = lunit2mp(1, inunit(1)), lunit2mp(2, inunit(2))
                            cmass_mp(imp) = rmass
                        end do

                    else if (char02(1:9) == 'fric_unit') then
                        do icol = 1, CARRAY_MXCOLM
                            if (char02(icol:icol) == ')') exit
                        end do
                        read (char02(11:icol - 1), *) char12
                        call util_unitstate(char12, inunit, instate)

                        cvalue = char02(1:icol)
                        call ukoto_rvalue2(lunout, csides(1, iequa), &
                                           fric, cvalue)
                        do imp = lunit2mp(1, inunit(1)), lunit2mp(2, inunit(2))
                            fric_mp(imp) = fric
                        end do

                    end if
                end do
            end do

            do iline = 1, nlines
                call ukoto_uiequa2(cwkinp(iline), nequat, csides)
                do iequa = 1, nequat
                    char02 = csides(1, iequa)

                    if (char02(1:8) == 'rmass_mp') then
!
                        write (6, *) ' setp_mass_fric:rmass_mp '
!
                        do icol = 1, CARRAY_MXCOLM
                            if (char02(icol:icol) == ')') exit
                        end do
                        read (char02(10:icol - 1), *) char12
                        call util_unitstate(char12, inunit, instate)

                        cvalue = char02(1:icol)
                        call ukoto_rvalue2(lunout, csides(1, iequa), &
                                           rmass, cvalue)
                        do imp = inunit(1), inunit(2)
                            cmass_mp(imp) = rmass
                        end do

                    else if (char02(1:7) == 'fric_mp') then
                        do icol = 1, CARRAY_MXCOLM
                            if (char02(icol:icol) == ')') exit
                        end do
                        read (char02(9:icol - 1), *) char12
                        call util_unitstate(char12, inunit, instate)

                        cvalue = char02(1:icol)
                        call ukoto_rvalue2(lunout, csides(1, iequa), &
                                           fric, cvalue)
                        do imp = inunit(1), inunit(2)
                            fric_mp(imp) = fric
                        end do

                    end if
                end do
            end do

#ifdef MPI_PAR
        end if
        call MPI_Bcast(cmass_mp, MXMP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(fric_mp, MXMP, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
#endif

    endif

#ifdef _DEBUG
    do imp = 1, MXMP
        write (6, *) imp, cmass_mp(imp), fric_mp(imp)
    enddo
#endif

!====================================================================
contains

    integer function imp2chemicaltype(imp)
        use const_index
        use var_struct, only: cmp2atom, cmp2seq, iclass_mp
        integer, intent(in) :: imp
        character(CARRAY_MSG_ERROR) :: error_message

        ! ------ DNA ------
        if (cmp2atom(imp) == ' O  ' .and. iclass_mp(imp) == CLASS%DNA) then
            imp2chemicaltype = CHEMICALTYPE%DP
        else if (cmp2atom(imp) == ' S' .and. iclass_mp(imp) == CLASS%DNA) then
            imp2chemicaltype = CHEMICALTYPE%DS
        else if (cmp2seq(imp) == ' DA' .and. iclass_mp(imp) == CLASS%DNA) then
            imp2chemicaltype = CHEMICALTYPE%DA
        else if (cmp2seq(imp) == ' DG' .and. iclass_mp(imp) == CLASS%DNA) then
            imp2chemicaltype = CHEMICALTYPE%DG
        else if (cmp2seq(imp) == ' DC' .and. iclass_mp(imp) == CLASS%DNA) then
            imp2chemicaltype = CHEMICALTYPE%DC
        else if (cmp2seq(imp) == ' DT' .and. iclass_mp(imp) == CLASS%DNA) then
            imp2chemicaltype = CHEMICALTYPE%DT
            ! ------ DNA2 ------
        else if (cmp2atom(imp) == 'DP  ') then
            imp2chemicaltype = CHEMICALTYPE%DP
        else if (cmp2atom(imp) == 'DS  ') then
            imp2chemicaltype = CHEMICALTYPE%DS
        else if (cmp2seq(imp) == 'DA ') then
            imp2chemicaltype = CHEMICALTYPE%DA
        else if (cmp2seq(imp) == 'DG ') then
            imp2chemicaltype = CHEMICALTYPE%DG
        else if (cmp2seq(imp) == 'DC ') then
            imp2chemicaltype = CHEMICALTYPE%DC
        else if (cmp2seq(imp) == 'DT ') then
            imp2chemicaltype = CHEMICALTYPE%DT
            ! ------ RNA ------
        else if (cmp2atom(imp) == ' P  ') then
            imp2chemicaltype = CHEMICALTYPE%P
        else if (cmp2atom(imp) == ' S  ') then
            imp2chemicaltype = CHEMICALTYPE%S
        else if (cmp2atom(imp) == ' Ab ') then
            imp2chemicaltype = CHEMICALTYPE%A
        else if (cmp2atom(imp) == ' Gb ') then
            imp2chemicaltype = CHEMICALTYPE%G
        else if (cmp2atom(imp) == ' Ub ') then
            imp2chemicaltype = CHEMICALTYPE%U
        else if (cmp2atom(imp) == ' Cb ') then
            imp2chemicaltype = CHEMICALTYPE%C
            ! ------ protein ------
        else if (cmp2seq(imp) == 'ALA') then
            imp2chemicaltype = CHEMICALTYPE%ALA
        else if (cmp2seq(imp) == 'ARG') then
            imp2chemicaltype = CHEMICALTYPE%ARG
        else if (cmp2seq(imp) == 'ASN') then
            imp2chemicaltype = CHEMICALTYPE%ASN
        else if (cmp2seq(imp) == 'ASP') then
            imp2chemicaltype = CHEMICALTYPE%ASP
        else if (cmp2seq(imp) == 'CYS') then
            imp2chemicaltype = CHEMICALTYPE%CYS
        else if (cmp2seq(imp) == 'GLN') then
            imp2chemicaltype = CHEMICALTYPE%GLN
        else if (cmp2seq(imp) == 'GLU') then
            imp2chemicaltype = CHEMICALTYPE%GLU
        else if (cmp2seq(imp) == 'GLY') then
            imp2chemicaltype = CHEMICALTYPE%GLY
        else if (cmp2seq(imp) == 'HIS') then
            imp2chemicaltype = CHEMICALTYPE%HIS
        else if (cmp2seq(imp) == 'ILE') then
            imp2chemicaltype = CHEMICALTYPE%ILE
        else if (cmp2seq(imp) == 'LEU') then
            imp2chemicaltype = CHEMICALTYPE%LEU
        else if (cmp2seq(imp) == 'LYS') then
            imp2chemicaltype = CHEMICALTYPE%LYS
        else if (cmp2seq(imp) == 'MET') then
            imp2chemicaltype = CHEMICALTYPE%MET
        else if (cmp2seq(imp) == 'PHE') then
            imp2chemicaltype = CHEMICALTYPE%PHE
        else if (cmp2seq(imp) == 'PRO') then
            imp2chemicaltype = CHEMICALTYPE%PRO
        else if (cmp2seq(imp) == 'SER') then
            imp2chemicaltype = CHEMICALTYPE%SER
        else if (cmp2seq(imp) == 'THR') then
            imp2chemicaltype = CHEMICALTYPE%THR
        else if (cmp2seq(imp) == 'TRP') then
            imp2chemicaltype = CHEMICALTYPE%TRP
        else if (cmp2seq(imp) == 'TYR') then
            imp2chemicaltype = CHEMICALTYPE%TYR
        else if (cmp2seq(imp) == 'VAL') then
            imp2chemicaltype = CHEMICALTYPE%VAL
        else if (cmp2seq(imp) == 'OTH') then
            imp2chemicaltype = CHEMICALTYPE%OTH
        else
            write (error_message, *) 'Error: unknown chemical type in setp_mass_fric; imp= ', imp
            call util_error(ERROR%STOP_ALL, error_message)
        endif
    endfunction imp2chemicaltype

end subroutine setp_mass_fric
