! setp_gb
!> @brief This subroutine is to read and set the parameters for calculating gb energy.

subroutine setp_gb()

    use const_maxsize
    use const_index
    use const_physical
    use var_inp, only: infile, outfile
    use var_struct, only: nmp_all, imp2type, cmp2seq
    use var_setp, only: ingb

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! -----------------------------------------------------------------------
    ! local variables
    integer :: imp, iline, nlines, iequa, nequat, luninp, lunout
    character(4) :: kfind
    character(CARRAY_MXCOLM) :: ctmp00
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)

    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    luninp = infile%inp
    lunout = outfile%data

    ! -----------------------------------------------------------------------
    ! read from input file
    rewind (luninp)
    call ukoto_uiread2(luninp, lunout, 'generalized_born', kfind, &
                       CARRAY_MXLINE, nlines, cwkinp)
    if (kfind /= 'FIND') then
        nlines = 0
    end if

    do iline = 1, nlines
        call ukoto_uiequa2(cwkinp(iline), nequat, csides)
        ctmp00 = cwkinp(iline)

        do iequa = 1, nequat
            cvalue = 'diele_mol'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                               ingb%diele_mol, cvalue)
            cvalue = 'cutoff_ele'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                               ingb%cutoff_ele, cvalue)
            cvalue = 'cutoff_gb'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                               ingb%cutoff_gb, cvalue)
            cvalue = 'pzpro'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                               ingb%pzpro, cvalue)
            cvalue = 'pznuc'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                               ingb%pznuc, cvalue)
        end do
    end do

    ingb%coef = JOUL2KCAL_MOL*1.0e10_PREC*ELE**2/(4.0e0_PREC*F_PI*EPSI_0)
    !write(*,*) ingb%coef

    ingb%pzgen = 2.260935
    ingb%para(1) = -2.771080
    ingb%para(2) = -8.808023*2/ingb%coef
    ingb%para(3) = 8.228746*2/ingb%coef
    ingb%para(4) = -63.50013*2/ingb%coef
    ingb%para(5) = 1.99

    do imp = 1, nmp_all
        if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'ALA') then
            ingb%rvdw(imp) = 3.201785
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'ARG') then
            ingb%rvdw(imp) = 4.241566
            ingb%charge(imp) = 1
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'ASN') then
            ingb%rvdw(imp) = 3.609638
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'ASP') then
            ingb%rvdw(imp) = 3.483978
            ingb%charge(imp) = -1
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'CYS') then
            ingb%rvdw(imp) = 3.331317
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'GLN') then
            ingb%rvdw(imp) = 3.833561
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'GLU') then
            ingb%rvdw(imp) = 3.722817
            ingb%charge(imp) = -1
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'GLY') then
            ingb%rvdw(imp) = 2.864972
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. (cmp2seq(imp) == 'HIS' &
                                           .or. cmp2seq(imp) == 'HIE' .or. cmp2seq(imp) == 'HID' &
                                           .or. cmp2seq(imp) == 'HIP' .or. cmp2seq(imp) == 'HSE' &
                                           .or. cmp2seq(imp) == 'HSD' .or. cmp2seq(imp) == 'HSP')) then
            ingb%rvdw(imp) = 3.949375
            ingb%charge(imp) = 1
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'ILE') then
            ingb%rvdw(imp) = 3.930980
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'LEU') then
            ingb%rvdw(imp) = 3.930980
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'LYS') then
            ingb%rvdw(imp) = 4.098825
            ingb%charge(imp) = 1
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'MET') then
            ingb%rvdw(imp) = 3.854543
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'PHE') then
            ingb%rvdw(imp) = 4.142545
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'PRO') then
            ingb%rvdw(imp) = 3.609923
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'SER') then
            ingb%rvdw(imp) = 3.255820
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'THR') then
            ingb%rvdw(imp) = 3.525524
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'TRP') then
            ingb%rvdw(imp) = 4.428257
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'TYR') then
            ingb%rvdw(imp) = 4.175115
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 1 .and. cmp2seq(imp) == 'VAL') then
            ingb%rvdw(imp) = 3.718989
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 2 .or. imp2type(imp) == 5 .or. imp2type(imp) == 8) then
            ingb%rvdw(imp) = 2.706779
            ingb%charge(imp) = -1
        else if (imp2type(imp) == 3 .and. cmp2seq(imp) == ' DA') then
            ingb%rvdw(imp) = 3.731850
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 3 .and. cmp2seq(imp) == ' DC') then
            ingb%rvdw(imp) = 3.503035
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 3 .and. cmp2seq(imp) == ' DG') then
            ingb%rvdw(imp) = 3.810942
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 3 .and. cmp2seq(imp) == ' DT') then
            ingb%rvdw(imp) = 3.677826
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 6 .and. cmp2seq(imp) == ' DA') then
            ingb%rvdw(imp) = 3.731850
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 6 .and. cmp2seq(imp) == ' DC') then
            ingb%rvdw(imp) = 3.503035
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 6 .and. cmp2seq(imp) == ' DG') then
            ingb%rvdw(imp) = 3.810942
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 6 .and. cmp2seq(imp) == ' DT') then
            ingb%rvdw(imp) = 3.432467
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 9 .and. cmp2seq(imp) == ' RA') then
            ingb%rvdw(imp) = 3.731850
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 9 .and. cmp2seq(imp) == ' RC') then
            ingb%rvdw(imp) = 3.503035
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 9 .and. cmp2seq(imp) == ' RG') then
            ingb%rvdw(imp) = 3.810942
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 9 .and. cmp2seq(imp) == ' RU') then
            ingb%rvdw(imp) = 3.432467
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 4 .or. imp2type(imp) == 7) then
            ingb%rvdw(imp) = 3.512040
            ingb%charge(imp) = 0
        else if (imp2type(imp) == 10) then
            ingb%rvdw(imp) = 3.557130
            ingb%charge(imp) = 0
        else
            ingb%rvdw(imp) = 0
            ingb%charge(imp) = 0
        end if
    end do

end subroutine setp_gb
