!write_top
!> @brief Outputs structure information with TOP style.

subroutine write_top()

    use, intrinsic :: iso_fortran_env
    use const_index
    use const_maxsize
    use const_physical
    use var_inp, only: outfile, path
    use var_setp, only: inmisc, inpara, inexv, indna2, inflp
    use var_struct, only: nunit_real, lunit2mp, ires_mp, cmp2seq, cmp2atom, &
        nba, iba2mp, ba_nat, coef_ba, ibd2mp, bd_nat, coef_bd, nbd, imp2unit, &
        iclass_unit, ndih, idih2mp, coef_dih_gauss, dih_nat, wid_dih_gauss, &
        cmass_mp, coef_aicg13_gauss, wid_aicg13_gauss, aicg13_nat, nfdih, &
        ifdih2mp, fdih_para, ncon, icon2mp, coef_go, go_nat
    use mod_unit
#ifdef MPI_PAR
    use mpiconst
#endif
    implicit none

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ofile
    integer :: imp, imp1, imp2, imp3, imp4, iunit, ibd, iba, idih, icon, iaa, &
               ipara, mp_offset, res_offset, ires, func, cgnr, itp_out
    character(len=3)   :: aa_name(20)
    character(len=8)   :: ptype = "       A"
    character(len=16)  :: mol_name(nunit_real)
    character(len=256) :: filename

    ! EXV sigma
    real(PREC) :: da_s ! DNA A sigma (nm)
    real(PREC) :: dg_s ! DNA G sigma (nm)
    real(PREC) :: dc_s ! DNA C sigma (nm)
    real(PREC) :: dt_s ! DNA T sigma (nm)
    real(PREC) :: dp_s ! DNA P sigma (nm)
    real(PREC) :: ds_s ! DNA S sigma (nm)
    ! EXV epsilon
    real(PREC) :: exv_e ! DNA EXV epsilon in kJ/mol
    ! EXV coefficients
    real(PREC) :: da_c6 ! DNA A EXV c6 coefficient
    real(PREC) :: dg_c6 ! DNA G EXV c6 coefficient
    real(PREC) :: dc_c6 ! DNA C EXV c6 coefficient
    real(PREC) :: dt_c6 ! DNA T EXV c6 coefficient
    real(PREC) :: dp_c6 ! DNA P EXV c6 coefficient
    real(PREC) :: ds_c6 ! DNA S EXV c6 coefficient
    real(PREC) :: da_c12 ! DNA A EXV c12 coefficient
    real(PREC) :: dg_c12 ! DNA G EXV c12 coefficient
    real(PREC) :: dc_c12 ! DNA C EXV c12 coefficient
    real(PREC) :: dt_c12 ! DNA T EXV c12 coefficient
    real(PREC) :: dp_c12 ! DNA P EXV c12 coefficient
    real(PREC) :: ds_c12 ! DNA S EXV c12 coefficient
    real(PREC) :: ala_c12 ! ALA EXV c12 coefficient
    real(PREC) :: arg_c12 ! ARG EXV c12 coefficient
    real(PREC) :: asn_c12 ! ASN EXV c12 coefficient
    real(PREC) :: asp_c12 ! ASP EXV c12 coefficient
    real(PREC) :: cys_c12 ! CYS EXV c12 coefficient
    real(PREC) :: gln_c12 ! GLN EXV c12 coefficient
    real(PREC) :: glu_c12 ! GLU EXV c12 coefficient
    real(PREC) :: gly_c12 ! GLY EXV c12 coefficient
    real(PREC) :: his_c12 ! HIS EXV c12 coefficient
    real(PREC) :: ile_c12 ! ILE EXV c12 coefficient
    real(PREC) :: leu_c12 ! LEU EXV c12 coefficient
    real(PREC) :: lys_c12 ! LYS EXV c12 coefficient
    real(PREC) :: met_c12 ! MET EXV c12 coefficient
    real(PREC) :: phe_c12 ! PHE EXV c12 coefficient
    real(PREC) :: pro_c12 ! PRO EXV c12 coefficient
    real(PREC) :: ser_c12 ! SER EXV c12 coefficient
    real(PREC) :: thr_c12 ! THR EXV c12 coefficient
    real(PREC) :: trp_c12 ! TRP EXV c12 coefficient
    real(PREC) :: tyr_c12 ! TYR EXV c12 coefficient
    real(PREC) :: val_c12 ! VAL EXV c12 coefficient
    real(PREC) :: ala_da_c12 ! ALA-DA EXV c12 coefficient
    real(PREC) :: ala_dg_c12 ! ALA-DG EXV c12 coefficient
    real(PREC) :: ala_dc_c12 ! ALA-DC EXV c12 coefficient
    real(PREC) :: ala_dt_c12 ! ALA-DT EXV c12 coefficient
    real(PREC) :: arg_da_c12 ! ARG-DA EXV c12 coefficient
    real(PREC) :: arg_dg_c12 ! ARG-DG EXV c12 coefficient
    real(PREC) :: arg_dc_c12 ! ARG-DC EXV c12 coefficient
    real(PREC) :: arg_dt_c12 ! ARG-DT EXV c12 coefficient
    real(PREC) :: asn_da_c12 ! ASN-DA EXV c12 coefficient
    real(PREC) :: asn_dg_c12 ! ASN-DG EXV c12 coefficient
    real(PREC) :: asn_dc_c12 ! ASN-DC EXV c12 coefficient
    real(PREC) :: asn_dt_c12 ! ASN-DT EXV c12 coefficient
    real(PREC) :: asp_da_c12 ! ASP-DA EXV c12 coefficient
    real(PREC) :: asp_dg_c12 ! ASP-DG EXV c12 coefficient
    real(PREC) :: asp_dc_c12 ! ASP-DC EXV c12 coefficient
    real(PREC) :: asp_dt_c12 ! ASP-DT EXV c12 coefficient
    real(PREC) :: cys_da_c12 ! CYS-DA EXV c12 coefficient
    real(PREC) :: cys_dg_c12 ! CYS-DG EXV c12 coefficient
    real(PREC) :: cys_dc_c12 ! CYS-DC EXV c12 coefficient
    real(PREC) :: cys_dt_c12 ! CYS-DT EXV c12 coefficient
    real(PREC) :: gln_da_c12 ! GLN-DA EXV c12 coefficient
    real(PREC) :: gln_dg_c12 ! GLN-DG EXV c12 coefficient
    real(PREC) :: gln_dc_c12 ! GLN-DC EXV c12 coefficient
    real(PREC) :: gln_dt_c12 ! GLN-DT EXV c12 coefficient
    real(PREC) :: glu_da_c12 ! GLU-DA EXV c12 coefficient
    real(PREC) :: glu_dg_c12 ! GLU-DG EXV c12 coefficient
    real(PREC) :: glu_dc_c12 ! GLU-DC EXV c12 coefficient
    real(PREC) :: glu_dt_c12 ! GLU-DT EXV c12 coefficient
    real(PREC) :: gly_da_c12 ! GLY-DA EXV c12 coefficient
    real(PREC) :: gly_dg_c12 ! GLY-DG EXV c12 coefficient
    real(PREC) :: gly_dc_c12 ! GLY-DC EXV c12 coefficient
    real(PREC) :: gly_dt_c12 ! GLY-DT EXV c12 coefficient
    real(PREC) :: his_da_c12 ! HIS-DA EXV c12 coefficient
    real(PREC) :: his_dg_c12 ! HIS-DG EXV c12 coefficient
    real(PREC) :: his_dc_c12 ! HIS-DC EXV c12 coefficient
    real(PREC) :: his_dt_c12 ! HIS-DT EXV c12 coefficient
    real(PREC) :: ile_da_c12 ! ILE-DA EXV c12 coefficient
    real(PREC) :: ile_dg_c12 ! ILE-DG EXV c12 coefficient
    real(PREC) :: ile_dc_c12 ! ILE-DC EXV c12 coefficient
    real(PREC) :: ile_dt_c12 ! ILE-DT EXV c12 coefficient
    real(PREC) :: leu_da_c12 ! LEU-DA EXV c12 coefficient
    real(PREC) :: leu_dg_c12 ! LEU-DG EXV c12 coefficient
    real(PREC) :: leu_dc_c12 ! LEU-DC EXV c12 coefficient
    real(PREC) :: leu_dt_c12 ! LEU-DT EXV c12 coefficient
    real(PREC) :: lys_da_c12 ! LYS-DA EXV c12 coefficient
    real(PREC) :: lys_dg_c12 ! LYS-DG EXV c12 coefficient
    real(PREC) :: lys_dc_c12 ! LYS-DC EXV c12 coefficient
    real(PREC) :: lys_dt_c12 ! LYS-DT EXV c12 coefficient
    real(PREC) :: met_da_c12 ! MET-DA EXV c12 coefficient
    real(PREC) :: met_dg_c12 ! MET-DG EXV c12 coefficient
    real(PREC) :: met_dc_c12 ! MET-DC EXV c12 coefficient
    real(PREC) :: met_dt_c12 ! MET-DT EXV c12 coefficient
    real(PREC) :: phe_da_c12 ! PHE-DA EXV c12 coefficient
    real(PREC) :: phe_dg_c12 ! PHE-DG EXV c12 coefficient
    real(PREC) :: phe_dc_c12 ! PHE-DC EXV c12 coefficient
    real(PREC) :: phe_dt_c12 ! PHE-DT EXV c12 coefficient
    real(PREC) :: pro_da_c12 ! PRO-DA EXV c12 coefficient
    real(PREC) :: pro_dg_c12 ! PRO-DG EXV c12 coefficient
    real(PREC) :: pro_dc_c12 ! PRO-DC EXV c12 coefficient
    real(PREC) :: pro_dt_c12 ! PRO-DT EXV c12 coefficient
    real(PREC) :: ser_da_c12 ! SER-DA EXV c12 coefficient
    real(PREC) :: ser_dg_c12 ! SER-DG EXV c12 coefficient
    real(PREC) :: ser_dc_c12 ! SER-DC EXV c12 coefficient
    real(PREC) :: ser_dt_c12 ! SER-DT EXV c12 coefficient
    real(PREC) :: thr_da_c12 ! THR-DA EXV c12 coefficient
    real(PREC) :: thr_dg_c12 ! THR-DG EXV c12 coefficient
    real(PREC) :: thr_dc_c12 ! THR-DC EXV c12 coefficient
    real(PREC) :: thr_dt_c12 ! THR-DT EXV c12 coefficient
    real(PREC) :: trp_da_c12 ! TRP-DA EXV c12 coefficient
    real(PREC) :: trp_dg_c12 ! TRP-DG EXV c12 coefficient
    real(PREC) :: trp_dc_c12 ! TRP-DC EXV c12 coefficient
    real(PREC) :: trp_dt_c12 ! TRP-DT EXV c12 coefficient
    real(PREC) :: tyr_da_c12 ! TYR-DA EXV c12 coefficient
    real(PREC) :: tyr_dg_c12 ! TYR-DG EXV c12 coefficient
    real(PREC) :: tyr_dc_c12 ! TYR-DC EXV c12 coefficient
    real(PREC) :: tyr_dt_c12 ! TYR-DT EXV c12 coefficient
    real(PREC) :: val_da_c12 ! VAL-DA EXV c12 coefficient
    real(PREC) :: val_dg_c12 ! VAL-DG EXV c12 coefficient
    real(PREC) :: val_dc_c12 ! VAL-DC EXV c12 coefficient
    real(PREC) :: val_dt_c12 ! VAL-DT EXV c12 coefficient
    ! Charge
    real(PREC) :: da_q ! DNA A Charge
    real(PREC) :: dg_q ! DNA G Charge
    real(PREC) :: dc_q ! DNA C Charge
    real(PREC) :: dt_q ! DNA T Charge
    real(PREC) :: dp_q ! DNA P Charge
    real(PREC) :: ds_q ! DNA S Charge
    ! Functions
    real(PREC) :: imp2charge

    ! Variable initialization
    da_s = indna2%sex(BASETYPE%A)/10.0_PREC
    dg_s = indna2%sex(BASETYPE%G)/10.0_PREC
    dc_s = indna2%sex(BASETYPE%C)/10.0_PREC
    dt_s = indna2%sex(BASETYPE%T)/10.0_PREC
    dp_s = indna2%sex(BASETYPE%P)/10.0_PREC
    ds_s = indna2%sex(BASETYPE%S)/10.0_PREC

    exv_e = 1.0_PREC

    da_c6 = exv_e*da_s**6
    dg_c6 = exv_e*dg_s**6
    dc_c6 = exv_e*dc_s**6
    dt_c6 = exv_e*dt_s**6
    dp_c6 = exv_e*dp_s**6
    ds_c6 = exv_e*ds_s**6

    da_c12 = exv_e*da_s**12
    dg_c12 = exv_e*dg_s**12
    dc_c12 = exv_e*dc_s**12
    dt_c12 = exv_e*dt_s**12
    dp_c12 = exv_e*dp_s**12
    ds_c12 = exv_e*ds_s**12

    da_q = 0.0_PREC
    dg_q = 0.0_PREC
    dc_q = 0.0_PREC
    dt_q = 0.0_PREC
    dp_q = -0.6_PREC
    ds_q = 0.0_PREC

    ala_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%ALA)/10.0_PREC)**12
    arg_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%ARG)/10.0_PREC)**12
    asn_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%ASN)/10.0_PREC)**12
    asp_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%ASP)/10.0_PREC)**12
    cys_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%CYS)/10.0_PREC)**12
    gln_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%GLN)/10.0_PREC)**12
    glu_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%GLU)/10.0_PREC)**12
    gly_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%GLY)/10.0_PREC)**12
    his_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%HIS)/10.0_PREC)**12
    ile_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%ILE)/10.0_PREC)**12
    leu_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%LEU)/10.0_PREC)**12
    lys_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%LYS)/10.0_PREC)**12
    met_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%MET)/10.0_PREC)**12
    phe_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%PHE)/10.0_PREC)**12
    pro_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%PRO)/10.0_PREC)**12
    ser_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%SER)/10.0_PREC)**12
    thr_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%THR)/10.0_PREC)**12
    trp_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%TRP)/10.0_PREC)**12
    tyr_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%TYR)/10.0_PREC)**12
    val_c12 = 4.184_PREC*inexv%exv_coef*(inexv%exv_sigma(CHEMICALTYPE%VAL)/10.0_PREC)**12

    ala_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%ALA))/20.0_PREC)**12
    ala_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%ALA))/20.0_PREC)**12
    ala_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%ALA))/20.0_PREC)**12
    ala_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%ALA))/20.0_PREC)**12
    arg_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%ARG))/20.0_PREC)**12
    arg_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%ARG))/20.0_PREC)**12
    arg_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%ARG))/20.0_PREC)**12
    arg_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%ARG))/20.0_PREC)**12
    asn_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%ASN))/20.0_PREC)**12
    asn_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%ASN))/20.0_PREC)**12
    asn_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%ASN))/20.0_PREC)**12
    asn_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%ASN))/20.0_PREC)**12
    asp_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%ASP))/20.0_PREC)**12
    asp_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%ASP))/20.0_PREC)**12
    asp_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%ASP))/20.0_PREC)**12
    asp_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%ASP))/20.0_PREC)**12
    cys_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%CYS))/20.0_PREC)**12
    cys_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%CYS))/20.0_PREC)**12
    cys_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%CYS))/20.0_PREC)**12
    cys_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%CYS))/20.0_PREC)**12
    gln_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%GLN))/20.0_PREC)**12
    gln_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%GLN))/20.0_PREC)**12
    gln_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%GLN))/20.0_PREC)**12
    gln_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%GLN))/20.0_PREC)**12
    glu_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%GLU))/20.0_PREC)**12
    glu_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%GLU))/20.0_PREC)**12
    glu_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%GLU))/20.0_PREC)**12
    glu_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%GLU))/20.0_PREC)**12
    gly_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%GLY))/20.0_PREC)**12
    gly_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%GLY))/20.0_PREC)**12
    gly_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%GLY))/20.0_PREC)**12
    gly_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%GLY))/20.0_PREC)**12
    his_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%HIS))/20.0_PREC)**12
    his_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%HIS))/20.0_PREC)**12
    his_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%HIS))/20.0_PREC)**12
    his_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%HIS))/20.0_PREC)**12
    ile_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%ILE))/20.0_PREC)**12
    ile_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%ILE))/20.0_PREC)**12
    ile_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%ILE))/20.0_PREC)**12
    ile_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%ILE))/20.0_PREC)**12
    leu_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%LEU))/20.0_PREC)**12
    leu_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%LEU))/20.0_PREC)**12
    leu_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%LEU))/20.0_PREC)**12
    leu_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%LEU))/20.0_PREC)**12
    lys_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%LYS))/20.0_PREC)**12
    lys_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%LYS))/20.0_PREC)**12
    lys_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%LYS))/20.0_PREC)**12
    lys_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%LYS))/20.0_PREC)**12
    met_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%MET))/20.0_PREC)**12
    met_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%MET))/20.0_PREC)**12
    met_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%MET))/20.0_PREC)**12
    met_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%MET))/20.0_PREC)**12
    phe_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%PHE))/20.0_PREC)**12
    phe_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%PHE))/20.0_PREC)**12
    phe_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%PHE))/20.0_PREC)**12
    phe_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%PHE))/20.0_PREC)**12
    pro_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%PRO))/20.0_PREC)**12
    pro_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%PRO))/20.0_PREC)**12
    pro_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%PRO))/20.0_PREC)**12
    pro_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%PRO))/20.0_PREC)**12
    ser_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%SER))/20.0_PREC)**12
    ser_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%SER))/20.0_PREC)**12
    ser_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%SER))/20.0_PREC)**12
    ser_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%SER))/20.0_PREC)**12
    thr_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%THR))/20.0_PREC)**12
    thr_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%THR))/20.0_PREC)**12
    thr_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%THR))/20.0_PREC)**12
    thr_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%THR))/20.0_PREC)**12
    trp_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%TRP))/20.0_PREC)**12
    trp_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%TRP))/20.0_PREC)**12
    trp_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%TRP))/20.0_PREC)**12
    trp_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%TRP))/20.0_PREC)**12
    tyr_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%TYR))/20.0_PREC)**12
    tyr_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%TYR))/20.0_PREC)**12
    tyr_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%TYR))/20.0_PREC)**12
    tyr_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%TYR))/20.0_PREC)**12
    val_da_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DA) + inexv%exv_sigma(CHEMICALTYPE%VAL))/20.0_PREC)**12
    val_dg_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DG) + inexv%exv_sigma(CHEMICALTYPE%VAL))/20.0_PREC)**12
    val_dc_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DC) + inexv%exv_sigma(CHEMICALTYPE%VAL))/20.0_PREC)**12
    val_dt_c12 = 4.184_PREC*inexv%exv_coef*((inexv%exv_sigma(CHEMICALTYPE%DT) + inexv%exv_sigma(CHEMICALTYPE%VAL))/20.0_PREC)**12

    aa_name = reshape((/"ALA", "ARG", "ASN", "ASP", "CYS", &
                        "GLN", "GLU", "GLY", "HIS", "ILE", &
                        "LEU", "LYS", "MET", "PHE", "PRO", &
                        "SER", "THR", "TRP", "TYR", "VAL"/), shape(aa_name))
    ires = -1

    !///////////////////////////////////////////////////////////////////////////////////////////////
    !-----------------------------------------------------------------------------------------------
    ! File: atom_types.itp
    !-----------------------------------------------------------------------------------------------
    ! Get itp file name
    if (path /= "") then
        filename = trim(adjustl(path))//"/atom_types.itp"
    else
        filename = "./atom_types.itp"
    end if

    ! Open itp file
    open (new_unit(itp_out), file=trim(adjustl(filename)), &
          status="replace", action="write")

    ! Write defaults section
    write (itp_out, "(A)") "[ defaults ]"
    write (itp_out, "(2(I5),A5,2(F5.1))") 1, 1, "   no", 1.0, 1.0
    write (itp_out, "(A)") ""

    ! Write atom type section
    write (itp_out, "(A)") "[ atomtypes ]"

    if (inmisc%class_flag(CLASS%DNA)) then
        write (output_unit, "(A)") "Error: Can't output TOP file for 3SPN model"
        stop
    end if

    if (inmisc%class_flag(CLASS%DNA2)) then
        write (itp_out, "(A)") ";DNA"
        write (itp_out, "(A)") &
            ";name  at.num    mass  charge   ptype         V(c6)        W(c12)"
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "DA", 1, inpara%cmass(CHEMICALTYPE%DA), da_q, ptype, da_c6, da_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "DG", 1, inpara%cmass(CHEMICALTYPE%DG), dg_q, ptype, dg_c6, dg_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "DC", 1, inpara%cmass(CHEMICALTYPE%DC), dc_q, ptype, dc_c6, dc_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "DT", 1, inpara%cmass(CHEMICALTYPE%DT), dt_q, ptype, dt_c6, dt_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "DP", 1, inpara%cmass(CHEMICALTYPE%DP), dp_q, ptype, dp_c6, dp_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "DS", 1, inpara%cmass(CHEMICALTYPE%DS), ds_q, ptype, ds_c6, ds_c12
    end if

    if (inmisc%class_flag(CLASS%PRO)) then
        write (itp_out, "(A)") ";PROTEIN"
        write (itp_out, "(A)") &
            ";name  at.num    mass  charge   ptype         V(c6)        W(c12)"
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  ALA", 6, inpara%cmass(CHEMICALTYPE%ALA), 0.0, ptype, 0.0, ala_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  ARG", 6, inpara%cmass(CHEMICALTYPE%ARG), 0.0, ptype, 0.0, arg_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  ASN", 6, inpara%cmass(CHEMICALTYPE%ASN), 0.0, ptype, 0.0, asn_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  ASP", 6, inpara%cmass(CHEMICALTYPE%ASP), 0.0, ptype, 0.0, asp_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  CYS", 6, inpara%cmass(CHEMICALTYPE%CYS), 0.0, ptype, 0.0, cys_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  GLN", 6, inpara%cmass(CHEMICALTYPE%GLN), 0.0, ptype, 0.0, gln_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  GLU", 6, inpara%cmass(CHEMICALTYPE%GLU), 0.0, ptype, 0.0, glu_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  GLY", 6, inpara%cmass(CHEMICALTYPE%GLY), 0.0, ptype, 0.0, gly_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  HIS", 6, inpara%cmass(CHEMICALTYPE%HIS), 0.0, ptype, 0.0, his_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  ILE", 6, inpara%cmass(CHEMICALTYPE%ILE), 0.0, ptype, 0.0, ile_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  LEU", 6, inpara%cmass(CHEMICALTYPE%LEU), 0.0, ptype, 0.0, leu_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  LYS", 6, inpara%cmass(CHEMICALTYPE%LYS), 0.0, ptype, 0.0, lys_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  MET", 7, inpara%cmass(CHEMICALTYPE%MET), 0.0, ptype, 0.0, met_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  PHE", 6, inpara%cmass(CHEMICALTYPE%PHE), 0.0, ptype, 0.0, phe_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  PRO", 6, inpara%cmass(CHEMICALTYPE%PRO), 0.0, ptype, 0.0, pro_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  SER", 6, inpara%cmass(CHEMICALTYPE%SER), 0.0, ptype, 0.0, ser_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  THR", 6, inpara%cmass(CHEMICALTYPE%THR), 0.0, ptype, 0.0, thr_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  TRP", 6, inpara%cmass(CHEMICALTYPE%TRP), 0.0, ptype, 0.0, trp_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  TYR", 6, inpara%cmass(CHEMICALTYPE%TYR), 0.0, ptype, 0.0, tyr_c12
        write (itp_out, "(A5,I8,2(F8.3),A8,2(ES14.4))") &
            "  VAL", 6, inpara%cmass(CHEMICALTYPE%VAL), 0.0, ptype, 0.0, val_c12
        write (itp_out, "(A)") ""
    end if
    write (itp_out, "(A)") ""

    ! Write pair type section
    if (inmisc%class_flag(CLASS%DNA)) then
        write (output_unit, "(A)") "Error: Can't output TOP file for 3SPN model"
        stop
    end if
    if (inmisc%class_flag(CLASS%DNA2)) then
        write (itp_out, "(A)") "[ pairtypes ]"
        write (itp_out, "(A)") ";DNA - Base Pairing"
        write (itp_out, "(A)") &
            ";   i    j func epsilon   alpha      eq     phi thetha1 thetha2       K"
        func = 1
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DT", func, indna2%ebp_AT*4.184_PREC, indna2%abp, indna2%sbp_AT/10.0_PREC, &
            indna2%pbp_AT, indna2%t1bp_AT, indna2%t2bp_AT, indna2%kbp
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DA", func, indna2%ebp_AT*4.184_PREC, indna2%abp, indna2%sbp_AT/10.0_PREC, &
            indna2%pbp_AT, indna2%t1bp_TA, indna2%t2bp_TA, indna2%kbp
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DC", func, indna2%ebp_CG*4.184_PREC, indna2%abp, indna2%sbp_CG/10.0_PREC, &
            indna2%pbp_CG, indna2%t1bp_GC, indna2%t2bp_GC, indna2%kbp
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DG", func, indna2%ebp_CG*4.184_PREC, indna2%abp, indna2%sbp_CG/10.0_PREC, &
            indna2%pbp_CG, indna2%t1bp_CG, indna2%t2bp_CG, indna2%kbp
        write (itp_out, "(A)") ";DNA - Base Stacking up5'-up3'"
        write (itp_out, "(A)") &
            ";   i    j func epsilon   alpha      eq  thetha       K"
        func = 2
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DA", "DA", func, indna2%ebstk(BPTYPE%AA)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%AA)/10.0_PREC, &
            indna2%tbstk(BPTYPE%AA), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DA", "DT", func, indna2%ebstk(BPTYPE%AT)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%AT)/10.0_PREC, &
            indna2%tbstk(BPTYPE%AT), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DA", "DG", func, indna2%ebstk(BPTYPE%AG)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%AG)/10.0_PREC, &
            indna2%tbstk(BPTYPE%AG), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DA", "DC", func, indna2%ebstk(BPTYPE%AC)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%AC)/10.0_PREC, &
            indna2%tbstk(BPTYPE%AC), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DT", "DA", func, indna2%ebstk(BPTYPE%TA)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%TA)/10.0_PREC, &
            indna2%tbstk(BPTYPE%TA), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DT", "DT", func, indna2%ebstk(BPTYPE%TT)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%TT)/10.0_PREC, &
            indna2%tbstk(BPTYPE%TT), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DT", "DG", func, indna2%ebstk(BPTYPE%TG)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%TG)/10.0_PREC, &
            indna2%tbstk(BPTYPE%TG), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DT", "DC", func, indna2%ebstk(BPTYPE%TC)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%TC)/10.0_PREC, &
            indna2%tbstk(BPTYPE%TC), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DG", "DA", func, indna2%ebstk(BPTYPE%GA)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%GA)/10.0_PREC, &
            indna2%tbstk(BPTYPE%GA), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DG", "DT", func, indna2%ebstk(BPTYPE%GT)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%GT)/10.0_PREC, &
            indna2%tbstk(BPTYPE%GT), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DG", "DG", func, indna2%ebstk(BPTYPE%GG)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%GG)/10.0_PREC, &
            indna2%tbstk(BPTYPE%GG), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DG", "DC", func, indna2%ebstk(BPTYPE%GC)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%GC)/10.0_PREC, &
            indna2%tbstk(BPTYPE%GC), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DC", "DA", func, indna2%ebstk(BPTYPE%CA)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%CA)/10.0_PREC, &
            indna2%tbstk(BPTYPE%CA), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DC", "DT", func, indna2%ebstk(BPTYPE%CT)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%CT)/10.0_PREC, &
            indna2%tbstk(BPTYPE%CT), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DC", "DG", func, indna2%ebstk(BPTYPE%CG)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%CG)/10.0_PREC, &
            indna2%tbstk(BPTYPE%CG), indna2%kbstk
        write (itp_out, "(2(A5),I5,5(F8.3))") &
            "DC", "DC", func, indna2%ebstk(BPTYPE%CC)*4.184_PREC, indna2%abstk, indna2%sbstk(BPTYPE%CC)/10.0_PREC, &
            indna2%tbstk(BPTYPE%CC), indna2%kbstk
        write (itp_out, "(A)") ";DNA - Base Cross-Stacking up5'-down5'"
        write (itp_out, "(A)") &
            ";   i    j func epsilon   alpha      eq thetha3     Kbp thetha4     Kcs"
        func = 3
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DA", func, indna2%ecstk1(BPTYPE%AA)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%AA)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%AA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DT", func, indna2%ecstk1(BPTYPE%AT)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%AT)/10.0_PREC, &
            indna2%t3cstk_AT, indna2%kbp, indna2%tcstk1(BPTYPE%AT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DG", func, indna2%ecstk1(BPTYPE%AG)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%AG)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%AG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DC", func, indna2%ecstk1(BPTYPE%AC)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%AC)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%AC), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DA", func, indna2%ecstk1(BPTYPE%TA)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%TA)/10.0_PREC, &
            indna2%t3cstk_AT, indna2%kbp, indna2%tcstk1(BPTYPE%TA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DT", func, indna2%ecstk1(BPTYPE%TT)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%TT)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%TT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DG", func, indna2%ecstk1(BPTYPE%TG)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%TG)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%TG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DC", func, indna2%ecstk1(BPTYPE%TC)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%TC)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%TC), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DA", func, indna2%ecstk1(BPTYPE%GA)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%GA)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%GA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DT", func, indna2%ecstk1(BPTYPE%GT)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%GT)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%GT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DG", func, indna2%ecstk1(BPTYPE%GG)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%GG)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%GG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DC", func, indna2%ecstk1(BPTYPE%GC)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%GC)/10.0_PREC, &
            indna2%t3cstk_CG, indna2%kbp, indna2%tcstk1(BPTYPE%GC), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DA", func, indna2%ecstk1(BPTYPE%CA)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%CA)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%CA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DT", func, indna2%ecstk1(BPTYPE%CT)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%CT)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%CT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DG", func, indna2%ecstk1(BPTYPE%CG)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%CG)/10.0_PREC, &
            indna2%t3cstk_CG, indna2%kbp, indna2%tcstk1(BPTYPE%CG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DC", func, indna2%ecstk1(BPTYPE%CC)*4.184_PREC, indna2%acstk, indna2%scstk1(BPTYPE%CC)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk1(BPTYPE%CC), indna2%kcstk
        write (itp_out, "(A)") ";DNA - Base Cross-Stacking down3'-up3'"
        write (itp_out, "(A)") &
            ";   i    j func epsilon   alpha      eq thetha3     Kbp thetha4     Kcs"
        func = 4
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DA", func, indna2%ecstk2(BPTYPE%AA)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%AA)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%AA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DT", func, indna2%ecstk2(BPTYPE%AT)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%AT)/10.0_PREC, &
            indna2%t3cstk_AT, indna2%kbp, indna2%tcstk2(BPTYPE%AT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DG", func, indna2%ecstk2(BPTYPE%AG)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%AG)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%AG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DA", "DC", func, indna2%ecstk2(BPTYPE%AC)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%AC)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%AC), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DA", func, indna2%ecstk2(BPTYPE%TA)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%TA)/10.0_PREC, &
            indna2%t3cstk_AT, indna2%kbp, indna2%tcstk2(BPTYPE%TA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DT", func, indna2%ecstk2(BPTYPE%TT)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%TT)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%TT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DG", func, indna2%ecstk2(BPTYPE%TG)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%TG)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%TG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DT", "DC", func, indna2%ecstk2(BPTYPE%TC)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%TC)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%TC), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DA", func, indna2%ecstk2(BPTYPE%GA)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%GA)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%GA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DT", func, indna2%ecstk2(BPTYPE%GT)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%GT)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%GT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DG", func, indna2%ecstk2(BPTYPE%GG)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%GG)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%GG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DG", "DC", func, indna2%ecstk2(BPTYPE%GC)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%GC)/10.0_PREC, &
            indna2%t3cstk_CG, indna2%kbp, indna2%tcstk2(BPTYPE%GC), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DA", func, indna2%ecstk2(BPTYPE%CA)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%CA)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%CA), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DT", func, indna2%ecstk2(BPTYPE%CT)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%CT)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%CT), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DG", func, indna2%ecstk2(BPTYPE%CG)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%CG)/10.0_PREC, &
            indna2%t3cstk_CG, indna2%kbp, indna2%tcstk2(BPTYPE%CG), indna2%kcstk
        write (itp_out, "(2(A5),I5,7(F8.3))") &
            "DC", "DC", func, indna2%ecstk2(BPTYPE%CC)*4.184_PREC, indna2%acstk, indna2%scstk2(BPTYPE%CC)/10.0_PREC, &
            0.0_PREC, indna2%kbp, indna2%tcstk2(BPTYPE%CC), indna2%kcstk
        write (itp_out, "(A)") ""
    end if

    if (inmisc%class_flag(CLASS%DNA2) .and. inmisc%class_flag(CLASS%PRO)) then
        write (itp_out, "(A)") "[ nonbond_params ] ; Protein - Non-Local - Go contacts"
        write (itp_out, "(A)") &
            ";   i    j    f          V(c6)         W(c12)"
        func = 1
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ALA", "DA", func, 0.0_PREC, ala_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ALA", "DG", func, 0.0_PREC, ala_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ALA", "DC", func, 0.0_PREC, ala_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ALA", "DT", func, 0.0_PREC, ala_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ARG", "DA", func, 0.0_PREC, arg_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ARG", "DG", func, 0.0_PREC, arg_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ARG", "DC", func, 0.0_PREC, arg_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ARG", "DT", func, 0.0_PREC, arg_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASN", "DA", func, 0.0_PREC, asn_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASN", "DG", func, 0.0_PREC, asn_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASN", "DC", func, 0.0_PREC, asn_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASN", "DT", func, 0.0_PREC, asn_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASP", "DA", func, 0.0_PREC, asp_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASP", "DG", func, 0.0_PREC, asp_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASP", "DC", func, 0.0_PREC, asp_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ASP", "DT", func, 0.0_PREC, asp_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "CYS", "DA", func, 0.0_PREC, cys_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "CYS", "DG", func, 0.0_PREC, cys_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "CYS", "DC", func, 0.0_PREC, cys_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "CYS", "DT", func, 0.0_PREC, cys_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLN", "DA", func, 0.0_PREC, gln_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLN", "DG", func, 0.0_PREC, gln_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLN", "DC", func, 0.0_PREC, gln_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLN", "DT", func, 0.0_PREC, gln_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLU", "DA", func, 0.0_PREC, glu_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLU", "DG", func, 0.0_PREC, glu_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLU", "DC", func, 0.0_PREC, glu_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLU", "DT", func, 0.0_PREC, glu_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLY", "DA", func, 0.0_PREC, gly_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLY", "DG", func, 0.0_PREC, gly_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLY", "DC", func, 0.0_PREC, gly_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "GLY", "DT", func, 0.0_PREC, gly_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "HIS", "DA", func, 0.0_PREC, his_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "HIS", "DG", func, 0.0_PREC, his_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "HIS", "DC", func, 0.0_PREC, his_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "HIS", "DT", func, 0.0_PREC, his_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ILE", "DA", func, 0.0_PREC, ile_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ILE", "DG", func, 0.0_PREC, ile_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ILE", "DC", func, 0.0_PREC, ile_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "ILE", "DT", func, 0.0_PREC, ile_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LEU", "DA", func, 0.0_PREC, leu_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LEU", "DG", func, 0.0_PREC, leu_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LEU", "DC", func, 0.0_PREC, leu_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LEU", "DT", func, 0.0_PREC, leu_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LYS", "DA", func, 0.0_PREC, lys_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LYS", "DG", func, 0.0_PREC, lys_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LYS", "DC", func, 0.0_PREC, lys_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "LYS", "DT", func, 0.0_PREC, lys_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "MET", "DA", func, 0.0_PREC, met_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "MET", "DG", func, 0.0_PREC, met_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "MET", "DC", func, 0.0_PREC, met_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "MET", "DT", func, 0.0_PREC, met_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PHE", "DA", func, 0.0_PREC, phe_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PHE", "DG", func, 0.0_PREC, phe_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PHE", "DC", func, 0.0_PREC, phe_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PHE", "DT", func, 0.0_PREC, phe_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PRO", "DA", func, 0.0_PREC, pro_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PRO", "DG", func, 0.0_PREC, pro_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PRO", "DC", func, 0.0_PREC, pro_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "PRO", "DT", func, 0.0_PREC, pro_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "SER", "DA", func, 0.0_PREC, ser_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "SER", "DG", func, 0.0_PREC, ser_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "SER", "DC", func, 0.0_PREC, ser_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "SER", "DT", func, 0.0_PREC, ser_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "THR", "DA", func, 0.0_PREC, thr_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "THR", "DG", func, 0.0_PREC, thr_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "THR", "DC", func, 0.0_PREC, thr_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "THR", "DT", func, 0.0_PREC, thr_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TRP", "DA", func, 0.0_PREC, trp_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TRP", "DG", func, 0.0_PREC, trp_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TRP", "DC", func, 0.0_PREC, trp_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TRP", "DT", func, 0.0_PREC, trp_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TYR", "DA", func, 0.0_PREC, tyr_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TYR", "DG", func, 0.0_PREC, tyr_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TYR", "DC", func, 0.0_PREC, tyr_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "TYR", "DT", func, 0.0_PREC, tyr_dt_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "VAL", "DA", func, 0.0_PREC, val_da_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "VAL", "DG", func, 0.0_PREC, val_dg_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "VAL", "DC", func, 0.0_PREC, val_dc_c12
        write (itp_out, "(A5,A5,I5,ES15.4,ES15.4)") "VAL", "DT", func, 0.0_PREC, val_dt_c12
        write (itp_out, "(A)") " "
    end if

    ! Close itp file
    close (itp_out)
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    !///////////////////////////////////////////////////////////////////////////////////////////////
    !-----------------------------------------------------------------------------------------------
    ! File: flexible_local_angle.itp
    !-----------------------------------------------------------------------------------------------
    ! Get itp file name
    if (path /= "") then
        filename = trim(adjustl(path))//"/flexible_local_angle.itp"
    else
        filename = "./flexible_local_angle.itp"
    end if

    ! Open itp file
    open (new_unit(itp_out), file=trim(adjustl(filename)), &
          status="replace", action="write")

    write (itp_out, "(A)") "[ flexible_local_angle ]"
    do iaa = 1, 20
        write (itp_out, "(A3,I2,I3)") aa_name(iaa), 1, 10
        write (itp_out, "(A)") ";             x              y             y2"
        do ipara = 1, 10
            write (itp_out, "(3ES15.4)") inflp%ang_para_x(ipara), &
                inflp%ang_para_y(iaa, ipara)*4.184_PREC, &
                inflp%ang_para_y2(iaa, ipara)*4.184_PREC
        end do
        write (itp_out, "(A)") ""
    end do

    ! Close itp file
    close (itp_out)
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    !///////////////////////////////////////////////////////////////////////////////////////////////
    !-----------------------------------------------------------------------------------------------
    ! File: unit_name.top
    !-----------------------------------------------------------------------------------------------
    ! Get output file unit number
    ofile = outfile%top

    ! Write include file
    write (ofile, "(A)") "; Atom types for coarse-grained models"
    write (ofile, "(A)") "#include atom_types.itp"
    if (inmisc%class_flag(CLASS%PRO)) then
        write (ofile, "(A)") "; Flexible local angle parameters"
        write (ofile, "(A)") "#include flexible_local_angle.itp"
    end if

    ! Init offset
    mp_offset = 0
    res_offset = 0

    do iunit = 1, nunit_real
        ! Get itp file name
        select case (iclass_unit(iunit))
        case (CLASS%DNA2)
            write (mol_name(iunit), "(A9,I0)") "dna_chain", iunit
        case (CLASS%PRO)
            write (mol_name(iunit), "(A7,I0)") "protein", iunit
        case default
            write (mol_name(iunit), "(A7,I0)") "unknown", iunit
        end select

        if (path /= "") then
            filename = trim(adjustl(path))//"/"// &
                       trim(adjustl(mol_name(iunit)))//".itp"
        else
            filename = "./"//trim(adjustl(mol_name(iunit)))//".itp"
        end if

        ! Open itp file
        open (new_unit(itp_out), file=trim(adjustl(filename)), &
              status="replace", action="write")

        ! Write molecule type section
        write (itp_out, "(A)") "[ moleculetype ]"
        write (itp_out, "(A)") ";name   nrexcl"
        select case (iclass_unit(iunit))
        case (CLASS%DNA2)
            write (itp_out, "(A,I0,I8)") "dna_chain", iunit, -1
        case (CLASS%PRO)
            write (itp_out, "(A,I0,I8)") "protein", iunit, -1
        case default
            write (itp_out, "(A,I0,I8)") "unknown", iunit, -1
        end select
        write (itp_out, "(A)") ""

        ! Write atom type section
        write (itp_out, "(A)") "[ atoms ]"
        write (itp_out, "(A)") &
            ";   nr    type   resnr  residu    atom    cgnr  charge    mass"
        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
            ires = ires_mp(imp)
            if (iclass_unit(iunit) == CLASS%DNA2 .and. &
                trim(adjustl(cmp2atom(imp))) == "DB") then
                ptype = cmp2seq(imp)
            else if (iclass_unit(iunit) == CLASS%PRO) then
                ptype = cmp2seq(imp)
            else
                ptype = cmp2atom(imp)
            end if
            cgnr = 1
            write (itp_out, "(I6,A8,I8,A8,A8,I8,F8.3,F8.3)") &
                imp - mp_offset, &
                trim(adjustl(ptype)), &
                ires - res_offset, &
                trim(adjustl(cmp2seq(imp))), &
                trim(adjustl(cmp2atom(imp))), &
                cgnr, &
                imp2charge(imp), &
                cmass_mp(imp)
        end do
        write (itp_out, "(A)") ""

        ! Write harmonic bonds
        write (itp_out, "(A)") "[ bonds ] ; Bond - Quadratic Harmonic"
        write (itp_out, "(A)") &
            ";   i    j    f             eq             k2"
        do ibd = 1, nbd
            imp1 = ibd2mp(1, ibd)
            imp2 = ibd2mp(2, ibd)
            func = 1
            if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit) then
                write (itp_out, "(I5,I5,I5,ES15.4,ES15.4)") &
                    imp1 - mp_offset, &
                    imp2 - mp_offset, &
                    func, &
                    bd_nat(ibd)/10.0_PREC, &
                    coef_bd(1, ibd)*418.4_PREC
            end if
        end do
        if (iclass_unit(iunit) == CLASS%DNA2) then
            write (itp_out, "(A)") "; DNA - Bond - Quartic Harmonic"
            write (itp_out, "(A)") &
                ";   i    j    f             eq             k4"
            do ibd = 1, nbd
                imp1 = ibd2mp(1, ibd)
                imp2 = ibd2mp(2, ibd)
                func = 21
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit) then
                    write (itp_out, "(I5,I5,I5,ES15.4,ES15.4)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        func, &
                        bd_nat(ibd)/10.0_PREC, &
                        coef_bd(2, ibd)*41840.0_PREC
                end if
            end do
        end if

        write (itp_out, "(A)") " "

        ! Write quadratic harmonic angles
        if (iclass_unit(iunit) == CLASS%DNA .or. &
            iclass_unit(iunit) == CLASS%DNA2) then
            write (itp_out, "(A)") "[ angles ] ; Angle - Quadratic Harmonic"
            write (itp_out, "(A)") &
                ";   i    j    k    f             eq              k"
            do iba = 1, nba
                imp1 = iba2mp(1, iba)
                imp2 = iba2mp(2, iba)
                imp3 = iba2mp(3, iba)
                func = 1
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit .and. &
                    imp2unit(imp3) == iunit) then
                    write (itp_out, "(I5,I5,I5,I5,ES15.4,ES15.4)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        imp3 - mp_offset, &
                        func, &
                        ba_nat(iba), &
                        coef_ba(1, iba)*4.184_PREC
                end if
            end do
            write (itp_out, "(A)") " "
        end if

        ! Write gaussian angles
        if (iclass_unit(iunit) == CLASS%PRO) then
            write (itp_out, "(A)") "[ angles ] ; Protein - Angle - Gaussian"
            write (itp_out, "(A)") &
                ";   i    j    k    f             eq              k              w"
            do iba = 1, nba
                imp1 = iba2mp(1, iba)
                imp2 = iba2mp(2, iba)
                imp3 = iba2mp(3, iba)
                func = 21
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit .and. &
                    imp2unit(imp3) == iunit) then
                    write (itp_out, "(I5,I5,I5,I5,ES15.4,ES15.4,ES15.4)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        imp3 - mp_offset, &
                        func, &
                        aicg13_nat(iba)/10.0_PREC, &
                        coef_aicg13_gauss(iba)*4.184_PREC, &
                        wid_aicg13_gauss(iba)/10.0_PREC
                end if
            end do
            write (itp_out, "(A)") " "
        end if

        ! Write flexible angles
        if (iclass_unit(iunit) == CLASS%PRO) then
            write (itp_out, "(A)") "[ angles ] ; Protein - Angle - Flexible"
            write (itp_out, "(A)") &
                ";   i    j    k    f"
            do iba = 1, nba
                imp1 = iba2mp(1, iba)
                imp2 = iba2mp(2, iba)
                imp3 = iba2mp(3, iba)
                func = 22
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit .and. &
                    imp2unit(imp3) == iunit) then
                    write (itp_out, "(I5,I5,I5,I5)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        imp3 - mp_offset, &
                        func
                end if
            end do
            write (itp_out, "(A)") " "
        end if

        ! Write gaussian dihedrals
        if (iclass_unit(iunit) == CLASS%DNA .or. &
            iclass_unit(iunit) == CLASS%DNA2) then
            write (itp_out, "(A)") "[ dihedrals ] ; DNA - Dihedral - Gaussian"
        else if (iclass_unit(iunit) == CLASS%PRO) then
            write (itp_out, "(A)") "[ dihedrals ] ; Protein - Dihedral - Gaussian"
        end if
        write (itp_out, "(A)") &
            ";   i    j    k    l    f             eq              k              w"
        do idih = 1, ndih
            imp1 = idih2mp(1, idih)
            imp2 = idih2mp(2, idih)
            imp3 = idih2mp(3, idih)
            imp4 = idih2mp(4, idih)
            func = 22
            if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit .and. &
                imp2unit(imp3) == iunit .and. imp2unit(imp4) == iunit) then
                write (itp_out, "(I5,I5,I5,I5,I5,ES15.4,ES15.4,ES15.4)") &
                    imp1 - mp_offset, &
                    imp2 - mp_offset, &
                    imp3 - mp_offset, &
                    imp4 - mp_offset, &
                    func, &
                    dih_nat(idih), &
                    coef_dih_gauss(idih)*4.184_PREC, &
                    wid_dih_gauss(idih)
            end if
        end do
        write (itp_out, "(A)") " "

        ! Write cosine dihedrals
        if (iclass_unit(iunit) == CLASS%DNA .or. &
            iclass_unit(iunit) == CLASS%DNA2) then
            write (itp_out, "(A)") "[ dihedrals ] ; DNA - Dihedral - Periodic"
            write (itp_out, "(A)") &
                ";   i    j    k    l    f             eq              k    n"
            do idih = 1, ndih
                imp1 = idih2mp(1, idih)
                imp2 = idih2mp(2, idih)
                imp3 = idih2mp(3, idih)
                imp4 = idih2mp(4, idih)
                func = 1
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit .and. &
                    imp2unit(imp3) == iunit .and. imp2unit(imp4) == iunit) then
                    write (itp_out, "(I5,I5,I5,I5,I5,ES15.4,ES15.4,I5)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        imp3 - mp_offset, &
                        imp4 - mp_offset, &
                        func, &
                        dih_nat(idih), &
                        2.0_PREC, &
                        1
                end if
            end do
            write (itp_out, "(A)") " "
        end if

        ! Write flexible dihedrals
        if (iclass_unit(iunit) == CLASS%PRO) then
            write (itp_out, "(A)") "[ dihedrals ] ; Protein - Flexible - Gaussian"
            write (itp_out, "(A)") &
                ";   i    j    k    l    f"// &
                "             c0             c1             c2"// &
                "             c3             c4             c5             c6"
            do idih = 1, nfdih
                imp1 = ifdih2mp(1, idih)
                imp2 = ifdih2mp(2, idih)
                imp3 = ifdih2mp(3, idih)
                imp4 = ifdih2mp(4, idih)
                func = 22
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit .and. &
                    imp2unit(imp3) == iunit .and. imp2unit(imp4) == iunit) then
                    write (itp_out, "(I5,I5,I5,I5,I5,7ES15.4)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        imp3 - mp_offset, &
                        imp4 - mp_offset, &
                        func, &
                        fdih_para(1, idih)*4.184_PREC, &
                        fdih_para(2, idih)*4.184_PREC, &
                        fdih_para(3, idih)*4.184_PREC, &
                        fdih_para(4, idih)*4.184_PREC, &
                        fdih_para(5, idih)*4.184_PREC, &
                        fdih_para(6, idih)*4.184_PREC, &
                        fdih_para(7, idih)*4.184_PREC
                end if
            end do
            write (itp_out, "(A)") " "
        end if

        ! Write non local go contacts
        if (iclass_unit(iunit) == CLASS%PRO) then
            write (itp_out, "(A)") "[ pairs ] ; Protein - Non-Local - Go contacts"
            write (itp_out, "(A)") &
                ";   i    j    f             eq              k"
            do icon = 1, ncon
                imp1 = icon2mp(1, icon)
                imp2 = icon2mp(2, icon)
                func = 1
                if (imp2unit(imp1) == iunit .and. imp2unit(imp2) == iunit) then
                    write (itp_out, "(I5,I5,I5,ES15.4,ES15.4)") &
                        imp1 - mp_offset, &
                        imp2 - mp_offset, &
                        func, &
                        go_nat(icon)/10.0_PREC, &
                        coef_go(icon)*4.184_PREC
                end if
            end do
            write (itp_out, "(A)") " "
        end if

        ! Close itp file
        close (itp_out)

        ! Update offset
        mp_offset = lunit2mp(2, iunit)
        res_offset = ires
    end do

    ! Write include file
    write (ofile, "(A)") "; Molecule topologies"
    do iunit = 1, nunit_real
        write (ofile, "(A)") "#include "//trim(adjustl(mol_name(iunit)))//".itp"
    end do
    write (ofile, "(A)") " "

    ! Write system name
    write (ofile, "(A)") "[ system ]"
    write (ofile, "(A)") "CafeMol course-grained system"
    write (ofile, "(A)") " "

    ! Write molecules
    write (ofile, "(A)") "[ molecules ]"
    do iunit = 1, nunit_real
        write (ofile, "(A16,I1)") mol_name(iunit), 1
    end do
    write (ofile, "(A)") " "

    ! Close top file
    close (ofile)
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

end subroutine write_top

function imp2charge(imp) result(charge)

    use var_setp, only: inele
    use var_struct, only: imp2type
    use const_index, only: MPTYPE, CHARGETYPE
    use const_maxsize, only: PREC

    ! Function arguments
    integer, intent(in) :: imp
    real(PREC) :: charge

    select case (imp2type(imp))
    case (MPTYPE%DNA2_PHOS)
        charge = inele%coef_charge_type(CHARGETYPE%P2)
    case default
        charge = 0.0
    end select

end function imp2charge
