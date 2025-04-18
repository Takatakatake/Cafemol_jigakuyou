! setpara
!> @brief Set parameters

! **********************************************************************
! set parameters
subroutine setpara(xyz_mp_init)

    use const_maxsize
    use const_index
    use const_physical
    use var_inp, only: infile, outfile, ifile_pdb, num_file, &
        i_run_mode, i_seq_read_style, &
        i_simulate_type, inperi, &
        i_aicg, i_go_native_read_style    ! aicg
    use var_setp, only: inpara, inmisc, ifix_mp, inmmc, inflp
    use var_struct, only: xyz_ref_mp, iontype_mp
    use var_mgo, only: inmgo
    use var_implig, only: inimplig  ! for implicit ligand

#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    real(PREC), intent(out) :: xyz_mp_init(SPACE_DIM, MXMP)

    ! -----------------------------------------------------------------------
    ! local variables
    integer :: ipdb
    integer :: lunout, lunpdb
    integer :: npdb, input_status
    integer, allocatable :: iatomnum(:)
    real(PREC), allocatable :: xyz(:, :, :)
    character(4), allocatable :: cname_ha(:, :)  ! aicg
    character(1), allocatable :: dssp(:)  ! aicg
    character(72) :: char72
    character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
    write (*, *) '#### start setpara'
#endif

    ! -----------------------------------------------------------------------
    lunout = outfile%data
    npdb = num_file%pdb

    ! -----------------------------------------------------------------------
    ! periodic boundary
    inperi%n_mirror_index = 0
    if (inperi%i_periodic == 1) call setp_periodic()

    ! -----------------------------------------------------------------------
    ! write the title on the filename_out
#ifdef MPI_PAR
    if (myrank == 0) then
#endif
        do ipdb = 1, npdb
            lunpdb = ifile_pdb(1, ipdb)

            if (ifile_pdb(5, ipdb) /= 1) cycle

            write (lunout, '(72("*"))')
            do
                read (lunpdb, '(a72)', iostat=input_status) char72
                if (input_status < 0) then
                    exit
                else if (input_status > 0) then
                    error_message = 'Error: input error in setpara'
                    call util_error(ERROR%STOP_ALL, error_message)
                end if

                if (char72(1:6) == 'HEADER' .or. char72(1:6) == 'COMPND' .or. &
                    char72(1:5) == 'TITLE') then
                    write (lunout, '(3x, a66)') char72(7:72)
                end if
            end do
        end do
#ifdef MPI_PAR
    end if
#endif

    ! -----------------------------------------------------------------------
    ! reading parameters from parameter files
#ifdef _DEBUG
    write (*, *) 'reading parameters from parameter files'
#endif
    call setp_mapara(infile%para_gen, outfile%data)
    call setp_mapara_pro(infile%para_pro, outfile%data)
    if (inmisc%class_flag(CLASS%DNA)) then
        call setp_mapara_dna()
    endif
    if (inmisc%force_flag(INTERACT%DNA2C)) then
        call setp_mapara_dna2(infile%para_dna2c, outfile%data)
    else if (inmisc%class_flag(CLASS%DNA2)) then
        call setp_mapara_dna2(infile%para_dna2, outfile%data)
    endif
    if (inmisc%class_flag(CLASS%ION)) then
        call setp_mapara_ion()
    endif
    if (inmisc%class_flag(CLASS%RNA)) then
        call setp_mapara_rna(infile%para_rna, outfile%data)
    endif
    if (inmisc%class_flag(CLASS%LIG)) then
        call setp_mapara_ligand()
    endif

    ! -----------------------------------------------------------------------
    ! reading parameters from parameter files
    if (inmisc%force_flag(INTERACT%HP)) then
        call setp_mapara_hp()
    endif

    if (inmisc%force_flag(INTERACT%ELE)) then
        call setp_mapara_ele()
    endif

    ! sasa
    if (inmisc%force_flag(INTERACT%SASA)) then
        call setp_mapara_sasa()
    endif

    ! gb
    if (inmisc%force_flag(INTERACT%GB)) then
        call setp_mapara_ele()
    endif

    ! excluded volume
    if (inmisc%force_flag(INTERACT%EXV)) then
        call setp_mapara_exv()
    endif

    ! protein-DNA sequence-specific integrations
    if (inmisc%force_flag(INTERACT%PRO_DNA_PWM)) then
        call setp_mapara_pro_dna_pwm()
        call setp_pro_dna_pwm()
    endif
    ! protein-DNA sequence-non-specific integrations
    if (inmisc%force_flag(INTERACT%PRO_DNA_NONSPEC)) then
        call setp_mapara_pro_dna_nonspec()
        call setp_pro_dna_nonspec()
    endif
    ! -----------------------------------------------------------------------
    ! reading simu parameters from input file
#ifdef _DEBUG
    write (*, *) 'reading simu parameters from input file'
#endif
    call setp_md_info()

    ! annealing
    if (i_run_mode == RUN%SA) call setp_anneal_para()

    ! searching Tf
    if (i_run_mode == RUN%SEARCH_TF) call setp_searchingtf_para()

    ! re-define parameters
    if (inmisc%i_redef_para == 1) call setp_redef_para()

    ! -----------------------------------------------------------------------
    ! Next two substitutions must be after both setp_mapara() and setp_redef_para().
    ! -----------------------------------------------------------------------

    ! flexible local potential parameter
    !if (inmisc%i_add_int == 1) call setp_mapara_flp()
    if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
        call setp_mapara_flp()
    end if

    if (inmisc%class_flag(CLASS%ION)) call setp_para_ion()

    ! -----------------------------------------------------------------------
    ! reading sequence
    allocate (iatomnum(MXMP))
    allocate (xyz(SPACE_DIM, MXATOM_MP, MXMP))
    allocate (cname_ha(MXATOM_MP, MXMP))    ! aicg
    allocate (dssp(MXMP))    ! aicg

    iontype_mp(1:MXMP) = 0
    if (i_seq_read_style == SEQREAD%PDB) then
        ! read xyz coordinate from PDB file
        call read_xyz(xyz_mp_init, iatomnum, xyz, cname_ha)   ! aicg
        xyz_ref_mp(1:SPACE_DIM, 1:MXMP) = xyz_mp_init(1:SPACE_DIM, 1:MXMP)

    else if (i_seq_read_style == SEQREAD%INPUT_SEQ) then
        ! read sequence from sequence field in input file
        call read_seq()
        call setp_make_dna(xyz_mp_init)

    else if (i_seq_read_style == SEQREAD%INPUT_LIPID) then
        error_message = 'Error: This read style is not implemented yet.'
        call util_error(ERROR%STOP_ALL, error_message)

    else if (i_seq_read_style == SEQREAD%CG) then
        ! read xyz coordinate from PDB file with CafeMol-CG style
        call read_xyz_cg(xyz_mp_init)
        xyz_ref_mp(1:SPACE_DIM, 1:MXMP) = xyz_mp_init(1:SPACE_DIM, 1:MXMP)

    end if

    ! Set residue-type-specific excluded volume radii
    if (inmisc%force_flag(INTERACT%EXV)) then
        if (inmisc%i_residue_exv_radii == 1) then
            call setp_exv()
        end if
    end if

    ! set constant exv radii for residues
    if (inmisc%i_set_const_radii == 1) call setp_set_const_radii()

    ! ----------------------------------------------------------------------
    ! Re-def mass & friction
    !if(inmisc%i_redef_mass_fric == 1)
    call setp_mass_fric()

    ! ----------------------------------------------------------------------
    ! Setup groups
    call setp_group()
    ! This must be called after setp_redef_mass_fric() to calculate mass of group

    ! ----------------------------------------------------------------------
    ! box
    if (inmisc%i_in_box == 1) call setp_box()

    ! ----------------------------------------------------------------------
    ! cap
    if (inmisc%i_in_cap == 1) call setp_cap()

    ! delete some interactions
    if (inmisc%i_del_int == 1) call setp_del_int()

    ! modified multi-canonical
    if (inmmc%i_modified_muca == 1) call setp_modified_muca()
    ! Wang-Landau MuCa
    if (inmmc%i_modified_muca == 2) call setp_modified_muca()

    ! try to read <<<< hydrophobic
    if (inmisc%force_flag(INTERACT%HP)) then
        call setp_hydrophobic()
    endif

    if (inmisc%force_flag(INTERACT%ELE)) call setp_electrostatic()

    ! try to read <<<< gb
    if (inmisc%force_flag(INTERACT%GB)) then
        call setp_gb()
        call setp_electrostatic()
    endif

    ! ----------------------------------------------------------------------
    ! default setting
    inmisc%factor_local_unit(1:MXUNIT, 1:MXUNIT) = 1.0e0_PREC
    inmisc%factor_go_unit(1:MXUNIT, 1:MXUNIT) = 1.0e0_PREC
    if (inmisc%i_energy_para == 1) call setp_energy_para()

    ! -----------------------------------------------------------------------
    inmisc%rneighbordist2_unit(1:MXUNIT, 1:MXUNIT) = inpara%rneighbor_dist**2
    if (inmisc%i_neigh_dist == 1) call setp_neigh_dist()

    ! -----------------------------------------------------------------------
    call setp_energy_unit()

    ! -----------------------------------------------------------------------
    ! bridge
    if (inmisc%i_bridge == 1) call setp_bridge_para()

    ! pulling
    if (inmisc%i_pulling == 1) call setp_pulling_para()

    ! anchor
    if (inmisc%i_anchor == 1) call setp_anchor_para()

    ! rest1d
    if (inmisc%i_rest1d == 1) call setp_rest1d_para()

    ! fix
    ifix_mp(1:MXMP) = 0
    if (inmisc%i_fix == 1) call setp_fix_para()

    ! cylinder
    if (inmisc%i_cylinder == 1) call setp_cylinder_para()

    !-----------------------------------------------------------------------
    ! parameter setting for AFM fitting interaction
    if (inmisc%i_afm_fitting == 1) call setp_afm_fitting()
    if (inmisc%i_stage == 1) call setp_stage()

    !-----------------------------------------------------------------------
    ! parameter setting for elastic network model (enm)
    if (inmisc%force_flag(INTERACT%ENM)) call setp_para_enm()

    ! -----------------------------------------------------------------------
    ! parameter setting for multiple-Go model (mgo)
    if (inmgo%i_multi_mgo == 1) call setp_para_mgo()

    ! -----------------------------------------------------------------------
    ! parameter setting for mpc
    if (i_simulate_type == SIM%MPC) call setp_para_mpc()

    ! =======================================================================
    ! constructing potential
    ! memory allocation (mostly for var_struct module)
    call allocate_nativestruct()

    ! constructing native(reference) structures
    call setp_nativestruct(xyz_mp_init, iatomnum, xyz, cname_ha)

    ! Atomic interaction based coarse-grained model (AICG)
    if (i_go_native_read_style == NATIVEREAD%PDB) then
    if (inmisc%force_flag(INTERACT%AICG1) .OR. inmisc%force_flag_local(LINTERACT%L_AICG1)) then
        if (i_aicg == AICG%AUTO) then
            call setp_mapara_aicg()
            call read_dssp(dssp)
            call setp_aicg(iatomnum, xyz, cname_ha, dssp)
        else if (i_aicg == AICG%USER) then
            call read_paraaicg()
        end if
    end if
    if (inmisc%force_flag(INTERACT%AICG2) .OR. &
        inmisc%force_flag(INTERACT%AICG2_PLUS) .OR. &
        inmisc%force_flag_local(LINTERACT%L_AICG2) .OR. &
        inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
        if (i_aicg == AICG%AUTO) then
            call setp_mapara_aicg2()
            call setp_aicg2(iatomnum, xyz, cname_ha)
        else if (i_aicg == AICG%USER) then
            call read_paraaicg2()
        end if
    end if
    end if
    ! =======================================================================

    !------------------------------------------------------------------------
    ! implicit ligand
    if (inmisc%i_implig == 1) then
        call setp_para_implig() !! read parameter from input_file (<<<< implicit_ligand )
        !! parameter setting for implicit ligand model
        !! set initial state (BOUND, UN_BOUND) for implicit ligand
        if (inimplig%iexe_implig == 1) then
            call setp_bindsite_implig() !! read binding site information from input_file
            !! (<<<< binding_site) [IMPLIGSITE]
            call setp_con_implig(xyz_mp_init, iatomnum, xyz) ! fix the ligand mediated contact pair
        endif
    endif
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! flexible local potential
    !if(inmisc%i_del_int == 1 .and. inmisc%i_add_int == 1) then
    !   call setp_flexible_local()
    !endif
    !-----------------------------------------------------------------------

    deallocate (iatomnum)
    deallocate (xyz)
    deallocate (dssp)  ! aicg
    deallocate (cname_ha)  ! aicg

    ! fluctuation matching
    if (i_run_mode == RUN%FMAT) call setp_fmat_para()

    ! Energy minimization
    if (i_run_mode == RUN%EMIN) call setp_minimize_para()

#ifdef _DEBUG
    write (*, *) '#### end setpara'
#endif

end subroutine setpara
