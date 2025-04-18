! setp_md_info
!> @brief Reading parameters of "md_information" field in input file

! ******************************************************************
subroutine setp_md_info()

    !use, intrinsic :: iso_fortran_env
    use const_maxsize
    use const_index
    use const_physical
    use var_inp, only: infile, outfile, flg_rst
    use var_setp, only: insimu, inmisc, irand, mts, &
        inmmc   !mcanonical
    use var_replica, only: exchange_step, flg_rep, &
        n_replica_all, n_replica_mpi, irep2grep
    use mt_stream
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    !------------------------------------------------------------------
    ! intent(out) :: n_step_sim, n_step_neighbor, n_step_save, &
    !     n_tstep, tstep_size, tempk

    !------------------------------------------------------------------
    ! local variables

    integer       :: i
    integer       :: icol, isim
    integer       :: luninp, lunout
    integer       :: iline, nlines, iequa, nequat
    integer       :: irep, grep, istream
    integer       :: clock

    character(4)  :: kfind
    character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM)  :: cvalue
    character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MXCOLM) :: ctmp02
    character(CARRAY_MSG_ERROR) :: error_message

    integer(INT32) :: idstream
#ifdef MPI_PAR
    integer :: tn
#endif

    ! ---------------------------------------------------------------------
    luninp = infile%inp
    lunout = outfile%data

    ! ---------------------------------------------------------------------
    ! default setting
    insimu%n_step_sim = -1
    insimu%i_step_sim_init = 1  ! default
    insimu%n_tstep(1:MXSIM) = -1
    insimu%i_tstep_init = 1  ! default
    insimu%tstep_size = -1.0
    insimu%n_step_save = -1
    insimu%n_step_rst = -1
    insimu%n_step_neighbor = -1
    insimu%i_com_zeroing_ini = 0
    insimu%i_com_zeroing = -1
    insimu%i_no_trans_rot = -1
    insimu%tempk = -1.0
    insimu%tempk_ref = 0.0
    insimu%i_rand_type = 0
    insimu%n_seed = 0

    inmisc%i_redef_para = 0
    inmisc%i_in_box = 0
    inmisc%i_in_cap = 0
    inmisc%i_cylinder = 0
    inmisc%i_del_int = 0
    inmisc%i_energy_para = 0
    inmisc%i_neigh_dist = 0
    inmisc%i_mass_unit_expert = 0
    inmisc%i_mass_unit = 1
    inmisc%i_mass = 1
    inmisc%i_fric = 1
    inmisc%i_redef_mass_fric = 0
    inmisc%i_bridge = 0
    if (.not. flg_rep(REPTYPE%PULL)) then
        inmisc%i_pulling = 0
    endif
    inmisc%i_anchor = 0
    inmisc%i_rest1d = 0
    inmisc%i_fix = 0
    inmisc%i_implig = 0
    inmisc%i_afm_fitting = 0
    inmisc%i_stage = 0

    inmmc%i_modified_muca = 0
    inmisc%i_reset_struct = 0

    ! ---------------------------------------------------------------------

    if (flg_rst) then
        call read_rst(RSTBLK%STEP)
        exchange_step(:) = insimu%i_tstep_init + exchange_step(:)
        insimu%i_tstep_init = insimu%i_tstep_init + 1
    endif

    ! ---------------------------------------------------------------------

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        rewind (luninp)
        call ukoto_uiread2(luninp, lunout, 'md_information  ', kfind, &
                           CARRAY_MXLINE, nlines, cwkinp)
        if (kfind /= 'FIND') then
            error_message = 'Error: cannot find "md_information" field in the input file'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
            call ukoto_uiequa2(cwkinp(iline), nequat, csides)

            do iequa = 1, nequat
                cvalue = 'n_step_sim'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%n_step_sim, cvalue)
            end do
        end do

        ! -----------------------------------------------------------------
        ! checking n_step_sim

        if (insimu%n_step_sim <= 0) then
            error_message = 'Error: invalid value for n_step_sim'
            call util_error(ERROR%STOP_ALL, error_message)

        else if (insimu%n_step_sim > MXSIM) then
            error_message = 'Error: should increase MXSIM'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        do iline = 1, nlines
            call ukoto_uiequa2(cwkinp(iline), nequat, csides)

            do iequa = 1, nequat
                ctmp02 = csides(1, iequa)

                if (ctmp02(1:7) == 'n_tstep') then
                    do icol = 9, CARRAY_MXCOLM
                        if (ctmp02(icol:icol) == ')') exit
                    end do
                    read (ctmp02(9:icol - 1), *) isim
                    cvalue = ctmp02(1:icol)
                    call ukoto_lvalue2(lunout, csides(1, iequa), &
                                       insimu%n_tstep(isim), cvalue)
                end if

                cvalue = 'tstep_size'
                call ukoto_rvalue2(lunout, csides(1, iequa), &
                                   insimu%tstep_size, cvalue)

                cvalue = 'n_step_save'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%n_step_save, cvalue)

                cvalue = 'n_step_rst'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%n_step_rst, cvalue)

                cvalue = 'n_step_neighbor'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%n_step_neighbor, cvalue)

                cvalue = 'i_com_zeroing_ini'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%i_com_zeroing_ini, cvalue)

                cvalue = 'i_com_zeroing'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%i_com_zeroing, cvalue)

                cvalue = 'i_no_trans_rot'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%i_no_trans_rot, cvalue)

                cvalue = 'i_in_box'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_in_box, cvalue)

                cvalue = 'i_in_cap'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_in_cap, cvalue)

                cvalue = 'i_cylinder'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_cylinder, cvalue)

                cvalue = 'i_afm_fitting'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_afm_fitting, cvalue)

                cvalue = 'i_stage'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_stage, cvalue)

                cvalue = 'i_del_int'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_del_int, cvalue)

                cvalue = 'i_redef_para'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_redef_para, cvalue)

                cvalue = 'i_energy_para'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_energy_para, cvalue)

                cvalue = 'i_neigh_dist'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_neigh_dist, cvalue)

                cvalue = 'i_mass_unit_expert'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_mass_unit_expert, cvalue)

                cvalue = 'i_mass_unit'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_mass_unit, cvalue)

                cvalue = 'i_mass'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_mass, cvalue)

                cvalue = 'i_fric'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_fric, cvalue)

                cvalue = 'i_mass_fric'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_redef_mass_fric, cvalue)

                cvalue = 'i_bridge'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_bridge, cvalue)

                if (.not. flg_rep(REPTYPE%PULL)) then
                    cvalue = 'i_pulling'
                    call ukoto_ivalue2(lunout, csides(1, iequa), &
                                       inmisc%i_pulling, cvalue)
                endif

                cvalue = 'i_anchor'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_anchor, cvalue)

                cvalue = 'i_rest1d'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_rest1d, cvalue)

                cvalue = 'i_fix'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_fix, cvalue)

                cvalue = 'i_implig'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_implig, cvalue)

                cvalue = 'i_modified_muca'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmmc%i_modified_muca, cvalue)

                cvalue = 'i_reset_struct'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   inmisc%i_reset_struct, cvalue)

                cvalue = 'tempk'
                call ukoto_rvalue2(lunout, csides(1, iequa), &
                                   insimu%tempk, cvalue)

                cvalue = 'tempk_ref'
                call ukoto_rvalue2(lunout, csides(1, iequa), &
                                   insimu%tempk_ref, cvalue)

                cvalue = 'i_rand_type'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%i_rand_type, cvalue)

                cvalue = 'n_seed'
                call ukoto_ivalue2(lunout, csides(1, iequa), &
                                   insimu%n_seed, cvalue)
            end do

        end do

        if (insimu%n_seed <= -2 .or. insimu%n_seed == 0) then
            error_message = 'Error: invalid n_seed number (should be n_seed >= 1 or n_seed =-1(system time))'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (insimu%n_seed == -1) then
            call system_clock(count=clock)
            irand = clock
            write (*, *) 'set random seed from system clocktime: using seed = ', irand
            write (lunout, *) 'set random seed from system clocktime: using seed =', irand

        else if (flg_rst) then
            ! The random seed is changed when restarted.
            i = insimu%i_tstep_init
            do while (i > 1000000000)
                i = i/1000000000
            enddo
            irand = insimu%n_seed + int(i)
            write (*, *) 'reset random seed from n_seed+i_tstep_init: using seed = ', irand
            write (lunout, *) 'reset random seed from n_seed+i_tstep_init: using seed = ', irand

        else
            irand = insimu%n_seed
            write (*, *) 'set random seed form n_seed: using seed = ', irand
            write (lunout, *) 'set random seed from n_seed: using seed = ', irand
        endif

        ! -----------------------------------------------------------------
        ! set ntstep_all
        insimu%n_tstep_all = 0
        do i = 1, insimu%n_step_sim
            insimu%n_tstep_all = insimu%n_tstep_all + insimu%n_tstep(i)
        enddo

        ! -----------------------------------------------------------------
        ! checking input variables

        do i = 1, insimu%n_step_sim
            if (insimu%n_tstep(i) <= 0) then
                write (error_message, *) 'Error: invalid value for n_tstep(', i, ') = ', insimu%n_tstep(i)
                call util_error(ERROR%STOP_ALL, error_message)
            else if (insimu%n_tstep(i) > MX_NTSTEP) then
                write (error_message, *) 'Error: n_tstep(', i, ') > MX_NTSTEP defined in const_maxsize.F90'
                call util_error(ERROR%STOP_ALL, error_message)
            end if
        end do

        if (insimu%tstep_size <= 0) then
            error_message = 'Error: invalid value for tstep_size'
            call util_error(ERROR%STOP_ALL, error_message)
        else if (insimu%n_step_save <= 0) then
            error_message = 'Error: invalid value for n_step_save'
            call util_error(ERROR%STOP_ALL, error_message)
        else if (insimu%n_step_rst <= 0) then
            insimu%n_step_rst = insimu%n_step_save
        else if (insimu%n_step_neighbor <= 0) then
            error_message = 'Error: invalid value for n_step_neighbor'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (insimu%i_com_zeroing_ini == 0) then
            write (lunout, *) 'Do not move center of mass to ZERO of initial structure: i_com_zeroing_ini = 0'
        else if (insimu%i_com_zeroing_ini == 1) then
            write (lunout, *) 'Move center of mass to ZERO of initial structure: i_com_zeroing_ini = 1'
        else if (insimu%i_com_zeroing_ini == 2) then
            write (lunout, *) 'Move center of mass to (4500, 4500, 4500) if initial structure: i_com_zeroing_ini = 2'
        else
            error_message = 'Error: invalid value for i_com_zeroing_ini'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (insimu%i_com_zeroing == 0) then
            write (lunout, *) 'Do not move center of mass to ZERO: i_com_zeroing = 0 (only for output)'
        else if (insimu%i_com_zeroing == 1) then
            write (lunout, *) 'Move center of mass to ZERO: i_com_zeroing = 1 (only for output)'
        else if (insimu%i_com_zeroing == 2) then
            write (lunout, *) 'Move center of mass to (4500, 4500, 4500): i_com_zeroing = 2 (only for output)'
        else
            error_message = 'Error: invalid value for i_com_zeroing'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (insimu%i_no_trans_rot == 0) then
            write (lunout, *) 'Unprohibit translation and rotation: i_no_trans_rot = 0'
        else if (insimu%i_no_trans_rot == 1) then
            write (lunout, *) 'Prohibit translation and rotation: i_no_trans_rot = 1'
        else
            error_message = 'Error: invalid value for i_no_trans_rot'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! using box, inmisc%i_in_box
        if (inmisc%i_in_box == 0) then
        else if (inmisc%i_in_box == 1) then
            write (lunout, *) 'using box: i_in_box = 1'
        else
            error_message = 'Error: invalid value for i_in_box'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! using cap, inmisc%i_in_cap
        if (inmisc%i_in_cap == 0) then
        else if (inmisc%i_in_cap == 1) then
            write (lunout, *) 'using cap: i_in_cap = 1'
        else
            error_message = 'Error: invalid value for i_in_cap'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! delete interaction, inmisc%i_del_int
        if (inmisc%i_del_int == 0) then
        else if (inmisc%i_del_int == 1) then
            write (lunout, *) 'delete interaction: i_del_int = 1'
        else
            error_message = 'Error: invalid value for i_del_int'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! redefine parameters, inmisc%i_redef_para
        if (inmisc%i_redef_para == 0) then
        else if (inmisc%i_redef_para == 1) then
            write (lunout, *) 'redefine parameters: i_redef_para = 1'
        else
            error_message = 'Error: invalid value for i_redef_para'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! changing energy coefficient, inmisc%i_energy_para
        if (inmisc%i_energy_para == 0) then
        else if (inmisc%i_energy_para == 1) then
            write (lunout, *) 'changing energy coefficient: i_energy_para = 1'
        else
            error_message = 'Error: invalid value for i_energy_para'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! changing rneighbor_dist, inmisc%i_neigh_dist
        if (inmisc%i_neigh_dist == 0) then
        else if (inmisc%i_neigh_dist == 1) then
            write (lunout, *) 'changing rneihbor_dist: i_neigh_dist = 1'
        else
            error_message = 'Error: invalid value for i_neigh_dist'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! Define cafemol mass unit expert
        if (inmisc%i_mass_unit_expert == 0) then
            ! Default
            continue
        else if (inmisc%i_mass_unit_expert == 1) then
            ! Expert mode about mass and friction
            write (lunout, *) "Expert mode of mass and friction"
        else
            error_message = 'Error: invalid value for i_mass_unit_expert'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! Define cafemol mass unit
        if (inmisc%i_mass_unit == 0) then
            ! Old style: 1 cafe-mu = 13.7 amu
            write (lunout, *) "Old style mass unit style: 1 cafe-mu = 13.7 amu"
        else if (inmisc%i_mass_unit == 1) then
            ! Default: 1 cafe-mu = 1 amu
            write (lunout, *) "Mass unit: 1 cafe-mu = 1 amu"
        else
            error_message = 'Error: invalid value for i_mass_unit'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! Mass
        if (inmisc%i_mass == 0) then
            ! cmass is used for mass of all particles
            write (lunout, *) "Each paritcle has constant mass defined by rmass in general.para."
        else if (inmisc%i_mass == 1) then
            ! Default: particle-type dependent mass written in para/general.para are used.
       write (lunout, *) "Each particle has different mass depending on the particle type define by MASS variables in general.para."
        else
            error_message = 'Error: invalid value for i_mass'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! Friction
        if (inmisc%i_fric == 0) then
            ! fric_const is used for friction of all particles
            write (lunout, *) "All particles are subject to a friction of fric_const in para/gerneral.para"
        else if (inmisc%i_fric == 1) then
            ! Default: Friction is depedenet of mass
            write (lunout, *) "Each particle is subject to a friction: fric_mp=fric_scale*fric_typical_coef/mass_mp."
        else if (inmisc%i_fric == 11) then
            ! Use Stokes' law to derive frictions for each particles
            write (lunout, *) "Friction coefficients are derived by Stokes' law based on particle's radius."
        else
            error_message = 'Error: invalid value for i_fric'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! Check combination of i_mass_unit, i_mass, and i_fric
        if (inmisc%i_mass_unit_expert == 0) then
            if (inmisc%i_mass_unit == 0 .and. inmisc%i_mass == 0 .and. inmisc%i_fric == 0) then
                continue
            else if (inmisc%i_mass_unit == 1 .and. inmisc%i_mass >= 1 .and. inmisc%i_fric >= 1) then
                continue
            else
                error_message = 'Error: invalid combination of for i_mass_unit, i_mass, and i_fric'
                call util_error(ERROR%STOP_ALL, error_message)
            end if
        end if

        ! -----------------------------------------------------------------
        ! changing mass and friction, inmisc%i_redef_mass_fric
        if (inmisc%i_redef_mass_fric == 0) then
        else if (inmisc%i_redef_mass_fric == 1) then
            write (lunout, *) 'changing mass and friction: i_mass_fric = 1'
        else
            error_message = 'Error: invalid value for i_mass_fric'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying bridge option, inmisc%i_bridge
        if (inmisc%i_bridge == 0) then
        else if (inmisc%i_bridge == 1) then
            write (lunout, *) 'applying bridge option: i_bridge = 1'
        else
            error_message = 'Error: invalid value for i_bridge'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying pulling option, inmisc%i_pulling
        if (inmisc%i_pulling == 0) then
        else if (inmisc%i_pulling == 1) then
            write (lunout, *) 'applying pulling option: i_pulling = 1'
        else
            error_message = 'Error: invalid value for i_pulling'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying anchor option, inmisc%i_anchor
        if (inmisc%i_anchor == 0) then
        else if (inmisc%i_anchor == 1) then
            write (lunout, *) 'applying anchor option: i_anchor = 1'
        else
            error_message = 'Error: invalid value for i_anchor'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying rest1d option, inmisc%i_rest1d
        if (inmisc%i_rest1d == 0) then
        else if (inmisc%i_rest1d == 1) then
            write (lunout, *) 'applying 1D-restraint option: i_rest1d = 1'
        else
            error_message = 'Error: invalid value for i_rest1d'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying fix option, inmisc%i_fix
        if (inmisc%i_fix == 0) then
        else if (inmisc%i_fix == 1) then
            write (lunout, *) 'applying fix option: i_fix = 1'
        else
            error_message = 'Error: invalid value for i_fix'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying implicit_ligand option, inmisc%i_implig
        if (inmisc%i_implig == 0) then
        else if (inmisc%i_implig == 1) then
            write (lunout, *) 'applying implicit_ligand (with MD-MC) option: i_implig = 1'
        else
            error_message = 'Error: invalid value for i_implig'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying modified-multicanonical sampling, inmmc%i_modified_muca
        if (inmmc%i_modified_muca == 0) then
            !write (lunout, *) 'not using modified multicanonical sampling'
        else if (inmmc%i_modified_muca == 1) then
            write (lunout, *) 'using modified multicanonical sampling'
        else if (inmmc%i_modified_muca == 2) then
            write (lunout, *) 'using Wang-Landau multicanonical sampling'
        else
            error_message = 'Error: invalid value for i_modified_muca'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! -----------------------------------------------------------------
        ! applying structure resetting option, inmisc%i_reset_struct
        if (inmisc%i_reset_struct == 0) then
        else if (inmisc%i_reset_struct == 1) then
            write (lunout, *) 'structure will be reset to the initial structure at each learning step'
        else
            error_message = 'Error: invalid value for i_reset_struct'
            call util_error(ERROR%STOP_ALL, error_message)
        endif

        if (insimu%tempk < 0) then
            error_message = 'Error: invalid value for tempk'
            call util_error(ERROR%STOP_ALL, error_message)
        endif
        if (insimu%tempk < ZERO_JUDGE) then
            error_message = 'Error: tempk is too small'
            call util_error(ERROR%STOP_ALL, error_message)
        endif

        if (insimu%i_rand_type == 0) then
            write (lunout, *) 'using serial mt_stream for random number (default)'
        else if (insimu%i_rand_type == 1) then
            write (lunout, *) 'using parallel mt_stream for random number'
        else if (insimu%i_rand_type == 2) then
            write (lunout, *) 'using parallel(only open MP) mt_stream for random number'
        else
            error_message = 'Error: invalid value for i_rand_type'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (insimu%n_seed == 0) then
            error_message = 'Error: invalid value for n_seed'
            call util_error(ERROR%STOP_ALL, error_message)
        end if

#ifdef MPI_PAR
    end if

    call MPI_Bcast(insimu, insimu%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(inmmc, inmmc%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(irand, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

!  call sgrnd(irand)

    if (insimu%i_rand_type == 0) then
        allocate (mts(0:n_replica_mpi, 0:0))
    else
#ifdef MPI_PAR
        allocate (mts(0:n_replica_mpi*npar_mpi, 0:nthreads - 1))
#else
        allocate (mts(0:n_replica_mpi, 0:0))
#endif
    end if

    call set_mt19937
    call new(mts(0, 0))
    call init(mts(0, 0), irand)  ! init by array

    do irep = 1, n_replica_mpi

        if (insimu%i_rand_type == 0) then
            istream = irep
            grep = irep2grep(irep)
            idstream = grep
            if (idstream /= 0) then
                call create_stream(mts(0, 0), mts(istream, 0), idstream)
            end if
        else
#ifdef MPI_PAR
!$omp parallel private(tn, grep, istream, idstream)
            tn = 0
!$          tn = omp_get_thread_num()
            grep = irep2grep(irep)
            istream = irep
            idstream = (grep - 1)*nthreads + tn + 1
            if (idstream /= 0) then
                call create_stream(mts(0, 0), mts(istream, tn), idstream)
            end if

            if (insimu%i_rand_type == 1) then
                istream = irep + local_rank_mpi*n_replica_mpi
                idstream = ((grep - 1) + local_rank_mpi*n_replica_all)*nthreads + tn + 1
                if (idstream /= 0) then
                    call create_stream(mts(0, 0), mts(istream, tn), idstream)
                end if
            end if

!$omp end parallel
#else
            istream = irep
            grep = irep2grep(irep)
            idstream = grep
            if (idstream /= 0) then
                call create_stream(mts(0, 0), mts(istream, 0), idstream)
            end if
#endif
        end if

    end do

end subroutine setp_md_info
