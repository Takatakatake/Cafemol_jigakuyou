! simu_force
!> @brief Subroutine for force calculation

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! ************************************************************************
! subroutine for the force
! ************************************************************************
subroutine simu_force(force_mp, &  ! [ o]
                      irep &  ! [i ]
                      )

    use const_maxsize
    use const_physical
    use const_index
    use var_setp, only: inmisc, inele, inflp
    use var_struct, only: nunit_all, nmp_all
    use var_mgo, only: inmgo
    use time
    use mpiconst

    implicit none

    ! ---------------------------------------------------------------------
    real(PREC), intent(out) :: force_mp(SPACE_DIM, nmp_all)
    integer, intent(in)  :: irep

    ! ---------------------------------------------------------------------
    ! local variables
    integer    :: isys, istat
    real(PREC) :: ene_unit(nunit_all, nunit_all)
    real(PREC) :: force_mp_l(SPACE_DIM, nmp_all, 0:nthreads - 1)
    real(PREC) :: ene_unit_l(nunit_all, nunit_all, 0:nthreads - 1)
    integer :: tn, n
    real(PREC), allocatable :: force_mp_mgo(:, :, :, :, :)

    if (inmgo%i_multi_mgo >= 1) then
        allocate (force_mp_mgo(SPACE_DIM, nmp_all, &
                               inmgo%nstate_max_mgo, inmgo%nsystem_mgo, 0:nthreads - 1))
    else
        allocate (force_mp_mgo(SPACE_DIM, nmp_all, 1, 1, 0:nthreads - 1))
    endif

#ifdef _DEBUG
    write (6, *) 'simu_force: START'
#endif

#ifdef MPI_PAR
    call mpi_barrier(mpi_comm_local, ierr)
#endif

    TIME_S(tm_force_local)

#ifdef MPI_DEBUG
    print *, "nmp_all      = ", nmp_all
#endif

!$omp parallel private(tn)
    tn = 0
!$  tn = omp_get_thread_num()

    force_mp_l(1:SPACE_DIM, 1:nmp_all, tn) = 0.0_PREC

    if (inmgo%i_multi_mgo >= 1) then
        ene_unit_l(1:nunit_all, 1:nunit_all, tn) = 0.0_PREC

        do isys = 1, inmgo%nsystem_mgo
            do istat = 1, inmgo%nstate_mgo(isys)
                force_mp_mgo(1:SPACE_DIM, 1:nmp_all, istat, isys, tn) = 0.0_PREC
            end do
        end do
    end if

    call simu_force_bond(irep, force_mp_l(1, 1, tn), &
                         force_mp_mgo(1, 1, 1, 1, tn), &
                         ene_unit_l(1, 1, tn))

    call simu_force_bangle(irep, force_mp_l(1, 1, tn), &
                           force_mp_mgo(1, 1, 1, 1, tn), &
                           ene_unit_l(1, 1, tn))
    call simu_force_dih(irep, force_mp_l(1, 1, tn), &
                        force_mp_mgo(1, 1, 1, 1, tn), &
                        ene_unit_l(1, 1, tn))

    if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
        inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
        call simu_force_aicg13_gauss(irep, force_mp_l(1, 1, tn), &
                                     force_mp_mgo(1, 1, 1, 1, tn), &
                                     ene_unit_l(1, 1, tn))
    end if

    if (inmisc%force_flag_local(LINTERACT%L_AICG2)) then
        call simu_force_aicg14_gauss(irep, force_mp_l(1, 1, tn), &
                                     force_mp_mgo(1, 1, 1, 1, tn), &
                                     ene_unit_l(1, 1, tn))

    else if (inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS) .or. &
             inmisc%force_flag_local(LINTERACT%L_DNA2)) then

        call simu_force_dih_gauss(irep, force_mp_l(1, 1, tn), &
                                  force_mp_mgo(1, 1, 1, 1, tn), &
                                  ene_unit_l(1, 1, tn))
    end if
    if (inmisc%force_flag_local(LINTERACT%L_DNA2C)) then
        call simu_force_dih_3spn2c(irep, force_mp_l(1, 1, tn), &
                                   force_mp_mgo(1, 1, 1, 1, tn), &
                                   ene_unit_l(1, 1, tn))
    end if

    call simu_force_dih_harmonic(irep, force_mp_l(1, 1, tn), &
                                 force_mp_mgo(1, 1, 1, 1, tn), &
                                 ene_unit_l(1, 1, tn))

    if (inmisc%class_flag(CLASS%RNA)) then
        call simu_force_rna_stack(irep, force_mp_l(1, 1, tn))

        call simu_force_dtrna_stack(irep, force_mp_l(1, 1, tn))
        call simu_force_dtrna_hbond(irep, force_mp_l(1, 1, tn))
    end if

    ! For dna2 (3SPN.2 model)

    if (inmisc%class_flag(CLASS%DNA) .or. inmisc%class_flag(CLASS%DNA2)) then
        call simu_force_dna2_stack(irep, force_mp_l(1, 1, tn))
    end if

    ! Calculate flexible local interactions

    !if (inmisc%i_add_int == 1) then
    if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
        call simu_force_fbangle(irep, force_mp_l(1, 1, tn), &
                                force_mp_mgo(1, 1, 1, 1, tn))
        call simu_force_fdih(irep, force_mp_l(1, 1, tn), &
                             force_mp_mgo(1, 1, 1, 1, tn))
    end if

!$omp master
    TIME_E(tm_force_local)

    TIME_S(tm_force_go)
!$omp end master

    if (inmgo%i_multi_mgo >= 1) then
        call simu_force_nlocal_mgo(irep, force_mp_l(1, 1, tn), &
                                   force_mp_mgo(1, 1, 1, 1, tn), &
                                   ene_unit_l(1, 1, tn))
        call simu_force_nlocal_rna_bp(irep, force_mp_l(1, 1, tn))
    else if (inmisc%force_flag(INTERACT%ENM)) then
        call simu_force_enm(irep, force_mp_l(1, 1, tn))
    else
        call simu_force_nlocal_go(irep, force_mp_l(1, 1, tn))
        call simu_force_nlocal_morse(irep, force_mp_l(1, 1, tn))
        call simu_force_nlocal_rna_bp(irep, force_mp_l(1, 1, tn))
    end if

!$omp master
    TIME_E(tm_force_go)

    TIME_S(tm_force_pnl)
!$omp end master

    if (inmisc%force_flag(INTERACT%EXV)) then
        if (inmisc%i_residue_exv_radii == 0) then
            call simu_force_pnl(irep, force_mp_l(1, 1, tn))
        else if (inmisc%i_residue_exv_radii == 1) then
            call simu_force_pnl_restype(irep, force_mp_l(1, 1, tn))
        endif
    end if
    if (inmisc%force_flag(INTERACT%EXV_WCA)) then
        call simu_force_exv_wca(irep, force_mp_l(1, 1, tn))
    endif

    if (inmisc%class_flag(CLASS%DNA) .or. inmisc%class_flag(CLASS%DNA2)) then
        call simu_force_pnl2(irep, force_mp_l(1, 1, tn))
    end if

    if (inmisc%class_flag(CLASS%ION)) then
        call simu_force_ion(irep, force_mp_l(1, 1, tn))
    endif

    if (inmisc%class_flag(CLASS%DNA) .or. inmisc%class_flag(CLASS%DNA2)) then
        call simu_force_dna2_bp(irep, force_mp_l(1, 1, tn))
        call simu_force_dna2_ex(irep, force_mp_l(1, 1, tn))
    end if

!$omp master
    TIME_E(tm_force_pnl)

    TIME_S(tm_force_ele)
!$omp end master

    if (inmisc%force_flag(INTERACT%ELE)) then
        if (inele%i_calc_method == 0) then
            call simu_force_ele(irep, force_mp_l(1, 1, tn))
        else if (inele%i_calc_method == 1) then
            call simu_force_ele2(irep, force_mp_l(1, 1, tn))
        else if (inele%i_calc_method == 2) then
            call simu_force_ele3(irep, force_mp_l(1, 1, tn))
        end if
    endif

!$omp master
    TIME_E(tm_force_ele)

    ! ------------------------------ PDPWM ------------------------------
    TIME_S(tm_force_pdpwm)
!$omp end master

    if (inmisc%force_flag(INTERACT%PRO_DNA_PWM)) then
        call simu_force_pro_dna_pwm(irep, force_mp_l(1, 1, tn))
    endif

!$omp master
    TIME_E(tm_force_pdpwm)

    ! ------------------------------ PDNS ------------------------------
    TIME_S(tm_force_pdns)
!$omp end master

    if (inmisc%force_flag(INTERACT%PRO_DNA_NONSPEC)) then
        call simu_force_pro_dna_nonspec(irep, force_mp_l(1, 1, tn))
    endif

    ! ----------------------------- STAGE ------------------------------

    if (inmisc%i_stage == 1) then
        call simu_force_stage(irep, force_mp_l(1, 1, tn))
    end if

!$omp master
    TIME_E(tm_force_pdns)

    TIME_S(tm_force_hp)
!$omp end master

    if (inmisc%force_flag(INTERACT%HP)) then
        call simu_force_hp(irep, force_mp_l(1, 1, tn))
    endif
!$omp master
    TIME_E(tm_force_hp)

    TIME_S(tm_force_sasa) !sasa
!$omp end master

    if (inmisc%force_flag(INTERACT%SASA)) then
        call simu_force_sasa(irep, force_mp_l(1, 1, tn))
    endif

!$omp master
    TIME_E(tm_force_sasa) !sasa

    TIME_S(tm_force_gb) !gb
!$omp end master

    if (inmisc%force_flag(INTERACT%GB)) then
        call simu_force_gb(irep, force_mp_l(1, 1, tn))
    endif

!$omp master
    TIME_E(tm_force_gb) !gb
!$omp end master

!$omp end parallel

    TIME_S(tmc_force)
    do n = 1, nthreads - 1
        force_mp_l(1:SPACE_DIM, 1:nmp_all, 0) = &
            force_mp_l(1:SPACE_DIM, 1:nmp_all, 0) + &
            force_mp_l(1:SPACE_DIM, 1:nmp_all, n)
    end do
    TIME_E(tmc_force)

    if (inmgo%i_multi_mgo >= 1) then
        TIME_S(tmc_force)
        do n = 1, nthreads - 1
            ene_unit_l(1:nunit_all, 1:nunit_all, 0) = &
                ene_unit_l(1:nunit_all, 1:nunit_all, 0) + &
                ene_unit_l(1:nunit_all, 1:nunit_all, n)

            do isys = 1, inmgo%nsystem_mgo
                do istat = 1, inmgo%nstate_mgo(isys)
                    force_mp_mgo(1:SPACE_DIM, 1:nmp_all, istat, isys, 0) = &
                        force_mp_mgo(1:SPACE_DIM, 1:nmp_all, istat, isys, 0) + &
                        force_mp_mgo(1:SPACE_DIM, 1:nmp_all, istat, isys, n)
                end do
            end do
        end do

#ifdef MPI_PAR
        call mpi_allreduce(ene_unit_l, ene_unit, nunit_all**2, PREC_MPI, &
                           MPI_SUM, mpi_comm_local, ierr)
#else
        ene_unit(1:nunit_all, 1:nunit_all) = ene_unit_l(1:nunit_all, 1:nunit_all, 0)
#endif
        TIME_E(tmc_force)

        TIME_S(tm_force_go)
        call simu_force_mgo(force_mp_l, force_mp_mgo, ene_unit)
        TIME_E(tm_force_go)
    end if

    TIME_S(tmc_force)
#ifdef MPI_PAR
!!  call mpi_allreduce(force_mp_l, force_mp, SPACE_DIM*nmp_real, PREC_MPI, &
!!                     MPI_SUM, mpi_comm_local, ierr)
    call allreduce(force_mp_l, force_mp)
#else
    force_mp(1:SPACE_DIM, 1:nmp_all) = force_mp_l(1:SPACE_DIM, 1:nmp_all, 0)
#endif
    TIME_E(tmc_force)

    ! ------------------------------- AFM ------------------------------
    if (inmisc%i_afm_fitting == 1) then
        call simu_force_afm_fitting(irep, force_mp)
    end if

    if (inmisc%i_bridge == 1) then
        call simu_force_bridge(irep, force_mp)
    end if

    if (inmisc%i_pulling == 1) then
        call simu_force_pulling(irep, force_mp)
    end if

    if (inmisc%i_anchor == 1) then
        call simu_force_anchor(irep, force_mp)
    end if

    if (inmisc%i_rest1d == 1) then
        call simu_force_rest1d(irep, force_mp)
    end if

    if (inmisc%i_in_box == 1) then
        call simu_force_box(irep, force_mp)
    end if

    if (inmisc%i_in_cap == 1) then
        call simu_force_cap(irep, force_mp)
    end if

    if (inmisc%i_cylinder == 1) then
        call simu_force_cylinder(irep, force_mp)
    end if

    if (inmisc%i_window == 1) then
        call simu_force_window(irep, force_mp)
    end if

    if (inmisc%i_winz == 1) then
        call simu_force_winz(irep, force_mp)
    end if

!! implicit ligand model
    !!if(inimplig%iexe_implig == 1) then
    if (inmisc%i_implig == 1) then
        call simu_force_implig(irep, force_mp)
    end if

#ifdef _DEBUG
    write (*, *) 'simu_force: END'
#endif

!---------------------------------------------------------------------------
#ifdef MPI_PAR
contains

    subroutine allreduce(force_mp_l, force_mp)
        use const_maxsize
        use var_struct, only: nmp_real
        use mpiconst

        implicit none

        real(PREC), intent(in)  :: force_mp_l(SPACE_DIM, nmp_all)
        real(PREC), intent(out) :: force_mp(SPACE_DIM, nmp_all)

        ! communication type
        integer, parameter :: comm_pack_sendrecv = 1
        integer, parameter :: comm_pack_gather = 2
        integer, parameter :: comm_reduce_bcast = 3
        integer, parameter :: comm_allreduce = 4
#ifdef RIKEN_TUNE2
        integer, parameter :: comm_reduce_scatter = 5
        integer, parameter :: comm = 5
#elif  RIKEN_TUNE3
        integer, parameter :: comm = 3
#elif  RIKEN_TUNE4
        integer, parameter :: comm_div_reduce_bcast = 6
        integer, parameter :: comm = 6
#else
        integer, parameter :: comm = 4
#endif

        ! pack
        real(PREC) :: force_mp_p(4, nmp_all)
        real(PREC) :: force_mp_pall(4, 5*nmp_all)
        integer :: nmp_p, nmp_pall(0:npar_mpi - 1)
        integer :: disp(0:npar_mpi - 1), count(0:npar_mpi - 1)
        integer, parameter :: send_tag = 1, recv_tag = 1
        integer :: status(MPI_STATUS_SIZE), n, imp_p
        real(8) :: st
        integer :: imp

#if RIKEN_TUNE2 || RIKEN_TUNE4
        integer :: irank, klen, ksta, kend, iidummy, ircnt, nnprocs
#endif

        select case (comm)
!--------------------------------------------------
        case (comm_pack_sendrecv)
!--------------------------------------------------
            st = mpi_wtime()

            nmp_p = 0
            if (myrank == 0) then
                do imp = 1, nmp_real
                    force_mp(1:SPACE_DIM, imp) = force_mp_l(1:SPACE_DIM, imp)
                end do
            else
                do imp = 1, nmp_real
                    if (any(force_mp_l(1:SPACE_DIM, imp) /= 0.0_PREC)) then
                        nmp_p = nmp_p + 1
                        force_mp_p(1:SPACE_DIM, nmp_p) = force_mp_l(1:SPACE_DIM, imp)
                        force_mp_p(4, nmp_p) = real(imp, PREC)
                    end if
                end do
            end if

            print *, "[time1] ", mpi_wtime() - st
            st = mpi_wtime()

            call mpi_gather(nmp_p, 1, MPI_INTEGER, &
                            nmp_pall, 1, MPI_INTEGER, &
                            0, mpi_comm_local, ierr)

            print *, "[time2] ", mpi_wtime() - st
            st = mpi_wtime()

            if (myrank == 0) then
                print *, "[debug] nmp_p = ", minval(nmp_pall(1:npar_mpi - 1)), &
                    maxval(nmp_pall(1:npar_mpi - 1))
                do n = 1, npar_mpi - 1
                    call mpi_recv(force_mp_p, 4*nmp_pall(n), PREC_MPI, n, recv_tag, mpi_comm_local, status, ierr)
                    do imp_p = 1, nmp_pall(n)
                        imp = int(force_mp_p(4, imp_p))
                        force_mp(1:SPACE_DIM, imp) = force_mp(1:SPACE_DIM, imp) + force_mp_p(1:SPACE_DIM, imp_p)
                    end do
                end do
            else
                call mpi_send(force_mp_p, 4*nmp_p, PREC_MPI, 0, send_tag, mpi_comm_local, status, ierr)
            end if

            print *, "[time3] ", mpi_wtime() - st
            st = mpi_wtime()

            call mpi_bcast(force_mp, SPACE_DIM*nmp_real, PREC_MPI, 0, mpi_comm_local, ierr)

            print *, "[time4] ", mpi_wtime() - st

!--------------------------------------------------
        case (comm_pack_gather)
!--------------------------------------------------
            st = mpi_wtime()

            nmp_p = 0
            do imp = 1, nmp_real
                force_mp(1:SPACE_DIM, imp) = 0.0_PREC
                if (any(force_mp_l(1:SPACE_DIM, imp) /= 0.0_PREC)) then
                    nmp_p = nmp_p + 1
                    force_mp_p(1:SPACE_DIM, nmp_p) = force_mp_l(1:SPACE_DIM, imp)
                    force_mp_p(4, nmp_p) = real(imp, PREC)
                end if
            end do

            print *, "[packing] ", mpi_wtime() - st
            st = mpi_wtime()

            call mpi_allgather(nmp_p, 1, MPI_INTEGER, &
                               nmp_pall, 1, MPI_INTEGER, &
                               mpi_comm_local, ierr)

            print *, "[gather ] ", mpi_wtime() - st

            if (sum(nmp_pall(0:npar_mpi - 1)) > 5*nmp_all) then
                print *, "[error] sum( nmp_pall(0:npar_mpi-1) ) > 5*nmp_all"
            end if

            st = mpi_wtime()

            disp(0) = 0
            count(0) = 4*nmp_pall(0)
            do n = 1, npar_mpi - 1
                disp(n) = disp(n - 1) + 4*nmp_pall(n - 1)
                count(n) = 4*nmp_pall(n)
            end do

            call mpi_allgatherv(force_mp_p, 4*nmp_p, PREC_MPI, &
                                force_mp_pall, count, disp, PREC_MPI, &
                                mpi_comm_local, ierr)

            print *, "[allgath] ", mpi_wtime() - st
            st = mpi_wtime()

            do imp_p = 1, sum(nmp_pall(0:npar_mpi - 1))
                imp = int(force_mp_pall(4, imp_p))
                force_mp(1:SPACE_DIM, imp) = force_mp(1:SPACE_DIM, imp) + force_mp_pall(1:SPACE_DIM, imp_p)
            end do

            print *, "[unpack ] ", mpi_wtime() - st
!    st = mpi_wtime()
!
!    call mpi_bcast(force_mp,SPACE_DIM*nmp_real,PREC_MPI,0,mpi_comm_local,ierr)
!
!    print *,"[time5] ", mpi_wtime()-st

!--------------------------------------------------
        case (comm_reduce_bcast)
!--------------------------------------------------
            call mpi_reduce(force_mp_l, force_mp, SPACE_DIM*nmp_real, PREC_MPI, &
                            MPI_SUM, 0, mpi_comm_local, ierr)

            call mpi_bcast(force_mp, SPACE_DIM*nmp_real, PREC_MPI, 0, mpi_comm_local, ierr)

#ifdef RIKEN_TUNE2
!--------------------------------------------------
        case (comm_reduce_scatter)
!--------------------------------------------------
            do irank = 0, npar_mpi - 1
                klen = (nmp_real - 1 + npar_mpi)/npar_mpi
                ksta = 1 + klen*irank
                kend = min(ksta + klen - 1, nmp_real)
                count(irank) = (kend - ksta + 1)*SPACE_DIM
                disp(irank) = (ksta - 1)*SPACE_DIM
            enddo
            klen = (nmp_real - 1 + npar_mpi)/npar_mpi
            ksta = 1 + klen*local_rank_mpi
            kend = min(ksta + klen - 1, nmp_real)
!
            call mpi_reduce_scatter(force_mp_l, force_mp(1, ksta), count, PREC_MPI, &
                                    MPI_SUM, mpi_comm_local, ierr)
            call mpi_allgatherv(mpi_in_place, iidummy, iidummy, &
                                force_mp(1, 1), count, disp, PREC_MPI, &
                                mpi_comm_local, ierr)
#elif RIKEN_TUNE4
!--------------------------------------------------
        case (comm_div_reduce_bcast)
!--------------------------------------------------
            nnprocs = 27
            do irank = 0, nnprocs - 1
                klen = (nmp_real - 1 + nnprocs)/nnprocs
                ksta = 1 + klen*irank
                kend = min(ksta + klen - 1, nmp_real)
                ircnt = (kend - ksta + 1)*SPACE_DIM
                call mpi_reduce(force_mp_l(1, ksta, 0), force_mp(1, ksta), ircnt, PREC_MPI, &
                                MPI_SUM, 0, mpi_comm_local, ierr)
            enddo
!
            call mpi_bcast(force_mp, SPACE_DIM*nmp_real, PREC_MPI, &
                           0, mpi_comm_local, ierr)
#endif

!--------------------------------------------------
        case default
!--------------------------------------------------
            call mpi_allreduce(force_mp_l, force_mp, SPACE_DIM*nmp_real, PREC_MPI, &
                               MPI_SUM, mpi_comm_local, ierr)
        end select

    end subroutine allreduce
#endif

end subroutine simu_force
