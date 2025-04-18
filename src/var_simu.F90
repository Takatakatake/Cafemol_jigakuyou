! var_simu
!> @brief Module defining the global variables used in mloop_simulator.F90

module var_simu

    use const_maxsize
    use const_physical
    use const_index
    implicit none

    ! time step
    integer, save        :: istep_sim, mstep_sim
    integer, save        :: imstep, mstep
    integer(L_INT), save :: istep, ntstep, nstep_opt_temp, ntstep_max
    integer(L_INT), save :: ibefore_time = 0

    !
    integer       :: n_exchange, max_exchange
    integer       :: iopt_stage

    ! physical variables
    real(PREC), save :: tempk
    real(PREC), save :: tstep, tsteph, tstep2
    real(PREC), save :: accelaf(SPACE_DIM)
    real(PREC), allocatable, save  :: velo_mp(:, :, :)  ! (SPACE_DIM, nmp_real, n_replica_mpi)
    real(PREC), allocatable, save  :: accel_mp(:, :, :) ! (SPACE_DIM, nmp_real, n_replica_mpi)
    real(PREC), allocatable, save :: force_mp(:, :)   ! (SPACE_DIM, nmp_all)
    real(PREC), allocatable, save :: rcmass_mp(:)    ! (nmp_all)

    ! mcanonical
    real(PREC), save :: e_md, fac_mmc, em_mid, em_depth, em_sigma
    real(PREC), save :: pnlet_muca(E_TYPE%MAX)
    real(PREC), allocatable, save :: pnle_unit_muca(:, :, :) ! (nunit_all, nunit_all, E_TYPE%MAX)

    ! Wang-Landau MuCa
    real(PREC), save :: muca_cv, muca_coef
    real(PREC), save :: muca_bias(1024), muca_bias_last(1024)
    real(PREC), save :: muca_bias_buf(1024) ! for MPI
    real(PREC), save :: muca_increment

    ! Langevin
    integer, parameter          :: nLAN_CONST = 4
    real(PREC), save              :: r_force(SPACE_DIM)
    real(PREC), save              :: tstep_fric_h, ulconst1, ulconst2
    real(PREC), allocatable, save :: rlan_const(:, :, :) ! (nLAN_CONST, mp, replica)

    ! Nose-Hoover
    ! Now, Replica is not available
    integer, save    :: ics, jcs, ncs
    real(PREC), save :: velo_yojou(MXCS), evcs(MXCS)
    real(PREC), save :: xyz_cs(MXCS), velo_cs(MXCS), cmass_cs(MXCS)

    ! energy
    real(PREC), allocatable, save :: pnlet(:, :)          ! (E_TYPE%MAX, replica)
    real(PREC), allocatable, save :: pnle_unit(:, :, :, :)  ! (unit, unit, E_TYPE%MAX, replica)
    real(PREC), allocatable, save :: qscore(:)           ! (replica)
    real(PREC), allocatable, save :: qscore_unit(:, :, :)  ! (unit, unit, replica)
    real(PREC), allocatable, save :: rg(:)               ! (replica)
    real(PREC), allocatable, save :: rg_unit(:, :)        ! (unit, replica)
    real(PREC), allocatable, save :: rmsd(:)             ! (replica)
    real(PREC), allocatable, save :: rmsd_unit(:, :)      ! (unit, replica)
    real(PREC), allocatable, save :: replica_energy(:, :) ! (2, replica)
#ifdef MPI_PAR
    real(PREC), allocatable, save :: replica_energy_l(:, :) ! (2, replica)
#endif

    ! sasa
    real(PREC), allocatable, save :: sasa(:)          ! (E_TYPE%MAX, replica)
    real(PREC), allocatable, save :: sasa_thread_local(:, :) ! (imp, threadid)

    logical, save :: flg_ppr_release(MXBRIDGE) = .false.

end module var_simu
