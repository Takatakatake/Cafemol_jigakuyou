! if_mloop
!> @brief Interface module for subroutines called by mloop_simulator

module if_mloop
    interface

        subroutine simu_boxmuller(r_boxmuller, nsize)
            use const_maxsize
            implicit none
            real(PREC), intent(inout) :: r_boxmuller(:)
            integer, intent(in)    :: nsize
        endsubroutine

        subroutine simu_energy_allrep(pnle_unit, pnlet, &
                                      velo_mp, replica_energy, flg_replica, tempk)
            use const_maxsize
            implicit none
            real(PREC), intent(out) :: pnle_unit(:, :, :, :)  ! (MXUNIT, MXUNIT, E_TYPE%MAX, replica)
            real(PREC), intent(out) :: pnlet(:, :)          ! (E_TYPE%MAX,replica)
            real(PREC), intent(in)  :: velo_mp(:, :, :)      ! (3, MXMP, replica)
            real(PREC), intent(out) :: replica_energy(:, :) ! (2, replica)
            real(PREC), intent(in)  :: tempk
            logical, intent(in)  :: flg_replica
        endsubroutine simu_energy_allrep

        subroutine simu_energy(irep, velo_mp, pnlet, pnle_unit)
            use const_maxsize
            implicit none
            integer, intent(in)  :: irep
            real(PREC), intent(in)  :: velo_mp(:, :)      ! (3, nmp_real)
            real(PREC), intent(out) :: pnlet(:)          ! (E_TYPE%MAX)
            real(PREC), intent(out) :: pnle_unit(:, :, :)  ! (nunit_all, nunit_all, E_TYPE%MAX)
        endsubroutine simu_energy

        subroutine write_traject_file(ibefore_time, istep, tempk, velo_mp)
            use const_maxsize
            implicit none
            integer(L_INT), intent(in) :: ibefore_time
            integer(L_INT), intent(in) :: istep
            real(PREC), intent(in) :: tempk
            real(PREC), intent(in) :: velo_mp(:, :, :)
        endsubroutine write_traject_file

        subroutine write_record_file(istep, velo_mp)
            use const_maxsize
            implicit none
            integer(L_INT), intent(in) :: istep
            real(PREC), intent(in) :: velo_mp(:, :, :)
        endsubroutine write_record_file

        subroutine simu_radiusg(rg_unit, rg)
            use const_maxsize
            implicit none
            real(PREC), intent(out) :: rg(:)               ! (replica)
            real(PREC), intent(out) :: rg_unit(:, :)        ! (unit, replica)
        endsubroutine simu_radiusg

        subroutine simu_rmsd(rmsd_unit, rmsd)
            use const_maxsize
            implicit none
            real(PREC), intent(out) :: rmsd(:)             ! (replica)
            real(PREC), intent(out) :: rmsd_unit(:, :)      ! (unit, replica)
        endsubroutine simu_rmsd

        subroutine simu_replica_exchange(velo_mp, replica_energy, tempk)
            use const_maxsize
            implicit none
            real(PREC), intent(inout)  :: velo_mp(:, :, :)      ! (SPACE_DIM, MXMP, replica)
            real(PREC), intent(in)     :: replica_energy(:, :) ! (2, replica)
            real(PREC), intent(in)     :: tempk
        endsubroutine simu_replica_exchange

        subroutine simu_velo_adjst(velo_mp, irep)
            use const_maxsize
            implicit none
            real(PREC), intent(inout) :: velo_mp(:, :, :)
            integer, intent(in)    :: irep
        endsubroutine simu_velo_adjst

!   subroutine simu_velo_mrand(velo_mp, tempk_in)
!      use const_maxsize
!      implicit none
!      real(PREC), intent(out) :: velo_mp(:,:,:)
!      real(PREC), intent(in)  :: tempk_in
!   endsubroutine simu_velo_mrand

        subroutine simu_velo_settemp(velo_mp, irep, tempk)
            use const_maxsize
            implicit none
            real(PREC), intent(in)    :: tempk
            real(PREC), intent(inout) :: velo_mp(:, :, :)
            integer, intent(in)    :: irep
        endsubroutine simu_velo_settemp

        subroutine simu_velo_nosehoover(velo_mp, irep, tempk, velo_yojou)
            use const_maxsize
            implicit none
            real(PREC), intent(in) :: tempk
            integer, intent(in) :: irep
            real(PREC), intent(out) :: velo_yojou
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_velo_nosehoover

        subroutine simu_mc_implig(irep, istep, tempk)
            use const_maxsize
            use var_implig, only: inimplig
            use var_replica, only: n_replica_all
            implicit none
            integer, intent(in) :: irep
            integer(L_INT), intent(in) :: istep
            real(PREC), intent(in) :: tempk
        endsubroutine simu_mc_implig

        subroutine simu_replica_opt_temp(i_current_stage)
            implicit none
            integer, intent(out), optional :: i_current_stage
        endsubroutine simu_replica_opt_temp

!   subroutine simu_initialset_mpc(velo_mp, irep)
!      use const_maxsize
!      implicit none
!      integer,    intent(in) :: irep
!      real(PREC), intent(inout) :: velo_mp(:,:,:)
!   endsubroutine simu_initialset_mpc

        subroutine simu_rotate_velo_mpc(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_rotate_velo_mpc

        subroutine simu_rotate_velo_mpc_rev(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_rotate_velo_mpc_rev

        subroutine simu_velo_correct_simp_mpc(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_velo_correct_simp_mpc

        subroutine simu_rotate_velo_mpc2(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_rotate_velo_mpc2

        subroutine simu_rotate_velo_mpc_rev2(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_rotate_velo_mpc_rev2

        subroutine simu_velo_correct_simp_mpc2(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_velo_correct_simp_mpc2

        subroutine simu_rotate_velo_mpc_rev3(velo_mp, irep, velo_solv_l, velo_mp_l)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
            real(PREC), intent(inout) :: velo_solv_l(:, :)
            real(PREC), intent(inout) :: velo_mp_l(:, :, :)
        endsubroutine simu_rotate_velo_mpc_rev3

        subroutine simu_velo_correct_simp_mpc3(velo_mp, irep, velo_solv_l, velo_mp_l)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
            real(PREC), intent(inout) :: velo_solv_l(:, :)
            real(PREC), intent(inout) :: velo_mp_l(:, :, :)
        endsubroutine simu_velo_correct_simp_mpc3

        subroutine simu_stream_mpc2(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_stream_mpc2

        subroutine simu_stream_mpc4(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_stream_mpc4

        subroutine simu_rotate_velo_mpc_rev4(velo_mp, irep)
            use const_maxsize
            implicit none
            integer, intent(in) :: irep
            real(PREC), intent(inout) :: velo_mp(:, :, :)
        endsubroutine simu_rotate_velo_mpc_rev4

    endinterface
endmodule if_mloop
