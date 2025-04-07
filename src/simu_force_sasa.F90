!simu_force_sasa
!> @brief Calculates the force related to solvent accessible surface area (sasa)
!>        The values are added into "pnlet(ENERGY%SASA)" and  &
!>        "pnle_unit(ENERGY%SASA)".

subroutine simu_force_sasa(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: sasa_list
    use var_setp, only: insasa
    use var_struct, only: nmp_all, xyz_mp_rep, pxyz_mp_rep, &
        nmp_real, para_sasa, rad_sasa, surf, connect, cmp2atom
    use var_simu, only: sasa, sasa_thread_local
!$  use omp_lib
    use mpiconst

    implicit none

    ! --------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: imp1, imp2, imp, isasa
    integer :: imirror
    integer :: thread_id, max_threads
    real(PREC) :: c4ij, c4ji, fij, fji
    real(PREC) :: c2, c3ij, c3ji, bij, bji
    real(PREC) :: forij, forji
    real(PREC) :: dist(sasa_list(irep)%num_pairs)
    real(PREC) :: dist2(sasa_list(irep)%num_pairs)
    real(PREC) :: v21(SPACE_DIM, sasa_list(irep)%num_pairs)
    real(PREC) :: radsum(sasa_list(irep)%num_pairs)
    real(PREC) :: radsum2(sasa_list(irep)%num_pairs)

    ksta = 1
    kend = sasa_list(irep)%num_pairs

    sasa_thread_local = 1.0_PREC
!$OMP do private(imp)
    do imp = 1, nmp_real
        if (cmp2atom(imp) == ' CA ' .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ') sasa(imp) = surf(imp)
    end do
!$OMP end do

!----------------------------------------------------------------------------

    thread_id = 1
    max_threads = 1
!$  max_threads = omp_get_max_threads()

!$OMP do private(isasa,imp1,imp2,&
!$OMP&           c2,c3ij,c3ji,bij,bji,fij,fji)
    do isasa = ksta, kend
!$      thread_id = omp_get_thread_num() + 1

        imp1 = sasa_list(irep)%pairs(1, isasa)
        imp2 = sasa_list(irep)%pairs(2, isasa)

        if (inperi%i_periodic == 0) then
            v21(1:3, isasa) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
            imirror = sasa_list(irep)%pairs(3, isasa)
            v21(1:3, isasa) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

        dist2(isasa) = v21(1, isasa)*v21(1, isasa) + v21(2, isasa)*v21(2, isasa) + v21(3, isasa)*v21(3, isasa)
        radsum(isasa) = rad_sasa(imp1) + rad_sasa(imp2)
        radsum2(isasa) = radsum(isasa)*radsum(isasa)

        if (dist2(isasa) .gt. radsum2(isasa)) cycle

        dist(isasa) = sqrt(dist2(isasa))
        c2 = radsum(isasa) - dist(isasa)
        c3ij = (1.0e0_PREC + (rad_sasa(imp2) - rad_sasa(imp1))/dist(isasa))
        c3ji = (1.0e0_PREC + (rad_sasa(imp1) - rad_sasa(imp2))/dist(isasa))

        bij = para_sasa(imp1)*c2*c3ij
        bji = para_sasa(imp2)*c2*c3ji

        fij = 1.0e0_PREC - connect(imp1 - imp2)*bij
        fji = 1.0e0_PREC - connect(imp2 - imp1)*bji

        sasa_thread_local(imp1, thread_id) = sasa_thread_local(imp1, thread_id)*fij
        sasa_thread_local(imp2, thread_id) = sasa_thread_local(imp2, thread_id)*fji
    end do
!$OMP end do

!$OMP do private(imp, thread_id)
    do imp = 1, nmp_real
        do thread_id = 1, max_threads
            sasa(imp) = sasa(imp) * sasa_thread_local(imp, thread_id)
        end do
    end do
!$OMP end do


!-------------------------------------------------------------------------------
!$OMP do private(isasa,imp1,imp2,&
!$OMP&           c2,c3ij,c3ji,bij,bji,fij,fji,c4ij,c4ji,forij,forji)

    do isasa = ksta, kend

        imp1 = sasa_list(irep)%pairs(1, isasa)
        imp2 = sasa_list(irep)%pairs(2, isasa)

        if (dist2(isasa) .gt. radsum2(isasa)) cycle

        c2 = radsum(isasa) - dist(isasa)
        c3ij = (1.0e0_PREC + (rad_sasa(imp2) - rad_sasa(imp1))/dist(isasa))
        c3ji = (1.0e0_PREC + (rad_sasa(imp1) - rad_sasa(imp2))/dist(isasa))

        bij = para_sasa(imp1)*c2*c3ij
        bji = para_sasa(imp2)*c2*c3ji

        fij = 1.0e0_PREC - connect(imp1 - imp2)*bij
        fji = 1.0e0_PREC - connect(imp2 - imp1)*bji

        c4ij = para_sasa(imp1)*connect(imp1 - imp2)* &
               (c3ij + c2*(rad_sasa(imp2) - rad_sasa(imp1))/dist2(isasa))
        c4ji = para_sasa(imp2)*connect(imp2 - imp1)* &
               (c3ji + c2*(rad_sasa(imp1) - rad_sasa(imp2))/dist2(isasa))

        forij = insasa%coef_surf*sasa(imp1)*c4ij/fij/dist(isasa)
        forji = insasa%coef_surf*sasa(imp2)*c4ji/fji/dist(isasa)

! diagonal term

        force_mp(1:3, imp1) = force_mp(1:3, imp1) + forij*v21(1:3, isasa)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) - forji*v21(1:3, isasa)

! off-diagonal term

        force_mp(1:3, imp1) = force_mp(1:3, imp1) + forji*v21(1:3, isasa)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) - forij*v21(1:3, isasa)

    end do

!$OMP end do

end subroutine simu_force_sasa
