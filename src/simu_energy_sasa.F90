!simu_energy_sasa
!> @brief Calculates the energy related to solvent accessible surface area (sasa)
!>        The values are added into "pnlet(ENERGY%SASA)" and  &
!>        "pnle_unit(ENERGY%SASA)".

subroutine simu_energy_sasa(irep, pnlet)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi, outfile
    use var_neighbor_list, only: sasa_list
    use var_setp, only: insasa
    use var_struct, only: xyz_mp_rep, pxyz_mp_rep, &
        nmp_real, para_sasa, rad_sasa, surf, connect, cmp2atom
    use var_simu, only: istep, sasa, sasa_thread_local
!$  use omp_lib
    use mpiconst

    implicit none

    ! ------------------------------------------------------------------------
    integer, intent(in)  :: irep
    real(PREC), intent(out) :: pnlet(:)         ! (E_TYPE%MAX)
!  real(PREC), intent(out) :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

    ! ------------------------------------------------------------------------
    ! local variables
    integer :: lunout
    integer :: ksta, kend
    integer :: imp1, imp2
    integer :: imirror
    integer :: imp, isasa
    integer :: thread_id, max_threads
    real(PREC) :: dist2(sasa_list(irep)%num_pairs)
    real(PREC) :: dist(sasa_list(irep)%num_pairs)
    real(PREC) :: radsum
    real(PREC) :: v21(SPACE_DIM)
    real(PREC) :: sasa_t, c2, c3ij, c3ji, bij, bji, fij, fji

    lunout = outfile%data

    ksta = 1
    kend = sasa_list(irep)%num_pairs

    sasa_thread_local = 1.0_PREC
!$OMP do private(imp)
    do imp = 1, nmp_real
        if (cmp2atom(imp) == ' CA ' .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ') then
            sasa(imp) = surf(imp)
        end if
    end do
!$OMP end do

!----------------------------------------------------------------------------

    thread_id = 1
    max_threads = 1
!$  max_threads = omp_get_max_threads()

!$OMP do private(isasa,imp1,imp2,v21,dist2,&
!$OMP&           dist,radsum,c2,c3ij,c3ji,bij,bji,fij,fji)
    do isasa = ksta, kend
!$      thread_id = omp_get_thread_num() + 1

        imp1 = sasa_list(irep)%pairs(1, isasa)
        imp2 = sasa_list(irep)%pairs(2, isasa)

        if (cmp2atom(imp1) /= ' CA ' .and. cmp2atom(imp1) /= ' P  ' .and. cmp2atom(imp1) /= ' O  ') then
            write (*, *) '*********************************'
            write (*, *) 'CG atoms other than "CA, O, and P" not supported in sasa'
            stop
        end if

        if (cmp2atom(imp2) /= ' CA ' .and. cmp2atom(imp2) /= ' P  ' .and. cmp2atom(imp2) /= ' O  ') then
            write (*, *) '*********************************'
            write (*, *) 'CG atoms other than "CA, O, and P" not supported in sasa'
            stop
        end if

!     iunit1 = imp2unit(imp1)
!     iunit2 = imp2unit(imp2)

        if (inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
            imirror = sasa_list(irep)%pairs(3, isasa)
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

        dist2(isasa) = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
        dist(isasa) = sqrt(dist2(isasa))
        radsum = rad_sasa(imp1) + rad_sasa(imp2)

        if (dist(isasa) .gt. radsum) cycle

        c2 = radsum - dist(isasa)
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

    if (istep == 0) then
        write (lunout, '(A45)') '**********************************************'
        write (lunout, '(A35)') 'residue No.,  SASA of ith residue'
    end if

!--------------------------------------------------------------------------
    sasa_t = 0.
!$OMP do private(imp)
    do imp = 1, nmp_real
        if (cmp2atom(imp) == ' CA ' .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ') then
            sasa_t = sasa_t + sasa(imp)
            if (istep == 0) then
                write (lunout, '(I10,f12.2)') imp, sasa(imp)
            end if
            pnlet(E_TYPE%SASA) = pnlet(E_TYPE%SASA) + insasa%coef_surf*sasa(imp)
        end if
    end do
!$OMP end do

! The sasa energy is not decomposed into subunit at this moment!!!!!
end subroutine simu_energy_sasa
