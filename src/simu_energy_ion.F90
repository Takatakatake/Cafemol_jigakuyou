! simu_energy_ion
!> @brief Calculate all the energy terms related to ion, &
!>        including LJ, hydration term, excluded volume term

subroutine simu_energy_ion(irep, pnle_unit, pnlet)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: exv_ion_list, hyd_ion_list
    use var_setp, only: inmisc, inion
    use var_struct, only: imp2unit, xyz_mp_rep, pxyz_mp_rep, iontype_mp
#ifdef MPI_PAR3
    use mpiconst
#endif

    implicit none

! ------------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout)   :: pnle_unit(:, :, :) ! (MXUNIT, MXUNIT, E_TYPE%MAX)
    real(PREC), intent(inout)   :: pnlet(:)         ! (E_TYPE%MAX)

    ! ------------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: imp1, imp2, iunit, junit, itype1, itype2
    integer :: ipnl, imirror
    real(PREC) :: dist1, dist2, cdist2, ddist1, ddist2
    real(PREC) :: roverdist2, roverdist6, roverdist12
    real(PREC) :: pnl, coef, cutoff2
    real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
    integer :: klen
#endif

    ! ------------------------------------------------------------------------

    if (inmisc%class_flag(CLASS%ION)) then

        ! ------------------------------------------------------------------------
        ! LJ and hydration interaction of ion
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
        klen = (hyd_ion_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, hyd_ion_list(irep)%num_pairs)
#else
        ksta = 1
        kend = hyd_ion_list(irep)%num_pairs
#endif
#else
        ksta = 1
        kend = hyd_ion_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,itype1,itype2,v21,dist2,dist1,&
!$OMP&           ddist1,ddist2,roverdist2, &
!$OMP&           roverdist6,roverdist12,pnl,&
!$OMP&           iunit,junit,imirror)
        do ipnl = ksta, kend

            imp1 = hyd_ion_list(irep)%pairs(1, ipnl)
            imp2 = hyd_ion_list(irep)%pairs(2, ipnl)
            itype1 = iontype_mp(imp1)
            itype2 = iontype_mp(imp2)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = hyd_ion_list(irep)%pairs(3, ipnl)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 > inion%cutofflj2(itype1, itype2)) cycle

            ! --------------------------------------------------------------------
            roverdist2 = inion%cdistlj2(itype1, itype2)/dist2
            roverdist6 = roverdist2*roverdist2*roverdist2
            roverdist12 = roverdist6*roverdist6

            dist1 = sqrt(dist2)
            ddist1 = dist1 - inion%cdistmh1(itype1, itype2)
            ddist2 = dist1 - inion%cdistmh2(itype1, itype2)

            pnl = inion%clj_energy(itype1, itype2)*(roverdist12 - roverdist6) &
                  + inion%cmh1_energy(itype1, itype2)*exp(-inion%rsigmamh1(itype1, itype2)*ddist1**2) &
                  + inion%cmh2_energy(itype1, itype2)*exp(-inion%rsigmamh2(itype1, itype2)*ddist2**2)

            ! --------------------------------------------------------------------
            ! sum of the energy
            pnlet(E_TYPE%HYD_ION) = pnlet(E_TYPE%HYD_ION) + pnl

            iunit = imp2unit(imp1)
            junit = imp2unit(imp2)
            pnle_unit(iunit, junit, E_TYPE%HYD_ION) = pnle_unit(iunit, junit, E_TYPE%HYD_ION) + pnl

        end do
!$OMP end do nowait

        ! ------------------------------------------------------------------------
        ! exvol ion
        ! for speed up
        cdist2 = inion%cdist_exv_ion**2
        cutoff2 = cdist2*inion%cutoff_exv_ion**2
        coef = 4.0e0_PREC*inion%cexv_ion

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
        klen = (exv_ion_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, exv_ion_list(irep)%num_pairs)
#else
        ksta = 1
        kend = exv_ion_list(irep)%num_pairs
#endif
#else
        ksta = 1
        kend = exv_ion_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,roverdist2, &
!$OMP&           roverdist6,roverdist12,pnl,iunit,junit,imirror)
        do ipnl = ksta, kend

            imp1 = exv_ion_list(irep)%pairs(1, ipnl)
            imp2 = exv_ion_list(irep)%pairs(2, ipnl)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = exv_ion_list(irep)%pairs(3, ipnl)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 > cutoff2) cycle

            ! --------------------------------------------------------------------
            roverdist2 = cdist2/dist2
            roverdist6 = roverdist2*roverdist2*roverdist2
            roverdist12 = roverdist6*roverdist6
            pnl = coef*(roverdist12 - roverdist6 + 0.25e0_PREC)

            ! --------------------------------------------------------------------
            ! sum of the energy
            pnlet(E_TYPE%EXV_ION) = pnlet(E_TYPE%EXV_ION) + pnl

            iunit = imp2unit(imp1)
            junit = imp2unit(imp2)
            pnle_unit(iunit, junit, E_TYPE%EXV_ION) = pnle_unit(iunit, junit, E_TYPE%EXV_ION) + pnl
        end do
!$OMP end do nowait

    endif  ! inmisc%class_flag(CLASS%ION)

end subroutine simu_energy_ion
