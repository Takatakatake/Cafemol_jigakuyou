! simu_force_pnl2
!> @brief Calculate forces generated by DNA-related terms, &
!>        including base stacking, base pairing, mismatched &
!>        base pairing, excluded volume term, electrostatic &
!>        interaction,  and solvation term

subroutine simu_force_pnl2(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_inp, only: inperi
    use var_neighbor_list, only: bp_at_list, bp_gc_list, bp_mismatch_list, &
        exv_dna_list, solv_list
    use var_setp, only: indna, inmisc
    use var_struct, only: xyz_mp_rep, pxyz_mp_rep, &
        nstack, istack2mp, stack_nat, nmp_all
    use var_replica, only: irep2grep
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! --------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ksta, kend
    integer :: imp1, imp2
    integer :: grep
    integer :: ipnl, istack, isolv, imirror
    real(PREC) :: dist1, dist2, cdist2, edist1, sdist
    real(PREC) :: roverdist2, roverdist4, roverdist8
    real(PREC) :: roverdist12, roverdist14
    real(PREC) :: dvdw_dr, coef, cutoff2, salpha
    real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
    integer :: klen
#endif

#ifdef _DEBUG
    write (*, *) '#### start simu_force_pnl2'
    write (*, *) 'nstack, ', nstack
#endif

    if (inmisc%class_flag(CLASS%DNA)) then
        ! --------------------------------------------------------------------
        ! base stacking DNA
        coef = 24.0e0_PREC*indna%cstack
#ifdef MPI_PAR
        klen = (nstack - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, nstack)
#ifdef MPI_DEBUG
        print *, "pnl2_1       = ", kend - ksta + 1
#endif
#else
        ksta = 1
        kend = nstack
#endif
!$OMP do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$OMP&           roverdist8,roverdist14,dvdw_dr,for,imirror)
        do istack = ksta, kend

            imp1 = istack2mp(1, istack)
            imp2 = istack2mp(2, istack)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
                call util_pbneighbor(v21, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            ! -----------------------------------------------------------------
            roverdist2 = stack_nat(istack)/dist2
            roverdist4 = roverdist2*roverdist2
            roverdist8 = roverdist4*roverdist4
            roverdist14 = roverdist2*roverdist4*roverdist8
            dvdw_dr = coef/stack_nat(istack)*(2.0e0_PREC*roverdist14 - roverdist8)
            if (dvdw_dr > DE_MAX) then
                dvdw_dr = DE_MAX
                !write (*, *) "base stacking DNA", imp1, imp2, dvdw_dr
            end if
            !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

            for(1:3) = dvdw_dr*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        end do
!$OMP end do nowait

        ! --------------------------------------------------------------------
        ! AT base pair DNA
        ! for speed up
        cutoff2 = (indna%cutoff_bp*indna%cdist_bp_at)**2
        cdist2 = indna%cdist_bp_at**2
        coef = 240.0e0_PREC*indna%cbp_at/cdist2
#ifdef _DEBUG
        write (*, *) 'AT base pair DNA'
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
        klen = (bp_at_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, bp_at_list(irep)%num_pairs)
#else
        ksta = 1
        kend = bp_at_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
        print *, "pnl2_2       = ", kend - ksta + 1
#endif
#else
        ksta = 1
        kend = bp_at_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$OMP&           roverdist8,roverdist12,roverdist14,dvdw_dr,for,imirror)
        do ipnl = ksta, kend

            imp1 = bp_at_list(irep)%pairs(1, ipnl)
            imp2 = bp_at_list(irep)%pairs(2, ipnl)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = bp_at_list(irep)%pairs(3, ipnl)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
            if (dist2 > cutoff2) cycle

            ! -----------------------------------------------------------------
            roverdist2 = cdist2/dist2
            roverdist4 = roverdist2*roverdist2
            roverdist8 = roverdist4*roverdist4
            roverdist12 = roverdist4*roverdist8
            roverdist14 = roverdist2*roverdist12
            dvdw_dr = coef*(roverdist14 - roverdist12)
            if (dvdw_dr > DE_MAX) then
!           write (*, *) "AT base pair DNA", imp1, imp2, dvdw_dr
                dvdw_dr = DE_MAX
            end if
            !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

            for(1:3) = dvdw_dr*v21(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
        end do
!$OMP end do nowait

        ! --------------------------------------------------------------------
        ! GC base pair DNA
        ! for speed up
        cutoff2 = (indna%cutoff_bp*indna%cdist_bp_gc)**2
        cdist2 = indna%cdist_bp_gc**2
        coef = 240.0e0_PREC*indna%cbp_gc/cdist2
#ifdef _DEBUG
        write (*, *) 'GC base pair DNA'
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
        klen = (bp_gc_list(irep)%num_pair - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, bp_gc_list(irep)%num_pairs)
#else
        ksta = 1
        kend = bp_gc_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
        print *, "pnl2_3       = ", kend - ksta + 1
#endif
#else
        ksta = 1
        kend = bp_gc_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$OMP&           roverdist8,roverdist12,roverdist14,dvdw_dr,for,imirror)
        do ipnl = ksta, kend

            imp1 = bp_gc_list(irep)%pairs(1, ipnl)
            imp2 = bp_gc_list(irep)%pairs(2, ipnl)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = bp_gc_list(irep)%pairs(3, ipnl)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 > cutoff2) cycle

            ! -----------------------------------------------------------------
            roverdist2 = cdist2/dist2
            roverdist4 = roverdist2*roverdist2
            roverdist8 = roverdist4*roverdist4
            roverdist12 = roverdist4*roverdist8
            roverdist14 = roverdist2*roverdist12
            dvdw_dr = coef*(roverdist14 - roverdist12)
            if (dvdw_dr > DE_MAX) then
!           write (*, *) "GC base pair DNA", imp1, imp2, dvdw_dr
                dvdw_dr = DE_MAX
            end if
            !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

            for(1:3) = dvdw_dr*v21(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
        end do
!$OMP end do nowait

        ! --------------------------------------------------------------------
        ! mismatch base pair DNA
        ! for speed up
        cutoff2 = (indna%cutoff_mbp*indna%cdist_mbp)**2
        cdist2 = indna%cdist_mbp**2
        coef = 24.0e0_PREC*indna%cmbp/cdist2
#ifdef _DEBUG
        write (*, *) 'mismatch base pair DNA'
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
        klen = (bp_mismatch_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, bp_mismatch_list(irep)%num_pairs)
#else
        ksta = 1
        kend = bp_mismatch_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
        print *, "pnl2_4       = ", kend - ksta + 1
#endif
#else
        ksta = 1
        kend = bp_mismatch_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$OMP&           roverdist8,roverdist14,dvdw_dr,for,imirror)
        do ipnl = ksta, kend

            imp1 = bp_mismatch_list(irep)%pairs(1, ipnl)
            imp2 = bp_mismatch_list(irep)%pairs(2, ipnl)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = bp_mismatch_list(irep)%pairs(3, ipnl)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 > cutoff2) cycle

            ! -----------------------------------------------------------------
            roverdist2 = cdist2/dist2
            roverdist4 = roverdist2*roverdist2
            roverdist8 = roverdist4*roverdist4
            roverdist14 = roverdist2*roverdist4*roverdist8
            dvdw_dr = coef*(2.0e0_PREC*roverdist14 - roverdist8)
            if (dvdw_dr > DE_MAX) then
                !write (*, *) "mismatch base pair DNA", imp1, imp2, dvdw_dr
                dvdw_dr = DE_MAX
            end if
            !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

            for(1:3) = dvdw_dr*v21(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
        end do
!$OMP end do nowait

        ! --------------------------------------------------------------------
        ! exvol DNA
        ! for speed up
        cutoff2 = (indna%cutoff_exv_dna*indna%cdist_exv_dna)**2
        cdist2 = indna%cdist_exv_dna**2
        coef = 24.0e0_PREC*indna%cexv_dna/cdist2
#ifdef _DEBUG
        write (*, *) 'exvol DNA'
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
        klen = (exv_dna_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, exv_dna_list(irep)%num_pairs)
#else
        ksta = 1
        kend = exv_dna_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
        print *, "pnl2_5       = ", kend - ksta + 1
#endif
#else
        ksta = 1
        kend = exv_dna_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$OMP&           roverdist8,roverdist14,dvdw_dr,for,imirror)
        do ipnl = ksta, kend

            imp1 = exv_dna_list(irep)%pairs(1, ipnl)
            imp2 = exv_dna_list(irep)%pairs(2, ipnl)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = exv_dna_list(irep)%pairs(3, ipnl)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 > cutoff2) cycle

            ! -----------------------------------------------------------------
            roverdist2 = cdist2/dist2
            roverdist4 = roverdist2*roverdist2
            roverdist8 = roverdist4*roverdist4
            roverdist14 = roverdist2*roverdist4*roverdist8
            dvdw_dr = coef*(2.0e0_PREC*roverdist14 - roverdist8)
            if (dvdw_dr > DE_MAX) then
                dvdw_dr = DE_MAX
                !write (*, *) "exvol DNA", imp1, imp2, dvdw_dr
            end if
            !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

            for(1:3) = dvdw_dr*v21(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
        end do
!$OMP end do nowait

    endif

    ! --------------------------------------------------------------------
    ! solvation DNA
    if (inmisc%force_flag(INTERACT%DNA)) then
        grep = irep2grep(irep)
        ! for speed up
        salpha = 1.0/indna%cralpha_solv_dna
        sdist = indna%cdist_solv_dna
        cutoff2 = (indna%cutoff_solv_dna/salpha + sdist)**2
        coef = 2.0e0_PREC*salpha*indna%coef_solv_dna(grep)
#ifdef _DEBUG
        write (*, *) 'solvation DNA'
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_SOLV
        klen = (solv_list(irep)%num_pairs - 1 + npar_mpi)/npar_mpi
        ksta = 1 + klen*local_rank_mpi
        kend = min(ksta + klen - 1, solv_list(irep)%num_pairs)
#else
        ksta = 1
        kend = solv_list(irep)%num_pairs
#endif
#ifdef MPI_DEBUG
        print *, "pnl2_7       = ", kend - ksta + 1
#endif
#else
        ksta = 1
        kend = solv_list(irep)%num_pairs
#endif
!$OMP do private(imp1,imp2,v21,dist2,dist1,edist1, &
!$OMP&           dvdw_dr,for,imirror)
        do isolv = ksta, kend
            imp1 = solv_list(irep)%pairs(1, isolv)
            imp2 = solv_list(irep)%pairs(2, isolv)

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            else
                imirror = solv_list(irep)%pairs(3, isolv)
                v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
            end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 > cutoff2) cycle

            ! -----------------------------------------------------------------
            dist1 = sqrt(dist2)
            edist1 = exp(-salpha*(dist1 - sdist))

            dvdw_dr = -coef/dist1*edist1*(1.0 - edist1)
            if (dvdw_dr > DE_MAX) then
                !write (*, *) "solvation DNA", imp1, imp2, dvdw_dr
                dvdw_dr = DE_MAX
            end if
            !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

            for(1:3) = dvdw_dr*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        end do
!$OMP end do nowait

    endif

end subroutine simu_force_pnl2
