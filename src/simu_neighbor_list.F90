!simu_neighbor_list
!> @brief Constructs the neighboring list for non-local interactions. &
!>        In the runup to caluculations, this subroutine calls        &
!>        "simu_neighbor_pre" to specify which unit pair(s) contain   &
!>        non-local contacts.

subroutine simu_neighbor_list(irep, num_pairs, pairs)

    use if_neighbor
    use const_maxsize
    use const_index
    use var_inp, only: inperi
    use var_setp, only: inmisc
    use var_struct, only: lunit2mp, xyz_mp_rep, pxyz_mp_rep, &
        imp2unit, nmp_real
    use mpiconst

    implicit none

    ! -------------------------------------------------------------------
    integer, intent(in)  :: irep
    integer, intent(out) :: num_pairs(0:nthreads - 1)
    integer, intent(out) :: pairs(2, MXMPNEIGHBOR*nmp_real/nthreads, 0:nthreads-1)

    ! -------------------------------------------------------------------
    ! local variables
    integer :: n
    integer :: klen, ksta, kend
    integer :: imp, jmp, iunit, junit
    integer :: imirror
    integer :: ineigh_unit(MXUNIT, MXUNIT)
    real(PREC) :: dist2, v21(3)
    character(CARRAY_MSG_ERROR) :: error_message
    integer :: imp_l
    integer :: all_pairs(2, nmp_real*(nmp_real - 1)/2)

    ! -------------------------------------------------------------------
    ! calc neigh_unit
    call simu_neighbor_pre(xyz_mp_rep(:, :, irep), ineigh_unit)

    n = 0
    do imp = 1, nmp_real-1
        do jmp = imp+1, nmp_real
            n = n + 1
            all_pairs(:, n) = [imp, jmp]
        end do
    end do

    num_pairs(:) = 0

!$OMP parallel
!$OMP do private(klen,ksta,kend,imp_l,imp,iunit,jmp,junit,dist2,v21,imirror)
    do n = 0, nthreads - 1
        klen = (nmp_real*(nmp_real - 1)/2 - 1 + nthreads)/nthreads
        ksta = 1 + klen*n
        kend = min(ksta + klen - 1, nmp_real*(nmp_real - 1)/2)

        do imp_l = ksta, kend
            imp = all_pairs(1, imp_l)
            jmp = all_pairs(2, imp_l)

            iunit = imp2unit(imp)
            junit = imp2unit(jmp)

            if (ineigh_unit(iunit, junit) /= 1) cycle

            if (inperi%i_periodic == 0) then
                v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
            else
                v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
                call util_pbneighbor(v21, imirror)
            end if

            dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

            if (dist2 < inmisc%rneighbordist2_unit(iunit, junit)) then
                num_pairs(n) = num_pairs(n) + 1
                pairs(1:2, num_pairs(n), n) = [imp, jmp]
            end if
        end do
    end do
!$OMP end do nowait
!$OMP end parallel

    if (any(num_pairs > MXMPNEIGHBOR*nmp_real/nthreads)) then
        error_message = 'Error: too many neighbor pairs in simu_neighbor_list'
        call util_error(ERROR%STOP_ALL, error_message)
    end if

end subroutine simu_neighbor_list
