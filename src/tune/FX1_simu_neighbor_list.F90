! *********************************************************************
subroutine simu_neighbor_list(irep, ineigh2mp, lmp2neigh)
  
  use if_neighbor
  use const_maxsize
  use const_index
  use var_setp, only : inmisc
  use var_struct, only : nmp_real, lunit2mp, xyz_mp_rep, imp2unit, nmp_all
  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer, intent(in)  :: irep
  integer, intent(out) :: lmp2neigh(MXMP/nthreads      ,0:nthreads-1)
  integer, intent(out) :: ineigh2mp(MXMPNEIGHBOR*nmp_all/nthreads,0:nthreads-1)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, iunit, junit
  integer :: ineighbor(0:nthreads-1)
  integer :: ineigh_unit(MXUNIT, MXUNIT)
  real(PREC) :: dist2
  character(CARRAY_MSG_ERROR) :: error_message
  real(PREC) :: xyz_tmp1,xyz_tmp2,xyz_tmp3

#ifdef MPI_PAR2
  integer :: imp_l
  integer :: klen, ksta, kend
#endif
  integer :: n

  ! -------------------------------------------------------------------
  ! calc neigh_unit
  call simu_neighbor_pre(xyz_mp_rep(:,:,irep), ineigh_unit)

  ! -------------------------------------------------------------------
  ! calc ineigh2mp
  ineighbor(0:nthreads-1) = 0

#ifdef MPI_PAR2
!$omp parallel 
!$omp do private(klen,ksta,kend,imp,iunit,jmp,junit,dist2,xyz_tmp1,xyz_tmp2,xyz_tmp3)
  do n = 0, nthreads-1
  klen=(nmp_l-1+nthreads)/nthreads
  ksta=1+klen*n
  kend=min(ksta+klen-1,nmp_l)

  do imp_l = ksta, kend
     imp = imp_l2g(imp_l)
#else
  n = 0
  do imp = 1, nmp_real - 1
#endif

     xyz_tmp1 = xyz_mp_rep(1, imp, irep)
     xyz_tmp2 = xyz_mp_rep(2, imp, irep)
     xyz_tmp3 = xyz_mp_rep(3, imp, irep)

     iunit = imp2unit(imp)
     jmp = imp + 1
     do while (jmp <= nmp_real)
        junit = imp2unit(jmp)

        if(ineigh_unit(iunit, junit) == 1) then
           dist2 = (xyz_mp_rep(1, jmp, irep) - xyz_tmp1)**2 + &
                   (xyz_mp_rep(2, jmp, irep) - xyz_tmp2)**2 + &
                   (xyz_mp_rep(3, jmp, irep) - xyz_tmp3)**2
           if(dist2 < inmisc%rneighbordist2_unit(iunit, junit)) then
              ineighbor(n) = ineighbor(n) + 1
              ineigh2mp(ineighbor(n),n) = jmp
           end if
        else
           ! jump to last point of 'junit'
           jmp = lunit2mp(2, junit)
        endif

        jmp = jmp + 1
     end do
#ifdef MPI_PAR2
     lmp2neigh(imp_l-ksta+1,n) = ineighbor(n)
  end do
  end do
!$omp end do nowait
!$omp end parallel
#else
     lmp2neigh(imp,n) = ineighbor(n)
  end do
#endif

  !if(any(ineighbor > MXNEIGHBOR/nthreads)) then
  if(any(ineighbor > MXMPNEIGHBOR*nmp_all/nthreads)) then
     error_message = 'Error: too big ineighbor in simu_neighbor_list'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end subroutine simu_neighbor_list
