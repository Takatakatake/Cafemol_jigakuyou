subroutine end_modeling(read_list, n_ref, xyz_4site, base, aamodel)

  use const_para
  use interface_modeling, only : util_bestfit
  implicit none

!-----------------------------------------------------------
  integer, intent(in) :: read_list(4), n_ref, base
  real(8),intent(in) :: xyz_4site(3,4)
  real(8),intent(out) :: aamodel(3,22)

!-----------------------------------------------------------
  real(8) :: distance(6,2), min_f, f, rmsd, pre_aa(3,28)
  integer :: i_ref, i, bestmodel_id

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate distance (CG model)
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  call calc_d(xyz_4site,distance(:,1))

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compare to each structure in 'ACGT_library.txt'
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  bestmodel_id = 0
  min_f = 100000.0d0

  do i_ref=1,n_ref
    call calc_d(dna_ref(:,read_list,i_ref,base),distance(:,2))

    f = 0.0d0
    do i=1,6
      f = f+(distance(i,1)-distance(i,2))**2
    end do
    if (f<min_f) then
      min_f = f
      bestmodel_id = i_ref
    end if
  end do
  !write(*,*) bestmodel_id
  !write(*,*) min_f
!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Suparpose the most suitable structure 
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  call util_bestfit(xyz_4site, 4, dna_ref(:,:,bestmodel_id,base), 28,&
                    4, pre_aa, (/1,2,3,4/), read_list, rmsd)
  aamodel=pre_aa(:,7:28)
  !write(*,*) 'rmsd = ',rmsd
  return 

!===========================================================
!===========================================================
contains

  subroutine calc_d(xyz,dist)

  !---------------------------------------------------------------
  real(8), intent(in) :: xyz(3,4)
  real(8), intent(out) :: dist(6)
  !---------------------------------------------------------------
  integer :: i,j,d_count
  !---------------------------------------------------------------
  d_count=0

  do i=1,3
    do j=i+1,4
      d_count=d_count+1
      dist(d_count)=dot_product(xyz(:,i)-xyz(:,j),&
                                  xyz(:,i)-xyz(:,j))
    end do
  end do

  do i=1,6
    dist(i)=sqrt(dist(i))
  end do
  
end subroutine calc_d
!===========================================================
!===========================================================

end subroutine end_modeling 
