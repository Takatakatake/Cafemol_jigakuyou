subroutine main_modeling(n_ref, xyz_5site, base, aamodel)

  use const_para
  use interface_modeling, only : util_bestfit
  implicit none

!-----------------------------------------------------------
  real(8),intent(in) :: xyz_5site(3,5) 
  integer, intent(in) :: n_ref, base
  real(8),intent(out) :: aamodel(3,22)

!-----------------------------------------------------------
  real(8) :: coord0(3,7),dp(3), ref0(3,3), after_fit(3,7), &
             rmsd, min_f, distance(5), f, bmodel(3,31), pre_aa(3,31) 
  integer :: bestmodel_id, i, j

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=
  ! Arrange CG model
  !=-=-=-=-=-=-=-=-=-=
   
  ! coord0: 3 points for superposition
  ! ref0: reference 3 points for superposition, superpose coord0 on ref0
  ! after_fit: 7 points after superposition
  ! xyz_5site(xyz,1:DS-, 2:DP, 3:DS, 4:DB, 5:DP+)
  ref0=0.0d0
  coord0=0.0d0
  coord0(:,1:5) = xyz_5site(:,:)
  coord0(:,6) = ( xyz_5site(:,2)+xyz_5site(:,3)+xyz_5site(:,5) )/3.0d0

  dp(1) = dot_product(xyz_5site(:,5)-xyz_5site(:,2),coord0(:,6)-xyz_5site(:,2))
  dp(2) = dot_product(xyz_5site(:,5)-xyz_5site(:,2),xyz_5site(:,5)-xyz_5site(:,2)) 
  dp(3) = dot_product(coord0(:,6)-xyz_5site(:,2),coord0(:,6)-xyz_5site(:,2))

  ref0(1,2) = -1.0d0*dp(1)/sqrt(dp(2))
  ref0(2,2) = -1.0d0*sqrt( dp(3) - ( ref0(1,2) )**2 )
  ref0(2,3) = ref0(2,2)

  coord0(:,7) = xyz_5site(:,2) + ( xyz_5site(:,5)-xyz_5site(:,2) ) * &
                dp(1)/dp(2)

  !superposition
  call util_bestfit(ref0(:,:), 3, coord0(:,:), 7, 3, after_fit(:,:), &
                   (/1,2,3/), (/6,2,7/), rmsd)
  !write(*,*) rmsd

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=
  ! Search the best model
  !=-=-=-=-=-=-=-=-=-=-=-=

  bestmodel_id = 0
  min_f = 1000000.0d0

  do i=1,n_ref
    do j=1,5
      distance(j)=dot_product(after_fit(:,j)-dna_ref(:,j,i,base),&
                              after_fit(:,j)-dna_ref(:,j,i,base))
    end do

    f = 0.0d0
    do j=1,5
      f = f+kpara(j)*distance(j)
    end do

    if (f<min_f)then
      min_f = f
      bestmodel_id = i
    end if
  end do
  !write(*,*) bestmodel_id,min_f
!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Suparpose the most suitable structure 
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  bmodel(:,1:28) = dna_ref(:,:,bestmodel_id,base)
  bmodel(:,29:31) = ref0(:,:)

  call util_bestfit(coord0(:,:), 7, bmodel, 31, 3, pre_aa, &
                   (/6,2,7/), (/29,30,31/), rmsd) 
  aamodel(:,:) = pre_aa(:,7:28)

  return 

end subroutine main_modeling 
  

  
