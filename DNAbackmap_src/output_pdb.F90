subroutine output_pdb(base, phos_type, res_id, atomid, atoms, chain_id)

  use const_para
  implicit none

!-----------------------------------------------------------
  integer, intent(in) :: base, phos_type(2), res_id
  integer, intent(inout) :: atomid
  real(8),intent(in) :: atoms(3,22)
  character(1), intent(in) :: chain_id

!-----------------------------------------------------------
  integer,parameter :: n_atoms(4) =(/20, 19, 21, 22/)
  integer :: i
  character(1),parameter :: base_type(4) = (/'T','C','A','G'/)

!-----------------------------------------------------------
  !=-=-=-=-=
  ! Output
  !=-=-=-=-=

   do i=1,n_atoms(base)
     if ((phos_type(1)==1) .or. (i>3)) then
       atomid=atomid+1
       write(outdev,'(a4,2x,i5,2x,a3,2x,a1,a1,1x,a1,i4,4x,3f8.3,a23,a1)')&
           'ATOM',atomid,particle_name(i,base),'D',base_type(base), &
           chain_id,res_id,atoms(:,i),'  1.00  0.00           ', &
           particle_name(i,base)(1:1)
     end if
   end do
   if (phos_type(2)==0) then
     write(outdev,'(a3)')'TER'
   end if

  return

end subroutine output_pdb    

