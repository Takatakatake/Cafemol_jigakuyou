subroutine in_out_file(in_line, in_pdb, out_pdb)

  use const_para
  implicit none

!-----------------------------------------------------------
  character(maxname), intent(in) :: in_line
  character(maxname), intent(out) :: in_pdb, out_pdb

!-----------------------------------------------------------
  integer :: char_term(4), current_point, i

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=
  ! Read filenames 
  !=-=-=-=-=-=-=-=-=

   char_term=0
   current_point=1
   
   do i=9,maxname
     if ( (in_line(i:i)/=' ') .and. (mod(current_point,2)==1) ) then
       char_term(current_point) = i
       current_point = current_point+1
     else if ( (in_line(i:i)==' ') .and. (mod(current_point,2)==0) ) then
       char_term(current_point) = i-1
       current_point = current_point+1
     end if
   end do

   in_pdb = in_line(char_term(1):char_term(2))
   out_pdb = in_line(char_term(3):char_term(4))

  return

end subroutine in_out_file  

