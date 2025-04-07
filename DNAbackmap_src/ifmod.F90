module interface_modeling


use const_para
interface

subroutine read_references()
  implicit none
end subroutine read_references

subroutine check_bp(xyz_pdb, seq, chain, n_3site, n_read_model)
  use const_para  
  implicit none
  real(8), intent(in) :: xyz_pdb(3,3,mx3site)
  integer, intent(in) :: seq(mx3site), n_3site
  character(1), intent(in) :: chain(mx3site)
  integer, intent(out) :: n_read_model(mx3site)
end subroutine check_bp

subroutine read_pdb(pdb_titl,xyz_pdb,n_3site,phos_type,res,seq,n_read_model,chain)
  use const_para
  implicit none
  character(maxname), intent(in) :: pdb_titl
  real(8), intent(out) :: xyz_pdb(3, 3, mx3site)
  integer, intent(out) :: n_3site
  integer, intent(out) :: phos_type(2, mx3site),res(mx3site), &
                          seq(mx3site), n_read_model(mx3site)
  character(1), intent(out) :: chain(mx3site)
end subroutine read_pdb

subroutine in_out_file(in_line, in_pdb, out_pdb)
  use const_para
  implicit none
  character(maxname), intent(in) :: in_line
  character(maxname), intent(out) :: in_pdb, out_pdb
end subroutine in_out_file


subroutine main_modeling(n_ref, xyz_5site, base, aamodel)
  implicit none
  real(8),intent(in) :: xyz_5site(3,5) 
  integer, intent(in) :: n_ref, base
  real(8),intent(out) :: aamodel(3,22)
end subroutine main_modeling

subroutine end_modeling(read_list, n_ref, xyz_4site, base, aamodel)
  implicit none
  integer, intent(in) :: read_list(4), n_ref, base
  real(8),intent(in) :: xyz_4site(3,4)
  real(8),intent(out) :: aamodel(3,22)
end subroutine end_modeling

subroutine output_pdb(base, phos_type, res_id, atomid, atoms, chain_id)
  implicit none
  integer, intent(in) :: base, phos_type(2), res_id
  integer, intent(inout) :: atomid
  real(8),intent(in) :: atoms(3,22)
  character(1), intent(in) :: chain_id
end subroutine output_pdb

subroutine util_bestfit(coord1, nat1, coord2, nat2, nat, coord3, &
     list1, list2, rmsd)
  implicit none
  integer, intent(in) :: nat1, nat2, nat
  integer, intent(in) :: list1(nat), list2(nat)
  real(8), intent(in) :: coord1(3, nat1), coord2(3, nat2)
  real(8), intent(out) :: coord3(3, nat2), rmsd
end subroutine util_bestfit


end interface
end module interface_modeling
