subroutine read_references()



  use const_para
  implicit none


!-----------------------------------------------------------

!-----------------------------------------------------------
  real(8) :: tmp_xyz(3)
  integer :: io_status, model_id, seq_id, particle_id
  character(43) ref_line
  character(7) ref_left
!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Read fragment library 
  !=-=-=-=-=-=-=-=-=-=-=-=-=

  allocate(dna_ref(3,28,n_ref_max,4))
  dna_ref=0.0d0
  ! dna_ref : (TCAG, fragment id, xyz, particle type)
  ! particle order T (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N1,C2,O2,N3,C4,O4,C5,C7,C6)
  ! particle order C (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N1,C2,O2,N3,C4,N4,C5,C6)
  ! particle order A (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N9,C8,N7,C5,C6,N6,N1,C2,N3,C4)
  ! particle order G (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4)

  open (refdev, file=ref_title, status='OLD', action='READ', iostat=io_status)
  !write(*,*) 'io_status',io_status
  if (io_status>0) then
    stop 'Cannot open nucleotide library (TCAG_library.txt).'
  end if

  do
    read(refdev,"(a43)",iostat=io_status) ref_line
    if (io_status<0) then
      exit
    end if
    read(ref_line,'(a7,i6,1x,i1,1x,i2,1x,3f8.3)') ref_left, model_id, & 
          seq_id, particle_id, tmp_xyz(1:3)
    dna_ref(:,particle_id,model_id,seq_id) = tmp_xyz(1:3)

  end do
  close(refdev)

  return

end subroutine read_references

