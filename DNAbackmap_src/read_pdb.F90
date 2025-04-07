subroutine read_pdb(pdb_titl,xyz_pdb,n_3site,phos_type,res,seq,n_read_model,chain)



  use const_para
  use interface_modeling, only : check_bp
  implicit none


!-----------------------------------------------------------
  character(maxname), intent(in) :: pdb_titl
  real(8), intent(out) :: xyz_pdb(3, 3, mx3site)
  integer, intent(out) :: n_3site
  integer, intent(out) :: phos_type(2, mx3site),res(mx3site), &
                          seq(mx3site), n_read_model(mx3site)
  character(1), intent(out) :: chain(mx3site)

!-----------------------------------------------------------
  integer :: io_status, i_3site
  integer :: sugar_id, base_id, phos_other_id
  character(54) :: pdb_line
  character(5) :: info_sugar, info_base
  character(1) :: phos_other_chain,tmp_seq

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=
  ! Open .pdb file 
  !=-=-=-=-=-=-=-=-=

  open (pdbdev, file=pdb_titl, status='OLD', action='READ', iostat=io_status)
  if (io_status>0) then
      stop 'Cannot open .pdb file.'
  end if
  write(*,*) pdb_titl
!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=
  ! Count residue number 
  !=-=-=-=-=-=-=-=-=-=-=-=
  ! Count residue number, check missing sugar or missing base.
  ! Residues must be arranged by residue ID.

  sugar_id = 0
  base_id = 0

  do
    read(pdbdev,'(a54)',iostat=io_status) pdb_line
    if (io_status<0) then
      exit
    end if
    
    if (pdb_line(1:4)=='ATOM') then
      ! read chain&residue ID of sugar/base
      if (pdb_line(13:14)=='DS') then
        sugar_id = sugar_id+1 
        info_sugar = pdb_line(22:26) ! chain and residue ID
      else if (pdb_line(13:14)=='DB') then
        base_id = base_id+1 
        info_base = pdb_line(22:26)
      end if
      ! check chain&residue ID
      if ( (pdb_line(13:14)/='DP') .and. &
           (sugar_id==base_id) .and. (info_base/=info_sugar) ) then
        stop 'A sugar or base might be missing. Please arrange beads according to the residue ID.'
      end if

    end if

  end do
  
  ! check number of sugar&base
  if (sugar_id /= base_id) then
     stop 'Missing sugar or base beads'
  end if

  ! record the number of nucleotides
  n_3site = sugar_id  
  !write(*,*)'Totally ',n_3site,' nucleotides.'

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=
  ! Read coordinates 
  !=-=-=-=-=-=-=-=-=-=
  ! read coordinates, check phosphate exist or not
  ! xyz_pdb (nucleotide id,  1:x 2:y 3:z, 1:phosphate 2:sugar 3:base)
  ! phos_type (5' side of sugar, 3' side of sugar)  <- phosphate exist(1) or not(0)
  rewind(pdbdev)
  xyz_pdb = 0.0d0
  sugar_id = 0
  base_id = 0
  phos_type = 0


  do
    read(pdbdev,'(a54)',iostat=io_status) pdb_line
    if (io_status<0) then
      exit
    end if

    if (pdb_line(1:4)=='ATOM') then

      !phosphate
      if (pdb_line(13:14)=='DP') then
        read(pdb_line(22:26),'(a1,i4)') phos_other_chain, phos_other_id
        if ((phos_other_id==res(sugar_id)).and.(phos_other_chain==chain(sugar_id))) then
          phos_type(1,sugar_id) = 1
          read(pdb_line(31:54),'(3f8.3)') xyz_pdb(1:3,1,sugar_id)
        else
          read(pdb_line(31:54),'(3f8.3)') xyz_pdb(1:3,1,sugar_id+1)
        end if

      ! sugar
      else if (pdb_line(13:14)=='DS') then
        sugar_id = sugar_id+1
        read(pdb_line(19:54),'(a1,2x,a1,i4,4x,3f8.3)') & 
             tmp_seq,chain(sugar_id),res(sugar_id),xyz_pdb(1:3,2,sugar_id)
        ! sequence T:1 C:2 A:3 G:4
        if (tmp_seq=='T') then
          seq(sugar_id) = 1
        else if (tmp_seq=='C') then
          seq(sugar_id) = 2
        else if (tmp_seq=='A') then
          seq(sugar_id) = 3
        else if (tmp_seq=='G') then
          seq(sugar_id) = 4
        end if
        ! check phosphate (5' side of sugar) 
        if ((res(sugar_id)==phos_other_id) .and. &
            (chain(sugar_id)==phos_other_chain)) then
          phos_type(1,sugar_id) = 1
        end if

      ! base
      else if (pdb_line(13:14)=='DB') then
        base_id = base_id+1
        read(pdb_line(31:54),'(3f8.3)') xyz_pdb(1:3,3,base_id)

      end if

    end if

  end do

  !check 3' side phosphate
  do i_3site=1,n_3site-1
    if ((res(i_3site+1)==res(i_3site)+1).and.(chain(i_3site+1)==chain(i_3site))&
          .and.(phos_type(1,i_3site+1)==1))then
        phos_type(2,i_3site) = 1
    end if
  end do

  close(pdbdev)

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Read Watson-Click pairing information
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  call check_bp(xyz_pdb, seq, chain, n_3site, n_read_model)


  return

end subroutine read_pdb


