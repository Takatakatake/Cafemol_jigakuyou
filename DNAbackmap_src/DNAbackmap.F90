program DNAbackmap

  use const_para
  use interface_modeling, only : read_references, read_pdb, in_out_file, &
                                 end_modeling, main_modeling, output_pdb
  implicit none

!-----------------------------------------------------------
  real(8) :: xyz_pdb(3,3,mx3site)
  real(8) :: xyz_4site(3,4), xyz_5site(3,5), aamodel(3,22)
  integer :: phos_type(2,mx3site),res(mx3site), &
             seq(mx3site), n_read_model(mx3site)
  integer :: n_3site, io_status, atomid, i_3site, read_list(4)
  character(1) :: chain(mx3site) 
  character(maxname) :: in_title,in_line,in_pdb,out_pdb
  ! for calculation time
  integer :: t0,tfin, count_rate, total_time 

!-----------------------------------------------------------
  call system_clock(t0)
!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=
  ! Set atom name
  !=-=-=-=-=-=-=-=-=-=-=

  ! dna_ref : (xyz, particle type, fragment id, TCAG)
  ! particle order T (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N1,C2,O2,N3,C4,O4,C5,C7,C6)
  ! particle order C (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N1,C2,O2,N3,C4,N4,C5,C6)
  ! particle order A (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N9,C8,N7,C5,C6,N6,N1,C2,N3,C4)
  ! particle order G (DS-,DP,DS,DB,DP+,DS+,P,OP1,OP2,O5',C5',C4',O4',C3',O3',C2',C1',N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4)

   particle_name = '   '

   particle_name(:,1) = (/"P  ","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'", &
                          "C2'","C1'","N1 ","C2 ","O2 ","N3 ","C4 ","O4 ","C5 ", &
                          "C7 ","C6 ","   ","   "/)

   particle_name(:,2) = (/"P  ","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'", &
                          "C2'","C1'","N1 ","C2 ","O2 ","N3 ","C4 ","N4 ","C5 ", &
                          "C6 ","   ","   ","   "/)

   particle_name(:,3) = (/"P  ","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'", &
                          "C2'","C1'","N9 ","C8 ","N7 ","C5 ","C6 ","N6 ","N1 ", &
                          "C2 ","N3 ","C4 ","   "/)

   particle_name(:,4) = (/"P  ","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'", &
                          "C2'","C1'","N9 ","C8 ","N7 ","C5 ","C6 ","O6 ","N1 ", &
                          "C2 ","N2 ","N3 ","C4 "/)


!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=
  ! Read fragment library
  !=-=-=-=-=-=-=-=-=-=-=-=
   call read_references()
!-----------------------------------------------------------
  call system_clock(tfin, count_rate)
  write(*,'(f10.5,a)') ( 1.0d0*(tfin-t0) )/(count_rate*1.0d0) ,'sec for reading fragment library.' 
!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=
  ! Open input file
  !=-=-=-=-=-=-=-=-=
  if (iargc()==1) then
    call getarg(1,in_title)
    write(*,*) 'Filemane:',in_title
  else
    stop 'Please specify one input file. [ ./DNAbackmap  inp.txt ]'
  end if 
  
  open (inpdev, file=in_title, status='OLD', action='READ', iostat=io_status)
  if (io_status>0) then
    stop 'Cannot open .pdb file.'
  end if

!----------------------------------
  !=-=-=-=-=-=-=-=
  ! Modeling loop
  !=-=-=-=-=-=-=-=

  do
    read(inpdev,'(a)',iostat=io_status) in_line
    if (io_status<0) then
      exit
    end if
!----------------------------------
  !!=-=-=-=-=-=-=-=-=
  !! Read .pdb file
  !!-=-=-=-=-=-=-=-=
    !write(*,*) in_line
    if (in_line(1:8) == 'FILENAME') then
      call in_out_file(in_line, in_pdb, out_pdb)

      write(*,*) 'Input File  : ',trim(adjustl(in_pdb))
      write(*,*) 'Output File : ',trim(adjustl(out_pdb))
      call read_pdb(in_pdb, xyz_pdb, n_3site, phos_type, res,&
                    seq, n_read_model, chain)

!-----------------------------------------------------------
  !!=-=-=-=-=-=-=-=-=-=
  !! Open output file
  !!=-=-=-=-=-=-=-=-=-=

      open (outdev, file=out_pdb, status='REPLACE', &
                    action='WRITE', iostat=io_status)
      if (io_status>0) then
        stop 'Cannot open output file.'
      end if
      atomid=0 ! atom ID in output file

!-----------------------------------------------------------
  !!=-=-=-=-=-=
  !! Modeling
  !!=-=-=-=-=-=


      do i_3site=1,n_3site
        if ((phos_type(1,i_3site)==0).and.(phos_type(2,i_3site)==0)) then
          write(*,*) "@ ", i_3site, 'th and',i_3site+1,"th nucleotide"
          stop 'phosphate missing'

        else if (phos_type(1,i_3site)==0)then ! 5' nucleotide use DS,DB,DP+,and DS+
          read_list=(/3, 4, 5, 6/)
          xyz_4site(:,1:2) = xyz_pdb(:,2:3,i_3site)
          xyz_4site(:,3:4) = xyz_pdb(:,1:2,i_3site+1)
          call end_modeling(read_list, n_read_model(i_3site), &
                            xyz_4site, seq(i_3site), aamodel)
     
        else if (phos_type(2,i_3site)==0)then ! 3' nucleotide use DS-,DP,DS,and DB
          read_list=(/1, 2, 3, 4/)
          xyz_4site(:,1) = xyz_pdb(:,2,i_3site-1)
          xyz_4site(:,2:4) = xyz_pdb(:,1:3,i_3site)
	  !write(*,*)xyz_4site
          call end_modeling(read_list, n_read_model(i_3site), &
                            xyz_4site, seq(i_3site), aamodel)

        else ! the others
          xyz_5site(:,1) = xyz_pdb(:,2,i_3site-1)
          xyz_5site(:,2:4) = xyz_pdb(:,1:3,i_3site)
          xyz_5site(:,5) = xyz_pdb(:,1,i_3site+1)
          call main_modeling(n_read_model(i_3site), xyz_5site,seq(i_3site), aamodel)
        end if

        call output_pdb(seq(i_3site), phos_type(:,i_3site), res(i_3site), &
                        atomid, aamodel, chain(i_3site))

      end do

  !---------------------------------------------------------
      close(outdev)

    end if
  
  end do
  close(inpdev)
!-----------------------------------------------------------
  call system_clock(tfin, count_rate)
  write(*,'(f10.5,a)') ( 1.0d0*(tfin-t0) )/(count_rate*1.0d0) ,'sec for modeling.' 
!-----------------------------------------------------------


  stop
end program DNAbackmap
