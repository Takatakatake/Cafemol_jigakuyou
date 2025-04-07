subroutine check_bp(xyz_pdb, seq, chain, n_3site, n_read_model)

  use const_para  

  implicit none
  !-----------------------------------------------------------
  real(8), intent(in) :: xyz_pdb(3,3,mx3site)
  integer, intent(in) :: seq(mx3site), n_3site
  character(1), intent(in) :: chain(mx3site)
  integer, intent(out) :: n_read_model(mx3site)

!-----------------------------------------------------------
  integer :: nbase(4), pair_pu(mx3site, 2), siteid_conv(mx3site, 4)
  real(8) :: xyz_by_seq(3,mx3site,2,4), pair_rms(mx3site, 2),dp, &
             superpo(3,4), after_align(3,4), rmsd 
  integer :: i_3site, tmp_i, i, j, ipu, ipy
 

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=
  ! Initialize variables
  !=-=-=-=-=-=-=-=-=-=-=-=

  do i_3site = 1, n_3site
    n_read_model(i_3site)=n_ref_list(seq(i_3site),1)  
  end do

  nbase=0
  xyz_by_seq = 0.0d0
  pair_pu = 0
  pair_rms = 2.0d0

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Divide the coordinates
  !=-=-=-=-=-=-=-=-=-=-=-=-=

  do i_3site = 1,n_3site
    tmp_i=seq(i_3site)
    nbase(tmp_i) = nbase(tmp_i) +1
    xyz_by_seq(:,nbase(tmp_i),1:2, tmp_i) = xyz_pdb(:,2:3,i_3site)
    siteid_conv(nbase(tmp_i), tmp_i) = i_3site
  end do

!-----------------------------------------------------------
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Check Watson-Click pairing
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  do i = 3,4 ! 1:T, 2: C, 3: A, 4: G
    tmp_i=i-2

    do ipu = 1, nbase(i)
      do ipy = 1, nbase(tmp_i)
	! avoid n th and n+1 th base pair
        if ( ( chain(siteid_conv(ipu,i)) /= chain(siteid_conv(ipy,tmp_i)) ).or. &
             ( abs( siteid_conv(ipu,i) - siteid_conv(ipy,tmp_i) ) > 1 ) ) then

          dp = dot_product(xyz_by_seq(:,ipu,2,i)-xyz_by_seq(:,ipy,2,tmp_i), &
                           xyz_by_seq(:,ipu,2,i)-xyz_by_seq(:,ipy,2,tmp_i) )
          ! square distance threshold 
          if ( dp < d2_threshold ) then
            ! calculate RMSD of 4 points
            superpo(:,1:2) = xyz_by_seq(:,ipy, 1:2,tmp_i)
            superpo(:,3:4) = xyz_by_seq(:, ipu, 1:2, i) 
            call util_bestfit(bpref(:,:,tmp_i), 4, superpo, 4, 4, after_align, &
                          (/1,2,3,4/), (/1,2,3,4/), rmsd) 
            if ( rmsd < pair_rms(ipy, tmp_i) ) then
               pair_rms(ipy, tmp_i) = rmsd
	       pair_pu(ipy ,tmp_i) = ipu
            end if
  
          end if

        end if
      end do      
    end do

  end do

  do i=1,2
    do ipy = 1, nbase(i)

      if ( pair_rms(ipy, i) < rms_threshold ) then
        !write(*,*) 'd',siteid_conv(ipy,i), siteid_conv(pair_pu(ipy, i), i+2)
        !write(*,*) 'r',siteid_conv(ipy,i), siteid_conv(pair_pu(ipy, i), i+2)
        !write(*,*) pair_rms(ipy, i)
        n_read_model(siteid_conv(ipy, i)) = n_ref_list(i,2)
        n_read_model(siteid_conv(pair_pu(ipy, i), i+2)) = n_ref_list(i+2,2)
      end if

    end do
  end do
  !do i=1,n_3site
    !write(*,*) i,n_read_model(i)
  !end do
  return

end subroutine check_bp

