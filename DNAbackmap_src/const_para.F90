module const_para

  implicit none

! fragment library
    character(72),parameter :: ref_title = 'TCAG_fragment.txt'
    real(8) ,allocatable, save :: dna_ref(:,:,:,:)

! device number
    integer, parameter :: inpdev  = 10
    integer, parameter :: outdev  = 11
    integer, parameter :: refdev =12
    integer, parameter :: pdbdev  = 13

! parameters for evaluation function
    real(8), parameter :: kpara(5) =  (/0.1d0, 0.75d0,1.0d0,1.0d0,0.1d0/)

! the number of flagments
    integer, parameter :: n_ref_max = 6309
    ! n_ref_list(base type, dataset type)
    ! base type 1:T 2:C 3:A 4:G
    ! dataset type 1:use all structures 2:use structures which form WC pair in PDB
    integer, parameter :: n_ref_list(4,2) = reshape( &
          (/5118,5940,4980,6309,4213,5126,4205,4983/),shape(n_ref_list))

! max characters
    integer, parameter :: maxname = 180
!max nucleotides
    integer, parameter :: mx3site = 2000
! particles' name
    character(3), save :: particle_name(22,4)

! reference coordinates for base-pairing check (xyz, S => B S => B, 1: T => A 2:C => G)
   real(8), parameter :: bpref(3,4,2) = reshape( &
          (/-5.822d0, -4.049d0, 14.707d0, -1.802d0, -2.888d0, 13.857d0, &
             2.053d0,  6.793d0, 12.341d0,  1.243d0,  2.197d0, 13.278d0, &
            -5.752d0, -4.067d0, 48.571d0, -2.035d0, -2.563d0, 47.606d0, &
             2.082d0,  6.729d0, 46.065d0,  0.886d0,  2.157d0, 47.130d0/), &
          shape(bpref))

! threshold of DB-DB distance for Watson-Click pairing 
  real(8), parameter :: d2_threshold = 39.69d0
  real(8), parameter :: rms_threshold  = 0.65d0


! precision for util_bestfit
  integer, parameter :: PREC = 8
end module const_para

