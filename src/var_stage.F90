!var_afmfit
!> @brief Contains parameters for L-J stage potential
! vim: tabstop=3 shiftwidth=3

module var_stage
    use const_maxsize
!$  use omp_lib
    implicit none

    type stage_parameter
        ! particles
        integer :: group_id !< group of particles that feels this potential
        real(PREC) :: height          !< position of stage along z coordinate
        real(PREC) :: cutoff_ratio    !< cutoff distance relative to the sigma value
        real(PREC) :: eps             !< epsilon
        real(PREC) :: distance_offset !< offset of distance between stage and particle
    end type stage_parameter

    type(stage_parameter), save :: stage

end module var_stage
