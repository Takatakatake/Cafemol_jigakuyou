!var_afmfit
!> @brief Contains parameters for AFM Fitting Potential
! vim: tabstop=3 shiftwidth=3

module var_afmfit
    use const_maxsize
    implicit none

    type afm_fitting_parameter
        ! image information
        integer    :: n_x
        integer    :: n_y
        integer    :: n_pixel ! == n_x * n_y
        real(PREC) :: pixel_size_x !< size of pixel along x axis
        real(PREC) :: pixel_size_y !< size of pixel along y axis
        real(PREC), allocatable :: x_at(:)     !< x coordinate at i-th pixel
        real(PREC), allocatable :: y_at(:)     !< y coordinate at i-th pixel
        real(PREC), allocatable :: z_ref_at(:) !< reference height information
        real(PREC), allocatable :: z_sim_at(:) !< height in the simulation

        ! forcefield parameters
        real(PREC) :: k            !< overall magnitude
        real(PREC) :: z0           !< parameter to supress numerical error
        real(PREC) :: beta         !< for softmax
        real(PREC) :: rbeta        !< r is for reciprocal, for speedup
        real(PREC) :: sigma_x      !< for gaussian along x coord
        real(PREC) :: sigma_y      !< for gaussian along y coord
        real(PREC) :: rsigma_x     !< 1 / sigma_x
        real(PREC) :: rsigma_y     !< 1 / sigma_y
        real(PREC) :: rsigma_x_sq  !< 1 / sigma_x**2
        real(PREC) :: rsigma_y_sq  !< 1 / sigma_y**2
        real(PREC) :: cutoff_sigma !< cutoff for gaussian
        real(PREC) :: cutoff_exp   !< cutoff for exp(z/beta)
        integer    :: n_max_pixel  !< number of pixels inside of the cutoff

        ! particles
        integer :: group_id !< group of particles that feels this potential

        ! forcefield calculation utilities
        integer :: num_max_threads
        integer, allocatable :: imp2pxl(:)  ! imp -> affecting pixels
        real(PREC), allocatable :: sumexp(:)   ! buffer
        real(PREC), allocatable :: r_sumexp(:) ! buffer
        real(PREC), allocatable :: sumexp_thread_local(:, :)
        real(PREC), allocatable :: z_ref_sim_sum_thloc(:)
        real(PREC), allocatable :: z_sim_sim_sum_thloc(:)
        real(PREC) :: z_ref_sq_sum, r_z_ref_sq_sum ! for speedup
    end type afm_fitting_parameter

    type(afm_fitting_parameter), save :: afm

end module var_afmfit
