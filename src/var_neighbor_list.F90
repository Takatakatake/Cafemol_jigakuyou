module var_neighbor_list
    use const_maxsize, only: PREC
    implicit none

    type neighbor_list_type
        integer :: capacity, num_coefs, num_pairs
        integer, allocatable :: pairs(:, :) ! (2 or 3, capacity)
        real(PREC), allocatable :: coefs(:, :) ! (num_coefs, capacity)
    end type

    !> parameters for neighboring (generalized born) list
    type(neighbor_list_type), allocatable :: ele_gb_list(:) ! (REPLICA)
    type(neighbor_list_type), allocatable :: pnl_gb_list(:) ! (REPLICA)

    !> parameters for electrostatic
    type(neighbor_list_type), allocatable :: ele_list(:)

    !> parameters for elctrostatic(K computer)
    integer, allocatable, save :: lele_k(:, :)        ! (ncharge, REPLICA)
    integer, allocatable, save :: iele2charge_k(:, :, :) ! (ncharge, ncharge, REPLICA)

    !> parameters for solvation (DNA) potential
    type(neighbor_list_type), allocatable :: solv_list(:)

    !> parameters for hydrophobic interaction
    type(neighbor_list_type), allocatable :: hp_list(:)

    !> parameters for neighboring (general) list
    type(neighbor_list_type), allocatable :: exv_list(:)
    type(neighbor_list_type), allocatable :: bp_at_list(:)
    type(neighbor_list_type), allocatable :: bp_gc_list(:)
    type(neighbor_list_type), allocatable :: bp_mismatch_list(:)
    type(neighbor_list_type), allocatable :: exv_dna_list(:)
    type(neighbor_list_type), allocatable :: exv_dna2_list(:)
    type(neighbor_list_type), allocatable :: exv_wca_list(:)
    type(neighbor_list_type), allocatable :: pro_dna_nonspec_list(:)
    type(neighbor_list_type), allocatable :: pro_dna_pwm_list(:)
    type(neighbor_list_type), allocatable :: exv_ion_list(:)
    type(neighbor_list_type), allocatable :: hyd_ion_list(:)
    type(neighbor_list_type), allocatable :: sasa_list(:)

    private alloc_neighbor_list

contains

    pure subroutine allocate_neighbor_list(list, capacity, num_coefs, stat)
        use var_replica, only: n_replica_mpi
        implicit none

        type(neighbor_list_type), allocatable, intent(inout) :: list(:)
        integer, intent(in) :: capacity
        integer, intent(in), optional :: num_coefs
        integer, intent(out) :: stat
        integer :: i

        allocate (list(n_replica_mpi), stat=stat)
        if (stat /= 0) return

        if (present(num_coefs)) then
            do i = 1, n_replica_mpi
                call alloc_neighbor_list(list(i), capacity, num_coefs, stat)
                if (stat /= 0) return
            end do
        else
            do i = 1, n_replica_mpi
                call alloc_neighbor_list(list(i), capacity, stat=stat)
                if (stat /= 0) return
            end do
        end if
    end subroutine

    pure subroutine alloc_neighbor_list(list, capacity, num_coefs, stat)
        use var_inp, only: inperi
        implicit none

        type(neighbor_list_type), intent(out) :: list
        integer, intent(in) :: capacity
        integer, intent(in), optional :: num_coefs
        integer, intent(out) :: stat

        list%capacity = capacity
        list%num_pairs = 0

        allocate (list%pairs(2 + inperi%n_mirror_index, capacity), stat=stat)
        if (stat /= 0) return
        list%pairs(:, :) = 0

        if (.not. present(num_coefs)) then
            list%num_coefs = 0
        else if (num_coefs == 0) then
            list%num_coefs = 0
        else
            allocate (list%coefs(num_coefs, capacity), stat=stat)
            if (stat /= 0) return
            list%coefs(:, :) = 0.0
            list%num_coefs = num_coefs
        end if
    end subroutine

    pure subroutine reserve_neighbor_list(list, template, stat)
        implicit none

        type(neighbor_list_type), intent(out) :: list
        type(neighbor_list_type), intent(in) :: template
        integer, intent(out) :: stat

        call alloc_neighbor_list(list, template%capacity, template%num_coefs, stat)
    end subroutine

    pure subroutine deallocate_neighbor_list(list)
        use var_replica, only: n_replica_mpi
        implicit none
        type(neighbor_list_type), allocatable, intent(inout) :: list(:)
        integer :: i

        do i = 1, n_replica_mpi
            call dealloc_neighbor_list(list(i))
        end do
        deallocate (list)
    end subroutine

    pure subroutine dealloc_neighbor_list(list)
        implicit none
        type(neighbor_list_type), intent(inout) :: list

        if (allocated(list%pairs)) deallocate (list%pairs)
        if (allocated(list%coefs)) deallocate (list%coefs)
    end subroutine

    pure subroutine push_pair(list, imp, jmp)
        implicit none
        type(neighbor_list_type), intent(inout) :: list
        integer, intent(in) :: imp, jmp
        list%num_pairs = list%num_pairs + 1
        list%pairs(1:2, list%num_pairs) = [imp, jmp]
    end subroutine

    pure subroutine clear_list(list)
        implicit none
        type(neighbor_list_type), intent(inout) :: list
        list%num_pairs = 0
    end subroutine

end module
